!> Implements sponge regions in SIS2 model
module SIS_sponge


!! Module to read in Time of the ice fields, ice concentration, ice thickness, relaxation time scale
! TODO: Clean not needed modules
use MOM_coms,          only : sum_across_PEs
use MOM_coms,          only : PE_here   !! debugging
use MOM_unit_scaling,  only : unit_scale_type
use ice_grid,          only : ice_grid_type

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : file_exists, MOM_read_data, slasher
                            
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_sum_output,    only : SIS_sum_out_CS, write_ice_statistics! , SIS_sum_output_init
use SIS_types,         only : ice_state_type, IST_chksum, IST_bounds_check, total_sfc_flux_type         
use SIS_types,         only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS2_ice_thm,      only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,      only : enth_from_TS, Temp_from_En_S, enthalpy_liquid, calculate_T_freeze

implicit none; private

#include <SIS2_memory.h>

public initialize_isponge, apply_isponge, set_up_isponge_field, SIS_sponge_end
public adjust_IOfluxes_isponge

!> A structure for creating arrays of pointers to 3D arrays
type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3D array [various]
  character(len=15)               :: fld_name    !< Name of the ice field being relaxed
end type p3d
!> A structure for creating arrays of pointers to 2D arrays
type, public :: p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array [various]
end type p2d
!> A structure for 2D arrays
type, public :: f2d
  real, allocatable, dimension(:,:) :: fld
end type f2d
 
!> This control structure holds memory and parameters for the SIS_sponge module
type, public :: isponge_CS ; private
  logical, public :: use_isponge = .false.  !< If true, ice tracer fields may be relaxed somewhere in the domain
  integer, public :: itest, jtest    ! debugging, output at idices on PE
  integer :: num_col    !< The number of sponge points within the computational domain.
  integer, public :: fldno = 0  !< The number of fields which have already been
                        !! registered by calls to set_up_sponge_field
  integer, pointer :: col_i(:) => NULL() !< Array of the i-indicies of each of the columns being damped.
  integer, pointer :: col_j(:) => NULL() !< Array of the j-indicies of each of the columns being damped.
  real, pointer :: Iresttime_col(:) => NULL() !< The inverse restoring time of each column [T-1 ~> s-1].
  type(p3d) :: var(MAX_FIELDS_RLX_)     !< Pointers to the fields that will be relaxed
  type(p2d) :: Ref_val(MAX_FIELDS_RLX_) !< Relaxation values - The values to which the fields are 
                                        ! relaxed (linear_index, ice_cat)
  type(f2d) :: Old_val(MAX_FIELDS_RLX_) !< Keep old values of relaxed fields prior to relaxation
                                        ! to estimate heat/ salt fluxes to compensate ice changes
  real, allocatable, dimension(:) :: &
      Enth_out_ocn_old, & ! Negative of the enthalpy extracted from ice by water fluxes to the ocean [Q R Z ~> J m-2]
      flux_salt_old  ! The flux of salt out of the ocean [1e3 S R Z T-1 ~> kgSalt m-2 s-1]
end type isponge_CS

contains

!> This subroutine determines the number of points which are within ice sponges in
!! this computational domain.  Only points that have positive values of
!! Iresttime and which mask2dT indicates are ocean points are included in the
!! sponges.  
!subroutine initialize_isponge(Iresttime, IST, G, IG, param_file, CS)
subroutine initialize_isponge(param_file, Iresttime, G, IG, CS, itest, jtest)
!  type(ice_state_type),    intent(in) :: IST        !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: Iresttime  !< The inverse of the restoring time [T-1 ~> s-1].
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  integer, intent(in) :: itest, jtest

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  logical :: use_isponge
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat
  integer :: col, total_isponge_cols

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
!  write(mesg,'("SIS_sponge: isc/iec=",I4,"/",I4," jsc/jec=",I4,"/",I4," icencat=",I3)') &
!       isc, iec, jsc, jec, ncat 
!  call SIS_mesg(mesg)

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_sponge: initialize_isponge called with "// &
                            "an associated control structure.")
    return
  endif

  call SIS_mesg("SIS_sponge: Checking use SIS sponge.")
! Set default, read and log parameters
! get_param (procedure --> get_param_logical) - checks if variable is set true in the param_file
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SIS_SPONGE", use_isponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (use_isponge) then
    call SIS_mesg("SIS_sponge: use_isponge: True")
  else
    call SIS_mesg("SIS_sponge: use_isponge: False")
  endif
  if (.not.use_isponge) return
  allocate(CS)

  CS%use_isponge = use_isponge
  CS%itest = itest
  CS%jtest = jtest
  CS%num_col = 0 ; CS%fldno = 0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo

  write(mesg,'("SIS_sponge: num_col=",I8)') CS%num_col
  call SIS_mesg(mesg)

  if (CS%num_col > 0) then
    allocate(CS%Iresttime_col(CS%num_col), source=0.0)
    allocate(CS%col_i(CS%num_col), source=0)
    allocate(CS%col_j(CS%num_col), source=0)
    allocate(CS%Enth_out_ocn_old(CS%num_col), source=0.0)
    allocate(CS%flux_salt_old(CS%num_col), source=0.0)

    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo

  endif

!  call SIS_mesg("Calling sum_across_PEs")
  total_isponge_cols = CS%num_col
  call sum_across_PEs(total_isponge_cols)

  write(mesg,'("SIS_sponge: total isponge cols=",I8)') total_isponge_cols
  call SIS_mesg(mesg)
  
  call log_param(param_file, mdl, "!Total isponge columns", total_isponge_cols, &
                 "The total number of ice columns where sponges are applied.")

end subroutine initialize_isponge

!> This subroutine stores the reference profile for the SIS variable whose
!! address is given by f_ptr. Reference profile = values towards which the
!! SIS field is being relaxed to. 
!! Current version assumes 2D fields only (+ 1D for categories) such as Hice
subroutine set_up_isponge_field(sp_val, f_ptr, G, IG, CS, ncat, fld_name)
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  integer,                 intent(in) :: ncat
  real, dimension(SZI_(G), SZJ_(G), ncat), &
                           intent(in) :: sp_val     !< The reference profiles of the quantity 
                                                    !! being registered [various]
  real, dimension(SZI_(G), SZJ_(G), ncat), &
                   target, intent(in) :: f_ptr      !< a pointer to the field which will be relaxed [various]
  character(*),             intent(in) :: fld_name   !< Name of the relaxed field
  type(isponge_CS),     pointer       :: CS         !< A pointer to the control structure for this module that
                                                    !! is set by a previous call to initialize_sponge.

  integer :: j, k, col
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1

  if (CS%fldno > MAX_FIELDS_RLX_) then
    write(mesg,'("Increase MAX_FIELDS_RLX_ to at least ",I3," in SIS_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call SIS_error(FATAL,"set_up_isponge_field: "//mesg)
  endif

  allocate(CS%Ref_val(CS%fldno)%p(CS%num_col,IG%CatIce), source=0.0)
  allocate(CS%Old_val(CS%fldno)%fld(CS%num_col,IG%CatIce), source=0.0)
  do col=1,CS%num_col
    do k=1,IG%CatIce
      CS%Ref_val(CS%fldno)%p(col,k) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr
  CS%var(CS%fldno)%fld_name = fld_name

end subroutine set_up_isponge_field

!> This subroutine applies relaxation ("damping") to ice thickness (by categories) and ice concentration
!! tracers for every column where the relaxation time scale > 0.
subroutine apply_isponge(dt_slow, CS, IG, IOF, IST, US, OSS)
  real,                      intent(in)  :: dt_slow   !< The amount of time covered by this call [T ~> s].
  type(isponge_CS),          pointer     :: CS   !< A pointer that is set to point to the ice sponge control
                                                 !! structure for this module
  type(ice_grid_type),       intent(in)  :: IG    !< The sea-ice specific grid type
  type(ice_ocean_flux_type), intent(inout)  :: IOF  !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(ice_state_type),   intent(inout)  :: IST  !< A type describing the state of the sea ice
  type(unit_scale_type),     intent(in)  :: US   !< A structure with unit conversion factors
  type(ocean_sfc_state_type), intent(in) :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.

  ! Local variables
!  real :: H_rescale_ice, H_rescale_snow
  real :: damp         ! The timestep times the local damping coefficient [nondim].
  real :: I1pdamp      ! I1pdamp is 1/(1 + damp). [nondim]
  real :: p_old        ! debugging  
  real :: dt           ! time step in s
  real :: s_ice_bulk   ! ice bulk S for filling S values in the newly ceated ice 
  real, allocatable :: sice(:), tfi(:)
  real :: I_Nk, enth_ice, Tfrz, coeff, enth_Tfrz
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  integer :: c, i, j, k, l, m, NkIce
  integer :: current_pe
  real    :: Idt_slow    ! The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new
  real    :: iconc_tot, dlt_ice_tot
  real    :: dlt_iconc, dlt_ithk, enthalpy_ocn, enthalpy_ocn_tfrz
  real    :: dlt_salt, dlt_heat, dlt_water, dlt_snow
  real    :: dlt_ice           ! total change of ice due to conc and thickness relaxation
  real    :: ice_salin         ! average ice column S gSalt kg-1 
  real    :: water_ice_ocn, heat_ice_ocn, salt_ice_ocn
  real    :: enthalpy_ocn0
  logical :: f_debug

  f_debug = .true.
  NkIce = IG%NkIce  
  I_Nk  = 1. / NkIce
  dt = dt_slow*US%T_to_s
  s_ice_bulk = 3.0*US%ppt_to_S

  if (CS%num_col == 0) return
    
  allocate(sice(NkIce), tfi(NkIce), source=-999.)

  current_pe = PE_here()
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)             
    damp = dt * CS%Iresttime_col(c); I1pdamp = 1.0 / (1.0 + damp)

    do k=1,IG%CatIce
      do m=1,CS%fldno
!        p_old = CS%var(m)%p(i,j,k)  ! debugging
        CS%Old_val(m)%fld(c,k) = CS%var(m)%p(i,j,k)  

        CS%var(m)%p(i,j,k) = I1pdamp * &
           (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(c,k)*damp)

        if (i==CS%itest .and. j==CS%jtest .and. f_debug) then
          if (CS%var(m)%fld_name(1:5) == 'mHice') then
            coeff = US%RZ_to_kg_m2
          else
            coeff = 1.0
          endif

          write(mesg,'(A8," k=",I2," old:=",D12.4," new=",D12.4,&
                      " 1/tau=",D12.3," refval=",D12.3)') &
            CS%var(m)%fld_name(1:8), k, CS%Old_val(m)%fld(c,k)*coeff, &
            CS%var(m)%p(i,j,k)*coeff, &
            CS%Iresttime_col(c), CS%Ref_val(m)%p(c,k)*coeff
          write(*,'(A)') trim(mesg)
        endif
      enddo
      ! Adjust enth and S in the newly formed ice if needed:
      ! Note ice enthalpy < 0
      do l=1,NkIce  
        if (IST%sal_ice(i,j,k,l) < s_ice_bulk) &
          IST%sal_ice(i,j,k,l) = s_ice_bulk
        ! Enth should be at least enth(T_freez)
        sice(l) = IST%sal_ice(i,j,k,l)
      enddo

!      write(mesg,'("SIS_sponge: calling calculate_T_Freeze")')
!      call SIS_mesg(mesg)
      call calculate_T_Freeze(sice, tfi, IST%ITV)

      do l=1,NkIce
        enth_ice = IST%enth_ice(i,j,k,l)
        enth_Tfrz = enth_from_TS(tfi(l), sice(l), IST%ITV)

        if (enth_ice > enth_Tfrz) &
          IST%enth_ice(i,j,k,l) = enth_Tfrz

!        if (i==CS%itest .and. j==CS%jtest .and. f_debug) then
!          write(mesg,'("apply_isp:  k=",I2," l=",I2," sice=",F8.3," Tfrz=",F9.3&
!                      " enth_orig=",D12.3," enth_final=",D12.3)') &
!                 k, l, sice(l)*US%S_to_ppt, tfi(l)*US%C_to_degC, enth_ice*US%Q_to_J_kg, &
!                 IST%enth_ice(i,j,k,l)*US%Q_to_J_kg
!          write(*,'(A)') trim(mesg)
!        endif

      enddo
    enddo  ! CatIce

  ! Keep ocean T near Tfrz to avoid rapid ice melt if dlt cice > 0
  ! and keep it slightly above if dlt cice <0 and new cice = 0
  ! need to know dz of the upper ocean layer - keep for debugging now but delte later
    iconc_tot = 0.0
    dlt_ice_tot = 0.0
    do k=1,IG%CatIce
      do m=1,CS%fldno
        fld_name = CS%var(m)%fld_name
        select case (trim(fld_name))
          case('mH_ice')
            ithk_old = CS%Old_val(m)%fld(c,k)
            ithk_new = CS%var(m)%p(i,j,k)
          case('part_size')
            iconc_old = CS%Old_val(m)%fld(c,k)
            iconc_new = CS%var(m)%p(i,j,k)
            iconc_tot = iconc_tot + iconc_new
        end select
      enddo

      dlt_iconc = iconc_new - iconc_old
      dlt_ithk  = ithk_new - ithk_old           ! ice mass change, kg m-2
      dlt_ice = ithk_new*iconc_new - ithk_old*iconc_old
      dlt_ice_tot = dlt_ice_tot + dlt_ice
    enddo

    enthalpy_ocn = enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
    enthalpy_ocn_tfrz = enthalpy_liquid_freeze(OSS%s_surf(i,j), IST%ITV) 
    enthalpy_ocn0 = enthalpy_liquid(0.0, OSS%s_surf(i,j), IST%ITV)
! need to know dz of the upper grid cell in MOM to do this:
!    if (iconc_tot > 1.e-30 .and. dlt_ithk > 0.0 .and. enthalpy_ocn > enthalpy_ocn_tfrz) then
!      dlt_heat = enthalpy_ocn - enthalpy_ocn_tfrz 
!      IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) - dlt_heat*rho_water*dz_ocean
!    endif

    if (i==CS%itest .and. j==CS%jtest .and. f_debug) then  
      write(mesg, '("enthalpy_ocn=",D12.4," enthalpy_tfrz=",D12.4," iconc=",D12.4," dltIce=",D12.4,&
            " enth0=",D12.4," sst=",F6.2)') &
           enthalpy_ocn*US%Q_to_J_kg, enthalpy_ocn_tfrz*US%Q_to_J_kg, &
           iconc_tot, dlt_ice_tot*US%RZ_to_kg_m2,&
           enthalpy_ocn0*US%Q_to_J_kg, OSS%SST_C(i,j)*US%C_to_degC
      write(*,'(A)') trim(mesg)
    endif
! Above is Probably not needed 

    CS%Enth_out_ocn_old(c) = IOF%Enth_Mass_out_ocn(i,j)
    CS%flux_salt_old(c) = IOF%flux_salt(i,j)
  
  enddo

  if (allocated(sice)) deallocate(sice)
  if (allocated(tfi)) deallocate(tfi)

end subroutine apply_isponge

!> Compute surplus S and heat fluxes due to changes in sea ice thickness and concentrations
!! Ice changes due to relaxation are not conservative ! i.e. heat and salt fluxes should be 0
!! To do this,Need to add/subtract surplus heat and salt fluxes with opposite signs to cancel 
!! fluxes estimated later in SIS slow thermodyn code
!! where ice-ocean fluxes ice are evaluated based on ice changes
! enthalpy of water - E required to warm water from Tfrz to SST or vice versa for freezing
! IOF%sal_ice:   The salinity of the sea ice, units="g/kg", conversion=US%S_to_ppt
! IOF%flux_salt: The flux of salt out of the ocean [1e3 S R Z T-1 ~> kgSalt m-2 s-1],
!                scale=US%S_to_ppt*US%RZ_T_to_kg_m2s
! IST%enth_ice: enthalpy of ice in each cat. and fractional thickn. layer, [Q ~> J kg-1]  scale=US%Q_to_J_kg
! IOF%Enth_Mass_in_ocn: The enthalpy introduced to the ice by water fluxes from the ocean [Q R Z ~> J m-2].
! IOF%Enth_Mass_out_ocn: Negative of the enthalpy extracted from the ice by water fluxes to
!                          the ocean [Q R Z ~> J m-2], scale=US%QRZ_T_to_W_m2*US%T_to_s
!> Compute surplus S and heat fluxes due to changes in sea ice thickness and concentrations
!! Ice changes due to relaxation are not conservative ! i.e. heat and salt fluxes should be 0
!! To do this,Need to add/subtract surplus heat and salt fluxes with opposite signs to cancel 
!! fluxes estimated later in SIS slow thermodyn code
!! where ice-ocean fluxes ice are evaluated based on ice changes
subroutine adjust_IOfluxes_isponge(dt_slow, CS, IG, IST, IOF, OSS, US)
  real,                       intent(in)  :: dt_slow !< The thermodynamic step [T ~> s].
  type(isponge_CS),           pointer     :: CS    !< A pointer that is set to point to the ice sponge control
                                                !! structure for this module
  type(ice_grid_type),        intent(in)  :: IG    !< The sea-ice specific grid type
  type(ice_state_type),       intent(in)  :: IST   !< A type describing the state of the sea ice
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model  
  type(ocean_sfc_state_type), intent(inout) :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  integer :: c, i, j, k, m, l, NkIce 
  integer :: current_pe                                  ! debugging
  real    :: Idt_slow    ! The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new
  real    :: dlt_iconc, dlt_ithk, I_Nk
  real    :: dlt_salt, &    ! g/kg * kg/m2 salt change due to relax over a cat
             dlt_heat, &    ! J/m2=enthalpy * kg/m2 change due to relaxn over a cat
             dlt_water, dlt_snow
  real    :: dlt_ice           ! change of ice mass due to conc and thickness relax changes, over a cat
  real    :: ice_salin         ! average ice column S gSalt kg-1 over a cat
  real    :: water_ice_ocn, heat_ice_ocn, salt_ice_ocn
  real    :: dlt_heat_tot, dlt_salt_tot, dlt_ice_tot !> heat, salt, iceM changes over all cat

!  dt_slow  = dt_slow*US%T_to_s
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0 / dt_slow
  NkIce = IG%NkIce
  I_Nk   = 1. / NkIce
  ! Estimate ice volume change by categories:
  if (CS%num_col == 0) return
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)
    dlt_heat_tot = 0.0 ; dlt_ice_tot = 0.0 ; dlt_salt_tot = 0.0
    do k=1,IG%CatIce
      do m=1,CS%fldno
        fld_name = CS%var(m)%fld_name
        select case (trim(fld_name))
          case('mH_ice')
            ithk_old = CS%Old_val(m)%fld(c,k)
            ithk_new = CS%var(m)%p(i,j,k)
          case('part_size')
            iconc_old = CS%Old_val(m)%fld(c,k)
            iconc_new = CS%var(m)%p(i,j,k)
        end select
      enddo

      ! Estimate ice change during the relaxation
      ! Negative change - "melt", positive - "growth"
      ! melt results in negative ice-ocean heat and salt fluxes
      ! To cancel relaxation-caused fluxes, need to add this amount
      ! that will be subtracted later, i.e. net flux = 0
      ! total ice change = dlt_conc*hice_old + dlt_hice*conc_new
      ! total snow change = dlt_conc*hsnow, because snow thickness is not changed 
      ! during relaxation
      ! enthalpy of water - E required to warm water from Tfrz to SST or vice versa for freezing

      dlt_iconc = iconc_new - iconc_old
      dlt_ithk  = ithk_new - ithk_old           ! ice mass change, kg m-2
      dlt_ice = ithk_new*iconc_new - ithk_old*iconc_old

      dlt_salt = 0.0 ; dlt_heat = 0.0 ; dlt_water = 0.0 
      ice_salin = 0.0
      ! enthalpy of ice/snow
      ! enthalpy to ocean = Negative of the enthalpy extracted from the ice 
      do l=1,NkIce 
        dlt_salt = dlt_salt + (dlt_ice * I_Nk) * IST%sal_ice(i,j,k,l) !> S change due to relax kg/m2*g/kg = g/m2
        dlt_heat = dlt_heat + (dlt_ice * I_Nk) * IST%enth_ice(i,j,k,l) ! enth<0 when dh>0
        ice_salin = ice_salin + IST%sal_ice(i,j,k,l) * I_Nk
      enddo 
      dlt_snow = dlt_iconc * IST%mH_snow(i,j,k)  ! kg m-2
      dlt_heat = dlt_heat + dlt_snow * IST%enth_snow(i,j,k,1)  ! kg m-2 * J kg-1 => J m-2

      dlt_ice_tot  = dlt_ice_tot + dlt_ice
      dlt_heat_tot = dlt_heat_tot + dlt_heat
      dlt_salt_tot = dlt_salt_tot + dlt_salt

      !> Total water dumped to / frozen from the ocean and associated S content per m2
      !! salt flux to ocean  <0 for ice melting, >0 for freezing
      water_ice_ocn = 0.0 ; heat_ice_ocn = 0.0 ; salt_ice_ocn = 0.0
      water_ice_ocn = dlt_ice + dlt_snow            ! kg m-2
      if (dlt_snow >= 0.0) then
      !> No salt flux from the ocean for new snow
        salt_ice_ocn = dlt_ice * (OSS%s_surf(i,j) - ice_salin)
      else
        salt_ice_ocn = dlt_ice * (OSS%s_surf(i,j) - ice_salin) + dlt_snow * OSS%s_surf(i,j) 
      endif

!> Debug
      if (i==CS%itest .and. j==CS%jtest) then
        write(mesg, '("cat=",I," iconc_new=",D12.3," iconc_old=",D12.3," ithk_new=",D12.3," ithk_old=",D12.3)') &
             k,iconc_new, iconc_old, ithk_new*US%RZ_to_kg_m2, ithk_old*US%RZ_to_kg_m2
        write(*,'(A)') trim(mesg)
        write(mesg, '("adjust_IO: dlt_iconc=",D12.4," dlt_ithk=",D12.4," dlt_ice=",D12.4, &
              " dlt_salt=",D12.4," dlt_heat=",D12.4," ice_salin=",D12.4," dlt_snow=",D12.4)') &
              dlt_iconc, dlt_ithk*US%RZ_to_kg_m2, dlt_ice*US%RZ_to_kg_m2, dlt_salt*US%S_to_ppt,&
              dlt_heat*US%Q_to_J_kg, ice_salin*US%S_to_ppt, dlt_snow*US%RZ_to_kg_m2
        write(*,'(A)') trim(mesg)
!        write(mesg, '("flux_salt=",D12.4," dltSFlux=",D12.4,&
!              " Enth_Mass=",D12.4," dltEnthFlux=",D12.4)') &
!          IOF%flux_salt(i,j)*US%S_to_ppt*US%RZ_T_to_kg_m2s, &
!          salt_ice_ocn*(0.001*Idt_slow)*US%RZ_T_to_kg_m2s,  &
!          IOF%Enth_Mass_out_ocn(i,j)*US%QRZ_T_to_W_m2*US%T_to_s, &
!          dlt_heat*US%QRZ_T_to_W_m2*US%T_to_s
        write(mesg, '("dlt_ice_tot=",D12.4," dlt_heat_tot=",D12.4," dlt_salt_tot=",D12.4)') &
             dlt_ice_tot*US%RZ_to_kg_m2, dlt_heat_tot*US%QRZ_T_to_W_m2*US%T_to_s, &
             dlt_salt_tot*0.001*US%RZ_to_kg_m2
        write(*,'(A)') trim(mesg)

      endif
!!  Debug ===
    enddo

!> Save fluxes with opposite sign to cancel out added fluxes later in the code
!! so that for positive S flux, it is subtracted here not added
!! Enth_Mass_out_ocn is negative of the enthalpy*mass extracted from the ice to the ocean
!! here, add same enthalpy*mass to offset the negative change added somewhere later
!! For flux_salt, Note the conversion here from g m-2 to kg m-2 s-1.
!      IOF%flux_salt(i,j) = IOF%flux_salt(i,j) - salt_ice_ocn * (0.001*Idt_slow) 
    IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) + dlt_heat_tot

  enddo

end subroutine adjust_IOfluxes_isponge

subroutine adjust_IOfluxes_isponge_OLD(dt_slow, CS, IG, IST, IOF, OSS, US)
  real,                       intent(in)  :: dt_slow !< The thermodynamic step [T ~> s].
  type(isponge_CS),           pointer     :: CS    !< A pointer that is set to point to the ice sponge control
                                                !! structure for this module
  type(ice_grid_type),        intent(in)  :: IG    !< The sea-ice specific grid type
  type(ice_state_type),       intent(in)  :: IST   !< A type describing the state of the sea ice
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model  
  type(ocean_sfc_state_type), intent(inout) :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  integer :: c, i, j, k, m, l, NkIce 
  integer :: current_pe                                  ! debugging
  real    :: Idt_slow    ! The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new
  real    :: heat_old, fsalt_old, dlt_fsalt, dlt_heat
  real    :: dlt_iconc, dlt_ithk, I_Nk
  real    :: dlt_salt, dlt_water, dlt_snow
  real    :: dlt_ice           ! total change of ice due to conc and thickness relaxation
  real    :: ice_salin         ! average ice column S gSalt kg-1 
  real    :: water_ice_ocn, heat_ice_ocn, salt_ice_ocn

  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0 / dt_slow
  NkIce = IG%NkIce
  I_Nk   = 1. / NkIce
  ! Estimate ice volume change by categories:
  if (CS%num_col == 0) return
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)
    heat_old  = CS%Enth_out_ocn_old(c)
    fsalt_old = CS%flux_salt_old(c) 
    dlt_heat  = IOF%Enth_Mass_out_ocn(i,j) - heat_old
    dlt_fsalt = IOF%flux_salt(i,j) - fsalt_old 

!> Debug
    do k=1,IG%CatIce
      do m=1,CS%fldno
        fld_name = CS%var(m)%fld_name
        select case (trim(fld_name))
          case('mH_ice')
            ithk_old = CS%Old_val(m)%fld(c,k)
            ithk_new = CS%var(m)%p(i,j,k)
          case('part_size')
            iconc_old = CS%Old_val(m)%fld(c,k)
            iconc_new = CS%var(m)%p(i,j,k)
        end select
      enddo
      if (i==CS%itest .and. j==CS%jtest) then
        write(mesg, '("adjust_IO: cat=",I," iconc_new=",D12.3," iconc_old=",D12.3,&
              " ithk_new=",D12.3," ithk_old=",D12.3)') &
             k,iconc_new, iconc_old, ithk_new*US%RZ_to_kg_m2, ithk_old*US%RZ_to_kg_m2
        write(*,'(A)') trim(mesg)
      endif
    enddo
    if (i==CS%itest .and. j==CS%jtest) then
      write(mesg, '(">>enth_old=",D12.3," enth_new=",D12.3," dlt_heat=",D12.3)') &
            heat_old*US%QRZ_T_to_W_m2*US%T_to_s, &
            IOF%Enth_Mass_out_ocn(i,j)*US%QRZ_T_to_W_m2*US%T_to_s, &
            dlt_heat*US%QRZ_T_to_W_m2*US%T_to_s 
      write(*,'(A)') trim(mesg)
      write(mesg, '(">>fsalt_old=",D12.3," fsalt_new=",D12.3," dlt_fsalt=",D12.3)')&
            fsalt_old*US%RZ_T_to_kg_m2s, IOF%flux_salt(i,j)*US%RZ_T_to_kg_m2s, &
            dlt_fsalt*US%RZ_T_to_kg_m2s
      write(*,'(A)') trim(mesg)
    endif

!> Save fluxes with opposite sign to cancel out added fluxes later in the code
!! so that for positive S flux, it is subtracted here not added
!! same for enthalpy
    IOF%flux_salt(i,j) = IOF%flux_salt(i,j) - dlt_fsalt
    IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) - dlt_heat

  enddo

end subroutine adjust_IOfluxes_isponge_OLD


!> Deallocate memory associated with the SIS_optics module
subroutine SIS_sponge_end(CS)
  type(isponge_CS), pointer :: CS !< The ice sponge control structure that is deallocated here
  
  deallocate(CS) 
                 
end subroutine SIS_sponge_end

end module SIS_sponge 

