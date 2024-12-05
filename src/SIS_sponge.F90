!> Implements relaxation regions in SIS2 model
!> mainly for sea ice concentration (partial area) and thickness by categories
!> This is a quick fix for the problem caused by SIS closed boundaries 
!> unrealistic ice concentrations and thicknesses are simulated when sea ice is pulled off / piled up along the 
!> closed boundaries
!>
!> Dmitry Dukhovskoy NOAA OAR PSL 2024
!> 
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
use SIS_optics,        only : VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF    ! debugging only, delete later

implicit none; private

#include <SIS2_memory.h>

public initialize_isponge, apply_isponge, set_up_isponge_field, SIS_sponge_end
public adjust_IOfluxes_isponge
public check_IOF, check_FIA

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
!> A structure for 3D arrays
type, public :: f3d
  real, allocatable, dimension(:,:,:) :: fld3
end type f3d
 
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
                                        ! to estimate heat/ salt fluxes to compensate ice changes <-- These are 
                                        ! not needed, will need to get rid off
  real, allocatable, dimension(:,:,:) :: & 
        dlt_enth_ice             ! Change of ice enthalpy during relaxation by ice thickn, cat
!  real, allocatable, dimension(:) :: &
!      Enth_out_ocn_old, & ! Negative of the enthalpy extracted from ice by water fluxes to the ocean [Q R Z ~> J m-2]
!      flux_salt_old  ! The flux of salt out of the ocean [1e3 S R Z T-1 ~> kgSalt m-2 s-1]
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
    allocate(CS%dlt_enth_ice(CS%num_col,IG%CatIce,IG%NkIce), source=0.0) ! <-- this is probably not needed
!    allocate(CS%Enth_out_ocn_old(CS%num_col), source=0.0)
!    allocate(CS%flux_salt_old(CS%num_col), source=0.0)

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
  real    :: dlt_enth
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
        sice(l) = IST%sal_ice(i,j,k,l)
      enddo

!      write(mesg,'("SIS_sponge: calling calculate_T_Freeze")')
!      call SIS_mesg(mesg)
! Enth should be at least enth(T_freez)
! Make ice T below T frz and/or keep at ocean SST if it is < ice Tfrz
! to prevent rapid ice melt in the relaxation zone
      call calculate_T_Freeze(sice, tfi, IST%ITV)
      tfi = min(tfi-0.1*US%degC_to_C, OSS%SST_C(i,j)*US%degC_to_C)

      do l=1,NkIce
        enth_ice = IST%enth_ice(i,j,k,l)
        enth_Tfrz = enth_from_TS(tfi(l), sice(l), IST%ITV)

        dlt_enth = 0.0
        if (enth_ice > enth_Tfrz) then
!          dlt_enth = enth_Tfrz - IST%enth_ice(i,j,k,l)
          IST%enth_ice(i,j,k,l) = enth_Tfrz
        endif
!        if (i==CS%itest .and. j==CS%jtest .and. f_debug) then
!          write(mesg,'("apply_isp:  k=",I2," l=",I2," sice=",F8.3," iceTfrz=",F9.3&
!                      " enth_orig=",D12.3," enth_final=",D12.3)') &
!                 k, l, sice(l)*US%S_to_ppt, tfi(l)*US%C_to_degC, enth_ice*US%Q_to_J_kg, &
!                 IST%enth_ice(i,j,k,l)*US%Q_to_J_kg
!          write(*,'(A)') trim(mesg)
!        endif
!        CS%dlt_enth_ice(c,k,l) = dlt_enth
      enddo
    enddo  ! CatIce

!    iconc_tot = 0.0
!    dlt_ice_tot = 0.0
!    do k=1,IG%CatIce
!      do m=1,CS%fldno
!        fld_name = CS%var(m)%fld_name
!        select case (trim(fld_name))
!          case('mH_ice')
!            ithk_old = CS%Old_val(m)%fld(c,k)
!            ithk_new = CS%var(m)%p(i,j,k)
!          case('part_size')
!            iconc_old = CS%Old_val(m)%fld(c,k)
!            iconc_new = CS%var(m)%p(i,j,k)
!            iconc_tot = iconc_tot + iconc_new
!        end select
!      enddo
!
!      dlt_iconc = iconc_new - iconc_old
!      dlt_ithk  = ithk_new - ithk_old           ! ice mass change, kg m-2
!      dlt_ice = ithk_new*iconc_new - ithk_old*iconc_old
!      dlt_ice_tot = dlt_ice_tot + dlt_ice
!    enddo
!
!    enthalpy_ocn = enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
!    enthalpy_ocn_tfrz = enthalpy_liquid_freeze(OSS%s_surf(i,j), IST%ITV) 
!    enthalpy_ocn0 = enthalpy_liquid(0.0, OSS%s_surf(i,j), IST%ITV)
!
!    if (i==CS%itest .and. j==CS%jtest .and. f_debug) then  
!      write(mesg, '("enthalpy_ocn=",D12.4," enthalpy_tfrz=",D12.4," iconc=",D12.4," dltIce=",D12.4,&
!            " enth0=",D12.4," sst=",F6.2)') &
!           enthalpy_ocn*US%Q_to_J_kg, enthalpy_ocn_tfrz*US%Q_to_J_kg, &
!           iconc_tot, dlt_ice_tot*US%RZ_to_kg_m2,&
!           enthalpy_ocn0*US%Q_to_J_kg, OSS%SST_C(i,j)*US%C_to_degC
!      write(*,'(A)') trim(mesg)
!    endif
!
!    CS%Enth_out_ocn_old(c) = IOF%Enth_Mass_out_ocn(i,j)
!    CS%flux_salt_old(c) = IOF%flux_salt(i,j)
  
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
! Ice heat fluxes assigned to Ice in set_ocean_top_fluxes subroutine:
! Ice%flux_t(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_sh_ocn_top(i,j)
subroutine adjust_IOfluxes_isponge(dt_slow, CS, IG, IST, IOF, OSS, FIA, US)
  real,                       intent(in)  :: dt_slow !< The thermodynamic step [T ~> s].
  type(isponge_CS),           pointer     :: CS    !< A pointer that is set to point to the ice sponge control
                                                !! structure for this module
  type(ice_grid_type),        intent(in)  :: IG    !< The sea-ice specific grid type
  type(ice_state_type),       intent(in)  :: IST   !< A type describing the state of the sea ice
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model  
  type(ocean_sfc_state_type), intent(inout) :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  integer :: c, i, j, k, m, l, NkIce 
  integer :: ncat
  integer :: current_pe                                  ! debugging
  real    :: Idt_slow    ! The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new
  real    :: dlt_iconc, dlt_ithk, I_Nk
  real    :: dlt_salt, &    ! g/kg * kg/m2 salt change due to relax over a cat
             dlt_heat, &    ! J/m2=enthalpy[J/kg] * kg/m2 change due to relaxn over a cat
             dlt_water, dlt_snow
  real    :: dlt_ice           ! change of ice mass due to conc and thickness relax changes, over a cat
  real    :: ice_salin         ! average ice column S gSalt kg-1 over a cat
  real    :: water_ice_ocn, heat_ice_ocn, salt_ice_ocn
  real    :: dlt_heat_tot, dlt_salt_tot, dlt_ice_tot !> heat, salt, iceM changes over all cat
  real    :: heat_water    ! heat of water added to / taken from the ocean during relaxation [J/m2]
  real    :: dmm ! debug

!  dt_slow  = dt_slow*US%T_to_s
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0 / dt_slow
  NkIce = IG%NkIce
  I_Nk  = 1. / NkIce
  ncat  = IG%CatIce 

  ! Estimate ice volume change by categories:
  if (CS%num_col == 0) return
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)
    dlt_heat_tot = 0.0 ; dlt_ice_tot = 0.0 ; dlt_salt_tot = 0.0
    do k=1,ncat
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
      ! To cancel these surplus heat fluxes cause by relaxation need to add this amount with opposite sign
      ! enthalpy of water - E required to warm water from Tfrz to SST or vice versa for freezing
      dlt_iconc = iconc_new - iconc_old
      dlt_ithk  = ithk_new - ithk_old           ! ice mass change, kg m-2
      dlt_ice = ithk_new*iconc_new - ithk_old*iconc_old

      dlt_salt = 0.0 ; dlt_heat = 0.0  
      ice_salin = 0.0  ! for debugging
      ! enthalpy of ice/snow
      ! enthalpy to ocean = Negative of the enthalpy extracted from the ice 
      do l=1,NkIce 
        dlt_salt = dlt_salt + (dlt_ice * I_Nk) * IST%sal_ice(i,j,k,l) !> S change due to relax kg/m2*g/kg = g/m2
        dlt_heat = dlt_heat + (dlt_ice * I_Nk) * IST%enth_ice(i,j,k,l) ! kg/m2*J/kg = J/m2 
        ice_salin = ice_salin + IST%sal_ice(i,j,k,l) * I_Nk
      enddo 
      dlt_snow = dlt_iconc * IST%mH_snow(i,j,k)  ! kg m-2
      dlt_heat = dlt_heat + dlt_snow * IST%enth_snow(i,j,k,1)  ! kg m-2 * J kg-1 => J m-2

      dlt_ice_tot  = dlt_ice_tot + dlt_ice      ! kg/m2 
      dlt_heat_tot = dlt_heat_tot + dlt_heat    ! J/m2
      dlt_salt_tot = dlt_salt_tot + dlt_salt

!> Debug
      if (i==CS%itest .and. j==CS%jtest) then
        write(mesg, '("adjust_IO: cat=",I," iconc_new=",D12.3," iconc_old=",D12.3," ithk_new=",D12.3," ithk_old=",D12.3)') &
             k,iconc_new, iconc_old, ithk_new*US%RZ_to_kg_m2, ithk_old*US%RZ_to_kg_m2
        write(*,'(A)') trim(mesg)
        write(mesg, '("adjust_IO: dlt_iconc=",D12.4," dlt_ithk=",D12.4," dlt_ice=",D12.4, &
              " dlt_salt=",D12.4," dlt_heat=",D12.4," ice_salin=",D12.4," dlt_snow=",D12.4)') &
              dlt_iconc, dlt_ithk*US%RZ_to_kg_m2, dlt_ice*US%RZ_to_kg_m2, dlt_salt*US%S_to_ppt,&
              dlt_heat*US%QRZ_T_to_W_m2*US%T_to_s, ice_salin*US%S_to_ppt, dlt_snow*US%RZ_to_kg_m2
        write(*,'(A)') trim(mesg)
      endif
!!  Debug ===
    enddo   ! ice cat

    if (i==CS%itest .and. j==CS%jtest) then
! 0.001*Idt_slow - conversions kg / kg --> g/(kg*sec)
      write(mesg, '("dlt_ice_tot=",D12.4," dlt_heat_tot[J/m2]=",D12.4," dlt_salt_tot [kg/m2]=",D12.4)') &
           dlt_ice_tot*US%RZ_to_kg_m2, dlt_heat_tot*US%QRZ_T_to_W_m2*US%T_to_s, &
           dlt_salt_tot*0.001*US%RZ_to_kg_m2
      write(*,'(A)') trim(mesg)
    endif
! Subtract water enthalpy that is added to / taken from the ocean
    heat_water = (dlt_ice_tot + dlt_snow) * enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
    dlt_heat_tot = dlt_heat_tot - heat_water
    if (i==CS%itest .and. j==CS%jtest) then
      write(mesg, '("water heat=",D12.4," added water heat: dlt_heat_tot[J/m2]=",D12.4)') &
           heat_water*US%QRZ_T_to_W_m2*US%T_to_s, dlt_heat_tot*US%QRZ_T_to_W_m2*US%T_to_s
      write(*,'(A)') trim(mesg)
    endif


!> Save fluxes with opposite sign to cancel out added fluxes later in the code
! No need to adjust salt fluxes, only heat flux - added to sens heat sent to the ocean W/m2
!      IOF%flux_salt(i,j) = IOF%flux_salt(i,j) - salt_ice_ocn * (0.001*Idt_slow) 
!    IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) + dlt_heat_tot
    dmm = IOF%flux_sh_ocn_top(i,j)
    IOF%flux_sh_ocn_top(i,j) = IOF%flux_sh_ocn_top(i,j) - dlt_heat_tot*Idt_slow  ! J/m2*1/sec = W/m2

    if (i==CS%itest .and. j==CS%jtest) then     
      write(mesg, '("==> Final step adj: old flux_sh=",D12.3," new flux_sh=",D12.3)') &
            dmm*US%QRZ_T_to_W_m2*US%T_to_s, IOF%flux_sh_ocn_top(i,j)*US%QRZ_T_to_W_m2*US%T_to_s
      write(*,'(A)') trim(mesg)
    endif  

  enddo

end subroutine adjust_IOfluxes_isponge

subroutine check_FIA(dt_slow, CS, FIA, US, IG, txtinfo)
! Check Fast Ice Averages of fields (mostly fluxes) updates
! frazil_left  - The frazil heat flux that has not yet been
!! consumed in making ice [Q R Z ~> J m-2]. This array is decremented by the ice
!! model as the heat flux is used up.
  real,                     intent(in)  :: dt_slow !< The amount of time covered by this call [T ~> s].
  type(isponge_CS),         pointer     :: CS      !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(fast_ice_avg_type),  intent(in)  :: FIA     !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(unit_scale_type),    intent(in)  :: US      !< A structure with unit conversion factors
  type(ice_grid_type),      intent(in)  :: IG      !< The sea-ice specific grid type
  character(len=*),  intent(in) :: txtinfo

  integer :: c, i, j, k
  real :: flux_fraz_left, flux_bmelt, Idt_slow
  character(len=256) :: mesg

  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0 / dt_slow

  ! Estimate ice volume change by categories:
  if (CS%num_col == 0) return
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)
    if (i==CS%itest .and. j==CS%jtest) then
      flux_fraz_left = US%QRZ_T_to_W_m2*FIA%frazil_left(i,j)*Idt_slow
      write(*,'(A)') trim(txtinfo)
      write(mesg, '("Flux Frazil heat left W/m2 =",D12.4)') flux_fraz_left
      write(*, '(A)') trim(mesg)
      do k=1,IG%CatIce
        flux_bmelt = US%QRZ_T_to_W_m2*FIA%bmelt(i,j,k)*Idt_slow
        write(mesg, '("cat ",I," Flux bottom melt W/m2 =",D12.4)') k, flux_bmelt
        write(*, '(A)') trim(mesg)
      enddo
    endif
  enddo

end subroutine check_FIA

subroutine check_IOF(CS, IOF, OSS, US, txtinfo)
! Check ice-ocen fluxes in IOF for test location
! debugging
! Note that IOF fluxes are weighted by ocean partial area in SIS_thermodynamics routine:
! IOF%flux_sh_ocn_top(i,j) = part_ocn * FIA%flux_sh_top(i,j,0)
! 
  type(isponge_CS),           pointer     :: CS    !< A pointer that is set to point to the ice sponge control
                                                !! structure for this module
  type(ice_ocean_flux_type),  intent(in)  :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model  
  type(ocean_sfc_state_type), intent(in)  :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(unit_scale_type),      intent(in)  :: US  !< A structure with unit conversion factors
  character(len=*),  intent(in) :: txtinfo
 
  integer :: c, i, j
  real :: flux_t, flux_q, flux_sw_visdir, flux_sw_visdif, flux_sw_nirdir, flux_sw_nirdif
  real :: flux_lw, flux_lh, flux_salt, sst_c
  character(len=256) :: mesg
 
  ! Estimate ice volume change by categories:
  if (CS%num_col == 0) return
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)  
    if (i==CS%itest .and. j==CS%jtest) then
      flux_t = US%QRZ_T_to_W_m2*IOF%flux_sh_ocn_top(i,j)
      flux_q = US%RZ_T_to_kg_m2s*IOF%evap_ocn_top(i,j)
      flux_sw_visdir = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,VIS_DIR)
      flux_sw_visdif = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,VIS_DIF)
      flux_sw_nirdir = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,NIR_DIR)
      flux_sw_nirdif = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,NIR_DIF)
      flux_lw = US%QRZ_T_to_W_m2*IOF%flux_lw_ocn_top(i,j)
      flux_lh = US%QRZ_T_to_W_m2*IOF%flux_lh_ocn_top(i,j)
      flux_salt = US%S_to_ppt*US%RZ_T_to_kg_m2s*IOF%flux_salt(i,j)
      sst_c = US%C_to_degC*OSS%SST_C(i,j)

      write(*,'(A)') trim(txtinfo)
      write(mesg, '("Fluxes W/m2: sens t=",D12.4," sw_visdir=",D12.4," sw_visdif=",D12.4," sw_nirdir=",D12.4," sw_nirdif=",D12.4)') &
            flux_t, flux_sw_visdir, flux_sw_visdif, flux_sw_nirdir, flux_sw_nirdif
      write(*, '(A)') trim(mesg)
      write(mesg, '("Fluxes W/m2: latent lh=",D12.4," longwv lw=",D12.4," salt kg/m2*s=",D12.4," sst_c=",F8.3)') &
           flux_lh, flux_lw, flux_salt, sst_c
      write(*, '(A)') trim(mesg)
    endif
  enddo

end subroutine check_IOF

!> Deallocate memory associated with the SIS_optics module
subroutine SIS_sponge_end(CS)
  type(isponge_CS), pointer :: CS !< The ice sponge control structure that is deallocated here
  
  deallocate(CS) 
                 
end subroutine SIS_sponge_end

end module SIS_sponge 

