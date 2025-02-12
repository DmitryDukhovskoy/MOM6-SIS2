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
use MOM_coms,          only : sum_across_PEs, max_across_PEs
use MOM_coms,          only : PE_here   !! debugging
use MOM_time_manager,  only : time_type, set_date, get_time, get_date
use MOM_unit_scaling,  only : unit_scale_type
use ice_grid,          only : ice_grid_type

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : file_exists, MOM_read_data, slasher
use MOM_io,            only : axis_info
use MOM_interpolate,   only : init_external_field, get_external_field_info, time_interp_external_init
use MOM_interpolate,   only : time_interp_external
use MOM_interpolate,   only : external_field     
                            
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

public initialize_icerelax_file, apply_isponge, set_up_isponge_field, SIS_sponge_end
public global_to_local_ij
!public adjust_IOfluxes_isponge
!public check_IOF, check_FIA

!> A structure for creating arrays of pointers to 3D arrays
type, public :: p3d; private
  !integer :: id !< id for FMS external time interpolator
  integer :: nz_data !< The number of vertical levels in the input field.
  integer :: num_tlevs !< The number of time records contained in the file
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3D array [various]
  character(len=15)               :: fld_name    !< Name of the ice field being relaxed
end type p3d
!> A structure for creating arrays of pointers to 2D arrays
type, public :: p2d; private
  type(external_field) :: field !< Time interpolator field handle
  !integer :: nz_data !< The number of vertical levels in the input field
  integer :: ncat_data !< The number of sea ice categories
  integer :: num_tlevs !< The number of time records contained in the file
  real :: scale = 1.0  !< A multiplicative factor by which to rescale input data [various]
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array [various]
  character(len=:), allocatable  :: name  !< The name of the input field
  character(len=:), allocatable  :: long_name !< The long name of the input field
  character(len=:), allocatable  :: unit !< The unit of the input field
  type(axis_info),  allocatable  :: axes_data(:) !< Axis types for the input field
                                                 !! name, longname, cartesian("X","Y",...) ax_size,...
end type p2d
!
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
  integer         :: num_col         !< The number of sponge points within the computational domain.
  integer, public :: fldno = 0       !< The number of fields which have already been
                                     !! registered by calls to set_up_sponge_field

  integer, pointer :: col_i(:) => NULL()         !< Array of the i-indicies of each of the columns being damped.
  integer, pointer :: col_j(:) => NULL()         !< Array of the j-indicies of each of the columns being damped.
  real, pointer    :: Iresttime_col(:) => NULL() !< The inverse restoring time of each column [T-1 ~> s-1].

  type(p3d) :: var(MAX_FIELDS_RLX_)     !< Pointers to the fields that will be relaxed
  type(p2d) :: Ref_val(MAX_FIELDS_RLX_) !< Relaxation values - The values to which the fields are 
                                        ! relaxed distributed by ice cats (linear_index, ice_cat)
  type(f2d) :: Old_val(MAX_FIELDS_RLX_) !< Keep old values of relaxed fields prior to relaxation
                                        !! debug only, will need to get rid off later
  type(f2d) :: Ref_orig(MAX_FIELDS_RLX_) !< Relaxation values original input fields, i.e.
                                         !! on 2d grid not distributed by ice cats.
  logical :: time_varying_sponges       !< True if using newer sponge code
  logical :: spongeDataOngrid           !< True if the sponge data are on the model horizontal grid

!  real, allocatable, dimension(:,:,:) :: & 
!        dlt_enth_ice             ! Change of ice enthalpy during relaxation by ice thickn, cat
!  real, allocatable, dimension(:) :: &
!      Enth_out_ocn_old, & ! Negative of the enthalpy extracted from ice by water fluxes to the ocean [Q R Z ~> J m-2]
!      flux_salt_old  ! The flux of salt out of the ocean [1e3 S R Z T-1 ~> kgSalt m-2 s-1]
end type isponge_CS

contains

!> This subroutine sets the inverse restoration time (Idamp) for sea ice fields and
!! the values towards which the interface heights and an arbitrary
!! number of tracers should be restored within the relaxation zone. 
subroutine initialize_icerelax_file(param_file, G, IG, CS, US, IST, Time, itest, jtest)
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  type(unit_scale_type),   intent(in) :: US         !< A structure with unit conversion factors
  type(ice_state_type),    intent(in) :: IST        !< A type describing the state of the sea ice
  type(time_type),         intent(in) :: Time       !< The sea-ice model's clock,
  integer, optional, intent(in) :: itest, jtest     !< Test grid indices for debugging, on data domain

  real, dimension(SZI_(G),SZJ_(G))  :: Irelax  !< The inverse of the restoring time [T-1 ~> s-1].
  real, allocatable, dimension(:,:) :: rlx_H ! A temporary array for reading relax target ice thickness
                                             ! mean grid cell value  kg/m [R Z L ~> kg m-1]
  real, allocatable, dimension(:,:) :: rlx_C ! A temporary array for reading relax target ice partial area

  integer :: i, j, k, is, ie, js, je, ncat
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: year !< The current model year
  integer :: day  !< The current model year-day
  integer :: second !< The second of the day
  integer :: mon, hr, minute, itick
  integer :: start_of_day, num_days
  real :: max_rlxrate, rho_ice
  integer, dimension(4) :: siz
  character(len=40) :: ithck_var, iarea_var, rlxrate_var, rlx_unit
  character(len=40) :: mdl = "initialize_icerelax_file"
  character(len=50) :: rlx_long_name
  character(len=200) :: relaxrate_file, state_file  ! relax filenames: inverse time, target fields
  character(len=200) :: filename, inputdir ! Strings for file/path and path.
  character(len=256) :: mesg

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  Irelax = 0.0

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "ISPONGE_RELAX_FILE", relaxrate_file, &
                 "The name of the file with the sponge relaxation rates.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "ISPONGE_STATE_FILE", state_file, &
                 "The name of the file with the state to relax toward.", &
                 fail_if_missing=.true.)                
  call get_param(param_file, mdl, "ISPONGE_ITHCK_VAR", ithck_var, & 
                 "The name of the ice thickness variable in "//&
                 "ISPONGE_STATE_FILE.", default="ithkn")
  call get_param(param_file, mdl, "ISPONGE_IAREA_VAR", iarea_var, & 
                 "The name of the ice partial area/concentration variable in "//&
                 "ISPONGE_STATE_FILE.", default="ithkn")
  call get_param(param_file, mdl, "ISPONGE_RLXRATE_VAR", rlxrate_var, &
                 "The name of the relaxation rate variable in "//&
                 "ISPONGE_RELAX_FILE.", default="relax_rate")

  ! Read in relaxation rate, s-1, for ice thickness and partial area
  filename = trim(inputdir)//trim(relaxrate_file)
  call log_param(param_file, mdl, "INPUTDIR/ISPONGE_RELAX_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call SIS_error(FATAL, " initialize_icerelax_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, rlxrate_var, Irelax(:,:), G%Domain, scale=US%s_to_T) 
  max_rlxrate =  maxval(Irelax*US%T_to_s)
  call max_across_PEs(max_rlxrate)
  call get_date(Time, year, mon, day, hr, minute, second, itick)
  write(mesg,'("SIS Time:",i6,2("/",i2.2),1x,3(":",i2.2),"; max(Irelax)=",D13.4," s-1")') &
        year, mon, day, hr, minute, second, max_rlxrate
  call SIS_mesg(mesg) 
  write(mesg,'("SIS tick=",I)') itick
  call SIS_mesg(mesg)
  call get_time(Time, start_of_day, num_days)
  write(mesg,'("SIS start day=",I," num_days=",I)') start_of_day, num_days
  call SIS_mesg(mesg)

  call SIS_mesg('initialize_icerelax_file: Calling initialize_isponge') 

  if (present(itest) .and. present(jtest)) then
    write(mesg,'(A," itest/jtest=",2(i4,1x)," calling initialize_isponge")') &
         trim(mdl), itest, jtest
    write(*,'(A)')
    call initialize_isponge(param_file, Irelax, G, IG, CS, itest=itest, jtest=jtest)

  else
    call initialize_isponge(param_file, Irelax, G, IG, CS)
  endif

  ! Now register all of the fields which are nudged in the relaxation region.
  filename = trim(inputdir)//trim(state_file)
  call log_param(param_file, mdl, "INPUTDIR/ISPONGE_STATE_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call SIS_error(FATAL, " initialize_icerelax_files: Unable to open "//trim(filename))
!
  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
! Read target value: ice thickness, ice concentration - mean grid cell values
! need to distribute by categories
! First read in a 2D temporary array
  call SIS_mesg('initialize_icerelax_file: Calling set_up_isponge_field: mH_ice') 
  call set_up_isponge_field(filename, ithck_var, Time, G, IG, US, IST%mH_ice, CS, &
       'mH_ice', rlx_long_name='ice_thickness', rlx_unit='kg m-2', scale=US%m_to_Z * rho_ice)
  call SIS_mesg('initialize_icerelax_file: Calling set_up_isponge_field: part_size') 
  call set_up_isponge_field(filename, iarea_var, Time, G, IG, US, IST%part_size, CS, &
         'part_size', rlx_long_name='partial_area', rlx_unit='none')

end subroutine initialize_icerelax_file
!
!> This subroutine determines the number of points which are within ice sponges in
!! this computational domain.  Only points that have positive values of
!! Iresttime and which mask2dT indicates are ocean points are included in the
!! sponges.  
!subroutine initialize_isponge(Iresttime, IST, G, IG, param_file, CS)
subroutine initialize_isponge(param_file, Iresttime, G, IG, CS, itest, jtest, time_var_rlx, sponge_ongrid)
!  type(ice_state_type),    intent(in) :: IST        !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: Iresttime  !< The inverse of the restoring time [T-1 ~> s-1].
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  integer, optional, intent(in) :: itest, jtest     !< test grid indices for debugging
  logical, optional, intent(in) :: time_var_rlx, sponge_ongrid !< place-holders, currently both true

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "initialize_isponge"  ! This module's name.
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
    call SIS_mesg("initialize_isponge: use_isponge=True")
  else
    call SIS_mesg("initialize_isponge: use_isponge=False")
  endif
  if (.not.use_isponge) return
  allocate(CS)

  CS%time_varying_sponges = .true.
  CS%spongeDataOngrid = .true.
  if (present(time_var_rlx)) CS%time_varying_sponges = time_var_rlx
  if (present(sponge_ongrid)) CS%spongeDataOngrid = sponge_ongrid

  CS%use_isponge = use_isponge
  if (present(itest) .and. present(jtest)) then
    write(mesg,'(A," itest/jtest =",2(i5,1x))') trim(mdl), itest, jtest
    write(*,'(A)') trim(mesg)
    CS%itest = itest
    CS%jtest = jtest
  else
    CS%itest = -999
    CS%jtest = -999
  endif

  CS%num_col = 0 ; CS%fldno = 0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo

!  write(mesg,'("SIS_sponge: num_col=",I8)') CS%num_col
!  call SIS_mesg(mesg)

  if (CS%num_col > 0) then
    allocate(CS%Iresttime_col(CS%num_col), source=0.0)
    allocate(CS%col_i(CS%num_col), source=0)
    allocate(CS%col_j(CS%num_col), source=0)
!    allocate(CS%dlt_enth_ice(CS%num_col,IG%CatIce,IG%NkIce), source=0.0) ! <-- this is probably not needed
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

  write(mesg,'(A,": total isponge cols=",i8)') trim(mdl), total_isponge_cols
  call SIS_mesg(mesg)
  
  call log_param(param_file, mdl, "!Total isponge columns", total_isponge_cols, &
                 "The total number of ice columns where relaxation is applied.")

end subroutine initialize_isponge

!> This subroutine stores the reference profile for the SIS variable whose
!! address is given by f_ptr. Reference profile = values towards which the
!! SIS field is being relaxed to. 
!! Current version assumes 2D fields only (+ 1D for categories) such as Hice
subroutine set_up_isponge_field(filename, fieldname, Time, G, IG, US, f_ptr, CS, &
                                rlxfld_name, rlx_long_name, rlx_unit, scale)
  character(len=*),        intent(in) :: filename   !< The name of the file with the
                                                    !! time varying field data
  character(len=*),        intent(in) :: fieldname  !< The name of the field in the file
                                                    !! with the time varying field data
  type(time_type),         intent(in) :: Time       !< The current model time
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in) :: US         !< A structure with unit conversion factors
!  integer,                 intent(in) :: ncat
!  real, dimension(SZI_(G), SZJ_(G), ncat), &
!                           intent(in) :: sp_val     !< The reference profiles of the quantity 
!                                                    !! being registered [various]
  real, dimension(SZI_(G), SZJ_(G), IG%CatIce), &
                   target, intent(in) :: f_ptr      !< a pointer to the field which will be relaxed [various]
  type(isponge_CS),     pointer       :: CS         !< A pointer to the control structure for this module that
                                                    !! is set by a previous call to initialize_sponge.
  character(*),            intent(in) :: rlxfld_name !< Name of the relaxed field
  character(len=*),        optional,  &
                           intent(in) :: rlx_long_name !< The long name of the tracer field
                                                      !! if not given, use the sp_name
  character(len=*),        optional,  &        
                           intent(in) :: rlx_unit !< The unit of the tracer field
                                                 !! if not given, use 'none'
  real,          optional, intent(in) :: scale !< A factor by which to rescale the input data, including any
                                               !! contributions due to dimensional rescaling [various ~> 1].


  ! Local variables
  integer :: isd, ied, jsd, jed
  integer, dimension(4) :: fld_sz
  integer :: j, k, col
  character(len=256) :: mesg ! String for error messages
  character(len=256) :: long_name ! The long name of the tracer field
  character(len=256) :: unit ! The unit of the tracer field

  long_name = rlxfld_name; if (present(rlx_long_name)) long_name = rlx_long_name
  unit = 'none'; if (present(rlx_unit)) unit = rlx_unit

  if (.not.associated(CS)) return
  ! initialize time interpolator module
  call time_interp_external_init()
  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
  CS%fldno = CS%fldno + 1
  write(mesg,'("set_up_isponge: fldno=",I)') CS%fldno
  call SIS_mesg(mesg)
  if (CS%fldno > MAX_FIELDS_RLX_) then
    write(mesg,'("Increase MAX_FIELDS_RLX_ to at least ",I3," in SIS_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call SIS_error(FATAL,"set_up_isponge_field: "//mesg)
  endif
  ! get a unique time interp id for this field. Ice relax target fields are on-grid
  if (CS%spongeDataOngrid) then
    call SIS_mesg("set_up_isponge_field: calling init_external_field")
    CS%Ref_val(CS%fldno)%field = init_external_field(filename, fieldname, MOM_domain=G%Domain, &
               verbose=.true.)
  else
    call SIS_error(FATAL,"set_up_isponge_field: SIS2 relaxation fields on a not-native grid not implemented")
  endif
  CS%Ref_val(CS%fldno)%name = rlxfld_name
  CS%Ref_val(CS%fldno)%long_name = long_name
  CS%Ref_val(CS%fldno)%unit = unit
  fld_sz(1:4) = -1
  call get_external_field_info(CS%Ref_val(CS%fldno)%field, size=fld_sz, axes=CS%Ref_val(CS%fldno)%axes_data)
  !nz_data = fld_sz(3)
  !ncat_data = IG%CatIce
  CS%Ref_val(CS%fldno)%ncat_data = IG%CatIce !< individual relax fields should have same # of categories
  CS%Ref_val(CS%fldno)%num_tlevs = fld_sz(4)
  CS%Ref_val(CS%fldno)%scale = 1.0 ; if (present(scale)) CS%Ref_val(CS%fldno)%scale = scale
  ! initializes the target profile array for this field
  ! for all columns which will be masked

  allocate(CS%Ref_val(CS%fldno)%p(CS%num_col,IG%CatIce), source=0.0)
  allocate(CS%Old_val(CS%fldno)%fld(CS%num_col,IG%CatIce), source=0.0)
  allocate(CS%Ref_orig(CS%fldno)%fld(isd:ied,jsd:jed), source=0.0)

  CS%var(CS%fldno)%p => f_ptr    ! points to the actual ice fields that will be relaxed
  CS%var(CS%fldno)%fld_name = rlxfld_name

  write(mesg,'("set_up_isponge_field: ",A," fld_sz(1:4)=",4(I5,1x)," scale=",f14.6)') &
               rlxfld_name, fld_sz(1:4), CS%Ref_val(CS%fldno)%scale
  call SIS_mesg(mesg)

end subroutine set_up_isponge_field

!> This subroutine applies relaxation ("damping") to ice thickness (by categories) and ice concentration
!! tracers for every column where the relaxation time scale > 0.
subroutine apply_isponge(dt_slow, CS, G, IG, IST, US, OSS, Time)
  real,                      intent(in)  :: dt_slow   !< The amount of time covered by this call [T ~> s].
  type(isponge_CS),          pointer     :: CS     !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(ice_grid_type),       intent(in)  :: IG     !< The sea-ice specific grid type
  type(SIS_hor_grid_type),   intent(in)  :: G          !< The horizontal grid type
!  type(ice_ocean_flux_type), intent(inout)  :: IOF  !< A structure containing fluxes from the ice to
!                                                   !! the ocean that are calculated by the ice model.
  type(ice_state_type),   intent(inout)  :: IST    !< A type describing the state of the sea ice
  type(unit_scale_type),     intent(in)  :: US     !< A structure with unit conversion factors
  type(ocean_sfc_state_type), intent(in) :: OSS    !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(time_type),           intent(in)    :: Time !< The current model date

  ! Local variables
!  real :: H_rescale_ice, H_rescale_snow
  real :: damp         ! The timestep times the local damping coefficient [nondim].
  real :: I1pdamp      ! I1pdamp is 1/(1 + damp). [nondim]
  real :: p_old        ! debugging  
  real :: dt           ! time step in s
  real :: s_ice_bulk   ! ice bulk S for filling S values in the newly ceated ice 
  real, allocatable :: sice(:), tfi(:)
  real :: enth_ice, Tfrz, coeff, enth_Tfrz
  character(len=40)  :: mdl = "apply_isponge"  ! This subroutine's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  real    :: Idt_slow    ! The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new
  real    :: iconc_tot, dlt_ice_tot, iconc_tot_old
  real    :: dlt_iconc, dlt_ithk, enthalpy_ocn, enthalpy_ocn_tfrz
  real    :: dlt_salt, dlt_heat, dlt_water, dlt_snow
  real    :: dlt_ice           ! total change of ice due to conc and thickness relaxation
  real    :: ice_salin         ! average ice column S gSalt kg-1 
  real    :: water_ice_ocn, heat_ice_ocn, salt_ice_ocn
  real    :: enthalpy_ocn0
  real    :: dlt_enth
  logical :: f_debug
!
  real, dimension(:,:), allocatable  :: data_in  !< A buffer for storing the full 2-d time-interpolated array
  real, dimension(:,:),  allocatable   :: mask_in    !< A 2-d mask for extended input grid [nondim]

  !real :: dgr2rad  ! A conversion factor from degrees to radians [radians degree-1]
  real :: I_Nk
  integer :: id, jd, kd, jdp ! Input dataset data sizes
  type(axis_info), dimension(4) :: axes_data
  integer :: i, j, k, l, m, col
  integer :: ii, jj, iiG, jjG
  integer :: current_pe
  integer :: CatIce             !< The number of sea ice categories.
  integer :: NkIce      !< The number of vertical partitions within the sea ice.
  integer :: is, ie, js, je     !< compute domain indices
  integer :: isg, ieg, jsg, jeg !< global extent
  integer :: isd, ied, jsd, jed !< data domain indices
  integer :: nid, njd, isdG, iedG, jsdG, jedG
  !integer :: id_clock_read
  !integer :: turns
  integer, dimension(4) :: fld_sz

  !turns = G%HI%turns
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  ! Include halo points, i.e. data domain:
  nid  = G%ied - G%isd + 1
  njd  = G%jed - G%jsd + 1
  isdG = G%isd_global; iedG = isdG + nid
  jsdG = G%jsd_global; jedG = jsdG + njd
!
  f_debug = .true.
  CatIce = IG%CatIce
  NkIce  = IG%NkIce
  I_Nk  = 1. / NkIce
  dt = dt_slow*US%T_to_s
  s_ice_bulk = 3.0*US%ppt_to_S
  !dgr2rad = atan(1.0)/45.  

  if (CS%num_col == 0) return
  
  current_pe = PE_here()
  
! First get relax fields and interp. in time:
  allocate(data_in(isd:ied,jsd:jed))  
  allocate(sice(NkIce), tfi(NkIce), source=-999.)
  ! Debug
  if (CS%itest .gt. 0 .or. CS%jtest .gt. 0) then
    write(mesg,'("apply_isponge: itest/jtest=",2(i5,1x)," isd/ied=",2(i4,1x),"jsd/jed=",2(i4,1x))') &
        CS%itest, CS%jtest, isd, ied, jsd, jed
    write(*,'(A)') trim(mesg)
  endif

  !call SIS_mesg(mesg)
  do m=1,CS%fldno
    call time_interp_external(CS%Ref_val(m)%field, Time, data_in, verbose=.true.)
    CS%Ref_orig(m)%fld(:,:) = data_in(:,:)
    ! Debug
    do col=1,CS%num_col
      i = CS%col_i(col) ; j = CS%col_j(col)
      if (i.eq.CS%itest .and. j.eq.CS%jtest) then
        iiG = isdG + (i-1)  ; jjG = jsdG + (j-1)
        write(mesg,'("apply_isponge: test i/j=",2(i5,1x)," time_iterp data_in=",f8.4)') &
        iiG, jjG, data_in(i,j)
        write(*,'(A)') trim(mesg)
!        do ii=isd,ied ; do jj=jsd,jed
!          iiG = isdG + (ii-1)
!          jjG = jsdG + (jj-1)
!          write(mesg,'(" == Check: Global i/j  data_in = ",2(i4,1x), f8.4)') &
!               iiG, jjG, data_in(ii,jj)
!          write(*,'(A)') trim(mesg)
!        enddo ; enddo 
      endif
    enddo

  enddo

  ! COnvert input 2D fields --> 3D ice thicknesses and concentration by categories
  ! scale ice thickness m --> kg m-2 and unscale US%m_to_Z
  ! Simplest approach - place all ice into 1 category based on hice input
  ! Register relaxation fields:
  call redistribute_ice2cats_simple(CS, IG, G)

  do col=1,CS%num_col
    i = CS%col_i(col) ; j = CS%col_j(col)             
    damp = dt * CS%Iresttime_col(col); I1pdamp = 1.0 / (1.0 + damp)

    do k=1,IG%CatIce
      do m=1,CS%fldno
!        p_old = CS%var(m)%p(i,j,k)  ! debugging
        CS%Old_val(m)%fld(col,k) = CS%var(m)%p(i,j,k)  

        CS%var(m)%p(i,j,k) = I1pdamp * &
           (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(col,k)*damp)

        if (i==CS%itest .and. j==CS%jtest) then
          select case (trim(CS%var(m)%fld_name))
            case('mH_ice')    ; coeff = US%RZ_to_kg_m2
            case('part_size') ; coeff = 1.0
            case default
              write(mesg,'("SIS_sponge: Unknown relaxation field: ",A)') trim(CS%var(m)%fld_name)
              call SIS_error(FATAL,"apply_isponge: "//mesg)          
          end select

          write(mesg,'(A8," k=",I2," old:=",D12.4," new=",D12.4,&
                      " 1/tau=",D12.4," refval=",D12.4)') &
            CS%var(m)%fld_name(1:8), k, CS%Old_val(m)%fld(col,k)*coeff, &
            CS%var(m)%p(i,j,k)*coeff, &
            CS%Iresttime_col(col), CS%Ref_val(m)%p(col,k)*coeff
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
!        CS%dlt_enth_ice(col,k,l) = dlt_enth
      enddo
    enddo  ! CatIce

    iconc_tot = 0.0
    iconc_tot_old = 0.0
!    dlt_ice_tot = 0.0
    do k=1,IG%CatIce
      do m=1,CS%fldno
        fld_name = CS%var(m)%fld_name
        select case (trim(fld_name))
!          case('mH_ice')
!            ithk_old = CS%Old_val(m)%fld(col,k)
!            ithk_new = CS%var(m)%p(i,j,k)
          case('part_size')
            iconc_old = CS%Old_val(m)%fld(col,k)
            iconc_new = CS%var(m)%p(i,j,k)
            iconc_tot_old = iconc_tot_old + iconc_old
            iconc_tot = iconc_tot + iconc_new
        end select
      enddo
!      dlt_iconc = iconc_new - iconc_old
!      dlt_ithk  = ithk_new - ithk_old           ! ice mass change, kg m-2
!      dlt_ice = ithk_new*iconc_new - ithk_old*iconc_old
!      dlt_ice_tot = dlt_ice_tot + dlt_ice
    enddo
!
!    enthalpy_ocn = enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
!    enthalpy_ocn_tfrz = enthalpy_liquid_freeze(OSS%s_surf(i,j), IST%ITV) 
!    enthalpy_ocn0 = enthalpy_liquid(0.0, OSS%s_surf(i,j), IST%ITV)
!
    if (i==CS%itest .and. j==CS%jtest .and. f_debug) then  
      write(mesg, '("old conc=",F6.4," new conc=",F6.4," ice enth J/kg=",D12.3)') &
            iconc_tot_old, iconc_tot, IST%enth_ice(i,j,k,l)*US%Q_to_J_kg
!      write(mesg, '("enthalpy_ocn=",D12.4," enthalpy_tfrz=",D12.4," iconc=",D12.4," dltIce=",D12.4,&
!            " enth0=",D12.4," sst=",F6.2)') &
!           enthalpy_ocn*US%Q_to_J_kg, enthalpy_ocn_tfrz*US%Q_to_J_kg, &
!           iconc_tot, dlt_ice_tot*US%RZ_to_kg_m2,&
!           enthalpy_ocn0*US%Q_to_J_kg, OSS%SST_C(i,j)*US%C_to_degC
      write(*,'(A)') trim(mesg)
    endif
!
!    CS%Enth_out_ocn_old(col) = IOF%Enth_Mass_out_ocn(i,j)
!    CS%flux_salt_old(col) = IOF%flux_salt(i,j)
  
  enddo

  if (allocated(sice)) deallocate(sice)
  if (allocated(tfi)) deallocate(tfi)

end subroutine apply_isponge

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
!
!> Convert global indices (itestG,jtestG) to indices on current tile 
subroutine global_to_local_ij(G, itestG, jtestG, itest, jtest)
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  integer, intent(in) :: itestG, jtestG
  integer, intent(out) :: itest, jtest

  integer :: current_pe, nihalo, njhalo, iscG, iecG, jscG, jecG
  integer :: isdG, jsdG, iedG, jedG
  integer :: nic, njc, nid, njd

  character(len=50) :: mdl  ! subroutine name
  character(len=256) :: mesg

  mdl = 'global_to_local_ij'
  nihalo = G%Domain%nihalo
  njhalo = G%Domain%njhalo

  current_pe = PE_here()

  ! Exclude halo points, computational domain:
  nic = G%iec - G%isc + 1
  njc = G%jec - G%jsc + 1
  iscG = G%isd_global + nihalo; iecG = iscG + nic
  jscG = G%jsd_global + njhalo; jecG = jscG + njc

  ! Include halo points, i.e. data domain:
  nid  = G%ied - G%isd + 1
  njd  = G%jed - G%jsd + 1
  isdG = G%isd_global; iedG = isdG + nid
  jsdG = G%jsd_global; jedG = jsdG + njd

  ! Find test point:
  itest = 0; jtest = 0
  if (iscG <= itestG .and. itestG <= iecG .and. jscG <= jtestG .and. jtestG <= jecG) then
    itest = itestG - iscG + 1; jtest = jtestG - jscG + 1
  endif
                 
  if (itest > 0 .and. jtest > 0) then
    write(mesg, '(A," current_pe=",i7," Global i, j=", 2(i5,1x)," local itest, jtest=", 2(i5,1x))') &
         trim(mdl), current_pe, itestG, jtestG, itest, jtest
    write(*, '(A)') trim(mesg)
  endif        

end subroutine global_to_local_ij
!
!> Simplest redistribution of 2D hice and iconc into ice thickness categories
!! place all ice into 1 category based on original ice thickness (hice)
subroutine redistribute_ice2cats_simple(CS, IG, G, rescaled)
  type(isponge_CS),        pointer     :: CS       !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(ice_grid_type),     intent(in)  :: IG       !< The sea-ice specific grid type
  type(SIS_hor_grid_type), intent(in)  :: G        !< The horizontal grid type
  logical, optional,       intent(in)  :: rescaled !< true if input hice (m) converted to kg/m2 and scaled
                                                   !! default = .false.

  ! local variables
  integer :: isd, ied, jsd, jed             ! data domain indices
  integer :: m, i, j, k, col
  integer :: CatIce     !< The number of sea ice categories.
!  integer :: NkIce     ! the number of ice thickn. categories
  integer :: icat0

  real, allocatable, dimension(:,:) :: cice2d, hice2d
  real, allocatable, dimension(:) :: hLim_vals
  real :: Iscale !< inverse scale to "unscale" the data
  real :: scale_cf 
  real :: hice, cice  !< relax ice thikn (cell mean) and conc at a grid pnt
  character(len=40)  :: mdl = "redistribute_ice2cats_simple"  ! This module's name.
  character(len=256) :: mesg
  logical :: rescale_hice

  rescale_hice = .false.
  if (present(rescaled)) rescale_hice = rescaled

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  CatIce = IG%CatIce

  allocate(cice2d(isd:ied,jsd:jed), hice2d(isd:ied,jsd:jed))
  allocate(hLim_vals(CatIce+1))
  
  hLim_vals(:) = IG%cat_thick_lim(:)

  do m=1,CS%fldno
    scale_cf = CS%Ref_val(m)%scale
    Iscale = 1.0
    if (abs(1.-scale_cf).gt.1.e-10 .and. scale_cf.gt.0.) &
      Iscale = 1.0/scale_cf
    select case (trim(CS%var(m)%fld_name))
      case('mH_ice')
        hice2d = CS%Ref_orig(m)%fld 
        if (rescale_hice .and. abs(1.-scale_cf).gt.1.e-10) &
            hice2d = CS%Ref_orig(m)%fld*Iscale
      case('part_size') 
        cice2d = CS%Ref_orig(m)%fld ! conc is not scaled
    end select 
  enddo 

  do col=1,CS%num_col
    i = CS%col_i(col) ; j = CS%col_j(col)
    hice = hice2d(i,j)
    cice = cice2d(i,j)
    ! ice categories: SIS_state_initialization.F90
! Debug:
    if (i.eq.CS%itest .and. j.eq.CS%jtest) then
      do k=1,CatIce
        write(mesg,  '(A,"test: cat=",i1," hlim=",f7.3)') trim(mdl), k, hLim_vals(k)
        write(*,'(A)') trim(mesg)
      enddo
      write(mesg,'(A,"test: hice=",f7.3," cice=",f6.3," hLim min/max=",2(f7.2,1x))') &
            trim(mdl), hice, cice, hLim_vals(1), hLim_vals(CatIce)
      write(*,'(A)') trim(mesg) 
    endif

    if (hice .ge. hLim_vals(CatIce)) then
      icat0 = CatIce
    elseif (hice .lt. hLim_vals(1)) then
      icat0=1
    else
      do k=1,CatIce
        if (hice .ge. hLim_vals(k) .and. hice .lt. hLim_vals(k+1)) then
          icat0 = k
          exit
        endif
      enddo
    endif

    do m=1,CS%fldno ; do k=1,CatIce
      if (k.eq.icat0) then
        select case (trim(CS%var(m)%fld_name))
          case('mH_ice')
            CS%Ref_val(m)%p(col,k) = hice*CS%Ref_val(m)%scale
! Debug:
          if (i==CS%itest .and. j==CS%jtest) then
            write(mesg,'(A," mH_ice scale=",f14.6)') trim(mdl), CS%Ref_val(m)%scale
            write(*,'(A)') trim(mesg)
          endif

          case('part_size')
            CS%Ref_val(m)%p(col,k) = cice
        end select
      endif
    enddo ; enddo

  enddo

end subroutine redistribute_ice2cats_simple


!> Deallocate memory associated with the SIS_optics module
subroutine SIS_sponge_end(CS)
  type(isponge_CS), pointer :: CS !< The ice sponge control structure that is deallocated here
  
  deallocate(CS) 
                 
end subroutine SIS_sponge_end

end module SIS_sponge 

