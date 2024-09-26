!> Implements sponge regions in SIS2 model
module SIS_sponge


!! Module to read in Time of the ice fields, ice concentration, ice thickness, relaxation time scale

! TODO: Figure out what modules to use
! This is the best guess:
use MOM_coms,          only : sum_across_PEs
use ice_grid,          only : ice_grid_type
use ice_spec_mod,      only : get_sea_surface

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE
use MOM_data_override, only : data_override
use MOM_EOS,           only : EOS_type, calculate_density_derivs
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : file_exists, MOM_read_data, slasher
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type
                            
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_framework,     only : coupler_type_spawn, coupler_type_initialized
use SIS_framework,     only : coupler_type_increment_data, coupler_type_rescale_data
use SIS_framework,     only : coupler_type_send_data
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_optics,        only : VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF
use SIS_sum_output,    only : SIS_sum_out_CS, write_ice_statistics! , SIS_sum_output_init
use SIS_sum_output,    only : accumulate_bottom_input, accumulate_input_1, accumulate_input_2                    
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS, SIS_call_tracer_column_fns
use SIS_tracer_registry, only : SIS_unpack_passive_ice_tr, SIS_repack_passive_ice_tr
use SIS_tracer_registry, only : SIS_count_passive_tracers
use SIS_transport,     only : adjust_ice_categories, SIS_transport_CS
use SIS_types,         only : ice_state_type, IST_chksum, IST_bounds_check, total_sfc_flux_type                  
use SIS_types,         only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_utils,         only : post_avg
use SIS2_ice_thm,      only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,      only : ice_resize_SIS2, add_frazil_SIS2, rebalance_ice_layers
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,      only : enth_from_TS, Temp_from_En_S, enthalpy_liquid, calculate_T_freeze

implicit none; private

#include <SIS2_memory.h>

public initialize_isponge, apply_isponge, set_up_isponge_field, SIS_sponge_end

!> A structure for creating arrays of pointers to 3D arrays
type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3D array [various]
end type p3d
!> A structure for creating arrays of pointers to 2D arrays
type, public :: p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2D array [various]
end type p2d
 
!> This control structure holds memory and parameters for the SIS_sponge module
type, public :: isponge_CS ; private
  logical, public :: use_isponge = .false.  !< If true, ice tracer fields may be relaxed somewhere in the domain
  integer :: num_col    !< The number of sponge points within the computational domain.
  integer :: fldno = 0  !< The number of fields which have already been
                        !! registered by calls to set_up_sponge_field
  integer, pointer :: col_i(:) => NULL() !< Array of the i-indicies of each of the columns being damped.
  integer, pointer :: col_j(:) => NULL() !< Array of the j-indicies of each of the columns being damped.
  real, pointer :: Iresttime_col(:) => NULL() !< The inverse restoring time of each column [T-1 ~> s-1].
  type(p3d) :: var(MAX_FIELDS_RLX_)     !< Pointers to the fields that will be relaxed
  type(p2d) :: Ref_val(MAX_FIELDS_RLX_) !< The values to which the fields are relaxed
end type isponge_CS

contains

!> This subroutine determines the number of points which are within ice sponges in
!! this computational domain.  Only points that have positive values of
!! Iresttime and which mask2dT indicates are ocean points are included in the
!! sponges.  
!subroutine initialize_isponge(Iresttime, IST, G, IG, param_file, CS)
subroutine initialize_isponge(param_file, Iresttime, G, IG, CS)
!  type(ice_state_type),    intent(in) :: IST        !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: Iresttime  !< The inverse of the restoring time [T-1 ~> s-1].
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  logical :: use_isponge
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat
  integer :: col, total_isponge_cols

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  write(mesg,'("SIS_sponge: isc/iec=",I4,"/",I4," jsc/jec=",I4,"/",I4," icencat=",I3)') &
       isc, iec, jsc, jec, ncat 
  call SIS_mesg(mesg)

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

    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo

  endif

  call SIS_mesg("Calling sum_across_PEs")
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
subroutine set_up_isponge_field(sp_val, f_ptr, G, IG, CS, ncat)
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  integer,                 intent(in) :: ncat
  real, dimension(ncat, SZI_(G),SZJ_(G)), &
                           intent(in) :: sp_val     !< The reference profiles of the quantity 
                                                    !! being registered [various]
  real, dimension(ncat, SZI_(G),SZJ_(G)), &
                   target, intent(in) :: f_ptr      !< a pointer to the field which will be relaxed [various]
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

  allocate(CS%Ref_val(CS%fldno)%p(IG%CatIce,CS%num_col), source=0.0)
  do col=1,CS%num_col
    do k=1,IG%CatIce
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

end subroutine set_up_isponge_field


!> This subroutine applies relaxation ("damping") to ice thickness (by categories) and ice concentration
!! tracers for every column where the relaxation time scale > 0.
subroutine apply_isponge(dt, CS, IG, IST)
  real,                    intent(in)  :: dt    !< The amount of time covered by this call [T ~> s].
  type(isponge_CS),        pointer     :: CS    !< A pointer that is set to point to the ice sponge control
                                                !! structure for this module
  type(ice_grid_type),     intent(in)  :: IG    !< The sea-ice specific grid type
  type(ice_state_type), intent(inout)  :: IST   !< A type describing the state of the sea ice
   

  ! Local variables
!  real :: H_rescale_ice, H_rescale_snow
  real :: damp         ! The timestep times the local damping coefficient [nondim].
  real :: I1pdamp      ! I1pdamp is 1/(1 + damp). [nondim]
  real :: damp_1pdamp  ! damp_1pdamp is damp/(1 + damp). [nondim]
  real :: Idt          ! The inverse of the timestep [T-1 ~> s-1]
  character(len=40)  :: mdl = "SIS_sponge"  ! This module's name.
  character(len=256) :: mesg
  integer :: c, i, j, k, m

  ! Determine the thickness rescaling factors that are needed.
!  H_rescale_ice = 1.0 ; H_rescale_snow = 1.0

  write(mesg,'("apply_isponge: dt=",F16.2," sec, num_col=",I5)') dt, CS%num_col
  call SIS_mesg(mesg)
  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)             
    damp = dt * CS%Iresttime_col(c); I1pdamp = 1.0 / (1.0 + damp)

    do k=1,IG%CatIce
      do m=1,CS%fldno
        CS%var(m)%p(i,j,k) = I1pdamp * &
           (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
      enddo
    enddo

  enddo

end subroutine apply_isponge

!> Deallocate memory associated with the SIS_optics module
subroutine SIS_sponge_end(CS)
  type(isponge_CS), pointer :: CS !< The ice sponge control structure that is deallocated here
  
  deallocate(CS) 
                 
end subroutine SIS_sponge_end

end module SIS_sponge 

