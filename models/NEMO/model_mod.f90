! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

use        types_mod,        only : r8, i8, MISSING_R8, MISSING_I, obstypelength

use time_manager_mod,        only : time_type, set_time, set_date, &
                                    print_time, print_date,        &
                                    set_calendar_type,             &
                                    operator(+), operator(-)

use     location_mod,        only : location_type, get_close_type,          &
                                    loc_get_close_obs => get_close_obs,     &
                                    loc_get_close_state => get_close_state, &
                                    set_location, set_location_missing,     &
                                    get_location, write_location, is_vertical, &
                                    VERTISHEIGHT, VERTISSURFACE

use    utilities_mod,        only : register_module, error_handler, &
                                    E_ERR, E_MSG, file_exist,  &
                                    nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                                    to_upper, &
                                    find_namelist_in_file, check_namelist_read

use netcdf_utilities_mod,    only : nc_add_global_attribute, nc_synchronize_file, &
                                    nc_add_global_creation_time, &
                                    nc_begin_define_mode, nc_end_define_mode, &
                                    nc_open_file_readonly, nc_close_file, nc_check


use state_structure_mod,     only : add_domain, get_domain_size, &
                                    state_structure_info,        &
                                    get_variable_name,           &
                                    get_varid_from_kind,         &
                                    get_model_variable_indices,  &
                                    get_dart_vector_index, get_kind_index

use distributed_state_mod, only : get_state

use obs_kind_mod,            only : get_index_for_quantity

use ensemble_manager_mod,    only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod,       only : pert_model_copies, &
                                    init_time => fail_init_time, &
                                    init_conditions => fail_init_conditions, &
                                    convert_vertical_obs, convert_vertical_state, adv_1step

use dart_time_io_mod,      only : write_model_time

use netcdf

    implicit none
    private

    ! routines required by DART code - will be called from filter and other
    ! DART executables.
    public :: get_model_size,         &
              get_state_meta_data,    &
              model_interpolate,      &
              end_model,              &
              static_init_model,      &
              nc_write_model_atts,    &
              get_close_obs,          &
              get_close_state,        &
              pert_model_copies,      &
              convert_vertical_obs,   &
              convert_vertical_state, &
              read_model_time,        &
              adv_1step,              &
              init_time,              &
              init_conditions,        &
              shortest_time_between_assimilations, &
              write_model_time

    character(len=32)  :: calendar = 'Gregorian'
    character(len=256) :: model_analysis_filename = '1d_grid_T.nc'
    character(len=256) :: grid_filename = 'mesh_mask.nc'

    character(len=256), parameter :: source   = "model_mod.f90"
    logical :: module_initialized = .false.
    integer :: dom_id ! used to access the state structure
    type(time_type) :: model_time, model_timestep
    type(time_type) :: assimilation_time_step

    ! Example Namelist
    ! Use the namelist for options to be set at runtime.
    character(len=256) :: template_file = 'model_restart.nc'
    integer  :: time_step_days      = 0
    integer  :: time_step_seconds   = 3600

    integer, parameter :: MAX_STATE_VARIABLES = 3
    integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5
    integer, parameter :: VARNAME_INDEX = 1
    integer, parameter ::    KIND_INDEX = 2
    integer, parameter ::  MINVAL_INDEX = 3
    integer, parameter ::  MAXVAL_INDEX = 4
    integer, parameter :: REPLACE_INDEX = 5

    character(len=obstypelength) :: var_names(max_state_variables)
    real(r8) :: var_ranges(max_state_variables,2)
    logical  :: var_update(max_state_variables)
    integer  :: var_qtys(  max_state_variables)


    character(len=NF90_MAX_NAME) :: nemo_variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ''

    namelist /model_nml/    &
             template_file, time_step_days, time_step_seconds, &
             model_analysis_filename, grid_filename, nemo_variables

    ! Everything needed to describe a variable

    integer :: FVAL=-999.0 !SIVA: The FVAL is the fill value used for input netcdf files.
    integer :: nfields

    integer :: domain_id ! global variable for state_structure_mod routines
    integer :: nav_lon_x, nav_lat_y, nav_lev_z, time_t
    real(r8), allocatable, target :: nav_lon(:,:), nav_lat(:,:), deptht(:)
    character(len=512) :: string1, string2, string3

    character(len=32 ), parameter :: revision = "$Revision$"
    character(len=128), parameter :: revdate  = "$Date$"

    integer            :: debug = 0   ! turn up for more and more debug messages


    contains

    !------------------------------------------------------------------
    !
    ! Called to do one time initialization of the model. As examples,
    ! might define information about the model size or model timestep.
    ! In models that require pre-computed static data, for instance
    ! spherical harmonic weights, these would also be computed here.

    subroutine static_init_model()

    integer  :: ncid
    integer  :: iunit, io
    integer  :: model_size

    module_initialized = .true.

    ! Print module information to log file and stdout.
    call register_module(source)

    call find_namelist_in_file("input.nml", "model_nml", iunit)
    read(iunit, nml = model_nml, iostat = io)
    call check_namelist_read(iunit, io, "model_nml")

    ! Record the namelist values used for the run
    if (do_nml_file()) write(nmlfileunit, nml=model_nml)
    if (do_nml_term()) write(     *     , nml=model_nml)

    ! This time is both the minimum time you can ask the model to advance
    ! (for models that can be advanced by filter) and it sets the assimilation
    ! window.  All observations within +/- 1/2 this interval from the current
    ! model time will be assimilated. If this is not settable at runtime
    ! feel free to hardcode it and remove from the namelist.
    assimilation_time_step = set_time(time_step_seconds, &
                                      time_step_days)


    ncid = nc_open_file_readonly(model_analysis_filename,'static_init_model')

    call parse_variable_input(nemo_variables, ncid, model_analysis_filename, nfields)

    call nc_close_file(ncid, 'static_init_model', model_analysis_filename)

    call get_grid_dimensions(model_analysis_filename)

    call get_grid(model_analysis_filename)

    ! Define which variables are in the model state
    dom_id = add_domain( model_analysis_filename, nfields, &
                var_names  = var_names( 1:nfields),        &
                kind_list  = var_qtys(  1:nfields),        &
                clamp_vals = var_ranges(1:nfields,:),      &
                update_list= var_update(1:nfields)        )

    if ( debug > 4 .and. do_output()) call state_structure_info(dom_id)

    model_size = get_domain_size(dom_id)


    call set_calendar_type( calendar )

    ! tell the location module how we want to localize in the vertical
!    call set_vertical_localization_coord(vert_localization_coord)

    model_time = read_model_time(model_analysis_filename)

    end subroutine static_init_model

    !------------------------------------------------------------------
    ! Returns the number of items in the state vector as an integer.

    function get_model_size()

    integer(i8) :: get_model_size

    if ( .not. module_initialized ) call static_init_model

    get_model_size = get_domain_size(dom_id)

    end function get_model_size

    !-----------------------------------------------------------------------
    !

    subroutine model_interpolate(state_handle, ens_size, location, quantity, interp_val, istatus)

    type(ensemble_type),   intent(in)  :: state_handle
    integer,               intent(in)  :: ens_size
    type(location_type),   intent(in)  :: location
    integer,               intent(in)  :: quantity
    real(r8),              intent(out) :: interp_val(ens_size)
    integer,               intent(out) :: istatus(ens_size)

    ! model_interpolate will interpolate any variable in the DART vector to the given location.
    ! The first variable matching the quantity of interest will be used for the interpolation.
    !
    ! istatus =  0 ... success
    ! istatus =  1 ... unknown model level
    ! istatus =  2 ... vertical coordinate unsupported
    ! istatus =  3 ... quantity not in the DART vector
    ! istatus =  4 ... vertical  interpolation failed
    ! istatus = 11 ... latitude  interpolation failed
    ! istatus = 12 ... longitude interpolation failed
    ! istatus = 13 ... corner a retrieval failed
    ! istatus = 14 ... corner b retrieval failed
    ! istatus = 15 ... corner c retrieval failed
    ! istatus = 16 ... corner d retrieval failed

    ! Local storage
    real(r8)    :: loc_array(3), llon, llat, lheight
    integer(i8) :: base_offset, offset
    integer     :: ind
    integer     :: hgt_bot, hgt_top
    real(r8)    :: hgt_fract
    real(r8)    :: top_val(ens_size), bot_val(ens_size)
    integer     :: hstatus
    integer     :: i, varid

    if ( .not. module_initialized ) call static_init_model

    ! Successful istatus is 0
    interp_val = MISSING_R8
    istatus    = 0

    ! Get the individual locations values
    loc_array = get_location(location)
    llon      = loc_array(1)
    llat      = loc_array(2)
    lheight   = loc_array(3)

    if( is_vertical(location,"HEIGHT") ) then
       ! Nothing to do
    elseif ( is_vertical(location,"SURFACE") ) then
       ! Nothing to do
    elseif (is_vertical(location,"LEVEL")) then
       ! convert the level index to an actual depth
       ind = nint(loc_array(3))
       if ( (ind < 1) .or. (ind > size(deptht)) ) then
          lheight = deptht(ind)
       else
          istatus = 1
          return
       endif
    else   ! if pressure or undefined, we don't know what to do
       istatus = 2
       return
    endif

    ! determine which variable is the desired QUANTITY
    varid = get_varid_from_kind(dom_id, quantity)

    if (varid < 1) then
       istatus = 3
       return
    endif

    ! Do horizontal interpolations for the appropriate levels

    ! For Sea Surface Height don't need the vertical coordinate
    if( is_vertical(location,"SURFACE") ) then
       call lat_lon_interpolate(state_handle, ens_size, llon, llat, 1, varid, quantity, interp_val, istatus)
       return
    endif

    ! Get the bounding vertical levels and the fraction between bottom and top
    call height_bounds(lheight, nav_lev_z, deptht, hgt_bot, hgt_top, hgt_fract, hstatus)
    if(hstatus /= 0) then
       istatus = 4
       return
    endif

    call lat_lon_interpolate(state_handle, ens_size, llon, llat, hgt_top, varid, quantity, top_val, istatus)
    ! Failed istatus from interpolate means give up
    do i =1,ens_size
       if(istatus(i) /= 0) return
    enddo

    call lat_lon_interpolate(state_handle, ens_size, llon, llat, hgt_bot, varid, quantity, bot_val, istatus)
    ! Failed istatus from interpolate means give up
    do i =1,ens_size
       if(istatus(i) /= 0) return
    enddo
    ! Then weight them by the fraction and return
    interp_val = bot_val + hgt_fract * (top_val - bot_val)

    end subroutine model_interpolate

!    subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)
!    !------------------------------------------------------------------
!    ! Given a state handle, a location, and a state quantity,
!    ! interpolates the state variable fields to that location and returns
!    ! the values in expected_obs. The istatus variables should be returned as
!    ! 0 unless there is some problem in computing the interpolation in
!    ! which case a positive istatus should be returned.
!    !
!    ! For applications in which only perfect model experiments
!    ! with identity observations (i.e. only the value of a particular
!    ! state variable is observed), this can be a NULL INTERFACE.
!
!! given a state vector, a location, and a QTY_xxx, return the
!! value at from the closest grid location and an error code.  0 is success,
!! anything positive is an error.  (negative reserved for system use)
!! TODO: possibly interpolate ...
!!
!!       ERROR codes:
!!
!!       ISTATUS = 99:  general error in case something terrible goes wrong...
!!       ISTATUS = 88:  this kind is not in the state vector
!!       ISTATUS = 11:  Could not find a triangle that contains this lat/lon
!!       ISTATUS = 12:  Height vertical coordinate out of model range.
!!       ISTATUS = 13:  Missing value in interpolation.
!!       ISTATUS = 16:  Don't know how to do vertical velocity for now
!!       ISTATUS = 17:  Unable to compute pressure values
!!       ISTATUS = 18:  altitude illegal
!!       ISTATUS = 19:  could not compute u using RBF code
!!       ISTATUS = 101: Internal error; reached end of subroutine without
!!                      finding an applicable case.
!!
!
!! passed variables
!
!    type(ensemble_type), intent(in) :: state_handle
!    integer,             intent(in) :: ens_size
!    type(location_type), intent(in) :: location
!    integer,             intent(in) :: qty
!    real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
!    integer,            intent(out) :: istatus(ens_size)
!
!
!    real(r8) :: llv(3)    ! lon/lat/vert
!    real(r8) :: lon, lat, vert
!    if ( .not. module_initialized ) call static_init_model
!
!    ! Decode the location into bits for error messages ...
!    llv  = get_location(location)
!    lon  = llv(1)    ! degrees East [0,360)
!    lat  = llv(2)    ! degrees North [-90,90]
!    vert = llv(3)    ! depth in meters ... even 2D fields have a value of 0.0
!    print*, "lon, lat, vert:", lon, lat, vert
!
!
!!    surface_index = find_closest_surface_location(location, obs_kind)
!!    if (surface_index < 1) then ! nothing close
!!       istatus = 11
!!       return
!!    endif
!
!    ! This should be the result of the interpolation of a
!    ! given kind (itype) of variable at the given location.
!    expected_obs(:) = MISSING_R8
!
!    ! istatus for successful return should be 0.
!    ! Any positive number is an error.
!    ! Negative values are reserved for use by the DART framework.
!    ! Using distinct positive values for different types of errors can be
!    ! useful in diagnosing problems.
!    istatus = 0
!
!    end subroutine model_interpolate



subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)
!=======================================================================
!
real(r8),             intent(in) :: lheight
integer,              intent(in) :: nheights
real(r8),             intent(in) :: hgt_array(nheights)
integer,             intent(out) :: bot, top
real(r8),            intent(out) :: fract
integer,             intent(out) :: istatus

! Local variables
integer   :: i

if ( .not. module_initialized ) call static_init_model

! Succesful istatus is 0
istatus = 0

! The deptht array contains the depths of the center of the vertical grid boxes

! It is assumed that the top box is shallow and any observations shallower
! than the depth of this boxes center are just given the value of the
! top box.
if(lheight > hgt_array(1)) then
   top = 1
   bot = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   ! ??? Make sure this is consistent with the interpolation
   fract = 1.0_r8
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is shallower than this entry, it must be in this box
   if(lheight >= hgt_array(i)) then
      top = i -1
      bot = i
      fract = (lheight - hgt_array(bot)) / (hgt_array(top) - hgt_array(bot))
      return
   endif
end do

! Falling off the end means the location is lower than the deepest height
! Fail with istatus 2 in this case
istatus = 2

end subroutine height_bounds


subroutine lat_lon_interpolate(state_handle, ens_size, llon, llat, level, var_id, qty, interp_val, istatus)
!=======================================================================
!

! Subroutine to interpolate to a lat lon location for a given level

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size

real(r8),            intent(in) :: llon, llat
integer,             intent(in) :: level
integer,             intent(in) :: var_id, qty
integer,            intent(out) :: istatus(ens_size)
real(r8),           intent(out) :: interp_val(ens_size)

! Local storage
real(r8) :: lat_array(nav_lat_y), lon_array(nav_lon_x)
integer  :: lat_bot, lat_top, lon_bot, lon_top
real(r8) :: lat_fract, lon_fract
real(r8),dimension(ens_size) :: pa, pb, pc, pd, xbot, xtop
integer  :: lat_status, lon_status
logical  :: masked

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
istatus = 0

! Find out what latitude box and fraction
! The latitude grid being used depends on the variable type
! V is on the YG latitude grid

lat_array = nav_lat(1,:)
!if(qty == QTY_V_CURRENT_COMPONENT) lat_array = yg

call lat_bounds(llat, nav_lat_y, lat_array, lat_bot, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then
   istatus = 11
   return
endif

! Find out what longitude box and fraction
lon_array = nav_lon(:,1)
!if(qty == QTY_U_CURRENT_COMPONENT) lon_array = xg

call lon_bounds(llon, nav_lon_x, lon_array, lon_bot, lon_top, lon_fract, lon_status)

! Check for error on the longitude interpolation
if(lon_status /= 0) then
   istatus = 12
   return
endif


! Vector is laid out with lat outermost loop, lon innermost loop
! Find the bounding points for the lat lon box
! NOTE: For now, it is assumed that a real(r8) value of exactly 0.0 indicates
! that a particular gridded quantity is masked and not available. This is not
! the most robust way to do this, but may be sufficient since exact 0's are
! expected to happen rarely. Jeff Anderson believes that the only implication
! will be that an observation whos forward operator requires interpolating
! from a point that has exactly 0.0 (but is not masked) will not be
! assimilated.

pa = get_val(lon_bot, lat_bot, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 13
   return
endif
pb = get_val(lon_top, lat_bot, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 14
   return
endif
pc = get_val(lon_bot, lat_top, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 15
   return
endif
pd = get_val(lon_top, lat_top, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 16
   return
endif

xbot = pa + lon_fract * (pb - pa)
xtop = pc + lon_fract * (pd - pc)
interp_val = xbot + lat_fract * (xtop - xbot)

end subroutine lat_lon_interpolate




subroutine lat_bounds(llat, nlats, lat_array, bot, top, fract, istatus)

!=======================================================================
!

! Given a latitude llat, the array of latitudes for grid boundaries, and the
! number of latitudes in the grid, returns the indices of the latitude
! below and above the location latitude and the fraction of the distance
! between. istatus is returned as 0 unless the location latitude is
! south of the southernmost grid point (1 returned) or north of the
! northernmost (2 returned). If one really had lots of polar obs would
! want to worry about interpolating around poles.

real(r8),          intent(in) :: llat
integer,           intent(in) :: nlats
real(r8),          intent(in) :: lat_array(nlats)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer    :: i

if ( .not. module_initialized ) call static_init_model

! Default is success
istatus = 0

! Check for too far south or north
if(llat < lat_array(1)) then
   istatus = 1
   return
else if(llat > lat_array(nlats)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, nlats
   if(llat <= lat_array(i)) then
      bot = i - 1
      top = i
      fract = (llat - lat_array(bot)) / (lat_array(top) - lat_array(bot))
      return
   endif
end do

end subroutine lat_bounds



subroutine lon_bounds(llon, nlons, lon_array, bot, top, fract, istatus)

!=======================================================================
!

! Given a longitude llon, the array of longitudes for grid boundaries, and the
! number of longitudes in the grid, returns the indices of the longitude
! below and above the location longitude and the fraction of the distance
! between. istatus is returned as 0 unless the location longitude is
! not between any of the longitude box boundaries. This should be modified
! for global wrap-around grids.
! Algorithm fails for a silly grid that
! has only two longitudes separated by 180 degrees.

real(r8),          intent(in) :: llon
integer,           intent(in) :: nlons
real(r8),          intent(in) :: lon_array(nlons)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

if ( .not. module_initialized ) call static_init_model

! Default is success
istatus = 0

! This is inefficient, someone could clean it up
! Plus, it doesn't work for a global model that wraps around
do i = 2, nlons
   dist_bot = lon_dist(llon, lon_array(i - 1))
   dist_top = lon_dist(llon, lon_array(i))
   if(dist_bot >= 0 .and. dist_top < 0) then
      bot = i - 1
      top = i
      fract = dist_bot / (dist_bot + abs(dist_top))
      ! orig: fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
      return
   endif
end do

! Falling off the end means its in between. Add the wraparound check.
! For now, return istatus 1
istatus = 1

end subroutine lon_bounds



function lon_dist(lon1, lon2)
!=======================================================================
!

! Returns the smallest signed distance between lon1 and lon2 on the sphere
! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

real(r8), intent(in) :: lon1, lon2
real(r8)             :: lon_dist

if ( .not. module_initialized ) call static_init_model

lon_dist = lon1 - lon2
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then
   return
else if(lon_dist < -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist


function get_val(lon_index, lat_index, level, var_id, state_handle,ens_size, masked)
!=======================================================================
!

! Returns the value from a single level array given the lat and lon indices
integer,             intent(in)  :: lon_index, lat_index, level
integer,             intent(in)  :: var_id ! state variable
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
logical,             intent(out) :: masked
real(r8)                         :: get_val(ens_size)

integer(i8) :: state_index
integer :: i

if ( .not. module_initialized ) call static_init_model

state_index = get_dart_vector_index(lon_index, lat_index, level, dom_id, var_id)
get_val = get_state(state_index,state_handle)

! Masked returns false if the value is masked
! A grid variable is assumed to be masked if its value is FVAL.
! Just to maintain legacy, we also assume that A grid variable is assumed
! to be masked if its value is exactly 0.
! See discussion in lat_lon_interpolate.

! MEG CAUTION: THE ABOVE STATEMENT IS INCORRECT
! trans_mitdart already looks for 0.0 and makes them FVAL
! So, in the condition below we don't need to check for zeros
! The only mask is FVAL
masked = .false.
do i=1,ens_size
!   if(get_val(i) == FVAL .or. get_val(i) == 0.0_r8 ) masked = .true.
    if(get_val(i) == FVAL) masked = .true.
enddo

end function get_val


    !------------------------------------------------------------------
    ! Returns the smallest increment in time that the model is capable
    ! of advancing the state in a given implementation, or the shortest
    ! time you want the model to advance between assimilations.

    function shortest_time_between_assimilations()

    type(time_type) :: shortest_time_between_assimilations

    if ( .not. module_initialized ) call static_init_model

    shortest_time_between_assimilations = assimilation_time_step

    end function shortest_time_between_assimilations



    !------------------------------------------------------------------
    ! Given an integer index into the state vector, returns the
    ! associated location and optionally the physical quantity.

    subroutine get_state_meta_data(index_in, location, qty)

    integer(i8),         intent(in)  :: index_in
    type(location_type), intent(out) :: location
    integer,             intent(out), optional :: qty

    integer  :: iloc, vloc, jloc
    integer  :: myvarid, mydomid, myqty


    if ( .not. module_initialized ) call static_init_model


    iloc = -1

    call get_model_variable_indices(index_in, iloc, jloc, vloc, &
        var_id=myvarid, dom_id=mydomid, kind_index=myqty)


    if( iloc == -1 ) then
         write(string1,*) 'Problem, cannot find base_offst, index_in is: ', index_in
         call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
    endif

    location = set_location(nav_lon(iloc, 1), nav_lat(1, jloc), deptht(vloc), VERTISHEIGHT)

    if (present(qty)) then
       qty = myqty
       print*, qty

       if( qty == MISSING_I ) then
          write(string1,*) 'Cannot find DART QTY for indx ', index_in
          write(string2,*) 'variable "'//trim(get_variable_name(mydomid, myvarid))//'"'
          call error_handler(E_ERR, 'get_state_meta_data', string1, &
                     source, revision, revdate, text2=string2)
       endif
    endif

    end subroutine get_state_meta_data

    subroutine get_grid_dimensions(filename)

    character(len=200), intent(in) :: filename
    integer  :: ncid

    ncid = nc_open_file_readonly(filename,'static_init_model')

    call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid), &
                   'get_grid_dimensions', 'open '//trim(filename))

    nav_lon_x   = get_dimension_length(ncid, 'x',   filename)
    nav_lat_y   = get_dimension_length(ncid, 'y',   filename)
    nav_lev_z   = get_dimension_length(ncid, 'deptht',   filename)
    time_t      = get_dimension_length(ncid, 'time_counter',   filename)

    print*, 'grid_dims, x, y, z, t: ', nav_lon_x, nav_lat_y, nav_lev_z, time_t
!    call parse_variable_input(nemo_variables, ncid, filename, nfields)

    call nc_close_file(ncid, 'get_grid_dimensions', filename)

    end subroutine get_grid_dimensions

    subroutine get_grid(filename)

    integer    :: ncid, VarID
    character(len=200), intent(in) :: filename

    call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid), &
        'get_grid', 'open '//trim(filename))

    if (.not. allocated(nav_lon)) allocate(nav_lon(nav_lon_x, nav_lat_y))
    call nc_check(nf90_inq_varid(ncid, 'nav_lon', VarID), &
        'get_grid', 'inq_varid nav_lon'//trim(filename))

    call nc_check(nf90_get_var( ncid, VarID, nav_lon), &
        'get_grid', 'get_var nav_lon '//trim(filename))

    if (.not. allocated(nav_lat)) allocate(nav_lat(nav_lon_x, nav_lat_y))

    call nc_check(nf90_inq_varid(ncid, 'nav_lat', VarID), &
        'get_grid', 'inq_varid nav_lat'//trim(filename))

    call nc_check(nf90_get_var( ncid, VarID, nav_lat), &
        'get_grid', 'get_var nav_lat '//trim(filename))

    if (.not. allocated(deptht)) allocate(deptht(nav_lev_z))

    call nc_check(nf90_inq_varid(ncid, 'deptht', VarID), &
        'get_grid', 'inq_varid deptht'//trim(filename))

    call nc_check(nf90_get_var( ncid, VarID, deptht), &
        'get_grid', 'get_var deptht '//trim(filename))

    ! nemo example file has longitude < 0
    ! DART uses [0,360]
    where(nav_lon < 0.0_r8 )   nav_lon   = nav_lon + 360.0_r8


    end subroutine get_grid

    !-----------------------------------------------------------------------
    !>
    !> gets the length of a netCDF dimension given the dimension name.
    !> This bundles the nf90_inq_dimid and nf90_inquire_dimension routines
    !> into a slightly easier-to-use function.
    !>
    !> @param dimlen the length of the netCDF dimension in question
    !> @param ncid the netCDF file handle
    !> @param dimension_name the character string of the dimension name
    !> @param filename the name of the netCDF file (for error message purposes)
    !>

    function get_dimension_length(ncid, dimension_name, filename) result(dimlen)

    integer                      :: dimlen
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimension_name
    character(len=*), intent(in) :: filename

    integer :: DimID

    write(string1,*)'inq_dimid '//trim(dimension_name)//' '//trim(filename)
    write(string2,*)'inquire_dimension '//trim(dimension_name)//' '//trim(filename)

    call nc_check(nf90_inq_dimid(ncid, trim(dimension_name), DimID), &
                  'get_dimension_length',string1)
    call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
                  'get_dimension_length', string2)

    end function get_dimension_length




    subroutine parse_variable_input( state_variables, ncid, filename, ngood )

    character(len=*), intent(in)  :: state_variables(:,:)
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: ngood

    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
    character(len=NF90_MAX_NAME) :: dimname
    integer :: i, j, VarID, dimlen, numdims
    logical :: failure

    character(len=NF90_MAX_NAME) :: varname       ! column 1
    character(len=NF90_MAX_NAME) :: dartstr       ! column 2
    character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
    character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
    character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5

    real(r8) :: minvalue, maxvalue
    integer  :: ios

    if ( .not. module_initialized ) call static_init_model

    var_names  = 'no_variable_specified'
    var_qtys   = MISSING_I
    var_ranges = MISSING_R8
    var_update = .false.

    failure = .FALSE. ! perhaps all with go well

    ngood = 0
    MyLoop : do i = 1, MAX_STATE_VARIABLES

       if ( nemo_variables(1,i) == ' ' .and. nemo_variables(2,i) == ' ' ) exit MyLoop ! Found end of list.

       if ( any(state_variables(:,i) == ' ') ) then
          string1 = '...  model_nml:"variables" not fully specified'
          write(string2,*)'failing on line ',i
          call error_handler(E_ERR, 'parse_variable_input', string1, &
                     source, revision, revdate, text2=string2)
       endif

       varname      = trim(state_variables(VARNAME_INDEX,i))
       dartstr      = trim(state_variables(   KIND_INDEX,i))
       minvalstring = trim(state_variables( MINVAL_INDEX,i))
       maxvalstring = trim(state_variables( MAXVAL_INDEX,i))
       state_or_aux = trim(state_variables(REPLACE_INDEX,i))
       call to_upper(state_or_aux)


       ! Make sure DART kind is valid

       if( get_index_for_quantity(dartstr) < 0 ) then
          write(string1,'(''there is no quantity <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
          call error_handler(E_ERR,'parse_variable_input:',string1,source,revision,revdate)
       endif

       var_names(i) = trim(varname)
       var_qtys(i)  = get_index_for_quantity(dartstr)
       print*, 'var_qtys(i): ', var_qtys(i)

       read(minvalstring,*,iostat=ios) minvalue
       if (ios == 0) var_ranges(i,1) = minvalue

       read(maxvalstring,*,iostat=ios) maxvalue
       if (ios == 0) var_ranges(i,1) = maxvalue

       if (state_or_aux == 'UPDATE' ) var_update(i) = .true.

       ! Make sure DART kind is valid

       if( var_qtys(i) < 0 ) then
          write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') &
                trim(dartstr)
          call error_handler(E_ERR,'parse_variable_input',string1,source,revision,revdate)
       endif

       ! Make sure variable exists in model analysis variable list

       write(string1,'(''variable '',a,'' in '',a)') trim(varname), trim(filename)
       write(string2,'(''there is no '',a)') trim(string1)
       call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                     'parse_variable_input', trim(string2))

       ! Make sure variable is defined by (Time,nCells) or (Time,nCells,vertical)
       ! unable to support Edges or Vertices at this time.

       call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
                     'parse_variable_input', 'inquire '//trim(string1))

       DimensionLoop : do j = 1,numdims

          write(string2,'(''inquire dimension'',i2,'' of '',a)') j,trim(string1)
          call nc_check(nf90_inquire_dimension(ncid, dimIDs(j), len=dimlen, name=dimname), &
                                              'parse_variable_input', trim(string2))
          select case ( trim(dimname) )
             case ('time_counter')
                ! supported - do nothing
             case ('x')
                ! supported - do nothing
             case ('y')
                ! supported - do nothing
             case ('deptht')
                ! supported - do nothing
             case ('axis_nbounds')
                ! supported - do nothing
             case default
                write(string2,'(''unsupported dimension '',a,'' in '',a)') trim(dimname),trim(string1)
                call error_handler(E_MSG,'parse_variable_input',string2,source,revision,revdate)
                failure = .TRUE.
          end select

       enddo DimensionLoop

       if (failure) then
           string2 = 'unsupported dimension(s) are fatal'
           call error_handler(E_ERR,'parse_variable_input',string2,source,revision,revdate)
       endif

       ! Record the contents of the DART state vector

       if (debug > 0) then
          write(string1,*)'variable ',i,' is ',trim(varname), ' ', trim(dartstr)
          call error_handler(E_MSG,'parse_variable_input',string1)
       endif

       ngood = ngood + 1
    enddo MyLoop

    if (ngood == MAX_STATE_VARIABLES) then
       string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
       write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
       call error_handler(E_MSG,'parse_variable_input',string1,text2=string2)
    endif

    end subroutine parse_variable_input

    subroutine get_time_information
    end subroutine get_time_information

    !------------------------------------------------------------------
    ! Any model specific distance calcualtion can be done here
    subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                             num_close, close_ind, dist, ens_handle)

    type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
    integer,                       intent(in)    :: base_type     ! observation TYPE
    type(location_type),           intent(inout) :: base_loc      ! location of interest
    type(location_type),           intent(inout) :: locs(:)       ! obs locations
    integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
    integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
    integer,                       intent(out)   :: num_close     ! how many are close
    integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
    real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
    type(ensemble_type), optional, intent(in)    :: ens_handle

    character(len=*), parameter :: routine = 'get_close_obs'

    call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                              num_close, close_ind, dist, ens_handle)

    end subroutine get_close_obs


    !------------------------------------------------------------------
    ! Any model specific distance calcualtion can be done here
    subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                               num_close, close_ind, dist, ens_handle)

    type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
    type(location_type),           intent(inout) :: base_loc     ! location of interest
    integer,                       intent(in)    :: base_type    ! observation TYPE
    type(location_type),           intent(inout) :: locs(:)      ! state locations
    integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
    integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
    integer,                       intent(out)   :: num_close    ! how many are close
    integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
    real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
    type(ensemble_type), optional, intent(in)    :: ens_handle

    character(len=*), parameter :: routine = 'get_close_state'


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)


end subroutine get_close_state


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

	if (allocated(nav_lon)) deallocate(nav_lon)
	if (allocated(nav_lat)) deallocate(nav_lat)
	if (allocated(deptht))  deallocate(deptht)


end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type) :: read_model_time
type(time_type) :: january1
real(r8)        :: time_sec, minutes, hours
integer         :: VarID, ncid, seconds, days

call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid), &
    'read_model_time', 'open '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'time_counter', VarID), &
    'read_model_time', 'inq_varid time_counter '//trim(filename))

call nc_check(nf90_get_var( ncid, VarID, time_sec), &
    'read_model_time', 'get_var time_counter '//trim(filename))

minutes = time_sec  / 60
hours   = minutes   / 60
days    = int(hours / 24)
seconds = mod(time_sec,real(days))
print*, 'seconds, days: ', seconds, days
january1 = set_date(1900,1,1)

read_model_time = january1 + set_time(seconds,days)

call print_time(read_model_time,'read_model_time:')
call print_date(read_model_time,'read_model_time:')


end function read_model_time

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

