!> @file nesting_offl_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Offline nesting in larger-scale models. Boundary conditions for the simulation are read from
!> NetCDF file and are prescribed onto the respective arrays.
!> Further, a mass-flux correction is performed to maintain the mass balance.
!--------------------------------------------------------------------------------------------------!
 MODULE nesting_offl_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  diss,                                                                               &
               drho_air_zw,                                                                        &
               dzw,                                                                                &
               e,                                                                                  &
               pt,                                                                                 &
               pt_init,                                                                            &
               q,                                                                                  &
               q_init,                                                                             &
               rdf,                                                                                &
               rdf_sc,                                                                             &
               rho_air,                                                                            &
               rho_air_zw,                                                                         &
               s,                                                                                  &
               u,                                                                                  &
               u_init,                                                                             &
               ug,                                                                                 &
               v,                                                                                  &
               v_init,                                                                             &
               vg,                                                                                 &
               w,                                                                                  &
               zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  g,                                                                                  &
               pi

    USE chem_modules,                                                                              &
        ONLY:  chem_species, nesting_offline_chem

    USE control_parameters,                                                                        &
        ONLY:  air_chemistry,                                                                      &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               coupling_char,                                                                      &
               constant_diffusion,                                                                 &
               child_domain,                                                                       &
               debug_output_timestep,                                                              &
               dt_3d,                                                                              &
               dz,                                                                                 &
               end_time,                                                                           &
               humidity,                                                                           &
               initializing_actions,                                                               &
               message_string,                                                                     &
               nesting_offline,                                                                    &
               neutral,                                                                            &
               passive_scalar,                                                                     &
               rans_mode,                                                                          &
               rans_tke_e,                                                                         &
               rayleigh_damping_factor,                                                            &
               rayleigh_damping_height,                                                            &
               salsa,                                                                              &
               spinup_time,                                                                        &
               time_since_reference_point,                                                         &
               volume_flow

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE grid_variables

    USE indices,                                                                                   &
        ONLY:  nbgp, nx, nxl, nxlg, nxlu, nxr, nxrg, ny, nys, nysv, nysg, nyn, nyng, nzb, nz, nzt, &
               topo_top_ind, topo_flags

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  char_fill,                                                                          &
               char_lod,                                                                           &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               get_variable_pr,                                                                    &
               input_pids_dynamic,                                                                 &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               input_file_dynamic,                                                                 &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               pids_id

    USE pegrid

    USE salsa_mod,                                                                                 &
        ONLY:  salsa_nesting_offl_bc,                                                              &
               salsa_nesting_offl_init,                                                            &
               salsa_nesting_offl_input

    IMPLICIT NONE

!
!-- Define data type for nesting in larger-scale models like COSMO.
!-- Data type comprises u, v, w, pt, and q at lateral and top boundaries.
    TYPE nest_offl_type

       CHARACTER(LEN=16) ::  char_l = 'ls_forcing_left_'   !< leading substring for variables at left boundary
       CHARACTER(LEN=17) ::  char_n = 'ls_forcing_north_'  !< leading substring for variables at north boundary
       CHARACTER(LEN=17) ::  char_r = 'ls_forcing_right_'  !< leading substring for variables at right boundary
       CHARACTER(LEN=17) ::  char_s = 'ls_forcing_south_'  !< leading substring for variables at south boundary
       CHARACTER(LEN=15) ::  char_t = 'ls_forcing_top_'    !< leading substring for variables at top boundary

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names         !< list of variable in dynamic input file
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_l  !< names of mesoscale nested chemistry variables at left
                                                                           !< boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_n  !< names of mesoscale nested chemistry variables at north
                                                                           !< boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_r  !< names of mesoscale nested chemistry variables at right
                                                                           !< boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_s  !< names of mesoscale nested chemistry variables at south
                                                                           !< boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_t  !< names of mesoscale nested chemistry variables at top
                                                                           !< boundary

       INTEGER(iwp) ::  lod_east_pt  = 2  !< level-of-detail of input data of potential temperature at the eastern boundary
       INTEGER(iwp) ::  lod_east_qc  = 2  !< level-of-detail of input data of cloud-water mixture fraction at the eastern boundary
       INTEGER(iwp) ::  lod_east_qv  = 2  !< level-of-detail of input data of specific humidity at the eastern boundary
       INTEGER(iwp) ::  lod_east_u   = 2  !< level-of-detail of input data of the u-component at the eastern boundary
       INTEGER(iwp) ::  lod_east_v   = 2  !< level-of-detail of input data of the v-component at the eastern boundary
       INTEGER(iwp) ::  lod_east_w   = 2  !< level-of-detail of input data of the w-component at the eastern boundary
       INTEGER(iwp) ::  lod_north_pt = 2  !< level-of-detail of input data of potential temperature at the northern boundary
       INTEGER(iwp) ::  lod_north_qc = 2  !< level-of-detail of input data of cloud-water mixture fraction at the northern boundary
       INTEGER(iwp) ::  lod_north_qv = 2  !< level-of-detail of input data of specific humidity at the northern boundary
       INTEGER(iwp) ::  lod_north_u  = 2  !< level-of-detail of input data of the u-component at the northern boundary
       INTEGER(iwp) ::  lod_north_v  = 2  !< level-of-detail of input data of the v-component at the northern boundary
       INTEGER(iwp) ::  lod_north_w  = 2  !< level-of-detail of input data of the w-component at the northern boundary
       INTEGER(iwp) ::  lod_south_pt = 2  !< level-of-detail of input data of potential temperature at the southern boundary
       INTEGER(iwp) ::  lod_south_qc = 2  !< level-of-detail of input data of cloud-water mixture fraction at the southern boundary
       INTEGER(iwp) ::  lod_south_qv = 2  !< level-of-detail of input data of specific humidity at the southern boundary
       INTEGER(iwp) ::  lod_south_u  = 2  !< level-of-detail of input data of the u-component at the southern boundary
       INTEGER(iwp) ::  lod_south_v  = 2  !< level-of-detail of input data of the v-component at the southern boundary
       INTEGER(iwp) ::  lod_south_w  = 2  !< level-of-detail of input data of the w-component at the southern boundary
       INTEGER(iwp) ::  lod_top_pt   = 2  !< level-of-detail of input data of potential temperature at the top boundary
       INTEGER(iwp) ::  lod_top_qc   = 2  !< level-of-detail of input data of cloud-water mixture fraction at the top boundary
       INTEGER(iwp) ::  lod_top_qv   = 2  !< level-of-detail of input data of specific humidity at the top boundary
       INTEGER(iwp) ::  lod_top_u    = 2  !< level-of-detail of input data of the u-component at the top boundary
       INTEGER(iwp) ::  lod_top_v    = 2  !< level-of-detail of input data of the v-component at the top boundary
       INTEGER(iwp) ::  lod_top_w    = 2  !< level-of-detail of input data of the w-component at the top boundary
       INTEGER(iwp) ::  lod_west_pt  = 2  !< level-of-detail of input data of potential temperature at the western boundary
       INTEGER(iwp) ::  lod_west_qc  = 2  !< level-of-detail of input data of cloud-water mixture fraction at the western boundary
       INTEGER(iwp) ::  lod_west_qv  = 2  !< level-of-detail of input data of specific humidity at the western boundary
       INTEGER(iwp) ::  lod_west_u   = 2  !< level-of-detail of input data of the u-component at the western boundary
       INTEGER(iwp) ::  lod_west_v   = 2  !< level-of-detail of input data of the v-component at the western boundary
       INTEGER(iwp) ::  lod_west_w   = 2  !< level-of-detail of input data of the w-component at the western boundary
       INTEGER(iwp) ::  nt                !< number of time levels in dynamic input file
       INTEGER(iwp) ::  nzu               !< number of vertical levels on scalar grid in dynamic input file
       INTEGER(iwp) ::  nzw               !< number of vertical levels on w grid in dynamic input file
       INTEGER(iwp) ::  tind = 0          !< time index for reference time in mesoscale-offline nesting
       INTEGER(iwp) ::  tind_p = 0        !< time index for following time in mesoscale-offline nesting

       LOGICAL      ::  init = .FALSE.    !< flag indicating that offline nesting is already initialized

       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_l  !< flags inidicating whether left boundary data for chemistry is in
                                                                !< dynamic input file
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_n  !< flags inidicating whether north boundary data for chemistry is in
                                                                !< dynamic input file
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_r  !< flags inidicating whether right boundary data for chemistry is in
                                                                !< dynamic input file
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_s  !< flags inidicating whether south boundary data for chemistry is in
                                                                !< dynamic input file
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_t  !< flags inidicating whether top boundary data for chemistry is in
                                                                !< dynamic input file

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surface_pressure  !< time dependent surface pressure
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time              !< time levels in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu_atmos          !< vertical levels at scalar grid in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw_atmos          !< vertical levels at w grid in dynamic input file

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_l     !< potentital temperautre at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_n     !< potentital temperautre at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_r     !< potentital temperautre at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_s     !< potentital temperautre at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_top   !< potentital temperautre at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_l      !< mixing ratio at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_n      !< mixing ratio at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_r      !< mixing ratio at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_s      !< mixing ratio at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_top    !< mixing ratio at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_l      !< u-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_n      !< u-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_r      !< u-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_s      !< u-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_top    !< u-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_l      !< v-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_n      !< v-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_r      !< v-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_s      !< v-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_top    !< v-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_l      !< w-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_n      !< w-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_r      !< w-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_s      !< w-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_top    !< w-component at top boundary

       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_l   !< chemical species at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_n   !< chemical species at north boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_r   !< chemical species at right boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_s   !< chemical species at south boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_top !< chemical species at left boundary

    END TYPE nest_offl_type

    INTEGER(iwp) ::  i_bound     !< boundary grid point in x-direction for scalars, v, and w
    INTEGER(iwp) ::  i_bound_u   !< boundary grid point in x-direction for u
    INTEGER(iwp) ::  i_end       !< end index for array allocation along x-direction at norther/southern boundary
    INTEGER(iwp) ::  i_start     !< start index for array allocation along x-direction at norther/southern boundary (scalars, v, w)
    INTEGER(iwp) ::  i_start_u   !< start index for array allocation along x-direction at norther/southern boundary (u)
    INTEGER(iwp) ::  j_bound     !< boundary grid point in y-direction for scalars, u, and w
    INTEGER(iwp) ::  j_bound_v   !< boundary grid point in y-direction for v
    INTEGER(iwp) ::  j_end       !< end index for array allocation along y-direction at eastern/western boundary
    INTEGER(iwp) ::  j_start     !< start index for array allocation along y-direction at eastern/western boundary (scalars, u, w)
    INTEGER(iwp) ::  j_start_v   !< start index for array allocation along y-direction at eastern/western boundary (v)
    INTEGER(iwp) ::  lod         !< level-of-detail of lateral input data

    REAL(wp) ::  fac_dt              !< interpolation factor
    REAL(wp) ::  zi_ribulk = 0.0_wp  !< boundary-layer depth according to bulk Richardson criterion, i.e. the height where Ri_bulk
                                     !< exceeds the critical bulk Richardson number of 0.2

    TYPE(nest_offl_type) ::  nest_offl  !< data structure for data input at lateral and top boundaries (provided by Inifor)

    SAVE
    PRIVATE
!
!-- Public subroutines
    PUBLIC nesting_offl_bc,                                                                        &
           nesting_offl_calc_zi,                                                                   &
           nesting_offl_check_parameters,                                                          &
           nesting_offl_header,                                                                    &
           nesting_offl_init,                                                                      &
           nesting_offl_input,                                                                     &
           nesting_offl_interpolation_factor,                                                      &
           nesting_offl_mass_conservation,                                                         &
           nesting_offl_parin, interpolate_in_time
!
!-- Public variables
    PUBLIC zi_ribulk

    INTERFACE nesting_offl_bc
       MODULE PROCEDURE nesting_offl_bc
    END INTERFACE nesting_offl_bc

    INTERFACE nesting_offl_calc_zi
       MODULE PROCEDURE nesting_offl_calc_zi
    END INTERFACE nesting_offl_calc_zi

    INTERFACE nesting_offl_check_parameters
       MODULE PROCEDURE nesting_offl_check_parameters
    END INTERFACE nesting_offl_check_parameters

    INTERFACE nesting_offl_header
       MODULE PROCEDURE nesting_offl_header
    END INTERFACE nesting_offl_header

    INTERFACE nesting_offl_init
       MODULE PROCEDURE nesting_offl_init
    END INTERFACE nesting_offl_init

    INTERFACE nesting_offl_input
       MODULE PROCEDURE nesting_offl_input
    END INTERFACE nesting_offl_input

    INTERFACE nesting_offl_interpolation_factor
       MODULE PROCEDURE nesting_offl_interpolation_factor
    END INTERFACE nesting_offl_interpolation_factor

    INTERFACE nesting_offl_mass_conservation
       MODULE PROCEDURE nesting_offl_mass_conservation
    END INTERFACE nesting_offl_mass_conservation

    INTERFACE nesting_offl_parin
       MODULE PROCEDURE nesting_offl_parin
    END INTERFACE nesting_offl_parin

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads data at lateral and top boundaries derived from larger-scale model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_input

    INTEGER(iwp) ::  n   !< running index for chemistry variables

!
!-- Initialize INIFOR forcing in first call.
    IF ( .NOT. nest_offl%init )  THEN
#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), pids_id )
!
!--    At first, inquire all variable names.
       CALL inquire_num_variables( pids_id, num_var_pids )
!
!--    Allocate memory to store variable names.
       ALLOCATE( nest_offl%var_names(1:num_var_pids) )
       CALL inquire_variable_names( pids_id, nest_offl%var_names )
!
!--    Read time dimension, allocate memory and finally read time array
       CALL get_dimension_length( pids_id, nest_offl%nt, 'time' )

       IF ( check_existence( nest_offl%var_names, 'time' ) )  THEN
          ALLOCATE( nest_offl%time(0:nest_offl%nt-1) )
          CALL get_variable( pids_id, 'time', nest_offl%time )
       ENDIF
!
!--    Read vertical dimension of scalar und w grid
       CALL get_dimension_length( pids_id, nest_offl%nzu, 'z' )
       CALL get_dimension_length( pids_id, nest_offl%nzw, 'zw' )

       IF ( check_existence( nest_offl%var_names, 'z' ) )  THEN
          ALLOCATE( nest_offl%zu_atmos(1:nest_offl%nzu) )
          CALL get_variable( pids_id, 'z', nest_offl%zu_atmos )
       ENDIF
       IF ( check_existence( nest_offl%var_names, 'zw' ) )  THEN
          ALLOCATE( nest_offl%zw_atmos(1:nest_offl%nzw) )
          CALL get_variable( pids_id, 'zw', nest_offl%zw_atmos )
       ENDIF
!
!--    Read surface pressure
       IF ( check_existence( nest_offl%var_names, 'surface_forcing_surface_pressure' ) )  THEN
          ALLOCATE( nest_offl%surface_pressure(0:nest_offl%nt-1) )
          CALL get_variable( pids_id, 'surface_forcing_surface_pressure',                          &
                             nest_offl%surface_pressure )
       ENDIF
!
!--    Close input file
       CALL close_input_file( pids_id )
#endif
    ENDIF
!
!-- Check if dynamic driver data input is required.
    IF ( nest_offl%time(nest_offl%tind_p) <= MAX( time_since_reference_point, 0.0_wp)  .OR.        &
         .NOT.  nest_offl%init )  THEN
       CONTINUE
!
!-- Return otherwise
    ELSE
       RETURN
    ENDIF
!
!-- Start of CPU measurement
    CALL cpu_log( log_point_s(86), 'NetCDF input forcing', 'start' )

!
!-- Obtain time index for current point in time. Note, the time coordinate in the input file is
!-- always relative to the initial time in UTC, i.e. the time coordinate always starts at 0.0 even
!-- if the initial UTC is e.g. 7200.0. Further, since time_since_reference_point is negativ here
!-- when spinup is applied, use MAX function to obtain correct time index.
    nest_offl%tind = MINLOC( ABS( nest_offl%time - MAX( time_since_reference_point, 0.0_wp ) ),    &
                             DIM = 1 ) - 1
!
!--    Note, in case of restart runs, the time index for the boundary data may indicate a time in
!--    the future. This needs to be checked and corrected.
       IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.                             &
            nest_offl%time(nest_offl%tind) > time_since_reference_point )  THEN
          nest_offl%tind = nest_offl%tind - 1
       ENDIF
    nest_offl%tind_p = nest_offl%tind + 1
!
!-- Open file in read-only mode
#if defined ( __netcdf )
    CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), pids_id )
!
!-- Read data at lateral and top boundaries. Please note, at left and right domain boundary,
!-- yz-layers are read for u, v, w, pt and q.
!-- For the v-component, the data starts at nysv, while for the other quantities the data starts at
!-- nys. This is equivalent at the north and south domain boundary for the u-component (nxlu).
!-- Note, lateral data is also accessed by parallel IO, which is the reason why different arguments
!-- are passed depending on the boundary control flags. Cores that do not belong to the respective
!-- boundary only do a dummy read with count = 0, just in order to participate the collective
!-- operation. This is because collective parallel access shows better performance than just a
!-- conditional access.
!-- Read data for LOD 2, i.e. time-dependent xz-, yz-, and xy-slices.
    IF ( lod == 2 )  THEN
       CALL get_variable( pids_id, 'ls_forcing_left_u',                                            &
                          nest_offl%u_l,                                                           & ! array to be read
                          MERGE( nys+1, 1, bc_dirichlet_l),                                        & ! start index y direction
                          MERGE( nzb+1, 1, bc_dirichlet_l),                                        & ! start index z direction
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                             & ! start index time dimension
                          MERGE( nyn-nys+1, 0, bc_dirichlet_l),                                    & ! number of elements along y
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                                & ! number of elements alogn z
                          MERGE( 2, 0, bc_dirichlet_l),                                            & ! number of time steps (2 or 0)
                          .TRUE. )                                                                   ! parallel IO when compiled accordingly

       CALL get_variable( pids_id, 'ls_forcing_left_v',                                            &
                          nest_offl%v_l,                                                           &
                          MERGE( nysv, 1, bc_dirichlet_l),                                         &
                          MERGE( nzb+1, 1, bc_dirichlet_l),                                        &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                             &
                          MERGE( nyn-nysv+1, 0, bc_dirichlet_l),                                   &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                                &
                          MERGE( 2, 0, bc_dirichlet_l),                                            &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_left_w',                                            &
                          nest_offl%w_l,                                                           &
                          MERGE( nys+1, 1, bc_dirichlet_l),                                        &
                          MERGE( nzb+1, 1, bc_dirichlet_l),                                        &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                             &
                          MERGE( nyn-nys+1, 0, bc_dirichlet_l),                                    &
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_l),                                &
                          MERGE( 2, 0, bc_dirichlet_l),                                            &
                          .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_left_pt',                                        &
                             nest_offl%pt_l,                                                       &
                             MERGE( nys+1, 1, bc_dirichlet_l),                                     &
                             MERGE( nzb+1, 1, bc_dirichlet_l),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                          &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_l),                                 &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                             &
                             MERGE( 2, 0, bc_dirichlet_l),                                         &
                             .TRUE. )
       ENDIF

       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_left_qv',                                        &
                             nest_offl%q_l,                                                        &
                             MERGE( nys+1, 1, bc_dirichlet_l),                                     &
                             MERGE( nzb+1, 1, bc_dirichlet_l),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                          &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_l),                                 &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                             &
                             MERGE( 2, 0, bc_dirichlet_l),                                         &
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( nest_offl%var_names_chem_l, 1 )
             IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_l(n) ) )  THEN
                CALL get_variable( pids_id,                                                        &
                                   TRIM( nest_offl%var_names_chem_l(n) ),                          &
                                   nest_offl%chem_l(:,:,:,n),                                      &
                                   MERGE( nys+1, 1, bc_dirichlet_l),                               &
                                   MERGE( nzb+1, 1, bc_dirichlet_l),                               &
                                   MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                    &
                                   MERGE( nyn-nys+1, 0, bc_dirichlet_l),                           &
                                   MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                       &
                                   MERGE( 2, 0, bc_dirichlet_l),                                   &
                                   .TRUE. )
                nest_offl%chem_from_file_l(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Read data for eastern boundary
       CALL get_variable( pids_id, 'ls_forcing_right_u',                                           &
                          nest_offl%u_r,                                                           &
                          MERGE( nys+1, 1, bc_dirichlet_r),                                        &
                          MERGE( nzb+1, 1, bc_dirichlet_r),                                        &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                             &
                          MERGE( nyn-nys+1, 0, bc_dirichlet_r),                                    &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                                &
                          MERGE( 2, 0, bc_dirichlet_r),                                            &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_right_v',                                           &
                          nest_offl%v_r,                                                           &
                          MERGE( nysv, 1, bc_dirichlet_r),                                         &
                          MERGE( nzb+1, 1, bc_dirichlet_r),                                        &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                             &
                          MERGE( nyn-nysv+1, 0, bc_dirichlet_r),                                   &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                                &
                          MERGE( 2, 0, bc_dirichlet_r),                                            &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_right_w',                                           &
                          nest_offl%w_r,                                                           &
                          MERGE( nys+1, 1, bc_dirichlet_r),                                        &
                          MERGE( nzb+1, 1, bc_dirichlet_r),                                        &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                             &
                          MERGE( nyn-nys+1, 0, bc_dirichlet_r),                                    &
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_r),                                &
                          MERGE( 2, 0, bc_dirichlet_r),                                            &
                          .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_right_pt',                                       &
                             nest_offl%pt_r,                                                       &
                             MERGE( nys+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                          &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_r),                                 &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                             &
                             MERGE( 2, 0, bc_dirichlet_r),                                         &
                             .TRUE. )
       ENDIF

       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_right_qv',                                       &
                             nest_offl%q_r,                                                        &
                             MERGE( nys+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                          &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_r),                                 &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                             &
                             MERGE( 2, 0, bc_dirichlet_r),                                         &
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( nest_offl%var_names_chem_r, 1 )
             IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_r(n) ) )  THEN
                CALL get_variable( pids_id,                                                        &
                                   TRIM( nest_offl%var_names_chem_r(n) ),                          &
                                   nest_offl%chem_r(:,:,:,n),                                      &
                                   MERGE( nys+1, 1, bc_dirichlet_r),                               &
                                   MERGE( nzb+1, 1, bc_dirichlet_r),                               &
                                   MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                    &
                                   MERGE( nyn-nys+1, 0, bc_dirichlet_r),                           &
                                   MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                       &
                                   MERGE( 2, 0, bc_dirichlet_r),                                   &
                                   .TRUE. )
                nest_offl%chem_from_file_r(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Read data for northern boundary
       CALL get_variable( pids_id, 'ls_forcing_north_u',                                           &
                          nest_offl%u_n,                                                           &
                          MERGE( nxlu, 1, bc_dirichlet_n ),                                        &
                          MERGE( nzb+1, 1, bc_dirichlet_n ),                                       &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),                            &
                          MERGE( nxr-nxlu+1, 0, bc_dirichlet_n ),                                  &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),                               &
                          MERGE( 2, 0, bc_dirichlet_n ),                                           &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_north_v',                                           &
                          nest_offl%v_n,                                                           &
                          MERGE( nxl+1, 1, bc_dirichlet_n ),                                       &
                          MERGE( nzb+1, 1, bc_dirichlet_n ),                                       &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),                            &
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),                                   &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),                               &
                          MERGE( 2, 0, bc_dirichlet_n ),                                           &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_north_w',                                           &
                          nest_offl%w_n,                                                           &
                          MERGE( nxl+1, 1, bc_dirichlet_n ),                                       &
                          MERGE( nzb+1, 1, bc_dirichlet_n ),                                       &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),                            &
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),                                   &
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_n ),                               &
                          MERGE( 2, 0, bc_dirichlet_n ),                                           &
                          .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_north_pt',                                       &
                             nest_offl%pt_n,                                                       &
                             MERGE( nxl+1, 1, bc_dirichlet_n ),                                    &
                             MERGE( nzb+1, 1, bc_dirichlet_n ),                                    &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),                         &
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),                                &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),                            &
                             MERGE( 2, 0, bc_dirichlet_n ),                                        &
                             .TRUE. )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_north_qv',                                       &
                             nest_offl%q_n,                                                        &
                             MERGE( nxl+1, 1, bc_dirichlet_n ),                                    &
                             MERGE( nzb+1, 1, bc_dirichlet_n ),                                    &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),                         &
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),                                &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),                            &
                             MERGE( 2, 0, bc_dirichlet_n ),                                        &
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( nest_offl%var_names_chem_n, 1 )
             IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_n(n) ) )  THEN
                CALL get_variable( pids_id,                                                        &
                                   TRIM( nest_offl%var_names_chem_n(n) ),                          &
                                   nest_offl%chem_n(:,:,:,n),                                      &
                                   MERGE( nxl+1, 1, bc_dirichlet_n ),                              &
                                   MERGE( nzb+1, 1, bc_dirichlet_n ),                              &
                                   MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),                   &
                                   MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),                          &
                                   MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),                      &
                                   MERGE( 2, 0, bc_dirichlet_n ),                                  &
                                   .TRUE. )
                nest_offl%chem_from_file_n(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Read data for southern boundary
       CALL get_variable( pids_id, 'ls_forcing_south_u',                                           &
                          nest_offl%u_s,                                                           &
                          MERGE( nxlu, 1, bc_dirichlet_s ),                                        &
                          MERGE( nzb+1, 1, bc_dirichlet_s ),                                       &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),                            &
                          MERGE( nxr-nxlu+1, 0, bc_dirichlet_s ),                                  &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),                               &
                          MERGE( 2, 0, bc_dirichlet_s ),                                           &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_south_v',                                           &
                          nest_offl%v_s,                                                           &
                          MERGE( nxl+1, 1, bc_dirichlet_s ),                                       &
                          MERGE( nzb+1, 1, bc_dirichlet_s ),                                       &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),                            &
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),                                   &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),                               &
                          MERGE( 2, 0, bc_dirichlet_s ),                                           &
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_south_w',                                           &
                          nest_offl%w_s,                                                           &
                          MERGE( nxl+1, 1, bc_dirichlet_s ),                                       &
                          MERGE( nzb+1, 1, bc_dirichlet_s ),                                       &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),                            &
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),                                   &
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_s ),                               &
                          MERGE( 2, 0, bc_dirichlet_s ),                                           &
                          .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_south_pt',                                       &
                             nest_offl%pt_s,                                                       &
                             MERGE( nxl+1, 1, bc_dirichlet_s ),                                    &
                             MERGE( nzb+1, 1, bc_dirichlet_s ),                                    &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),                         &
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),                                &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),                            &
                             MERGE( 2, 0, bc_dirichlet_s ),                                        &
                             .TRUE. )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_south_qv',                                       &
                             nest_offl%q_s,                                                        &
                             MERGE( nxl+1, 1, bc_dirichlet_s ),                                    &
                             MERGE( nzb+1, 1, bc_dirichlet_s ),                                    &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),                         &
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),                                &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),                            &
                             MERGE( 2, 0, bc_dirichlet_s ),                                        &
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( nest_offl%var_names_chem_s, 1 )
             IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_s(n) ) )  THEN
                CALL get_variable( pids_id,                                                        &
                                   TRIM( nest_offl%var_names_chem_s(n) ),                          &
                                   nest_offl%chem_s(:,:,:,n),                                      &
                                   MERGE( nxl+1, 1, bc_dirichlet_s ),                              &
                                   MERGE( nzb+1, 1, bc_dirichlet_s ),                              &
                                   MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),                   &
                                   MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),                          &
                                   MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),                      &
                                   MERGE( 2, 0, bc_dirichlet_s ),                                  &
                                   .TRUE. )
                nest_offl%chem_from_file_s(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Top boundary
       CALL get_variable( pids_id, 'ls_forcing_top_u',                                             &
                          nest_offl%u_top(0:1,nys:nyn,nxlu:nxr),                                   &
                          nxlu, nys+1, nest_offl%tind+1,                                           &
                          nxr-nxlu+1, nyn-nys+1, 2, .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_top_v',                                             &
                          nest_offl%v_top(0:1,nysv:nyn,nxl:nxr),                                   &
                          nxl+1, nysv, nest_offl%tind+1,                                           &
                          nxr-nxl+1, nyn-nysv+1, 2, .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_top_w',                                             &
                          nest_offl%w_top(0:1,nys:nyn,nxl:nxr),                                    &
                          nxl+1, nys+1, nest_offl%tind+1,                                          &
                          nxr-nxl+1, nyn-nys+1, 2, .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_top_pt',                                         &
                             nest_offl%pt_top(0:1,nys:nyn,nxl:nxr),                                &
                             nxl+1, nys+1, nest_offl%tind+1,                                       &
                             nxr-nxl+1, nyn-nys+1, 2, .TRUE. )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_top_qv',                                         &
                             nest_offl%q_top(0:1,nys:nyn,nxl:nxr),                                 &
                             nxl+1, nys+1, nest_offl%tind+1,                                       &
                             nxr-nxl+1, nyn-nys+1, 2, .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( nest_offl%var_names_chem_t, 1 )
             IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_t(n) ) )  THEN
                CALL get_variable( pids_id,                                                        &
                                   TRIM( nest_offl%var_names_chem_t(n) ),                          &
                                   nest_offl%chem_top(0:1,nys:nyn,nxl:nxr,n),                      &
                                   nxl+1, nys+1, nest_offl%tind+1,                                 &
                                   nxr-nxl+1, nyn-nys+1, 2, .TRUE. )
                nest_offl%chem_from_file_t(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!-- Read data for LOD 1, i.e. time-dependent profiles. In constrast to LOD 2 where the amount of IO
!-- is larger, only the respective boundary processes read the data.
    ELSE
       IF ( bc_dirichlet_l )  THEN
          CALL get_variable( pids_id, 'ls_forcing_left_u',                                         &
                             nest_offl%u_l(0:1,:,1:1),                                             & ! array to be read
                             MERGE( nzb+1, 1, bc_dirichlet_l),                                     & ! start index z direction
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                          & ! start index time dimension
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                             & ! number of elements along z
                             MERGE( 2, 0, bc_dirichlet_l) )                                          ! number of time steps (2 or 0)
          CALL get_variable( pids_id, 'ls_forcing_left_v',                                         &
                             nest_offl%v_l(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_l),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                             &
                             MERGE( 2, 0, bc_dirichlet_l) )
          CALL get_variable( pids_id, 'ls_forcing_left_w',                                         &
                             nest_offl%w_l(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_l),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                          &
                             MERGE( nest_offl%nzw, 0, bc_dirichlet_l),                             &
                             MERGE( 2, 0, bc_dirichlet_l) )
          IF ( .NOT. neutral )  THEN
             CALL get_variable( pids_id, 'ls_forcing_left_pt',                                     &
                                nest_offl%pt_l(0:1,:,1:1),                                         &
                                MERGE( nzb+1, 1, bc_dirichlet_l),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                          &
                                MERGE( 2, 0, bc_dirichlet_l) )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( pids_id, 'ls_forcing_left_qv',                                     &
                                nest_offl%q_l(0:1,:,1:1),                                          &
                                MERGE( nzb+1, 1, bc_dirichlet_l),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                          &
                                MERGE( 2, 0, bc_dirichlet_l) )
          ENDIF
          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( nest_offl%var_names_chem_t, 1 )
                IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_t(n) ) )  THEN
                   CALL get_variable( pids_id, TRIM( nest_offl%var_names_chem_t(n) ),              &
                                      nest_offl%chem_l(0:1,:,1:1,n),                               &
                                      MERGE( nzb+1, 1, bc_dirichlet_l),                            &
                                      MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),                 &
                                      MERGE( nest_offl%nzu, 0, bc_dirichlet_l),                    &
                                      MERGE( 2, 0, bc_dirichlet_l) )
                   nest_offl%chem_from_file_l(n) = .TRUE.
                ENDIF
             ENDDO
          ENDIF
       ENDIF
       IF ( bc_dirichlet_r )  THEN
          CALL get_variable( pids_id, 'ls_forcing_right_u',                                        &
                             nest_offl%u_r(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                             &
                             MERGE( 2, 0, bc_dirichlet_r) )
          CALL get_variable( pids_id, 'ls_forcing_right_v',                                        &
                             nest_offl%v_r(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                             &
                             MERGE( 2, 0, bc_dirichlet_r) )
          CALL get_variable( pids_id, 'ls_forcing_right_w',                                        &
                             nest_offl%w_r(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                          &
                             MERGE( nest_offl%nzw, 0, bc_dirichlet_r),                             &
                             MERGE( 2, 0, bc_dirichlet_r) )
          IF ( .NOT. neutral )  THEN
             CALL get_variable( pids_id, 'ls_forcing_right_pt',                                    &
                                nest_offl%pt_r(0:1,:,1:1),                                         &
                                MERGE( nzb+1, 1, bc_dirichlet_r),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                          &
                                MERGE( 2, 0, bc_dirichlet_r) )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( pids_id, 'ls_forcing_right_qv',                                    &
                                nest_offl%q_r(0:1,:,1:1),                                          &
                                MERGE( nzb+1, 1, bc_dirichlet_r),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                          &
                                MERGE( 2, 0, bc_dirichlet_r) )
          ENDIF
          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( nest_offl%var_names_chem_t, 1 )
                IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_t(n) ) )  THEN
                   CALL get_variable( pids_id, TRIM( nest_offl%var_names_chem_t(n) ),              &
                                      nest_offl%chem_r(0:1,:,1:1,n),                               &
                                      MERGE( nzb+1, 1, bc_dirichlet_r),                            &
                                      MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),                 &
                                      MERGE( nest_offl%nzu, 0, bc_dirichlet_r),                    &
                                      MERGE( 2, 0, bc_dirichlet_r) )
                   nest_offl%chem_from_file_r(n) = .TRUE.
                ENDIF
             ENDDO
          ENDIF
       ENDIF
       IF ( bc_dirichlet_n )  THEN
          CALL get_variable( pids_id, 'ls_forcing_north_u',                                        &
                             nest_offl%u_n(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_n),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_n),                             &
                             MERGE( 2, 0, bc_dirichlet_n) )
          CALL get_variable( pids_id, 'ls_forcing_north_v',                                        &
                             nest_offl%v_n(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_n),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_n),                             &
                             MERGE( 2, 0, bc_dirichlet_n) )
          CALL get_variable( pids_id, 'ls_forcing_north_w',                                        &
                             nest_offl%w_n(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_n),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n),                          &
                             MERGE( nest_offl%nzw, 0, bc_dirichlet_n),                             &
                             MERGE( 2, 0, bc_dirichlet_n) )
          IF ( .NOT. neutral )  THEN
             CALL get_variable( pids_id, 'ls_forcing_north_pt',                                    &
                                nest_offl%pt_n(0:1,:,1:1),                                         &
                                MERGE( nzb+1, 1, bc_dirichlet_n),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_n),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_n),                          &
                                MERGE( 2, 0, bc_dirichlet_n) )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( pids_id, 'ls_forcing_north_qv',                                    &
                                nest_offl%q_n(0:1,:,1:1),                                          &
                                MERGE( nzb+1, 1, bc_dirichlet_n),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_n),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_n),                          &
                                MERGE( 2, 0, bc_dirichlet_n) )
          ENDIF
          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( nest_offl%var_names_chem_t, 1 )
                IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_t(n) ) )  THEN
                   CALL get_variable( pids_id, TRIM( nest_offl%var_names_chem_t(n) ),              &
                                      nest_offl%chem_n(0:1,:,1:1,n),                               &
                                      MERGE( nzb+1, 1, bc_dirichlet_n),                            &
                                      MERGE( nest_offl%tind+1, 1, bc_dirichlet_n),                 &
                                      MERGE( nest_offl%nzu, 0, bc_dirichlet_n),                    &
                                      MERGE( 2, 0, bc_dirichlet_n) )
                   nest_offl%chem_from_file_n(n) = .TRUE.
                ENDIF
             ENDDO
          ENDIF
       ENDIF
       IF ( bc_dirichlet_s )  THEN
          CALL get_variable( pids_id, 'ls_forcing_south_u',                                        &
                             nest_offl%u_s(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_s),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_s),                             &
                             MERGE( 2, 0, bc_dirichlet_s) )
          CALL get_variable( pids_id, 'ls_forcing_south_v',                                        &
                             nest_offl%v_s(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_s),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s),                          &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_s),                             &
                             MERGE( 2, 0, bc_dirichlet_s) )
          CALL get_variable( pids_id, 'ls_forcing_south_w',                                        &
                             nest_offl%w_s(0:1,:,1:1),                                             &
                             MERGE( nzb+1, 1, bc_dirichlet_s),                                     &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s),                          &
                             MERGE( nest_offl%nzw, 0, bc_dirichlet_s),                             &
                             MERGE( 2, 0, bc_dirichlet_s) )
          IF ( .NOT. neutral )  THEN
             CALL get_variable( pids_id, 'ls_forcing_south_pt',                                    &
                                nest_offl%pt_s(0:1,:,1:1),                                         &
                                MERGE( nzb+1, 1, bc_dirichlet_s),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_s),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_s),                          &
                                MERGE( 2, 0, bc_dirichlet_s) )
          ENDIF
          IF ( humidity )  THEN
             CALL get_variable( pids_id, 'ls_forcing_south_qv',                                    &
                                nest_offl%q_s(0:1,:,1:1),                                          &
                                MERGE( nzb+1, 1, bc_dirichlet_s),                                  &
                                MERGE( nest_offl%tind+1, 1, bc_dirichlet_s),                       &
                                MERGE( nest_offl%nzu, 0, bc_dirichlet_s),                          &
                                MERGE( 2, 0, bc_dirichlet_s) )
          ENDIF
          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( nest_offl%var_names_chem_t, 1 )
                IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_t(n) ) )  THEN
                   CALL get_variable( pids_id, TRIM( nest_offl%var_names_chem_t(n) ),              &
                                      nest_offl%chem_s(0:1,:,1:1,n),                               &
                                      MERGE( nzb+1, 1, bc_dirichlet_s),                            &
                                      MERGE( nest_offl%tind+1, 1, bc_dirichlet_s),                 &
                                      MERGE( nest_offl%nzu, 0, bc_dirichlet_s),                    &
                                      MERGE( 2, 0, bc_dirichlet_s) )
                   nest_offl%chem_from_file_s(n) = .TRUE.
                ENDIF
             ENDDO
          ENDIF
       ENDIF
!
!--    Read top boundary data, which is actually only a scalar value in the LOD 1 case.
       CALL get_variable( pids_id, 'ls_forcing_top_u',                                             &
                          nest_offl%u_top(0:1,1,1),                                                & ! array to be read
                          nest_offl%tind+1,                                                        & ! start index in time
                          2 )                                                                        ! number of elements to be read
       CALL get_variable( pids_id, 'ls_forcing_top_v',                                             &
                          nest_offl%v_top(0:1,1,1),                                                &
                          nest_offl%tind+1,                                                        &
                          2 )
       CALL get_variable( pids_id, 'ls_forcing_top_w',                                             &
                          nest_offl%w_top(0:1,1,1),                                                &
                          nest_offl%tind+1,                                                        &
                          2 )
       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_top_pt',                                         &
                             nest_offl%pt_top(0:1,1,1),                                            &
                             nest_offl%tind+1,                                                     &
                             2 )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_top_qv',                                         &
                             nest_offl%q_top(0:1,1,1),                                             &
                             nest_offl%tind+1,                                                     &
                             2 )
       ENDIF
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( nest_offl%var_names_chem_t, 1 )
             IF ( check_existence( nest_offl%var_names, nest_offl%var_names_chem_t(n) ) )  THEN
                CALL get_variable( pids_id, TRIM( nest_offl%var_names_chem_t(n) ),                 &
                                   nest_offl%chem_top(0:1,1,1,n),                                  &
                                   nest_offl%tind+1,                                               &
                                   2 )
                nest_offl%chem_from_file_t(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
    ENDIF


!
!-- Close input file
    CALL close_input_file( pids_id )
#endif
!
!-- Set control flag to indicate that boundary data has been initially input.
    nest_offl%init = .TRUE.
!
!-- Call offline nesting for salsa
    IF ( salsa )  CALL salsa_nesting_offl_input
!
!-- End of CPU measurement
    CALL cpu_log( log_point_s(86), 'NetCDF input forcing', 'stop' )

 END SUBROUTINE nesting_offl_input


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In this subroutine a constant mass within the model domain is guaranteed.
!> Larger-scale models may be based on a compressible equation system, which is not consistent with
!> PALMs incompressible equation system. In order to avoid a decrease or increase of mass during the
!> simulation, non-divergent flow through the lateral and top boundaries is compensated by the
!> vertical wind component at the top boundary.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_mass_conservation

    INTEGER(iwp) ::  i  !< grid index in x-direction
    INTEGER(iwp) ::  j  !< grid index in y-direction
    INTEGER(iwp) ::  k  !< grid index in z-direction

    REAL(wp) ::  d_area_t   !< inverse of the total area of the horizontal model domain
    REAL(wp) ::  w_correct  !< vertical velocity increment required to compensate non-divergent flow through the boundaries

    REAL(wp), DIMENSION(1:3) ::  volume_flow_l   !< local volume flow


    IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_mass_conservation', 'start' )

    CALL  cpu_log( log_point(58), 'offline nesting', 'start' )

    volume_flow   = 0.0_wp
    volume_flow_l = 0.0_wp

    d_area_t = 1.0_wp / ( ( nx + 1 ) * dx * ( ny + 1 ) * dy )

    IF ( bc_dirichlet_l )  THEN
       i = nxl
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1) + u(k,j,i) * dzw(k) * dy * rho_air(k)             &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r )  THEN
       i = nxr+1
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1) - u(k,j,i) * dzw(k) * dy * rho_air(k)             &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s )  THEN
       j = nys
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2) + v(k,j,i) * dzw(k) * dx * rho_air(k)             &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n )  THEN
       j = nyn+1
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2) - v(k,j,i) * dzw(k) * dx * rho_air(k)             &
                                * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
       ENDDO
    ENDIF
!
!-- Top boundary
    k = nzt
    DO  i = nxl, nxr
       DO  j = nys, nyn
          volume_flow_l(3) = volume_flow_l(3) - rho_air_zw(k) * w(k,j,i) * dx * dy
       ENDDO
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flow_l, volume_flow, 3, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flow = volume_flow_l
#endif

    w_correct = SUM( volume_flow ) * d_area_t * drho_air_zw(nzt)

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzt, nzt + 1
             w(k,j,i) = w(k,j,i) + w_correct                                                       &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
          ENDDO
       ENDDO
    ENDDO

    CALL  cpu_log( log_point(58), 'offline nesting', 'stop' )

    IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_mass_conservation', 'end' )

 END SUBROUTINE nesting_offl_mass_conservation


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the lateral and top boundary conditions in case the PALM domain is nested offline in a
!> mesoscale model. Further, average boundary data and determine mean profiles, further used for
!> correct damping in the sponge layer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_bc

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    INTEGER(iwp) ::  i  !< running index x-direction
    INTEGER(iwp) ::  j  !< running index y-direction
    INTEGER(iwp) ::  k  !< running index z-direction
    INTEGER(iwp) ::  n  !< running index for chemical species

    REAL(wp), DIMENSION(nzb:nzt+1) ::  pt_ref    !< reference profile for potential temperature
    REAL(wp), DIMENSION(nzb:nzt+1) ::  pt_ref_l  !< reference profile for potential temperature on subdomain
    REAL(wp), DIMENSION(nzb:nzt+1) ::  q_ref     !< reference profile for mixing ratio
    REAL(wp), DIMENSION(nzb:nzt+1) ::  q_ref_l   !< reference profile for mixing ratio on subdomain
    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_ref     !< reference profile for u-component
    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_ref_l   !< reference profile for u-component on subdomain
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_ref     !< reference profile for v-component
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_ref_l   !< reference profile for v-component on subdomain
    REAL(wp), DIMENSION(nzb:nzt+1) ::  var_1d    !< pre-interpolated profile for LOD1 mode

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_chem    !< reference profile for chemical species
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_chem_l  !< reference profile for chemical species on subdomain

    IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_bc', 'start' )

    CALL  cpu_log( log_point(58), 'offline nesting', 'start' )
!
!-- Initialize mean profiles, derived from boundary data, to zero.
    pt_ref   = 0.0_wp
    q_ref    = 0.0_wp
    u_ref    = 0.0_wp
    v_ref    = 0.0_wp

    pt_ref_l = 0.0_wp
    q_ref_l  = 0.0_wp
    u_ref_l  = 0.0_wp
    v_ref_l  = 0.0_wp
!
!-- If required, allocate temporary arrays to compute chemistry mean profiles
    IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
       ALLOCATE( ref_chem(nzb:nzt+1,1:UBOUND( chem_species, 1 ) )   )
       ALLOCATE( ref_chem_l(nzb:nzt+1,1:UBOUND( chem_species, 1 ) ) )
       ref_chem   = 0.0_wp
       ref_chem_l = 0.0_wp
    ENDIF
!
!-- Set boundary conditions of u-, v-, w-component, as well as q, and pt.
!-- Note, boundary values at the left boundary: i=-1 (v,w,pt,q) and i=0 (u), at the right boundary:
!-- i=nxr+1 (all), at the south boundary: j=-1 (u,w,pt,q) and j=0 (v), at the north boundary:
!-- j=nyn+1 (all).
!-- Please note, at the left (for u) and south (for v) boundary, values for u and v are set also at
!-- i/j=-1, since these values are used in boundary_conditions() to restore prognostic values.
!-- Further, sum up data to calculate mean profiles from boundary data, used for Rayleigh damping.
    IF ( bc_dirichlet_l  )  THEN
!
!--    u-component
       IF ( lod == 2 )  THEN
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                u(k,j,i_bound_u) = interpolate_in_time( nest_offl%u_l(0,k,j),                      &
                                                        nest_offl%u_l(1,k,j),                      &
                                                        fac_dt ) *                                 &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i_bound_u), 1 ) )
             ENDDO
             u(:,j,i_bound_u-1) = u(:,j,i_bound_u)
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j,i_bound_u)
          ENDDO
       ELSE
!
!--       Pre-interpolate profile before mapping onto the boundaries.
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%u_l(0,k,1),                                &
                                              nest_offl%u_l(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  j = nys, nyn
             u(nzb+1:nzt,j,i_bound_u) = var_1d(nzb+1:nzt) *                                        &
                                     MERGE( 1.0_wp, 0.0_wp,                                        &
                                            BTEST( topo_flags(nzb+1:nzt,j,i_bound_u), 1 ) )
             u(:,j,i_bound_u-1) = u(:,j,i_bound_u)
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j,i_bound_u)
          ENDDO
       ENDIF
!
!--    w-component
       IF ( lod == 2 )  THEN
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
                w(k,j,i_bound) = interpolate_in_time( nest_offl%w_l(0,k,j),                        &
                                                      nest_offl%w_l(1,k,j),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i_bound), 3 ) )
             ENDDO
             w(nzt,j,i_bound) = w(nzt-1,j,i_bound)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt-1
             var_1d(k) = interpolate_in_time( nest_offl%w_l(0,k,1),                                &
                                              nest_offl%w_l(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  j = nys, nyn
             w(nzb+1:nzt-1,j,i_bound) = var_1d(nzb+1:nzt-1) *                                      &
                                      MERGE( 1.0_wp, 0.0_wp,                                       &
                                             BTEST( topo_flags(nzb+1:nzt-1,j,i_bound), 3 ) )
             w(nzt,j,i_bound) = w(nzt-1,j,i_bound)
          ENDDO
       ENDIF
!
!--    v-component
       IF ( lod == 2 )  THEN
          DO  j = nysv, nyn
             DO  k = nzb+1, nzt
                v(k,j,i_bound) = interpolate_in_time( nest_offl%v_l(0,k,j),                        &
                                                      nest_offl%v_l(1,k,j),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i_bound), 2 ) )
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j,i_bound)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%v_l(0,k,1),                                &
                                              nest_offl%v_l(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  j = nysv, nyn
             v(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt) *                                          &
                                      MERGE( 1.0_wp, 0.0_wp,                                       &
                                             BTEST( topo_flags(nzb+1:nzt,j,i_bound), 2 ) )
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j,i_bound)
          ENDDO
       ENDIF
!
!--    Potential temperature
       IF ( .NOT. neutral )  THEN
          IF ( lod == 2 )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   pt(k,j,i_bound) = interpolate_in_time( nest_offl%pt_l(0,k,j),                   &
                                                          nest_offl%pt_l(1,k,j),                   &
                                                          fac_dt )
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j,i_bound)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%pt_l(0,k,1),                            &
                                                 nest_offl%pt_l(1,k,1),                            &
                                                 fac_dt )
             ENDDO
             DO  j = nys, nyn
                pt(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt)
                pt_ref_l(nzb+1:nzt)     = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j,i_bound)
             ENDDO
          ENDIF
       ENDIF
!
!--    Humidity
       IF ( humidity )  THEN
          IF ( lod == 2 )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   q(k,j,i_bound) = interpolate_in_time( nest_offl%q_l(0,k,j),                     &
                                                         nest_offl%q_l(1,k,j),                     &
                                                         fac_dt )
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j,i_bound)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%q_l(0,k,1),                             &
                                                 nest_offl%q_l(1,k,1),                             &
                                                 fac_dt )
             ENDDO
             DO  j = nys, nyn
                q(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt)
                q_ref_l(nzb+1:nzt)     = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j,i_bound)
             ENDDO
          ENDIF
       ENDIF
!
!--    Chemistry
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
             IF ( nest_offl%chem_from_file_l(n) )  THEN
                IF ( lod == 2 )  THEN
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,j,i_bound) = interpolate_in_time(                  &
                                                                        nest_offl%chem_l(0,k,j,n), &
                                                                        nest_offl%chem_l(1,k,j,n), &
                                                                        fac_dt                   )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)                            &
                                                + chem_species(n)%conc(nzb+1:nzt,j,i_bound)
                   ENDDO
                ELSE
                   DO  k = nzb+1, nzt
                      var_1d(k) = interpolate_in_time( nest_offl%chem_l(0,k,1,n),                  &
                                                       nest_offl%chem_l(1,k,1,n),                  &
                                                       fac_dt )
                   ENDDO
                   DO  j = nys, nyn
                      chem_species(n)%conc(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt)
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)                            &
                                                + chem_species(n)%conc(nzb+1:nzt,j,i_bound)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF

    ENDIF

    IF ( bc_dirichlet_r  )  THEN
!
!--    u-component
       IF ( lod == 2 )  THEN
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                u(k,j,i_bound_u) = interpolate_in_time( nest_offl%u_r(0,k,j),                      &
                                                        nest_offl%u_r(1,k,j),                      &
                                                        fac_dt ) *                                 &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i_bound_u), 1 ) )
             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j,i_bound_u)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%u_r(0,k,1),                                &
                                              nest_offl%u_r(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  j = nys, nyn
             u(nzb+1:nzt,j,i_bound_u) = var_1d(nzb+1:nzt) *                                        &
                                      MERGE( 1.0_wp, 0.0_wp,                                       &
                                             BTEST( topo_flags(nzb+1:nzt,j,i_bound_u), 1 ) )
             u_ref_l(nzb+1:nzt)       = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j,i_bound_u)
          ENDDO
       ENDIF
!
!--    w-component
       IF ( lod == 2 )  THEN
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
                w(k,j,i_bound) = interpolate_in_time( nest_offl%w_r(0,k,j),                        &
                                                      nest_offl%w_r(1,k,j),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i_bound), 3 ) )
             ENDDO
             w(nzt,j,i_bound) = w(nzt-1,j,i_bound)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt-1
             var_1d(k) = interpolate_in_time( nest_offl%w_r(0,k,1),                                &
                                              nest_offl%w_r(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  j = nys, nyn
             w(nzb+1:nzt-1,j,i_bound) = var_1d(nzb+1:nzt-1) *                                      &
                                  MERGE( 1.0_wp, 0.0_wp,                                           &
                                         BTEST( topo_flags(nzb+1:nzt-1,j,i_bound), 3 ) )
             w(nzt,j,i_bound) = w(nzt-1,j,i_bound)
          ENDDO
       ENDIF
!
!--    v-component
       IF ( lod == 2 )  THEN
          DO  j = nysv, nyn
             DO  k = nzb+1, nzt
                v(k,j,i_bound) = interpolate_in_time( nest_offl%v_r(0,k,j),                        &
                                                      nest_offl%v_r(1,k,j),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i_bound), 2 ) )
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j,i_bound)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%v_r(0,k,1),                                &
                                              nest_offl%v_r(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  j = nysv, nyn
             v(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt) *                                          &
                                    MERGE( 1.0_wp, 0.0_wp,                                         &
                                           BTEST( topo_flags(nzb+1:nzt,j,i_bound), 2 ) )
             v_ref_l(nzb+1:nzt)     = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j,i_bound)
          ENDDO
       ENDIF
!
!--    Potential temperature
       IF ( .NOT. neutral )  THEN
          IF ( lod == 2 )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   pt(k,j,i_bound) = interpolate_in_time( nest_offl%pt_r(0,k,j),                   &
                                                          nest_offl%pt_r(1,k,j),                   &
                                                          fac_dt )
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j,i_bound)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%pt_r(0,k,1),                            &
                                                 nest_offl%pt_r(1,k,1),                            &
                                                 fac_dt )
             ENDDO
             DO  j = nys, nyn
                pt(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt)
                pt_ref_l(nzb+1:nzt)     = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j,i_bound)
             ENDDO
          ENDIF
       ENDIF
!
!--    Humidity
       IF ( humidity )  THEN
          IF ( lod == 2 )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   q(k,j,i_bound) = interpolate_in_time( nest_offl%q_r(0,k,j),                     &
                                                         nest_offl%q_r(1,k,j),                     &
                                                         fac_dt )
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j,i_bound)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%q_r(0,k,1),                             &
                                                 nest_offl%q_r(1,k,1),                             &
                                                 fac_dt )
             ENDDO
             DO  j = nys, nyn
                q(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt)
                q_ref_l(nzb+1:nzt)     = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j,i_bound)
             ENDDO
          ENDIF
       ENDIF
!
!--    Chemistry
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
             IF ( nest_offl%chem_from_file_r(n) )  THEN
                IF ( lod == 2 )  THEN
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,j,i_bound) = interpolate_in_time(                  &
                                                                        nest_offl%chem_r(0,k,j,n), &
                                                                        nest_offl%chem_r(1,k,j,n), &
                                                                        fac_dt                   )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)                            &
                                                + chem_species(n)%conc(nzb+1:nzt,j,i_bound)
                   ENDDO
                ELSE
                   DO  k = nzb+1, nzt
                      var_1d(k) = interpolate_in_time( nest_offl%chem_r(0,k,1,n),                  &
                                                       nest_offl%chem_r(1,k,1,n),                  &
                                                       fac_dt )
                   ENDDO
                   DO  j = nys, nyn
                      chem_species(n)%conc(nzb+1:nzt,j,i_bound) = var_1d(nzb+1:nzt)
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)                            &
                                                + chem_species(n)%conc(nzb+1:nzt,j,i_bound)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF

    ENDIF

    IF ( bc_dirichlet_n )  THEN
!
!--    v-component
       IF ( lod == 2 )  THEN
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                v(k,j_bound_v,i) = interpolate_in_time( nest_offl%v_n(0,k,i),                      &
                                                        nest_offl%v_n(1,k,i),                      &
                                                        fac_dt ) *                                 &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j_bound_v,i), 2 ) )
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j_bound_v,i)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%v_n(0,k,1),                                &
                                              nest_offl%v_n(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  i = nxl, nxr
             v(nzb+1:nzt,j_bound_v,i) = var_1d(nzb+1:nzt) *                                        &
                                  MERGE( 1.0_wp, 0.0_wp,                                           &
                                         BTEST( topo_flags(nzb+1:nzt,j_bound_v,i), 2 ) )
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j_bound_v,i)
          ENDDO
       ENDIF
!
!--    w-component
       IF ( lod == 2 )  THEN
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt-1
                w(k,j_bound,i) = interpolate_in_time( nest_offl%w_n(0,k,i),                        &
                                                      nest_offl%w_n(1,k,i),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j_bound,i), 3 ) )
             ENDDO
             w(nzt,j_bound,i) = w(nzt-1,j_bound,i)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt-1
             var_1d(k) = interpolate_in_time( nest_offl%w_n(0,k,1),                                &
                                              nest_offl%w_n(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  i = nxl, nxr
             w(nzb+1:nzt-1,j_bound,i) = var_1d(nzb+1:nzt-1) *                                      &
                                  MERGE( 1.0_wp, 0.0_wp,                                           &
                                         BTEST( topo_flags(nzb+1:nzt-1,j_bound,i), 3 ) )
             w(nzt,j_bound,i) = w(nzt-1,j_bound,i)
          ENDDO
       ENDIF
!
!--    u-component
       IF ( lod == 2 )  THEN
          DO  i = nxlu, nxr
             DO  k = nzb+1, nzt
                u(k,j_bound,i) = interpolate_in_time( nest_offl%u_n(0,k,i),                        &
                                                      nest_offl%u_n(1,k,i),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j_bound,i), 1 ) )
             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j_bound,i)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%u_n(0,k,1),                                &
                                              nest_offl%u_n(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  i = nxlu, nxr
             u(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt) *                                          &
                                    MERGE( 1.0_wp, 0.0_wp,                                         &
                                           BTEST( topo_flags(nzb+1:nzt,j_bound,i), 1 ) )
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j_bound,i)
          ENDDO
       ENDIF
!
!--    Potential temperature
       IF ( .NOT. neutral )  THEN
          IF ( lod == 2 )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   pt(k,j_bound,i) = interpolate_in_time( nest_offl%pt_n(0,k,i),                   &
                                                          nest_offl%pt_n(1,k,i),                   &
                                                          fac_dt )
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j_bound,i)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%pt_n(0,k,1),                            &
                                                 nest_offl%pt_n(1,k,1),                            &
                                                 fac_dt )
             ENDDO
             DO  i = nxl, nxr
                pt(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt)
                pt_ref_l(nzb+1:nzt)     = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j_bound,i)
             ENDDO
          ENDIF
       ENDIF
!
!--    Humidity
       IF ( humidity )  THEN
          IF ( lod == 2 )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   q(k,j_bound,i) = interpolate_in_time( nest_offl%q_n(0,k,i),                     &
                                                         nest_offl%q_n(1,k,i),                     &
                                                         fac_dt )
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j_bound,i)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%q_n(0,k,1),                             &
                                                 nest_offl%q_n(1,k,1),                             &
                                                 fac_dt )
             ENDDO
             DO  i = nxl, nxr
                q(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt)
                q_ref_l(nzb+1:nzt)     = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j_bound,i)
             ENDDO
          ENDIF
       ENDIF
!
!--    Chemistry
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
             IF ( nest_offl%chem_from_file_n(n) )  THEN
                IF ( lod == 2 )  THEN
                   DO  i = nxl, nxr
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,j_bound,i) = interpolate_in_time(                  &
                                                                     nest_offl%chem_n(0,k,i,n),    &
                                                                     nest_offl%chem_n(1,k,i,n),    &
                                                                     fac_dt                    )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)                            &
                                                + chem_species(n)%conc(nzb+1:nzt,j_bound,i)
                   ENDDO
                ELSE
                   DO  k = nzb+1, nzt
                      var_1d(k) = interpolate_in_time( nest_offl%chem_n(0,k,1,n),                  &
                                                       nest_offl%chem_n(1,k,1,n),                  &
                                                       fac_dt )
                   ENDDO
                   DO  i = nxl, nxr
                      chem_species(n)%conc(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt)
                      ref_chem_l(nzb+1:nzt,n)                   = ref_chem_l(nzb+1:nzt,n) +        &
                                                                  chem_species(n)                  &
                                                                  %conc(nzb+1:nzt,j_bound,i)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    IF ( bc_dirichlet_s )  THEN
!
!--    v-component
       IF ( lod == 2 )  THEN
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                v(k,j_bound_v,i) = interpolate_in_time( nest_offl%v_s(0,k,i),                      &
                                                        nest_offl%v_s(1,k,i),                      &
                                                        fac_dt ) *                                 &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j_bound_v,i), 2 ) )
             ENDDO
             v(:,j_bound_v-1,i) = v(:,j_bound_v,i)
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j_bound_v,i)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%v_s(0,k,1),                                &
                                              nest_offl%v_s(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  i = nxl, nxr
             v(nzb+1:nzt,j_bound_v,i) = var_1d(nzb+1:nzt) *                                        &
                                      MERGE( 1.0_wp, 0.0_wp,                                       &
                                             BTEST( topo_flags(nzb+1:nzt,j_bound_v,i), 2 ) )
             v(:,j_bound_v-1,i) = v(:,j_bound_v,i)
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j_bound_v,i)
          ENDDO
       ENDIF
!
!--    w-component
       IF ( lod == 2 )  THEN
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt-1
                w(k,j_bound,i) = interpolate_in_time( nest_offl%w_s(0,k,i),                        &
                                                      nest_offl%w_s(1,k,i),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j_bound,i), 3 ) )
             ENDDO
             w(nzt,j_bound,i) = w(nzt-1,j_bound,i)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt-1
             var_1d(k) = interpolate_in_time( nest_offl%w_s(0,k,1),                                &
                                              nest_offl%w_s(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  i = nxl, nxr
             w(nzb+1:nzt-1,j_bound,i) = var_1d(nzb+1:nzt-1) *                                      &
                                      MERGE( 1.0_wp, 0.0_wp,                                       &
                                             BTEST( topo_flags(nzb+1:nzt-1,j_bound,i), 3 ) )
             w(nzt,j_bound,i) = w(nzt-1,j_bound,i)
          ENDDO
       ENDIF
!
!--    u-component
       IF ( lod == 2 )  THEN
          DO  i = nxlu, nxr
             DO  k = nzb+1, nzt
                u(k,j_bound,i) = interpolate_in_time( nest_offl%u_s(0,k,i),                        &
                                                      nest_offl%u_s(1,k,i),                        &
                                                      fac_dt ) *                                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j_bound,i), 1 ) )
             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j_bound,i)
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             var_1d(k) = interpolate_in_time( nest_offl%u_s(0,k,1),                                &
                                              nest_offl%u_s(1,k,1),                                &
                                              fac_dt )
          ENDDO
          DO  i = nxlu, nxr
             u(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt) *                                          &
                                      MERGE( 1.0_wp, 0.0_wp,                                       &
                                             BTEST( topo_flags(nzb+1:nzt,j_bound,i), 1 ) )
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j_bound,i)
          ENDDO
       ENDIF
!
!--    Potential temperature
       IF ( .NOT. neutral )  THEN
          IF ( lod == 2 )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   pt(k,j_bound,i) = interpolate_in_time( nest_offl%pt_s(0,k,i),                   &
                                                          nest_offl%pt_s(1,k,i),                   &
                                                          fac_dt )
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j_bound,i)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%pt_s(0,k,1),                            &
                                                 nest_offl%pt_s(1,k,1),                            &
                                                 fac_dt )
             ENDDO
             DO  i = nxl, nxr
                pt(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt)
                pt_ref_l(nzb+1:nzt)     = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j_bound,i)
             ENDDO
          ENDIF
       ENDIF
!
!--    Humidity
       IF ( humidity )  THEN
          IF ( lod == 2 )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   q(k,j_bound,i) = interpolate_in_time( nest_offl%q_s(0,k,i),                     &
                                                         nest_offl%q_s(1,k,i),                     &
                                                         fac_dt )
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j_bound,i)
             ENDDO
          ELSE
             DO  k = nzb+1, nzt
                var_1d(k) = interpolate_in_time( nest_offl%q_s(0,k,1),                             &
                                                 nest_offl%q_s(1,k,1),                             &
                                                 fac_dt )
             ENDDO
             DO  i = nxl, nxr
                q(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt)
                q_ref_l(nzb+1:nzt)     = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j_bound,i)
             ENDDO
          ENDIF
       ENDIF
!
!--    Chemistry
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
             IF ( nest_offl%chem_from_file_s(n) )  THEN
                IF ( lod == 2 )  THEN
                   DO  i = nxl, nxr
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,j_bound,i) = interpolate_in_time(                  &
                                                                        nest_offl%chem_s(0,k,i,n), &
                                                                        nest_offl%chem_s(1,k,i,n), &
                                                                        fac_dt  )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)                            &
                                                + chem_species(n)%conc(nzb+1:nzt,j_bound,i)
                   ENDDO
                ELSE
                   DO  k = nzb+1, nzt
                      var_1d(k) = interpolate_in_time( nest_offl%chem_s(0,k,1,n),                  &
                                                       nest_offl%chem_s(1,k,1,n),                  &
                                                       fac_dt )
                   ENDDO
                   DO  i = nxl, nxr
                      chem_species(n)%conc(nzb+1:nzt,j_bound,i) = var_1d(nzb+1:nzt)
                      ref_chem_l(nzb+1:nzt,n)                   = ref_chem_l(nzb+1:nzt,n) +        &
                                                                  chem_species(n)                  &
                                                                  %conc(nzb+1:nzt,j_bound,i)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDIF
!
!-- Top boundary
!-- u-component
    IF ( lod == 2 )  THEN
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             u(nzt+1,j,i) = interpolate_in_time( nest_offl%u_top(0,j,i),                           &
                                                 nest_offl%u_top(1,j,i),                           &
                                                 fac_dt ) *                                        &
                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzt+1,j,i), 1 ) )
             u_ref_l(nzt+1) = u_ref_l(nzt+1) + u(nzt+1,j,i)
          ENDDO
       ENDDO
    ELSE
       var_1d(nzt+1) = interpolate_in_time( nest_offl%u_top(0,1,1),                                &
                                            nest_offl%u_top(1,1,1),                                &
                                            fac_dt )
       u(nzt+1,nys:nyn,nxlu:nxr) = var_1d(nzt+1) *                                                 &
                                   MERGE( 1.0_wp, 0.0_wp,                                          &
                                          BTEST( topo_flags(nzt+1,nys:nyn,nxlu:nxr), 1 ) )
       u_ref_l(nzt+1) = u_ref_l(nzt+1) + SUM( u(nzt+1,nys:nyn,nxlu:nxr) )
    ENDIF
!
!--    For left boundary set boundary condition for u-component also at top grid point.
!--    Note, this has no effect on the numeric solution, only for data output.
    IF ( bc_dirichlet_l )  u(nzt+1,:,nxl) = u(nzt+1,:,nxlu)
!
!-- v-component
    IF ( lod == 2 )  THEN
       DO  i = nxl, nxr
          DO  j = nysv, nyn
             v(nzt+1,j,i) = interpolate_in_time( nest_offl%v_top(0,j,i),                           &
                                                 nest_offl%v_top(1,j,i),                           &
                                                 fac_dt ) *                                        &
                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzt+1,j,i), 2 ) )
             v_ref_l(nzt+1) = v_ref_l(nzt+1) + v(nzt+1,j,i)
          ENDDO
       ENDDO
    ELSE
       var_1d(nzt+1) = interpolate_in_time( nest_offl%v_top(0,1,1),                                &
                                            nest_offl%v_top(1,1,1),                                &
                                            fac_dt )
       v(nzt+1,nysv:nyn,nxl:nxr) = var_1d(nzt+1) *                                                 &
                                   MERGE( 1.0_wp, 0.0_wp,                                          &
                                          BTEST( topo_flags(nzt+1,nysv:nyn,nxl:nxr), 2 ) )
       v_ref_l(nzt+1) = v_ref_l(nzt+1) + SUM( v(nzt+1,nysv:nyn,nxl:nxr) )
    ENDIF
!
!-- For south boundary set boundary condition for v-component also at top grid point.
!-- Note, this has no effect on the numeric solution, only for data output.
    IF ( bc_dirichlet_s )  v(nzt+1,nys,:) = v(nzt+1,nysv,:)
!
!-- w-component
    IF ( lod == 2 )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             w(nzt,j,i) = interpolate_in_time( nest_offl%w_top(0,j,i),                             &
                                               nest_offl%w_top(1,j,i),                             &
                                               fac_dt ) *                                          &
                          MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzt,j,i), 3 ) )
             w(nzt+1,j,i) = w(nzt,j,i)
          ENDDO
       ENDDO
    ELSE
       var_1d(nzt) = interpolate_in_time( nest_offl%w_top(0,1,1),                                  &
                                          nest_offl%w_top(1,1,1),                                  &
                                          fac_dt )
       w(nzt,nys:nyn,nxl:nxr) = var_1d(nzt) *                                                      &
                                MERGE( 1.0_wp, 0.0_wp,                                             &
                                       BTEST( topo_flags(nzt,nys:nyn,nxl:nxr), 3 ) )
       w(nzt+1,nys:nyn,nxl:nxr) = w(nzt,nys:nyn,nxl:nxr)
    ENDIF
!
!-- Potential temperture
    IF ( .NOT. neutral )  THEN
       IF ( lod == 2 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                pt(nzt+1,j,i) = interpolate_in_time( nest_offl%pt_top(0,j,i),                      &
                                                     nest_offl%pt_top(1,j,i),                      &
                                                     fac_dt )
                pt_ref_l(nzt+1) = pt_ref_l(nzt+1) + pt(nzt+1,j,i)
             ENDDO
          ENDDO
       ELSE
          var_1d(nzt+1) = interpolate_in_time( nest_offl%pt_top(0,1,1),                            &
                                               nest_offl%pt_top(1,1,1),                            &
                                               fac_dt )
          pt(nzt+1,nys:nyn,nxl:nxr) = var_1d(nzt+1)
          pt_ref_l(nzt+1) = pt_ref_l(nzt+1) + SUM( pt(nzt+1,nys:nyn,nxl:nxr) )
       ENDIF
    ENDIF
!
!--    humidity
    IF ( humidity )  THEN
       IF ( lod == 2 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                q(nzt+1,j,i) = interpolate_in_time( nest_offl%q_top(0,j,i),                        &
                                                    nest_offl%q_top(1,j,i),                        &
                                                    fac_dt )
                q_ref_l(nzt+1) = q_ref_l(nzt+1) + q(nzt+1,j,i)
             ENDDO
          ENDDO
       ELSE
          var_1d(nzt+1) = interpolate_in_time( nest_offl%q_top(0,1,1),                             &
                                               nest_offl%q_top(1,1,1),                             &
                                               fac_dt )
          q(nzt+1,nys:nyn,nxl:nxr) = var_1d(nzt+1)
          q_ref_l(nzt+1) = q_ref_l(nzt+1) + SUM( q(nzt+1,nys:nyn,nxl:nxr) )
       ENDIF
    ENDIF

    IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
       DO  n = 1, UBOUND( chem_species, 1 )
          IF ( nest_offl%chem_from_file_t(n) )  THEN
             IF ( lod == 2 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      chem_species(n)%conc(nzt+1,j,i) = interpolate_in_time(                       &
                                                           nest_offl%chem_top(0,j,i,n),            &
                                                           nest_offl%chem_top(1,j,i,n),            &
                                                           fac_dt          )
                      ref_chem_l(nzt+1,n) = ref_chem_l(nzt+1,n) + chem_species(n)%conc(nzt+1,j,i)
                   ENDDO
                ENDDO
             ELSE
                var_1d(nzt+1) = interpolate_in_time( nest_offl%chem_top(0,1,1,n),                  &
                                                     nest_offl%chem_top(1,1,1,n),                  &
                                                     fac_dt )
                chem_species(n)%conc(nzt+1,nys:nyn,nxl:nxr) = var_1d(nzt+1)
                ref_chem_l(nzt+1,n) = ref_chem_l(nzt+1,n) +                                        &
                                      SUM( chem_species(n)%conc(nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
          ENDIF
       ENDDO
    ENDIF
!
!-- Moreover, set Neumann boundary condition for subgrid-scale TKE, passive scalar, dissipation, and
!-- chemical species if required.
    IF ( rans_mode  .AND.  rans_tke_e )  THEN
       IF (  bc_dirichlet_l )  diss(:,:,nxl-1) = diss(:,:,nxl)
       IF (  bc_dirichlet_r )  diss(:,:,nxr+1) = diss(:,:,nxr)
       IF (  bc_dirichlet_s )  diss(:,nys-1,:) = diss(:,nys,:)
       IF (  bc_dirichlet_n )  diss(:,nyn+1,:) = diss(:,nyn,:)
    ENDIF
!        IF ( .NOT. constant_diffusion )  THEN
!           IF (  bc_dirichlet_l )  e(:,:,nxl-1) = e(:,:,nxl)
!           IF (  bc_dirichlet_r )  e(:,:,nxr+1) = e(:,:,nxr)
!           IF (  bc_dirichlet_s )  e(:,nys-1,:) = e(:,nys,:)
!           IF (  bc_dirichlet_n )  e(:,nyn+1,:) = e(:,nyn,:)
!           e(nzt+1,:,:) = e(nzt,:,:)
!        ENDIF
!        IF ( passive_scalar )  THEN
!           IF (  bc_dirichlet_l )  s(:,:,nxl-1) = s(:,:,nxl)
!           IF (  bc_dirichlet_r )  s(:,:,nxr+1) = s(:,:,nxr)
!           IF (  bc_dirichlet_s )  s(:,nys-1,:) = s(:,nys,:)
!           IF (  bc_dirichlet_n )  s(:,nyn+1,:) = s(:,nyn,:)
!        ENDIF

    CALL exchange_horiz( u, nbgp )
    CALL exchange_horiz( v, nbgp )
    CALL exchange_horiz( w, nbgp )
    IF ( .NOT. neutral )  CALL exchange_horiz( pt, nbgp )
    IF ( humidity      )  CALL exchange_horiz( q,  nbgp )
    IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
       DO  n = 1, UBOUND( chem_species, 1 )
!
!--       Do local exchange only when necessary, i.e. when data is coming from dynamic file.
          IF ( nest_offl%chem_from_file_t(n) )  CALL exchange_horiz( chem_species(n)%conc, nbgp )
       ENDDO
    ENDIF
!
!-- Set top boundary condition at all horizontal grid points, also at the lateral boundary grid
!-- points.
    w(nzt+1,:,:) = w(nzt,:,:)
!
!-- Offline nesting for salsa
    IF ( salsa )  CALL salsa_nesting_offl_bc
!
!-- Calculate the mean profiles. These are later stored on u_init, v_init, etc., in order to adjust
!-- the Rayleigh damping under time-evolving atmospheric conditions accordingly - damping against
!-- the representative mean profiles, not against the initial profiles. Note, in LOD = 1 case no
!-- averaging is required.
#if defined( __parallel )
    CALL MPI_ALLREDUCE( u_ref_l, u_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( v_ref_l, v_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
    IF ( humidity )  THEN
       CALL MPI_ALLREDUCE( q_ref_l, q_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
    ENDIF
    IF ( .NOT. neutral )  THEN
       CALL MPI_ALLREDUCE( pt_ref_l, pt_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
    ENDIF
    IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
       CALL MPI_ALLREDUCE( ref_chem_l, ref_chem, ( nzt+1-nzb+1 ) * SIZE( ref_chem(nzb,:) ),        &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
    ENDIF
#else
    u_ref  = u_ref_l
    v_ref  = v_ref_l
    IF ( humidity )       q_ref    = q_ref_l
    IF ( .NOT. neutral )  pt_ref   = pt_ref_l
    IF ( air_chemistry  .AND.  nesting_offline_chem )  ref_chem = ref_chem_l
#endif
!
!-- Average data. Note, reference profiles up to nzt are derived from lateral boundaries, at the
!-- model top it is derived from the top boundary. Thus, number of input data is different from
!-- nzb:nzt compared to nzt+1.
!-- Derived from lateral boundaries.
    u_ref(nzb:nzt) = u_ref(nzb:nzt) / REAL( 2.0_wp * ( ny + 1 + nx     ), KIND = wp )
    v_ref(nzb:nzt) = v_ref(nzb:nzt) / REAL( 2.0_wp * ( ny   + nx + 1   ), KIND = wp )
    IF ( humidity )                                                                                &
       q_ref(nzb:nzt) = q_ref(nzb:nzt)   / REAL( 2.0_wp * ( ny + 1 + nx + 1 ), KIND = wp )
    IF ( .NOT. neutral )                                                                           &
       pt_ref(nzb:nzt) = pt_ref(nzb:nzt) / REAL( 2.0_wp * ( ny + 1 + nx + 1 ), KIND = wp )
    IF ( air_chemistry  .AND.  nesting_offline_chem )                                              &
       ref_chem(nzb:nzt,:) = ref_chem(nzb:nzt,:) / REAL( 2.0_wp * ( ny + 1 + nx + 1 ), KIND = wp )
!
!-- Derived from top boundary.
    u_ref(nzt+1) = u_ref(nzt+1) / REAL( ( ny + 1 ) * ( nx     ), KIND = wp )
    v_ref(nzt+1) = v_ref(nzt+1) / REAL( ( ny     ) * ( nx + 1 ), KIND = wp )
    IF ( humidity )                                                                                &
       q_ref(nzt+1) = q_ref(nzt+1)   / REAL( ( ny + 1 ) * ( nx + 1 ), KIND = wp )
    IF ( .NOT. neutral )                                                                           &
       pt_ref(nzt+1) = pt_ref(nzt+1) / REAL( ( ny + 1 ) * ( nx + 1 ), KIND = wp )
    IF ( air_chemistry  .AND.  nesting_offline_chem )                                              &
       ref_chem(nzt+1,:) = ref_chem(nzt+1,:) / REAL( ( ny + 1 ) * ( nx + 1 ),KIND = wp )
!
!-- Write onto init profiles, which are used for damping. Also set lower boundary condition for
!-- scalars (not required for u and v as these are zero at k=nzb).
    u_init = u_ref
    v_init = v_ref
    IF ( humidity      )  THEN
       q_init      = q_ref
       q_init(nzb) = q_init(nzb+1)
    ENDIF
    IF ( .NOT. neutral )  THEN
       pt_init      = pt_ref
       pt_init(nzb) = pt_init(nzb+1)
    ENDIF

    IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
       DO  n = 1, UBOUND( chem_species, 1 )
          IF ( nest_offl%chem_from_file_t(n) )  THEN
             chem_species(n)%conc_pr_init(:)   = ref_chem(:,n)
             chem_species(n)%conc_pr_init(nzb) = chem_species(n)%conc_pr_init(nzb+1)
          ENDIF
       ENDDO
    ENDIF
    IF ( ALLOCATED( ref_chem   ) )  DEALLOCATE( ref_chem   )
    IF ( ALLOCATED( ref_chem_l ) )  DEALLOCATE( ref_chem_l )
!
!-- Further, adjust Rayleigh damping height in case of time-changing conditions.
!-- Therefore, calculate boundary-layer depth first.
    CALL nesting_offl_calc_zi
    CALL adjust_sponge_layer

    CALL  cpu_log( log_point(58), 'offline nesting', 'stop' )

    IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_bc', 'end' )


 END SUBROUTINE nesting_offl_bc

!--------------------------------------------------------------------------------------------------!
! Description:
!--------------------------------------------------------------------------------------------------!
!> Determine the interpolation constant for time interpolation. The calculation is separated from
!> nesting_offl_bc in order to be independent on the order of calls.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_interpolation_factor
!
!-- Determine interpolation factor and limit it to 1. This is because t+dt can slightly exceed
!-- time(tind_p) before boundary data is updated again.
    fac_dt = ( time_since_reference_point - nest_offl%time(nest_offl%tind) + dt_3d ) /             &
             ( nest_offl%time(nest_offl%tind_p) - nest_offl%time(nest_offl%tind) )

    fac_dt = MIN( 1.0_wp, fac_dt )

 END SUBROUTINE nesting_offl_interpolation_factor

!--------------------------------------------------------------------------------------------------!
! Description:
!--------------------------------------------------------------------------------------------------!
!> Calculates the boundary-layer depth from the boundary data, according to bulk-Richardson
!> criterion.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_calc_zi

    INTEGER(iwp) :: i                             !< loop index in x-direction
    INTEGER(iwp) :: j                             !< loop index in y-direction
    INTEGER(iwp) :: k                             !< loop index in z-direction
    INTEGER(iwp) :: k_max_loc                     !< index of maximum wind speed along z-direction
    INTEGER(iwp) :: k_surface                     !< topography top index in z-direction
    INTEGER(iwp) :: num_boundary_gp_non_cyclic    !< number of non-cyclic boundaries, used for averaging ABL depth
    INTEGER(iwp) :: num_boundary_gp_non_cyclic_l  !< number of non-cyclic boundaries, used for averaging ABL depth

    REAL(wp) ::  ri_bulk                 !< bulk Richardson number
    REAL(wp) ::  ri_bulk_crit = 0.25_wp  !< critical bulk Richardson number
    REAL(wp) ::  topo_max                !< maximum topography level in model domain
    REAL(wp) ::  topo_max_l              !< maximum topography level in subdomain
    REAL(wp) ::  vpt_surface             !< near-surface virtual potential temperature
    REAL(wp) ::  zi_l                    !< mean boundary-layer depth on subdomain
    REAL(wp) ::  zi_local                !< local boundary-layer depth

    REAL(wp), DIMENSION(nzb:nzt+1) ::  vpt_col  !< vertical profile of virtual potential temperature at (j,i)-grid point
    REAL(wp), DIMENSION(nzb:nzt+1) ::  uv_abs   !< vertical profile of horizontal wind speed at (j,i)-grid point


!
!-- Calculate mean boundary-layer height from boundary data.
!-- Start with the left and right boundaries.
    zi_l      = 0.0_wp
    num_boundary_gp_non_cyclic_l = 0
    IF ( bc_dirichlet_l  .OR.  bc_dirichlet_r )  THEN
!
!--    Sum-up and store number of boundary grid points used for averaging ABL depth
       num_boundary_gp_non_cyclic_l = num_boundary_gp_non_cyclic_l + nxr - nxl + 1
!
!--    Determine index along x. Please note, index indicates boundary grid point for scalars.
       i = MERGE( -1, nxr + 1, bc_dirichlet_l )

       DO  j = nys, nyn
!
!--       Determine topography top index at current (j,i) index
          k_surface = topo_top_ind(j,i,0)
!
!--       Pre-compute surface virtual temperature. Therefore, use 2nd prognostic level according to
!--       Heinze et al. (2017).
          IF ( humidity )  THEN
             vpt_surface = pt(k_surface+2,j,i) * ( 1.0_wp + 0.61_wp * q(k_surface+2,j,i) )
             vpt_col     = pt(:,j,i) * ( 1.0_wp + 0.61_wp * q(:,j,i) )
          ELSE
             vpt_surface = pt(k_surface+2,j,i)
             vpt_col     = pt(:,j,i)
          ENDIF
!
!--       Calculate local boundary layer height from bulk Richardson number, i.e. the height where
!--       the bulk Richardson number exceeds its critical value of 0.25
!--       (according to Heinze et al., 2017).
!--       Note, no interpolation of u- and v-component is made, as both are mainly mean inflow
!--       profiles with very small spatial variation.
!--       Add a safety factor in case the velocity term becomes zero. This may happen if overhanging
!--       3D structures are directly located at the boundary, where velocity inside the building is
!--       zero (k_surface is the index of the lowest upward-facing surface).
          uv_abs(:) = SQRT( MERGE( u(:,j,i+1), u(:,j,i), bc_dirichlet_l )**2 + v(:,j,i)**2 )
!
!--       Determine index of the maximum wind speed
          k_max_loc = MAXLOC( uv_abs(:), DIM = 1 ) - 1

          zi_local = 0.0_wp
          DO  k = k_surface+1, nzt
             ri_bulk = zu(k) * g / vpt_surface *                                                   &
                       ( vpt_col(k) - vpt_surface ) / ( uv_abs(k) + 1E-5_wp )
!
!--          Check if critical Richardson number is exceeded. Further, check if there is a maxium in
!--          the wind profile in order to detect also ABL heights in the stable boundary layer.
             IF ( zi_local == 0.0_wp  .AND.  ( ri_bulk > ri_bulk_crit .OR. k == k_max_loc ) )      &
                zi_local = zu(k)
          ENDDO
!
!--       Assure that the minimum local boundary-layer depth is at least at the second vertical grid
!--       level.
          zi_l = zi_l + MAX( zi_local, zu(k_surface+2) )

       ENDDO

    ENDIF
!
!-- Do the same at the north and south boundaries.
    IF ( bc_dirichlet_s  .OR.  bc_dirichlet_n )  THEN

       num_boundary_gp_non_cyclic_l = num_boundary_gp_non_cyclic_l + nxr - nxl + 1

       j = MERGE( -1, nyn + 1, bc_dirichlet_s )

       DO  i = nxl, nxr
          k_surface = topo_top_ind(j,i,0)

          IF ( humidity )  THEN
             vpt_surface = pt(k_surface+2,j,i) * ( 1.0_wp + 0.61_wp * q(k_surface+2,j,i) )
             vpt_col     = pt(:,j,i) * ( 1.0_wp + 0.61_wp * q(:,j,i) )
          ELSE
             vpt_surface = pt(k_surface+2,j,i)
             vpt_col  = pt(:,j,i)
          ENDIF

          uv_abs(:) = SQRT( u(:,j,i)**2 + MERGE( v(:,j+1,i), v(:,j,i), bc_dirichlet_s )**2 )
!
!--       Determine index of the maximum wind speed
          k_max_loc = MAXLOC( uv_abs(:), DIM = 1 ) - 1

          zi_local = 0.0_wp
          DO  k = k_surface+1, nzt
             ri_bulk = zu(k) * g / vpt_surface *                                                   &
                       ( vpt_col(k) - vpt_surface ) / ( uv_abs(k) + 1E-5_wp )
!
!--          Check if critical Richardson number is exceeded. Further, check if there is a maxium in
!--          the wind profile in order to detect also ABL heights in the stable boundary layer.
             IF ( zi_local == 0.0_wp  .AND.  ( ri_bulk > ri_bulk_crit .OR. k == k_max_loc ) )      &
                zi_local = zu(k)
          ENDDO
          zi_l = zi_l + MAX( zi_local, zu(k_surface+2) )

       ENDDO

    ENDIF

#if defined( __parallel )
    CALL MPI_ALLREDUCE( zi_l, zi_ribulk, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( num_boundary_gp_non_cyclic_l, num_boundary_gp_non_cyclic,                  &
                        1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    zi_ribulk = zi_l
    num_boundary_gp_non_cyclic = num_boundary_gp_non_cyclic_l
#endif
    zi_ribulk = zi_ribulk / REAL( num_boundary_gp_non_cyclic, KIND = wp )
!
!-- Finally, check if boundary layer depth is not below the any topography.
!-- zi_ribulk will be used to adjust rayleigh damping height, i.e. the lower level of the sponge
!-- layer, as well as to adjust the synthetic turbulence generator accordingly. If Rayleigh damping
!-- would be applied near buildings, etc., this would spoil the simulation results.
    topo_max_l = zw(MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) ))

#if defined( __parallel )
    CALL MPI_ALLREDUCE( topo_max_l, topo_max, 1, MPI_REAL, MPI_MAX, comm2d, ierr )
#else
    topo_max     = topo_max_l
#endif
!        zi_ribulk = MAX( zi_ribulk, topo_max )

 END SUBROUTINE nesting_offl_calc_zi


!--------------------------------------------------------------------------------------------------!
! Description:
!--------------------------------------------------------------------------------------------------!
!> Adjust the height where the rayleigh damping starts, i.e. the lower level of the sponge layer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE adjust_sponge_layer

    INTEGER(iwp) :: k   !< loop index in z-direction

    REAL(wp) ::  rdh    !< updated Rayleigh damping height


    IF ( rayleigh_damping_height > 0.0_wp  .AND.  rayleigh_damping_factor > 0.0_wp )  THEN
!
!--    Update Rayleigh-damping height and re-calculate height-depending damping coefficients.
!--    Assure that rayleigh damping starts well above the boundary layer.
       rdh = MIN( MAX( zi_ribulk * 1.3_wp, 10.0_wp * dz(1) ),                                      &
                  0.8_wp * zu(nzt), rayleigh_damping_height )
!
!--       Update Rayleigh damping factor
       DO  k = nzb+1, nzt
          IF ( zu(k) >= rdh )  THEN
             rdf(k) = rayleigh_damping_factor *                                                    &
                      ( SIN( pi * 0.5_wp * ( zu(k) - rdh ) / ( zu(nzt) - rdh ) ) )**2
          ENDIF
       ENDDO
       rdf_sc = rdf

    ENDIF

 END SUBROUTINE adjust_sponge_layer

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Performs consistency checks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_check_parameters
!
!-- Check if offline nesting is applied in nested child domain.
    IF ( nesting_offline  .AND.  child_domain )  THEN
       message_string = 'Offline nesting is only applicable in root model.'
       CALL message( 'offline_nesting_check_parameters', 'PA0622', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE nesting_offl_check_parameters

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the parameter list nesting_offl_parameters
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /nesting_offl_parameters/  switch_off_module


!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, nesting_offl_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    nesting_offl_parameters namelist was found and read correctly. Enable the
!--    offline nesting.
       IF ( .NOT. switch_off_module )  nesting_offline = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    nesting_offl_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'nesting_offl_parameters', line )

    ENDIF

 END SUBROUTINE nesting_offl_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes information about offline nesting into HEADER file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_header ( io )

    INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file

    WRITE ( io, 1 )
    IF ( nesting_offline )  THEN
       WRITE ( io, 3 )
    ELSE
       WRITE ( io, 2 )
    ENDIF

1 FORMAT (//' Offline nesting in COSMO model:'/' -------------------------------'/)
2 FORMAT (' --> No offlince nesting is used (default) ')
3 FORMAT (' --> Offlince nesting is used. Boundary data is read from dynamic input file ')

 END SUBROUTINE nesting_offl_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays used to read boundary data from NetCDF file and initialize boundary data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_init

    INTEGER(iwp) ::  i   !< loop index for x-direction
    INTEGER(iwp) ::  j   !< loop index for y-direction
    INTEGER(iwp) ::  n   !< running index for chemical species

!
!-- Before arrays for the boundary data are allocated, the LOD of the dynamic input data at the
!-- boundaries is read.
#if defined ( __netcdf )
!
!-- Open file in read-only mode
    CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), pids_id )
!
!-- Read attributes for LOD. In order to gurantee that also older drivers, where attribute is not
!-- given, are working, do not abort the run but assume LOD2 forcing.
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_east_pt,  .FALSE., 'ls_forcing_left_pt', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_east_qv,  .FALSE., 'ls_forcing_left_qv', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_east_u,   .FALSE., 'ls_forcing_left_u',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_east_v,   .FALSE., 'ls_forcing_left_v',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_east_w,   .FALSE., 'ls_forcing_left_w',  .FALSE. )

    CALL get_attribute( pids_id, char_lod, nest_offl%lod_north_pt, .FALSE., 'ls_forcing_north_pt', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_north_qv, .FALSE., 'ls_forcing_north_qv', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_north_u,  .FALSE., 'ls_forcing_north_u',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_north_v,  .FALSE., 'ls_forcing_north_v',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_north_w,  .FALSE., 'ls_forcing_north_w',  .FALSE. )

    CALL get_attribute( pids_id, char_lod, nest_offl%lod_south_pt, .FALSE., 'ls_forcing_south_pt', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_south_qv, .FALSE., 'ls_forcing_south_qv', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_south_u,  .FALSE., 'ls_forcing_south_u',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_south_v,  .FALSE., 'ls_forcing_south_v',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_south_w,  .FALSE., 'ls_forcing_south_w',  .FALSE. )

    CALL get_attribute( pids_id, char_lod, nest_offl%lod_west_pt,  .FALSE., 'ls_forcing_right_pt', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_west_qv,  .FALSE., 'ls_forcing_right_qv', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_west_u,   .FALSE., 'ls_forcing_right_u',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_west_v,   .FALSE., 'ls_forcing_right_v',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_west_w,   .FALSE., 'ls_forcing_right_w',  .FALSE. )

    CALL get_attribute( pids_id, char_lod, nest_offl%lod_top_pt,   .FALSE., 'ls_forcing_top_pt', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_top_qv,   .FALSE., 'ls_forcing_top_qv', .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_top_u,    .FALSE., 'ls_forcing_top_u',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_top_v,    .FALSE., 'ls_forcing_top_v',  .FALSE. )
    CALL get_attribute( pids_id, char_lod, nest_offl%lod_top_w,    .FALSE., 'ls_forcing_top_w',  .FALSE. )

    CALL close_input_file( pids_id )
#endif
!
!-- Temporary workaround until most of the dynamic drivers contain a LOD attribute. So far INIFOR
!-- did not provide the LOD attribute. In order to still use these older dynamic drivers, provide a
!-- temporary workaround. If the LOD is not given, a NetCDF interal error will occur but the
!-- simulation will not be aborted since the no_abort flag is passed. However, the respective
!-- attribute value might be given an arbitrary number. Hence, check for valid LOD's and manually
!-- set them to LOD 2 (as assumed so far). Note, this workaround should be removed later (date of
!-- reference: 6. Oct. 2020).
    IF ( nest_offl%lod_east_pt /= 1  .AND.  nest_offl%lod_east_pt /= 2 )  nest_offl%lod_east_pt = 2
    IF ( nest_offl%lod_east_qv /= 1  .AND.  nest_offl%lod_east_qv /= 2 )  nest_offl%lod_east_qv = 2
    IF ( nest_offl%lod_east_u  /= 1  .AND.  nest_offl%lod_east_u  /= 2 )  nest_offl%lod_east_u  = 2
    IF ( nest_offl%lod_east_v  /= 1  .AND.  nest_offl%lod_east_v  /= 2 )  nest_offl%lod_east_v  = 2
    IF ( nest_offl%lod_east_w  /= 1  .AND.  nest_offl%lod_east_w  /= 2 )  nest_offl%lod_east_w  = 2

    IF ( nest_offl%lod_north_pt /= 1  .AND.  nest_offl%lod_north_pt /= 2 )  nest_offl%lod_north_pt = 2
    IF ( nest_offl%lod_north_qv /= 1  .AND.  nest_offl%lod_north_qv /= 2 )  nest_offl%lod_north_qv = 2
    IF ( nest_offl%lod_north_u  /= 1  .AND.  nest_offl%lod_north_u  /= 2 )  nest_offl%lod_north_u  = 2
    IF ( nest_offl%lod_north_v  /= 1  .AND.  nest_offl%lod_north_v  /= 2 )  nest_offl%lod_north_v  = 2
    IF ( nest_offl%lod_north_w  /= 1  .AND.  nest_offl%lod_north_w  /= 2 )  nest_offl%lod_north_w  = 2

    IF ( nest_offl%lod_south_pt /= 1  .AND.  nest_offl%lod_south_pt /= 2 )  nest_offl%lod_south_pt = 2
    IF ( nest_offl%lod_south_qv /= 1  .AND.  nest_offl%lod_south_qv /= 2 )  nest_offl%lod_south_qv = 2
    IF ( nest_offl%lod_south_u  /= 1  .AND.  nest_offl%lod_south_u  /= 2 )  nest_offl%lod_south_u  = 2
    IF ( nest_offl%lod_south_v  /= 1  .AND.  nest_offl%lod_south_v  /= 2 )  nest_offl%lod_south_v  = 2
    IF ( nest_offl%lod_south_w  /= 1  .AND.  nest_offl%lod_south_w  /= 2 )  nest_offl%lod_south_w  = 2

    IF ( nest_offl%lod_west_pt /= 1  .AND.  nest_offl%lod_west_pt /= 2 )  nest_offl%lod_west_pt = 2
    IF ( nest_offl%lod_west_qv /= 1  .AND.  nest_offl%lod_west_qv /= 2 )  nest_offl%lod_west_qv = 2
    IF ( nest_offl%lod_west_u  /= 1  .AND.  nest_offl%lod_west_u  /= 2 )  nest_offl%lod_west_u  = 2
    IF ( nest_offl%lod_west_v  /= 1  .AND.  nest_offl%lod_west_v  /= 2 )  nest_offl%lod_west_v  = 2
    IF ( nest_offl%lod_west_w  /= 1  .AND.  nest_offl%lod_west_w  /= 2 )  nest_offl%lod_west_w  = 2

    IF ( nest_offl%lod_top_pt /= 1  .AND.  nest_offl%lod_top_pt /= 2 )  nest_offl%lod_top_pt = 2
    IF ( nest_offl%lod_top_qv /= 1  .AND.  nest_offl%lod_top_qv /= 2 )  nest_offl%lod_top_qv = 2
    IF ( nest_offl%lod_top_u  /= 1  .AND.  nest_offl%lod_top_u  /= 2 )  nest_offl%lod_top_u  = 2
    IF ( nest_offl%lod_top_v  /= 1  .AND.  nest_offl%lod_top_v  /= 2 )  nest_offl%lod_top_v  = 2
    IF ( nest_offl%lod_top_w  /= 1  .AND.  nest_offl%lod_top_w  /= 2 )  nest_offl%lod_top_w  = 2
!
!-- For consistency, check if all boundary input variables have the same LOD.
    IF ( MAX( nest_offl%lod_east_pt,  nest_offl%lod_east_qv,  nest_offl%lod_east_u,                &
              nest_offl%lod_east_v,   nest_offl%lod_east_w,                                        &
              nest_offl%lod_north_pt, nest_offl%lod_north_qv, nest_offl%lod_north_u,               &
              nest_offl%lod_north_v,  nest_offl%lod_north_w,                                       &
              nest_offl%lod_south_pt, nest_offl%lod_south_qv, nest_offl%lod_south_u,               &
              nest_offl%lod_south_v,  nest_offl%lod_south_w,                                       &
              nest_offl%lod_north_pt, nest_offl%lod_north_qv, nest_offl%lod_north_u,               &
              nest_offl%lod_north_v,  nest_offl%lod_north_w,                                       &
              nest_offl%lod_top_pt,   nest_offl%lod_top_qv,   nest_offl%lod_top_u,                 &
              nest_offl%lod_top_v,    nest_offl%lod_top_w )                                        &
            /=                                                                                     &
         MIN( nest_offl%lod_east_pt,  nest_offl%lod_east_qv,  nest_offl%lod_east_u,                &
              nest_offl%lod_east_v,   nest_offl%lod_east_w,                                        &
              nest_offl%lod_north_pt, nest_offl%lod_north_qv, nest_offl%lod_north_u,               &
              nest_offl%lod_north_v,  nest_offl%lod_north_w,                                       &
              nest_offl%lod_south_pt, nest_offl%lod_south_qv, nest_offl%lod_south_u,               &
              nest_offl%lod_south_v,  nest_offl%lod_south_w,                                       &
              nest_offl%lod_north_pt, nest_offl%lod_north_qv, nest_offl%lod_north_u,               &
              nest_offl%lod_north_v,  nest_offl%lod_north_w,                                       &
              nest_offl%lod_top_pt,   nest_offl%lod_top_qv,   nest_offl%lod_top_u,                 &
              nest_offl%lod_top_v,    nest_offl%lod_top_w ) )  THEN
       message_string = 'A mixture of different LOD for the provided boundary data is not ' //     &
                        'possible.'
       CALL message( 'nesting_offl_init', 'PA0504', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- As all LODs are the same, store it.
    lod = nest_offl%lod_east_u
!
!-- Set index range according to the given LOD in order to allocate the input arrays.
    IF ( bc_dirichlet_l  .OR.  bc_dirichlet_r  )  THEN
       IF ( lod == 2 )  THEN
          j_start   = nys
          j_start_v = nysv
          j_end     = nyn
       ELSE
          j_start   = 1
          j_start_v = 1
          j_end     = 1
       ENDIF
    ENDIF

    IF ( bc_dirichlet_n  .OR.  bc_dirichlet_s )  THEN
       IF( lod == 2 )  THEN
          i_start   = nxl
          i_start_u = nxlu
          i_end     = nxr
       ELSE
          i_start   = 1
          i_start_u = 1
          i_end     = 1
       ENDIF
    ENDIF
!
!-- Allocate arrays for reading left/right boundary values. Arrays will incorporate 2 time levels in
!-- order to interpolate in between. Depending on the given LOD, the x-, or y-dimension will be
!-- either nxl:nxr, or nys:nyn (for LOD=2), or it reduces to one element for LOD=1. If the core has
!-- no lateral boundary, allocate a dummy array as well, in order to enable netcdf parallel access.
!-- Dummy arrays will be allocated with dimension length zero.
    IF ( bc_dirichlet_l )  THEN
       ALLOCATE( nest_offl%u_l(0:1,nzb+1:nzt,j_start:j_end)  )
       ALLOCATE( nest_offl%v_l(0:1,nzb+1:nzt,j_start_v:j_end) )
       ALLOCATE( nest_offl%w_l(0:1,nzb+1:nzt-1,j_start:j_end) )
       IF ( humidity )       ALLOCATE( nest_offl%q_l(0:1,nzb+1:nzt,j_start:j_end)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_l(0:1,nzb+1:nzt,j_start:j_end) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_l(0:1,nzb+1:nzt,j_start:j_end,1:UBOUND( chem_species, 1 )) )
    ELSE
       ALLOCATE( nest_offl%u_l(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%v_l(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%w_l(1:1,1:1,1:1)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_l(1:1,1:1,1:1)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_l(1:1,1:1,1:1)  )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_l(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
    ENDIF
    IF ( bc_dirichlet_r )  THEN
       ALLOCATE( nest_offl%u_r(0:1,nzb+1:nzt,j_start:j_end)  )
       ALLOCATE( nest_offl%v_r(0:1,nzb+1:nzt,j_start_v:j_end) )
       ALLOCATE( nest_offl%w_r(0:1,nzb+1:nzt-1,j_start:j_end) )
       IF ( humidity )       ALLOCATE( nest_offl%q_r(0:1,nzb+1:nzt,j_start:j_end)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_r(0:1,nzb+1:nzt,j_start:j_end) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_r(0:1,nzb+1:nzt,j_start:j_end,1:UBOUND( chem_species, 1 )) )
    ELSE
       ALLOCATE( nest_offl%u_r(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%v_r(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%w_r(1:1,1:1,1:1)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_r(1:1,1:1,1:1)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_r(1:1,1:1,1:1)  )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_r(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
    ENDIF
!
!-- Allocate arrays for reading north/south boundary values. Arrays will incorporate 2 time levels
!-- in order to interpolate in between. If the core has no boundary, allocate a dummy array, in
!-- order to enable netcdf parallel access. Dummy arrays will be allocated with dimension length
!-- zero.
    IF ( bc_dirichlet_n )  THEN
       ALLOCATE( nest_offl%u_n(0:1,nzb+1:nzt,i_start_u:i_end) )
       ALLOCATE( nest_offl%v_n(0:1,nzb+1:nzt,i_start:i_end)  )
       ALLOCATE( nest_offl%w_n(0:1,nzb+1:nzt-1,i_start:i_end) )
       IF ( humidity )       ALLOCATE( nest_offl%q_n(0:1,nzb+1:nzt,i_start:i_end)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_n(0:1,nzb+1:nzt,i_start:i_end) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_n(0:1,nzb+1:nzt,i_start:i_end,1:UBOUND( chem_species, 1 )) )
    ELSE
       ALLOCATE( nest_offl%u_n(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%v_n(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%w_n(1:1,1:1,1:1)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_n(1:1,1:1,1:1)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_n(1:1,1:1,1:1)  )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_n(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
    ENDIF
    IF ( bc_dirichlet_s )  THEN
       ALLOCATE( nest_offl%u_s(0:1,nzb+1:nzt,i_start_u:i_end) )
       ALLOCATE( nest_offl%v_s(0:1,nzb+1:nzt,i_start:i_end)  )
       ALLOCATE( nest_offl%w_s(0:1,nzb+1:nzt-1,i_start:i_end) )
       IF ( humidity )       ALLOCATE( nest_offl%q_s(0:1,nzb+1:nzt,i_start:i_end)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_s(0:1,nzb+1:nzt,i_start:i_end) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_s(0:1,nzb+1:nzt,i_start:i_end,1:UBOUND( chem_species, 1 )) )
    ELSE
       ALLOCATE( nest_offl%u_s(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%v_s(1:1,1:1,1:1)  )
       ALLOCATE( nest_offl%w_s(1:1,1:1,1:1)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_s(1:1,1:1,1:1)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_s(1:1,1:1,1:1)  )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_s(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
    ENDIF
!
!-- Allocate arrays for reading data at the top boundary. In contrast to the lateral boundaries,
!-- each core reads these data so that no dummy arrays need to be allocated.
    IF ( lod == 2 )  THEN
       ALLOCATE( nest_offl%u_top(0:1,nys:nyn,nxlu:nxr) )
       ALLOCATE( nest_offl%v_top(0:1,nysv:nyn,nxl:nxr) )
       ALLOCATE( nest_offl%w_top(0:1,nys:nyn,nxl:nxr)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_top(0:1,nys:nyn,nxl:nxr)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_top(0:1,nys:nyn,nxl:nxr) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_top(0:1,nys:nyn,nxl:nxr,1:UBOUND( chem_species, 1 )) )
    ELSE
       ALLOCATE( nest_offl%u_top(0:1,1:1,1:1) )
       ALLOCATE( nest_offl%v_top(0:1,1:1,1:1) )
       ALLOCATE( nest_offl%w_top(0:1,1:1,1:1)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_top(0:1,1:1,1:1)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_top(0:1,1:1,1:1) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_top(0:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
    ENDIF
!
!-- For chemical species, create the names of the variables. This is necessary to identify the
!-- respective variable and write it onto the correct array in the chem_species datatype.
    IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
       ALLOCATE( nest_offl%chem_from_file_l(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%chem_from_file_n(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%chem_from_file_r(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%chem_from_file_s(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%chem_from_file_t(1:UBOUND( chem_species, 1 )) )

       ALLOCATE( nest_offl%var_names_chem_l(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%var_names_chem_n(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%var_names_chem_r(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%var_names_chem_s(1:UBOUND( chem_species, 1 )) )
       ALLOCATE( nest_offl%var_names_chem_t(1:UBOUND( chem_species, 1 )) )
!
!--    Initialize flags that indicate whether the variable is on file or not. Please note, this is
!--    only necessary for chemistry variables.
       nest_offl%chem_from_file_l(:) = .FALSE.
       nest_offl%chem_from_file_n(:) = .FALSE.
       nest_offl%chem_from_file_r(:) = .FALSE.
       nest_offl%chem_from_file_s(:) = .FALSE.
       nest_offl%chem_from_file_t(:) = .FALSE.

       DO  n = 1, UBOUND( chem_species, 1 )
          nest_offl%var_names_chem_l(n) = nest_offl%char_l // TRIM(chem_species(n)%name)
          nest_offl%var_names_chem_n(n) = nest_offl%char_n // TRIM(chem_species(n)%name)
          nest_offl%var_names_chem_r(n) = nest_offl%char_r // TRIM(chem_species(n)%name)
          nest_offl%var_names_chem_s(n) = nest_offl%char_s // TRIM(chem_species(n)%name)
          nest_offl%var_names_chem_t(n) = nest_offl%char_t // TRIM(chem_species(n)%name)
       ENDDO
    ENDIF
!
!-- Offline nesting for salsa
    IF ( salsa )  CALL salsa_nesting_offl_init
!
!-- Before initial data input is initiated, check if dynamic input file is present.
    IF ( .NOT. input_pids_dynamic )  THEN
       message_string = 'nesting_offline = .TRUE. requires dynamic '  //                           &
                         'input file ' // TRIM( input_file_dynamic ) // TRIM( coupling_char )
       CALL message( 'nesting_offl_init', 'PA0546', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read COSMO data at lateral and top boundaries
    CALL nesting_offl_input
!
!-- Check if sufficient time steps are provided to cover the entire simulation. Note, dynamic input
!-- is only required for the 3D simulation, not for the soil/wall spinup. However, as the spinup
!-- time is added to the end_time, this must be considered here.
    IF ( end_time - spinup_time > nest_offl%time(nest_offl%nt-1) )  THEN
       message_string = 'end_time of the simulation exceeds the ' //                               &
                        'time dimension in the dynamic input file.'
       CALL message( 'nesting_offl_init', 'PA0183', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set indicies for boundary grid points
    IF ( bc_dirichlet_l  .OR.  bc_dirichlet_r )  THEN
       i_bound   = MERGE( nxl  - 1, nxr + 1, bc_dirichlet_l )
       i_bound_u = MERGE( nxlu - 1, nxr + 1, bc_dirichlet_l )
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_dirichlet_s )  THEN
       j_bound   = MERGE( nys  - 1, nyn + 1, bc_dirichlet_s )
       j_bound_v = MERGE( nysv - 1, nyn + 1, bc_dirichlet_s )
    ENDIF
!
!-- Initialize boundary data. Please note, do not initialize boundaries in case of restart runs.
!-- This case the boundaries are already initialized and the boundary data from file would be on the
!-- wrong time level.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
!
!--    Distinguish between LOD = 1 and LOD = 2 inititialization
       IF ( lod == 2 )  THEN
          IF ( bc_dirichlet_l )  THEN
             u(nzb+1:nzt,nys:nyn,i_bound_u) = nest_offl%u_l(0,nzb+1:nzt,nys:nyn)
             v(nzb+1:nzt,nysv:nyn,i_bound)  = nest_offl%v_l(0,nzb+1:nzt,nysv:nyn)
             w(nzb+1:nzt-1,nys:nyn,i_bound) = nest_offl%w_l(0,nzb+1:nzt-1,nys:nyn)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,nys:nyn,i_bound) = nest_offl%pt_l(0,nzb+1:nzt,nys:nyn)
             IF ( humidity      )  q(nzb+1:nzt,nys:nyn,i_bound)  = nest_offl%q_l(0,nzb+1:nzt,nys:nyn)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_l(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,nys:nyn,i_bound) =                            &
                                                             nest_offl%chem_l(0,nzb+1:nzt,nys:nyn,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_r )  THEN
             u(nzb+1:nzt,nys:nyn,i_bound_u) = nest_offl%u_r(0,nzb+1:nzt,nys:nyn)
             v(nzb+1:nzt,nysv:nyn,i_bound)  = nest_offl%v_r(0,nzb+1:nzt,nysv:nyn)
             w(nzb+1:nzt-1,nys:nyn,i_bound) = nest_offl%w_r(0,nzb+1:nzt-1,nys:nyn)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,nys:nyn,i_bound) = nest_offl%pt_r(0,nzb+1:nzt,nys:nyn)
             IF ( humidity      )  q(nzb+1:nzt,nys:nyn,i_bound)  = nest_offl%q_r(0,nzb+1:nzt,nys:nyn)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_r(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,nys:nyn,i_bound) =                            &
                                                             nest_offl%chem_r(0,nzb+1:nzt,nys:nyn,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF

          IF ( bc_dirichlet_n)  THEN
             u(nzb+1:nzt,j_bound,nxlu:nxr)  = nest_offl%u_n(0,nzb+1:nzt,nxlu:nxr)
             v(nzb+1:nzt,j_bound_v,nxl:nxr) = nest_offl%v_n(0,nzb+1:nzt,nxl:nxr)
             w(nzb+1:nzt-1,j_bound,nxl:nxr) = nest_offl%w_n(0,nzb+1:nzt-1,nxl:nxr)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,j_bound,nxl:nxr) = nest_offl%pt_n(0,nzb+1:nzt,nxl:nxr)
             IF ( humidity      )  q(nzb+1:nzt,j_bound,nxl:nxr)  = nest_offl%q_n(0,nzb+1:nzt,nxl:nxr)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_n(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,j_bound,nxl:nxr) =                            &
                                                             nest_offl%chem_n(0,nzb+1:nzt,nxl:nxr,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_s)  THEN
             u(nzb+1:nzt,j_bound,nxlu:nxr)  = nest_offl%u_s(0,nzb+1:nzt,nxlu:nxr)
             v(nzb+1:nzt,j_bound_v,nxl:nxr) = nest_offl%v_s(0,nzb+1:nzt,nxl:nxr)
             w(nzb+1:nzt-1,j_bound,nxl:nxr) = nest_offl%w_s(0,nzb+1:nzt-1,nxl:nxr)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,j_bound,nxl:nxr) = nest_offl%pt_s(0,nzb+1:nzt,nxl:nxr)
             IF ( humidity      )  q(nzb+1:nzt,j_bound,nxl:nxr)  = nest_offl%q_s(0,nzb+1:nzt,nxl:nxr)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_s(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,j_bound,nxl:nxr) =                            &
                                                             nest_offl%chem_s(0,nzb+1:nzt,nxl:nxr,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF

          u(nzt+1,nys:nyn,nxlu:nxr) = nest_offl%u_top(0,nys:nyn,nxlu:nxr)
          v(nzt+1,nysv:nyn,nxl:nxr) = nest_offl%v_top(0,nysv:nyn,nxl:nxr)
          w(nzt,nys:nyn,nxl:nxr)    = nest_offl%w_top(0,nys:nyn,nxl:nxr)
          w(nzt+1,nys:nyn,nxl:nxr)  = nest_offl%w_top(0,nys:nyn,nxl:nxr)
          IF ( .NOT. neutral )  pt(nzt+1,nys:nyn,nxl:nxr) = nest_offl%pt_top(0,nys:nyn,nxl:nxr)
          IF ( humidity )       q(nzt+1,nys:nyn,nxl:nxr)  = nest_offl%q_top(0,nys:nyn,nxl:nxr)
          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( chem_species, 1 )
                IF( nest_offl%chem_from_file_t(n) )  THEN
                   chem_species(n)%conc(nzt+1,nys:nyn,nxl:nxr) =                                   &
                                                             nest_offl%chem_top(0,nys:nyn,nxl:nxr,n)
                ENDIF
             ENDDO
          ENDIF
!
!--    LOD 1
       ELSE
          IF ( bc_dirichlet_l )  THEN
             DO  j = nys, nyn
                u(nzb+1:nzt,j,i_bound_u) = nest_offl%u_l(0,nzb+1:nzt,1)
                w(nzb+1:nzt-1,j,i_bound) = nest_offl%w_l(0,nzb+1:nzt-1,1)
             ENDDO
             DO  j = nysv, nyn
                v(nzb+1:nzt,j,i_bound)  = nest_offl%v_l(0,nzb+1:nzt,1)
             ENDDO
             IF ( .NOT. neutral )  THEN
                DO  j = nys, nyn
                   pt(nzb+1:nzt,j,i_bound) = nest_offl%pt_l(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( humidity      )  THEN
                DO  j = nys, nyn
                   q(nzb+1:nzt,j,i_bound)  = nest_offl%q_l(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_l(n) )  THEN
                      DO  j = nys, nyn
                         chem_species(n)%conc(nzb+1:nzt,j,i_bound) =                               &
                                                                   nest_offl%chem_l(0,nzb+1:nzt,1,n)
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_r )  THEN
             DO  j = nys, nyn
                u(nzb+1:nzt,j,i_bound_u) = nest_offl%u_r(0,nzb+1:nzt,1)
                w(nzb+1:nzt-1,j,i_bound) = nest_offl%w_r(0,nzb+1:nzt-1,1)
             ENDDO
             DO  j = nysv, nyn
                v(nzb+1:nzt,j,i_bound)  = nest_offl%v_r(0,nzb+1:nzt,1)
             ENDDO
             IF ( .NOT. neutral )  THEN
                DO  j = nys, nyn
                   pt(nzb+1:nzt,j,i_bound) = nest_offl%pt_r(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( humidity      )  THEN
                DO  j = nys, nyn
                   q(nzb+1:nzt,j,i_bound)  = nest_offl%q_r(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_r(n) )  THEN
                      DO  j = nys, nyn
                         chem_species(n)%conc(nzb+1:nzt,j,i_bound) =                               &
                                                                   nest_offl%chem_r(0,nzb+1:nzt,1,n)
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_n )  THEN
             DO  i = nxlu, nxr
                u(nzb+1:nzt,j_bound,i)  = nest_offl%u_n(0,nzb+1:nzt,1)
             ENDDO
             DO  i = nxl, nxr
                v(nzb+1:nzt,j_bound_v,i) = nest_offl%v_n(0,nzb+1:nzt,1)
                w(nzb+1:nzt-1,j_bound,i) = nest_offl%w_n(0,nzb+1:nzt-1,1)
             ENDDO
             IF ( .NOT. neutral )  THEN
                DO  i = nxl, nxr
                   pt(nzb+1:nzt,j_bound,i) = nest_offl%pt_n(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( humidity      )  THEN
                DO  i = nxl, nxr
                   q(nzb+1:nzt,j_bound,i)  = nest_offl%q_n(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_n(n) )  THEN
                      DO  i = nxl, nxr
                         chem_species(n)%conc(nzb+1:nzt,j_bound,i) =                               &
                                                                   nest_offl%chem_n(0,nzb+1:nzt,1,n)
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_s )  THEN
             DO  i = nxlu, nxr
                u(nzb+1:nzt,j_bound,i)  = nest_offl%u_s(0,nzb+1:nzt,1)
             ENDDO
             DO  i = nxl, nxr
                v(nzb+1:nzt,j_bound_v,i) = nest_offl%v_s(0,nzb+1:nzt,1)
                w(nzb+1:nzt-1,j_bound,i) = nest_offl%w_s(0,nzb+1:nzt-1,1)
             ENDDO
             IF ( .NOT. neutral )  THEN
                DO  i = nxl, nxr
                   pt(nzb+1:nzt,j_bound,i) = nest_offl%pt_s(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( humidity      )  THEN
                DO  i = nxl, nxr
                   q(nzb+1:nzt,j_bound,i)  = nest_offl%q_s(0,nzb+1:nzt,1)
                ENDDO
             ENDIF
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_s(n) )  THEN
                      DO  i = nxl, nxr
                         chem_species(n)%conc(nzb+1:nzt,j_bound,i) =                               &
                                                                   nest_offl%chem_s(0,nzb+1:nzt,1,n)
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF
          ENDIF

          u(nzt+1,nys:nyn,nxlu:nxr) = nest_offl%u_top(0,1,1)
          v(nzt+1,nysv:nyn,nxl:nxr) = nest_offl%v_top(0,1,1)
          w(nzt,nys:nyn,nxl:nxr)    = nest_offl%w_top(0,1,1)
          w(nzt+1,nys:nyn,nxl:nxr)  = nest_offl%w_top(0,1,1)
          IF ( .NOT. neutral )  pt(nzt+1,nys:nyn,nxl:nxr) = nest_offl%pt_top(0,1,1)
          IF ( humidity )       q(nzt+1,nys:nyn,nxl:nxr)  = nest_offl%q_top(0,1,1)
          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( chem_species, 1 )
                IF( nest_offl%chem_from_file_t(n) )  THEN
                   chem_species(n)%conc(nzt+1,nys:nyn,nxl:nxr) = nest_offl%chem_top(0,1,1,n)
                ENDIF
             ENDDO
          ENDIF
       ENDIF
!
!--    In case of offline nesting the pressure forms itself based on the prescribed lateral
!--    boundary conditions. Hence, explicit forcing by pressure gradients via geostrophic wind
!--    components is not necessary and would be canceled out by the perturbation pressure otherwise.
!--    For this reason, set geostrophic wind components to zero.
       ug(nzb+1:nzt) = 0.0_wp
       vg(nzb+1:nzt) = 0.0_wp

    ENDIF
!
!-- After boundary data is initialized, mask topography at the boundaries for the velocity
!-- components.
    u = MERGE( u, 0.0_wp, BTEST( topo_flags, 1 ) )
    v = MERGE( v, 0.0_wp, BTEST( topo_flags, 2 ) )
    w = MERGE( w, 0.0_wp, BTEST( topo_flags, 3 ) )
!
!-- Initial calculation of the boundary layer depth from the prescribed boundary data. This is
!-- required for initialize the synthetic turbulence generator correctly.
    CALL nesting_offl_calc_zi
!
!-- After boundary data is initialized, ensure mass conservation. Not necessary in restart runs.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
       CALL nesting_offl_mass_conservation
    ENDIF

 END SUBROUTINE nesting_offl_init

!--------------------------------------------------------------------------------------------------!
! Description:
!--------------------------------------------------------------------------------------------------!
!> Interpolation function, used to interpolate boundary data in time.
!--------------------------------------------------------------------------------------------------!
 FUNCTION interpolate_in_time( var_t1, var_t2, fac  )

    REAL(wp)            :: fac                  !< interpolation factor
    REAL(wp)            :: interpolate_in_time  !< time-interpolated boundary value
    REAL(wp)            :: var_t1               !< boundary value at t1
    REAL(wp)            :: var_t2               !< boundary value at t2

    interpolate_in_time = ( 1.0_wp - fac ) * var_t1 + fac * var_t2

 END FUNCTION interpolate_in_time



 END MODULE nesting_offl_mod
