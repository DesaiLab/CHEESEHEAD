!> @synthetic_turbulence_generator_mod.f90
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
!
! Authors:
! --------
! @author Tobias Gronemeier, Matthias Suehring, Atsushi Inagaki, Micha Gryschka, Christoph Knigge
!
!
! Description:
! ------------
!> The module generates turbulence at the inflow boundary based on a method by Xie and Castro (2008)
!> utilizing a Lund rotation (Lund, 1998) and a mass-flux correction by Kim et al. (2013).
!> The turbulence is correlated based on length scales in y- and z-direction and a time scale for
!> each velocity component. The profiles of length and time scales, mean u, v, w, e and pt, and all
!> components of the Reynolds stress tensor can be either read from file STG_PROFILES, or will be
!> parametrized within the boundary layer.
!>
!> @todo Enable cyclic_fill
!>       Implement turbulence generation for e and pt
!> @note <Enter notes on the module>
!> @bug  Height information from input file is not used. Profiles from input must match with current
!>       PALM grid.
!>       In case of restart, velocity seeds differ from precursor run if a11, a22, or a33 are zero.
!--------------------------------------------------------------------------------------------------!
 MODULE synthetic_turbulence_generator_mod


    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               ddzw,                                                                               &
               drho_air,                                                                           &
               mean_inflow_profiles,                                                               &
               pt,                                                                                 &
               pt_init,                                                                            &
               q,                                                                                  &
               q_init,                                                                             &
               u,                                                                                  &
               u_init,                                                                             &
               v,                                                                                  &
               v_init,                                                                             &
               w,                                                                                  &
               zu,                                                                                 &
               zw,                                                                                 &
               heatflux_output_conversion

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  g,                                                                                  &
               kappa,                                                                              &
               pi

    USE control_parameters,                                                                        &
        ONLY:  bc_lr,                                                                              &
               bc_ns,                                                                              &
               child_domain,                                                                       &
               coupling_char,                                                                      &
               debug_output_timestep,                                                              &
               dt_3d,                                                                              &
               e_init,                                                                             &
               humidity,                                                                           &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               length,                                                                             &
               message_string,                                                                     &
               nesting_offline,                                                                    &
               neutral,                                                                            &
               num_mean_inflow_profiles,                                                           &
               random_generator,                                                                   &
               rans_mode,                                                                          &
               restart_data_format_output,                                                         &
               restart_string,                                                                     &
               syn_turb_gen,                                                                       &
               time_since_reference_point,                                                         &
               turbulent_inflow

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlu,                                                                               &
               nxr,                                                                                &
               ny,                                                                                 &
               nys,                                                                                &
               nysv,                                                                               &
               nyn,                                                                                &
               topo_flags,                                                                         &
               topo_min_level

    USE kinds

#if defined( __parallel )
    USE MPI
#endif

    USE nesting_offl_mod,                                                                          &
        ONLY:  nesting_offl_calc_zi,                                                               &
               zi_ribulk

    USE pegrid,                                                                                    &
        ONLY:  comm1dx,                                                                            &
               comm1dy,                                                                            &
               comm2d,                                                                             &
               ierr,                                                                               &
               myidx,                                                                              &
               myidy,                                                                              &
               npex,                                                                               &
               npey

    USE pmc_interface,                                                                             &
        ONLY : rans_mode_parent

    USE random_generator_parallel,                                                                 &
        ONLY:  init_parallel_random_generator,                                                     &
               random_dummy,                                                                       &
               random_number_parallel,                                                             &
               random_seed_parallel

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rrd_mpi_io,                                                                         &
               wrd_mpi_io

    USE transpose_mod,                                                                             &
        ONLY:  nzb_x,                                                                              &
               nzt_x

    USE statistics

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h,                                                                         &
               surf_lsm_h,                                                                         &
               surf_usm_h

    IMPLICIT NONE

    INTEGER(iwp) ::  id_stg_left        !< left lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_north       !< north lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_right       !< right lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_south       !< south lateral boundary core id in case of turbulence generator
    INTEGER(iwp) ::  k_zi               !< vertical index of boundary-layer top
    INTEGER(iwp) ::  mergp_limit = 1000 !< maximum length scale (in gp)
    INTEGER(iwp) ::  mergp_x            !< maximum length scale in x (in gp)
    INTEGER(iwp) ::  mergp_xy           !< maximum horizontal length scale (in gp)
    INTEGER(iwp) ::  mergp_y            !< maximum length scale in y (in gp)
    INTEGER(iwp) ::  mergp_z            !< maximum length scale in z (in gp)
    INTEGER(iwp) ::  nzb_x_stg          !< lower bound of z coordinate (required for transposing z on PEs along x)
    INTEGER(iwp) ::  nzt_x_stg          !< upper bound of z coordinate (required for transposing z on PEs along x)
    INTEGER(iwp) ::  nzb_y_stg          !< lower bound of z coordinate (required for transposing z on PEs along y)
    INTEGER(iwp) ::  nzt_y_stg          !< upper bound of z coordinate (required for transposing z on PEs along y)
#if defined( __parallel )
    INTEGER(iwp) ::  stg_type_xz        !< MPI type for full z range
    INTEGER(iwp) ::  stg_type_xz_small  !< MPI type for small z range
    INTEGER(iwp) ::  stg_type_yz        !< MPI type for full z range
    INTEGER(iwp) ::  stg_type_yz_small  !< MPI type for small z range
#endif

    INTEGER(iwp), DIMENSION(3) ::  nr_non_topo_xz = 0  !< number of non-topography grid points at xz cross-sections,
                                                       !< required for bias correction of imposed perturbations
    INTEGER(iwp), DIMENSION(3) ::  nr_non_topo_yz = 0  !< number of non-topography grid points at yz cross-sections,
                                                       !< required for bias correction of imposed perturbations

#if defined( __parallel )
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  displs_xz      !< displacement for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  recv_count_xz  !< receive count for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  displs_yz      !< displacement for MPI_GATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  recv_count_yz  !< receive count for MPI_GATHERV
#endif
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nux            !< length scale of u in x direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nuy            !< length scale of u in y direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nuz            !< length scale of u in z direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nvx            !< length scale of v in x direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nvy            !< length scale of v in y direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nvz            !< length scale of v in z direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nwx            !< length scale of w in x direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nwy            !< length scale of w in y direction (in gp)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nwz            !< length scale of w in z direction (in gp)

    INTEGER(isp), DIMENSION(:), ALLOCATABLE   ::  id_rand_xz   !< initial random IDs at xz inflow boundary
    INTEGER(isp), DIMENSION(:), ALLOCATABLE   ::  id_rand_yz   !< initial random IDs at yz inflow boundary
    INTEGER(isp), DIMENSION(:,:), ALLOCATABLE ::  seq_rand_xz  !< initial random seeds at xz inflow boundary
    INTEGER(isp), DIMENSION(:,:), ALLOCATABLE ::  seq_rand_yz  !< initial random seeds at yz inflow boundary


    LOGICAL ::  adjustment_step               = .FALSE.  !< control flag indicating that time and lenght scales have been updated and
                                                         !< no time correlation to the timestep before should be considered
    LOGICAL ::  compute_velocity_seeds_local  = .FALSE.   !< switch to decide whether velocity seeds are computed locally
                                                         !< or if computation is distributed over several processes
    LOGICAL ::  file_stg_exist                = .FALSE.  !< flag indicating whether parameter file for Reynolds stress and length scales exist
    LOGICAL ::  parametrize_inflow_turbulence = .FALSE.  !< flag indicating that inflow turbulence is either read from file
                                                         !< (.FALSE.) or if it parametrized
    LOGICAL ::  velocity_seed_initialized     = .FALSE.  !< true after first call of stg_main


    REAL(wp) ::  blend                     !< value to create gradually and smooth blending of Reynolds stress and length
                                           !< scales above the boundary layer
    REAL(wp) ::  blend_coeff = -9.3_wp     !< coefficient used to ensure that blending functions decreases to 1/10 after
                                           !< one length scale above ABL top
    REAL(wp) ::  d_l                       !< blend_coeff/length_scale
    REAL(wp) ::  d_nxy                     !< inverse of the total number of xy-grid points
    REAL(wp) ::  dt_stg_adjust = 1800.0_wp !< time interval for adjusting turbulence statistics
    REAL(wp) ::  dt_stg_call = 0.0_wp      !< time interval for calling synthetic turbulence generator
    REAL(wp) ::  scale_l                   !< scaling parameter used for turbulence parametrization - Obukhov length
    REAL(wp) ::  scale_us                  !< scaling parameter used for turbulence parametrization - friction velocity
    REAL(wp) ::  scale_wm                  !< scaling parameter used for turbulence parametrization - momentum scale
    REAL(wp) ::  time_stg_adjust = 0.0_wp  !< time counter for adjusting turbulence information
    REAL(wp) ::  time_stg_call = 0.0_wp    !< time counter for calling generator

    REAL(wp), DIMENSION(3) ::  mc_factor = 1.0_wp  !< correction factor for the u,v,w-components to maintain original mass flux


    REAL(wp),DIMENSION(:), ALLOCATABLE ::  r11  !< Reynolds parameter
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  r21  !< Reynolds parameter
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  r22  !< Reynolds parameter
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  r31  !< Reynolds parameter
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  r32  !< Reynolds parameter
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  r33  !< Reynolds parameter

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  a11  !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  a21  !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  a22  !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  a31  !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  a32  !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  a33  !< coefficient for Lund rotation
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tu   !< Lagrangian time scale of u
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tv   !< Lagrangian time scale of v
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tw   !< Lagrangian time scale of w

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bux     !< filter function for u in x direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  buy     !< filter function for u in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  buz     !< filter function for u in z direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bvx     !< filter function for v in x direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bvy     !< filter function for v in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bvz     !< filter function for v in z direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bwx     !< filter function for w in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bwy     !< filter function for w in y direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bwz     !< filter function for w in z direction
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fu_xz   !< velocity seed for u at xz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fuo_xz  !< velocity seed for u at xz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fu_yz   !< velocity seed for u at yz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fuo_yz  !< velocity seed for u at yz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fv_xz   !< velocity seed for v at xz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fvo_xz  !< velocity seed for v at xz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fv_yz   !< velocity seed for v at yz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fvo_yz  !< velocity seed for v at yz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fw_xz   !< velocity seed for w at xz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fwo_xz  !< velocity seed for w at xz plane with new random number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fw_yz   !< velocity seed for w at yz plane
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fwo_yz  !< velocity seed for w at yz plane with new random number

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dist_xz !< disturbances for parallel/crosswind/vertical component at north/south boundary
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dist_yz !< disturbances for parallel/crosswind/vertical component  at north/south boundary

!
!-- PALM interfaces:
!-- Adjust time and lenght scales, Reynolds stress, and filter functions
    INTERFACE stg_adjust
       MODULE PROCEDURE stg_adjust
    END INTERFACE stg_adjust
!
!-- Input parameter checks to be done in check_parameters
    INTERFACE stg_check_parameters
       MODULE PROCEDURE stg_check_parameters
    END INTERFACE stg_check_parameters

!
!-- Calculate filter functions
    INTERFACE stg_filter_func
       MODULE PROCEDURE stg_filter_func
    END INTERFACE stg_filter_func

!
!-- Generate velocity seeds at south and north domain boundary
    INTERFACE stg_generate_seed_xz
       MODULE PROCEDURE stg_generate_seed_xz
    END INTERFACE stg_generate_seed_xz
!
!-- Generate velocity seeds at left and/or right domain boundary
    INTERFACE stg_generate_seed_yz
       MODULE PROCEDURE stg_generate_seed_yz
    END INTERFACE stg_generate_seed_yz

!
!-- Output of information to the header file
    INTERFACE stg_header
       MODULE PROCEDURE stg_header
    END INTERFACE stg_header

!
!-- Initialization actions
    INTERFACE stg_init
       MODULE PROCEDURE stg_init
    END INTERFACE stg_init

!
!-- Main procedure of synth. turb. gen.
    INTERFACE stg_main
       MODULE PROCEDURE stg_main
    END INTERFACE stg_main

!
!-- Reading of NAMELIST parameters
    INTERFACE stg_parin
       MODULE PROCEDURE stg_parin
    END INTERFACE stg_parin

!
!-- Reading of parameters for restart runs
    INTERFACE stg_rrd_global
       MODULE PROCEDURE stg_rrd_global_ftn
       MODULE PROCEDURE stg_rrd_global_mpi
    END INTERFACE stg_rrd_global

!
!-- Writing of binary output for restart runs
    INTERFACE stg_wrd_global
       MODULE PROCEDURE stg_wrd_global
    END INTERFACE stg_wrd_global

    SAVE

    PRIVATE

!
!-- Public interfaces
    PUBLIC  stg_adjust,                                                                            &
            stg_check_parameters,                                                                  &
            stg_header,                                                                            &
            stg_init,                                                                              &
            stg_main,                                                                              &
            stg_parin,                                                                             &
            stg_rrd_global,                                                                        &
            stg_wrd_global

!
!-- Public variables
    PUBLIC  dt_stg_call,                                                                           &
            dt_stg_adjust,                                                                         &
            id_stg_left,                                                                           &
            id_stg_north,                                                                          &
            id_stg_right,                                                                          &
            id_stg_south,                                                                          &
            parametrize_inflow_turbulence,                                                         &
            time_stg_adjust,                                                                       &
            time_stg_call


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for synthetic turbulence generator
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_check_parameters

    IF ( .NOT. syn_turb_gen  .AND.  .NOT. rans_mode  .AND.                                         &
          nesting_offline )  THEN
       message_string = 'Synthetic turbulence generator is required ' //                           &
                        'if offline nesting is applied and PALM operates in LES mode.'
       CALL message( 'stg_check_parameters', 'PA0520', 0, 0, 0, 6, 0 )
    ENDIF

    IF ( .NOT. syn_turb_gen  .AND.  child_domain                                                   &
         .AND. rans_mode_parent  .AND.  .NOT. rans_mode )  THEN
       message_string = 'Synthetic turbulence generator is required when nesting is applied ' //   &
                        'and parent operates in RANS-mode but current child in LES mode.'
       CALL message( 'stg_check_parameters', 'PA0524', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( syn_turb_gen )  THEN

       IF ( child_domain  .AND.  .NOT. rans_mode  .AND.  .NOT. rans_mode_parent )  THEN
          message_string = 'Using synthetic turbulence generator is not allowed in LES-LES nesting.'
          CALL message( 'stg_check_parameters', 'PA0620', 1, 2, 0, 6, 0 )

       ENDIF

       IF ( child_domain  .AND.  rans_mode  .AND.  rans_mode_parent )  THEN
          message_string = 'Using synthetic turbulence generator is not allowed in RANS-RANS nesting.'
          CALL message( 'stg_check_parameters', 'PA0621', 1, 2, 0, 6, 0 )

       ENDIF

       IF ( .NOT. nesting_offline  .AND.  .NOT. child_domain )  THEN

          IF ( INDEX( initializing_actions, 'set_constant_profiles' ) == 0                         &
          .AND.  INDEX( initializing_actions, 'read_restart_data' ) == 0 )  THEN
             message_string = 'Using synthetic turbulence generator requires ' //                  &
                              '%initializing_actions = "set_constant_profiles" or ' //             &
                              ' "read_restart_data", if not offline nesting is applied.'
             CALL message( 'stg_check_parameters', 'PA0015', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( bc_lr /= 'dirichlet/radiation' )  THEN
             message_string = 'Using synthetic turbulence generator requires &bc_lr = ' //         &
                              ' "dirichlet/radiation", if not offline nesting is applied.'
             CALL message( 'stg_check_parameters', 'PA0035', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( bc_ns /= 'cyclic' )  THEN
             message_string = 'Using synthetic turbulence generator requires &bc_ns = ' //         &
                              ' "cyclic", if not offline nesting is applied.'
             CALL message( 'stg_check_parameters', 'PA0037', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       In case of idealized setups (i.e. no mesoscale nesting), a STG_PROFILES file must be
!--       available. Though turbulence parameters could also be parametrized in this case, this
!--       has not been tested thoroughly.
          INQUIRE( FILE = 'STG_PROFILES' // TRIM( coupling_char ), EXIST = file_stg_exist  )
          IF ( .NOT. file_stg_exist )  THEN
             message_string = 'Using the synthetic turbulence generator in idealized setups ' //   &
                              '(no mesoscale nesting or self nesting), requires a ' //             &
                              'STG_PROFILES file to have proper information about turbulent ' //   &
                              'length and time scales, as well as Reynolds stress components.'
             CALL message( 'stg_check_parameters', 'PA0753', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

       IF ( turbulent_inflow )  THEN
          message_string = 'Using synthetic turbulence generator in combination ' //               &
                           '&with turbulent_inflow = .T. is not allowed'
          CALL message( 'stg_check_parameters', 'PA0039', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Synthetic turbulence generator requires the parallel random generator
       IF ( random_generator /= 'random-parallel' )  THEN
          message_string = 'Using synthetic turbulence generator requires random_generator = ' //  &
                           'random-parallel.'
          CALL message( 'stg_check_parameters', 'PA0421', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

 END SUBROUTINE stg_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for synthetic turbulence generator
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_header ( io )

    INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file

!
!-- Write synthetic turbulence generator Header
    IF ( syn_turb_gen )  THEN
       WRITE( io, 1 )
       IF ( parametrize_inflow_turbulence )  THEN
          WRITE( io, 4 ) dt_stg_adjust
       ELSE
          WRITE( io, 5 )
       ENDIF
    ENDIF

1   FORMAT (//' Synthetic turbulence generator information:'/                                      &
              ' ------------------------------------------'/)
4   FORMAT ('    imposed turbulence statistics are parametrized and ajdusted to boundary-layer development each ', F8.2, ' s' )
5   FORMAT ('    imposed turbulence is read from file' )

 END SUBROUTINE stg_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the synthetic turbulence generator
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_init

#if defined( __parallel )
    INTEGER(KIND=MPI_ADDRESS_KIND) :: extent   !< extent of new MPI type
    INTEGER(KIND=MPI_ADDRESS_KIND) :: tob = 0  !< dummy variable
#endif

    INTEGER(iwp) :: i         !> grid index in x-direction
    INTEGER(iwp) :: j         !> loop index
    INTEGER(iwp) :: k         !< index
#if defined( __parallel )
    INTEGER(iwp) :: newtype   !< dummy MPI type
    INTEGER(iwp) :: realsize  !< size of REAL variables
#endif

    INTEGER(iwp), DIMENSION(3) ::  nr_non_topo_xz_l = 0 !< number of non-topography grid points at xz-cross-section on subdomain
    INTEGER(iwp), DIMENSION(3) ::  nr_non_topo_yz_l = 0 !< number of non-topography grid points at yz-cross-section on subdomain


!
!-- Dummy variables used for reading profiles
    REAL(wp) :: d1   !< u profile
    REAL(wp) :: d2   !< v profile
    REAL(wp) :: d3   !< w profile
    REAL(wp) :: d5   !< e profile
    REAL(wp) :: luy  !< length scale for u in y direction
    REAL(wp) :: luz  !< length scale for u in z direction
    REAL(wp) :: lvy  !< length scale for v in y direction
    REAL(wp) :: lvz  !< length scale for v in z direction
    REAL(wp) :: lwy  !< length scale for w in y direction
    REAL(wp) :: lwz  !< length scale for w in z direction
#if defined( __parallel )
    REAL(wp) :: nnz  !< increment used to determine processor decomposition of z-axis along x and y direction
#endif
    REAL(wp) :: zz   !< height


#if defined( __parallel )
    CALL MPI_BARRIER( comm2d, ierr )
#endif
!
!-- Create mpi-datatypes for exchange in case of non-local but distributed computation of the
!-- velocity seeds. This option is useful in case large turbulent length scales are present, where
!-- the computational effort becomes large and needs to be parallelized. For parameterized turbulence
!-- the length scales are small and computing the velocity seeds locally is faster (no overhead by
!-- communication).
    IF ( .NOT. compute_velocity_seeds_local )  THEN
#if defined( __parallel )
!
!--    Determine processor decomposition of z-axis along x- and y-direction
       nnz = nz / npex
       nzb_x_stg = 1 + myidx * INT( nnz )
       nzt_x_stg = ( myidx + 1 ) * INT( nnz )

       IF ( MOD( nz , npex ) /= 0  .AND.  myidx == id_stg_right )                              &
          nzt_x_stg = nzt_x_stg + myidx * ( nnz - INT( nnz ) )

       IF ( nesting_offline   .OR.  ( child_domain  .AND.  rans_mode_parent                        &
                               .AND.  .NOT.  rans_mode ) )  THEN
          nnz = nz / npey
          nzb_y_stg = 1 + myidy * INT( nnz )
          nzt_y_stg = ( myidy + 1 ) * INT( nnz )

          IF ( MOD( nz , npey ) /= 0  .AND.  myidy == id_stg_north )                           &
             nzt_y_stg = nzt_y_stg + myidy * ( nnz - INT( nnz ) )
       ENDIF

!
!--    Define MPI type used in stg_generate_seed_yz to gather vertical splitted velocity seeds
       CALL MPI_TYPE_SIZE( MPI_REAL, realsize, ierr )
       extent = 1 * realsize
!
!--    Set-up MPI datatyp to involve all cores for turbulence generation at yz layer
!--    stg_type_yz: yz-slice with vertical bounds nzb:nzt+1
       CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt-nzb+2,nyn-nys+1],                                    &
               [1,nyn-nys+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
       CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_yz, ierr )
       CALL MPI_TYPE_COMMIT( stg_type_yz, ierr )
       CALL MPI_TYPE_FREE( newtype, ierr )

       ! stg_type_yz_small: yz-slice with vertical bounds nzb_x_stg:nzt_x_stg+1
       CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt_x_stg-nzb_x_stg+2,nyn-nys+1],                        &
               [1,nyn-nys+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
       CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_yz_small, ierr )
       CALL MPI_TYPE_COMMIT( stg_type_yz_small, ierr )
       CALL MPI_TYPE_FREE( newtype, ierr )

       ! Receive count and displacement for MPI_GATHERV in stg_generate_seed_yz
       ALLOCATE( recv_count_yz(npex), displs_yz(npex) )

       recv_count_yz       = nzt_x_stg-nzb_x_stg + 1
       recv_count_yz(npex) = recv_count_yz(npex) + 1

       DO  j = 1, npex
          displs_yz(j) = 0 + (nzt_x_stg-nzb_x_stg+1) * (j-1)
       ENDDO
!
!--    Set-up MPI datatyp to involve all cores for turbulence generation at xz layer
!--    stg_type_xz: xz-slice with vertical bounds nzb:nzt+1
       IF ( nesting_offline  .OR.  ( child_domain .AND.  rans_mode_parent                          &
                                     .AND.  .NOT.  rans_mode ) )  THEN
          CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt-nzb+2,nxr-nxl+1],                                 &
                  [1,nxr-nxl+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
          CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_xz, ierr )
          CALL MPI_TYPE_COMMIT( stg_type_xz, ierr )
          CALL MPI_TYPE_FREE( newtype, ierr )

          ! stg_type_yz_small: xz-slice with vertical bounds nzb_x_stg:nzt_x_stg+1
          CALL MPI_TYPE_CREATE_SUBARRAY( 2, [nzt_y_stg-nzb_y_stg+2,nxr-nxl+1],                     &
                  [1,nxr-nxl+1], [0,0], MPI_ORDER_FORTRAN, MPI_REAL, newtype, ierr )
          CALL MPI_TYPE_CREATE_RESIZED( newtype, tob, extent, stg_type_xz_small, ierr )
          CALL MPI_TYPE_COMMIT( stg_type_xz_small, ierr )
          CALL MPI_TYPE_FREE( newtype, ierr )

          ! Receive count and displacement for MPI_GATHERV in stg_generate_seed_yz
          ALLOCATE( recv_count_xz(npey), displs_xz(npey) )

          recv_count_xz       = nzt_y_stg-nzb_y_stg + 1
          recv_count_xz(npey) = recv_count_xz(npey) + 1

          DO  j = 1, npey
             displs_xz(j) = 0 + (nzt_y_stg-nzb_y_stg+1) * (j-1)
          ENDDO

       ENDIF
#endif
    ENDIF
!
!-- Allocate required arrays.
!-- In case no offline nesting or self nesting is used, the arrary mean_inflow profiles is required.
!-- Check if it is already allocated, else allocate and initialize it appropriately. Note, in case
!-- turbulence and inflow information is read from file, mean_inflow_profiles is already allocated
!-- and initialized appropriately.
    IF ( .NOT. nesting_offline  .AND.  .NOT.  child_domain )  THEN
       IF ( .NOT. ALLOCATED( mean_inflow_profiles ) )  THEN
          ALLOCATE( mean_inflow_profiles(nzb:nzt+1,1:num_mean_inflow_profiles) )
          mean_inflow_profiles = 0.0_wp
          mean_inflow_profiles(:,1) = u_init
          mean_inflow_profiles(:,2) = v_init
!
!--       Even though potential temperature and humidity are not perturbed, they need to be
!--       initialized appropriately.
          IF ( .NOT. neutral )                                                                     &
             mean_inflow_profiles(:,4) = pt_init
          IF ( humidity )                                                                          &
             mean_inflow_profiles(:,6) = q_init
       ENDIF
    ENDIF
!
!-- Assign initial profiles. Note, this is only required if turbulent inflow from the left is
!-- desired, not in case of any of the nesting (offline or self nesting) approaches.
    IF ( .NOT. nesting_offline  .AND.  .NOT.  child_domain )  THEN
       u_init = mean_inflow_profiles(:,1)
       v_init = mean_inflow_profiles(:,2)
       e_init = MAXVAL( mean_inflow_profiles(:,5) )
    ENDIF

    ALLOCATE ( a11(nzb:nzt+1), a21(nzb:nzt+1), a22(nzb:nzt+1),                                     &
               a31(nzb:nzt+1), a32(nzb:nzt+1), a33(nzb:nzt+1),                                     &
               nux(nzb:nzt+1), nuy(nzb:nzt+1), nuz(nzb:nzt+1),                                     &
               nvx(nzb:nzt+1), nvy(nzb:nzt+1), nvz(nzb:nzt+1),                                     &
               nwx(nzb:nzt+1), nwy(nzb:nzt+1), nwz(nzb:nzt+1),                                     &
               r11(nzb:nzt+1), r21(nzb:nzt+1), r22(nzb:nzt+1),                                     &
               r31(nzb:nzt+1), r32(nzb:nzt+1), r33(nzb:nzt+1),                                     &
               tu(nzb:nzt+1),  tv(nzb:nzt+1),  tw(nzb:nzt+1) )

    r11 = 0.0_wp
    r21 = 0.0_wp
    r22 = 0.0_wp
    r31 = 0.0_wp
    r32 = 0.0_wp
    r33 = 0.0_wp
    tu  = 0.0_wp
    tv  = 0.0_wp
    tw  = 0.0_wp

    ALLOCATE ( dist_xz(nzb:nzt+1,nxl:nxr,3) )
    ALLOCATE ( dist_yz(nzb:nzt+1,nys:nyn,3) )
    dist_xz = 0.0_wp
    dist_yz = 0.0_wp
!
!-- Check if a STG_PROFILES file exist.
    INQUIRE( FILE = 'STG_PROFILES' // TRIM( coupling_char ), EXIST = file_stg_exist  )

    IF ( file_stg_exist )  THEN

       OPEN( 90, FILE = 'STG_PROFILES' // TRIM( coupling_char ), STATUS = 'OLD', FORM = 'FORMATTED' )
!
!--    Skip header
       READ( 90, * )
!
!--    Read inflow profile
!--    Height levels of profiles in input profile is as follows:
!--    zu: luy, luz, tu, lvy, lvz, tv, r11, r21, r22, d1, d2, d5 zw: lwy, lwz, tw, r31, r32, r33, d3
!--    WARNING: zz is not used at the moment
       DO  k = nzb, nzt
          READ( 90, * ) zz, luy, luz, tu(k), lvy, lvz, tv(k), lwy, lwz, tw(k), r11(k), r21(k),     &
                        r22(k), r31(k), r32(k), r33(k), d1, d2, d3, d5

!
!--       Convert length scales from meter to number of grid points.
          IF ( k /= nzb )  THEN
             nuz(k) = INT( luz * ddzw(k) )
             nvz(k) = INT( lvz * ddzw(k) )
             nwz(k) = INT( lwz * ddzw(k) )
          ELSE
             nuz(k) = INT( luz * ddzw(k+1) )
             nvz(k) = INT( lvz * ddzw(k+1) )
             nwz(k) = INT( lwz * ddzw(k+1) )
          ENDIF

          nuy(k) = INT( luy * ddy )
          nvy(k) = INT( lvy * ddy )
          nwy(k) = INT( lwy * ddy )
!
!--       Set length scales in x-direction. As a workaround assume isotropic turbulence in x- and
!--       y-direction.
          nwx(k) = nwy(k)
          nvx(k) = nvy(k)
          nux(k) = nuy(k)
!
!--       Save Mean inflow profiles
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' ) THEN
             mean_inflow_profiles(k,1) = d1
             mean_inflow_profiles(k,2) = d2
            !  mean_inflow_profiles(k,4) = d4
             mean_inflow_profiles(k,5) = d5
          ENDIF
       ENDDO
!
!--    Set length scales at the surface and top boundary. At the surface the lengths scales are
!--    simply overwritten.
       nwx(nzb) = nwy(nzb+1)
       nvx(nzb) = nvy(nzb+1)
       nux(nzb) = nuy(nzb+1)
       nuy(nzb) = nuy(nzb+1)
       nuz(nzb) = nuz(nzb+1)
       nvy(nzb) = nvy(nzb+1)
       nvz(nzb) = nvz(nzb+1)
       nwy(nzb) = nwy(nzb+1)
       nwz(nzb) = nwz(nzb+1)

       nwx(nzt+1) = nwy(nzt)
       nvx(nzt+1) = nvy(nzt)
       nux(nzt+1) = nuy(nzt)
       nuy(nzt+1) = nuy(nzt)
       nuz(nzt+1) = nuz(nzt)
       nvy(nzt+1) = nvy(nzt)
       nvz(nzt+1) = nvz(nzt)
       nwy(nzt+1) = nwy(nzt)
       nwz(nzt+1) = nwz(nzt)

       CLOSE( 90 )
!
!--    Calculate coefficient matrix from Reynolds stress tensor (Lund rotation)
       CALL calc_coeff_matrix
!
!-- No information about turbulence and its length scales are available. Instead, parametrize
!-- turbulence which is imposed at the boundaries. Set flag which indicates that turbulence is
!-- parametrized, which is done later when energy-balance models are already initialized. This is
!-- because the STG needs information about surface properties such as roughness to build
!-- 'realistic' turbulence profiles.
    ELSE
!
!--    Precalulate the inversion number of xy-grid points - later used for mass conservation
       d_nxy = 1.0_wp / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )
!
!--    Set flag indicating that turbulence is parametrized
       parametrize_inflow_turbulence = .TRUE.
!
!--    In case of dirichlet inflow boundary conditions only at one lateral boundary, i.e. in the
!--    case no offline or self nesting is applied but synthetic turbulence shall be parametrized
!--    nevertheless, the boundary-layer depth needs to be determined first.
       IF ( .NOT. nesting_offline  .AND.  .NOT.  child_domain )                                    &
          CALL nesting_offl_calc_zi
!
!--    Determine boundary-layer depth, which is used to initialize lenght scales
       CALL calc_scaling_variables
!
!--    Parametrize Reynolds-stress tensor, diagonal elements as well as r21 (v'u'), r31 (w'u'),
!--    r32 (w'v'). Parametrization follows Rotach et al. (1996) and is based on boundary-layer depth,
!--    friction velocity and velocity scale.
       CALL parametrize_turbulence
!
!--    Calculate coefficient matrix from Reynolds stress tensor (Lund rotation)
       CALL calc_coeff_matrix

    ENDIF
!
!-- Initial filter functions
    CALL stg_setup_filter_function
!
!-- Allocate velocity seeds for turbulence at xz-layer
    ALLOCATE ( fu_xz(nzb:nzt+1,nxl:nxr), fuo_xz(nzb:nzt+1,nxl:nxr),                                &
               fv_xz(nzb:nzt+1,nxl:nxr), fvo_xz(nzb:nzt+1,nxl:nxr),                                &
               fw_xz(nzb:nzt+1,nxl:nxr), fwo_xz(nzb:nzt+1,nxl:nxr) )

!
!-- Allocate velocity seeds for turbulence at yz-layer
    ALLOCATE ( fu_yz(nzb:nzt+1,nys:nyn), fuo_yz(nzb:nzt+1,nys:nyn),                                &
               fv_yz(nzb:nzt+1,nys:nyn), fvo_yz(nzb:nzt+1,nys:nyn),                                &
               fw_yz(nzb:nzt+1,nys:nyn), fwo_yz(nzb:nzt+1,nys:nyn) )

    fu_xz  = 0.0_wp
    fuo_xz = 0.0_wp
    fv_xz  = 0.0_wp
    fvo_xz = 0.0_wp
    fw_xz  = 0.0_wp
    fwo_xz = 0.0_wp

    fu_yz  = 0.0_wp
    fuo_yz = 0.0_wp
    fv_yz  = 0.0_wp
    fvo_yz = 0.0_wp
    fw_yz  = 0.0_wp
    fwo_yz = 0.0_wp

#if defined( __parallel )
    CALL MPI_BARRIER( comm2d, ierr )
#endif

!
!-- In case of restart, calculate velocity seeds fu, fv, fw from former time step.
!   Bug: fu, fv, fw are different in those heights where a11, a22, a33 are 0 compared to the prerun.
!-- This is mostly for k=nzt+1.
    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
       IF ( myidx == id_stg_left  .OR.  myidx == id_stg_right )  THEN

          IF ( myidx == id_stg_left  )  i = -1
          IF ( myidx == id_stg_right )  i = nxr+1

          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                IF  ( a11(k) > 10E-8_wp )  THEN
                   fu_yz(k,j) = ( u(k,j,i) - u_init(k) ) / a11(k)
                ELSE
                   fu_yz(k,j) = 10E-8_wp
                ENDIF

                IF  ( a22(k) > 10E-8_wp )  THEN
                   fv_yz(k,j) = ( v(k,j,i) - a21(k) * fu_yz(k,j) - v_init(k) ) / a22(k)
                ELSE
                   fv_yz(k,j) = 10E-8_wp
                ENDIF

                IF  ( a33(k) > 10E-8_wp )  THEN
                   fw_yz(k,j) = ( w(k,j,i) - a31(k) * fu_yz(k,j) - a32(k) * fv_yz(k,j) ) / a33(k)
                ELSE
                   fw_yz(k,j) = 10E-8_wp
                ENDIF
             ENDDO
          ENDDO
       ENDIF

       IF ( myidy == id_stg_south  .OR.  myidy == id_stg_north )  THEN

          IF ( myidy == id_stg_south )  j = -1
          IF ( myidy == id_stg_north )  j = nyn+1

          DO  i = nxl, nxr
             DO  k = nzb, nzt+1
!
!--             In case the correlation coefficients are very small, the velocity seeds may become
!--             very large finally creating numerical instabilities in the synthetic turbulence
!--             generator. Empirically, a value of 10E-8 seems to be sufficient.
                IF  ( a11(k) > 10E-8_wp )  THEN
                   fu_xz(k,i) = ( u(k,j,i) - u_init(k) ) / a11(k)
                ELSE
                   fu_xz(k,i) = 10E-8_wp
                ENDIF

                IF  ( a22(k) > 10E-8_wp )  THEN
                   fv_xz(k,i) = ( v(k,j,i) - a21(k) * fu_xz(k,i) - v_init(k) ) / a22(k)
                ELSE
                   fv_xz(k,i) = 10E-8_wp
                ENDIF

                IF  ( a33(k) > 10E-8_wp )  THEN
                   fw_xz(k,i) = ( w(k,j,i) - a31(k) * fu_xz(k,i) - a32(k) * fv_xz(k,i) ) / a33(k)
                ELSE
                   fw_xz(k,i) = 10E-8_wp
                ENDIF

             ENDDO
          ENDDO
       ENDIF
    ENDIF
!
!-- Count the number of non-topography grid points at the boundaries where perturbations are imposed.
!-- This number is later used for bias corrections of the perturbations, i.e. to force their mean to
!-- be zero. Please note, due to the asymetry of u and v along x and y direction, respectively,
!-- different cases must be distinguished.
    IF ( myidx == id_stg_left  .OR.  myidx == id_stg_right )  THEN
!
!--    Number of grid points where perturbations are imposed on u
       IF ( myidx == id_stg_left  )  i = nxl
       IF ( myidx == id_stg_right )  i = nxr+1

       nr_non_topo_yz_l(1) = SUM( MERGE( 1, 0, BTEST( topo_flags(nzb:nzt,nys:nyn,i), 1 ) ) )
!
!--    Number of grid points where perturbations are imposed on v and w
       IF ( myidx == id_stg_left  )  i = nxl-1
       IF ( myidx == id_stg_right )  i = nxr+1

       nr_non_topo_yz_l(2) = SUM( MERGE( 1, 0, BTEST( topo_flags(nzb:nzt,nysv:nyn,i), 2 ) ) )
       nr_non_topo_yz_l(3) = SUM( MERGE( 1, 0, BTEST( topo_flags(nzb:nzt,nys:nyn,i), 3 ) ) )

#if defined( __parallel )
       CALL MPI_ALLREDUCE( nr_non_topo_yz_l, nr_non_topo_yz, 3, MPI_INTEGER, MPI_SUM, comm1dy, ierr )
#else
       nr_non_topo_yz = nr_non_topo_yz_l
#endif
    ENDIF

    IF ( myidy == id_stg_south  .OR.  myidy == id_stg_north )  THEN
!
!--    Number of grid points where perturbations are imposed on v
       IF ( myidy == id_stg_south )  j = nys
       IF ( myidy == id_stg_north )  j = nyn+1

       nr_non_topo_xz_l(2) = SUM( MERGE( 1, 0, BTEST( topo_flags(nzb:nzt,j,nxl:nxr), 2 ) ) )
!
!--    Number of grid points where perturbations are imposed on u and w
       IF ( myidy == id_stg_south )  j = nys-1
       IF ( myidy == id_stg_north )  j = nyn+1

       nr_non_topo_xz_l(1) = SUM( MERGE( 1, 0, BTEST( topo_flags(nzb:nzt,j,nxlu:nxr), 1 ) ) )
       nr_non_topo_xz_l(3) = SUM( MERGE( 1, 0, BTEST( topo_flags(nzb:nzt,j,nxl:nxr), 3 ) ) )

#if defined( __parallel )
       CALL MPI_ALLREDUCE( nr_non_topo_xz_l, nr_non_topo_xz, 3, MPI_INTEGER, MPI_SUM, comm1dx, ierr )
#else
       nr_non_topo_xz = nr_non_topo_xz_l
#endif
    ENDIF
!
!-- Initialize random number generator at xz- and yz-layers. Random numbers are initialized at each
!-- core. In case there is only inflow from the left, it is sufficient to generate only random
!-- numbers for the yz-layer, else random numbers for the xz-layer are also required.
    ALLOCATE ( id_rand_yz(-mergp_limit+nys:nyn+mergp_limit) )
    ALLOCATE ( seq_rand_yz(5,-mergp_limit+nys:nyn+mergp_limit) )
    id_rand_yz  = 0
    seq_rand_yz = 0

    CALL init_parallel_random_generator( ny, -mergp_limit+nys, nyn+mergp_limit, id_rand_yz, seq_rand_yz )

    IF ( nesting_offline  .OR.  ( child_domain .AND.  rans_mode_parent                             &
         .AND.  .NOT.  rans_mode ) )  THEN
       ALLOCATE ( id_rand_xz(-mergp_limit+nxl:nxr+mergp_limit) )
       ALLOCATE ( seq_rand_xz(5,-mergp_limit+nxl:nxr+mergp_limit) )
       id_rand_xz  = 0
       seq_rand_xz = 0

       CALL init_parallel_random_generator( nx, -mergp_limit+nxl, nxr+mergp_limit, id_rand_xz, seq_rand_xz )
    ENDIF

 END SUBROUTINE stg_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate filter function bxx from length scale nxx following Eg.9 and 10 (Xie and Castro, 2008)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_filter_func( nxx, bxx, mergp )

    INTEGER(iwp) ::  k     !< loop index
    INTEGER(iwp) ::  mergp !< passed length scale in grid points
    INTEGER(iwp) ::  n_k   !< length scale nXX in height k
    INTEGER(iwp) ::  nf    !< index for length scales

    INTEGER(iwp), DIMENSION(nzb:nzt+1) ::  nxx  !< length scale (in gp)

    REAL(wp) :: bdenom        !< denominator for filter functions bXX
    REAL(wp) :: qsi = 1.0_wp  !< minimization factor

    REAL(wp), DIMENSION(-mergp:mergp,nzb:nzt+1) ::  bxx  !< filter function


    bxx = 0.0_wp

    DO  k = nzb, nzt+1
       bdenom = 0.0_wp
       n_k    = nxx(k)
       IF ( n_k /= 0 )  THEN

!
!--       ( Eq.10 )^2
          DO  nf = -n_k, n_k
             bdenom = bdenom + EXP( -qsi * pi * ABS( nf ) / n_k )**2
          ENDDO

!
!--       ( Eq.9 )
          bdenom = SQRT( bdenom )
          DO  nf = -n_k, n_k
             bxx(nf,k) = EXP( -qsi * pi * ABS( nf ) / n_k ) / bdenom
          ENDDO
       ENDIF
    ENDDO

 END SUBROUTINE stg_filter_func


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate filter function bxx from length scale nxx following Eg.9 and 10
!> (Xie and Castro, 2008)
!------------------------------------------------------------------------------!
 SUBROUTINE stg_setup_filter_function

    INTEGER(iwp) :: j  !< dummy value to calculate maximum length scale index
    INTEGER(iwp) :: k  !< loop index along vertical direction
!
!-- Define the size of the filter functions and allocate them.
    mergp_x = 0
    mergp_y = 0
    mergp_z = 0

    DO  k = nzb, nzt+1
       j = MAX( ABS(nux(k)), ABS(nvx(k)), ABS(nwx(k)) )
       IF ( j > mergp_x )  mergp_x = j
    ENDDO
    DO  k = nzb, nzt+1
       j = MAX( ABS(nuy(k)), ABS(nvy(k)), ABS(nwy(k)) )
       IF ( j > mergp_y )  mergp_y = j
    ENDDO
    DO  k = nzb, nzt+1
       j = MAX( ABS(nuz(k)), ABS(nvz(k)), ABS(nwz(k)) )
       IF ( j > mergp_z )  mergp_z = j
    ENDDO

    mergp_xy = MAX( mergp_x, mergp_y )

    IF ( ALLOCATED( bux ) )  DEALLOCATE( bux )
    IF ( ALLOCATED( bvx ) )  DEALLOCATE( bvx )
    IF ( ALLOCATED( bwx ) )  DEALLOCATE( bwx )
    IF ( ALLOCATED( buy ) )  DEALLOCATE( buy )
    IF ( ALLOCATED( bvy ) )  DEALLOCATE( bvy )
    IF ( ALLOCATED( bwy ) )  DEALLOCATE( bwy )
    IF ( ALLOCATED( buz ) )  DEALLOCATE( buz )
    IF ( ALLOCATED( bvz ) )  DEALLOCATE( bvz )
    IF ( ALLOCATED( bwz ) )  DEALLOCATE( bwz )

    ALLOCATE ( bux(-mergp_x:mergp_x,nzb:nzt+1),                                                    &
               buy(-mergp_y:mergp_y,nzb:nzt+1),                                                    &
               buz(-mergp_z:mergp_z,nzb:nzt+1),                                                    &
               bvx(-mergp_x:mergp_x,nzb:nzt+1),                                                    &
               bvy(-mergp_y:mergp_y,nzb:nzt+1),                                                    &
               bvz(-mergp_z:mergp_z,nzb:nzt+1),                                                    &
               bwx(-mergp_x:mergp_x,nzb:nzt+1),                                                    &
               bwy(-mergp_y:mergp_y,nzb:nzt+1),                                                    &
               bwz(-mergp_z:mergp_z,nzb:nzt+1) )
!
!-- Create filter functions
    CALL stg_filter_func( nux, bux, mergp_x ) !filter ux
    CALL stg_filter_func( nuy, buy, mergp_y ) !filter uy
    CALL stg_filter_func( nuz, buz, mergp_z ) !filter uz
    CALL stg_filter_func( nvx, bvx, mergp_x ) !filter vx
    CALL stg_filter_func( nvy, bvy, mergp_y ) !filter vy
    CALL stg_filter_func( nvz, bvz, mergp_z ) !filter vz
    CALL stg_filter_func( nwx, bwx, mergp_x ) !filter wx
    CALL stg_filter_func( nwy, bwy, mergp_y ) !filter wy
    CALL stg_filter_func( nwz, bwz, mergp_z ) !filter wz

 END SUBROUTINE stg_setup_filter_function

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &stg_par for synthetic turbulence generator
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /stg_par/  compute_velocity_seeds_local,                                              &
                        dt_stg_adjust,                                                             &
                        dt_stg_call,                                                               &
                        switch_off_module


!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, stg_par, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    stg_par namelist was found and read correctly. Set flag that indicates that the synthetic
!--    turbulence generator is switched on.
       IF ( .NOT. switch_off_module )  syn_turb_gen = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    stg_par namelist was found but contained errors. Print an error message including the line
!--    that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'stg_par', line )

    ENDIF

 END SUBROUTINE stg_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_rrd_global_ftn( found )

    LOGICAL, INTENT(OUT) ::  found !< flag indicating if variable was found

    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'time_stg_adjust' )
          READ ( 13 )  time_stg_adjust

       CASE ( 'time_stg_call' )
          READ ( 13 )  time_stg_call

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE stg_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_rrd_global_mpi

    CALL rrd_mpi_io( 'time_stg_adjust', time_stg_adjust )
    CALL rrd_mpi_io( 'time_stg_call', time_stg_call )

 END SUBROUTINE stg_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_wrd_global

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

      CALL wrd_write_string( 'time_stg_adjust' )
       WRITE ( 14 )  time_stg_adjust

       CALL wrd_write_string( 'time_stg_call' )
       WRITE ( 14 )  time_stg_call

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       CALL wrd_mpi_io( 'time_stg_adjust', time_stg_adjust )
       CALL wrd_mpi_io( 'time_stg_call', time_stg_call )

    ENDIF

 END SUBROUTINE stg_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create turbulent inflow fields for u, v, w with prescribed length scales and Reynolds stress
!> tensor after a method of Xie and Castro (2008), modified following suggestions of Kim et al.
!> (2013), and using a Lund rotation (Lund, 1998).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_main

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    INTEGER(iwp) ::  i  !< grid index in x-direction
    INTEGER(iwp) ::  j  !< loop index in y-direction
    INTEGER(iwp) ::  k  !< loop index in z-direction

    LOGICAL  ::  stg_call = .FALSE.  !< control flag indicating whether turbulence was updated or only restored from last call

    REAL(wp) ::  dt_stg  !< time interval the STG is called

    REAL(wp), DIMENSION(3) ::  mc_factor_l  !< local mass flux correction factor

    IF ( debug_output_timestep )  CALL debug_message( 'stg_main', 'start' )
!
!-- Calculate time step which is needed for filter functions
    dt_stg = MAX( dt_3d, dt_stg_call )
!
!-- Check if synthetic turbulence generator needs to be called and new perturbations need to be
!-- created or if old disturbances can be imposed again.
    IF ( time_stg_call >= dt_stg_call  .AND.                                                       &
         intermediate_timestep_count == intermediate_timestep_count_max  )  THEN
       stg_call = .TRUE.
    ELSE
       stg_call = .FALSE.
    ENDIF
!
!-- Initial value of fu, fv, fw
    IF ( time_since_reference_point == 0.0_wp .AND. .NOT. velocity_seed_initialized )  THEN
!
!--    Generate turbulence at the left and right boundary. Random numbers for the yz-planes at the
!--    left/right boundary are generated by the left-sided mpi ranks only. After random numbers are
!--    calculated, they are distributed to all other mpi ranks in the model, so that the velocity
!--    seed matrices are available on all mpi ranks (i.e. also on the right-sided boundary mpi ranks).
!--    In case of offline nesting, this implies, that the left- and the right-sided lateral boundary
!--    have the same initial seeds.
!--    Note, in case of inflow from the right only, only turbulence at the left boundary is required.
       IF ( .NOT. ( nesting_offline  .OR.  ( child_domain .AND.  rans_mode_parent                  &
            .AND.  .NOT.  rans_mode ) ) )  THEN
          CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fu_yz, id_stg_left )
          CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fv_yz, id_stg_left )
          CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fw_yz, id_stg_left )
       ELSE
          CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fu_yz, id_stg_left, id_stg_right )
          CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fv_yz, id_stg_left, id_stg_right )
          CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fw_yz, id_stg_left, id_stg_right )
!
!--       Generate turbulence at the south and north boundary. Random numbers for the xz-planes at
!--       the south/north boundary are generated by the south-sided mpi ranks only. Please see also
!--       comment above.
          CALL stg_generate_seed_xz( nux, nuz, bux, buz, fu_xz, id_stg_south, id_stg_north )
          CALL stg_generate_seed_xz( nvx, nvz, bvx, bvz, fv_xz, id_stg_south, id_stg_north )
          CALL stg_generate_seed_xz( nwx, nwz, bwx, bwz, fw_xz, id_stg_south, id_stg_north )
       ENDIF
       velocity_seed_initialized = .TRUE.
    ENDIF
!
!-- New set of fu, fv, fw. Note, for inflow from the left side only, velocity seeds are only
!-- required at the left boundary, while in case of offline nesting or RANS-LES nesting velocity
!-- seeds are required also at the right, south and north boundaries.
    IF ( stg_call )  THEN
       IF ( .NOT. ( nesting_offline  .OR.                                                          &
                  ( child_domain .AND.  rans_mode_parent                                           &
                                 .AND.  .NOT.  rans_mode ) ) )  THEN
          CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fuo_yz, id_stg_left )
          CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fvo_yz, id_stg_left )
          CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fwo_yz, id_stg_left )

       ELSE
          CALL stg_generate_seed_yz( nuy, nuz, buy, buz, fuo_yz, id_stg_left, id_stg_right )
          CALL stg_generate_seed_yz( nvy, nvz, bvy, bvz, fvo_yz, id_stg_left, id_stg_right )
          CALL stg_generate_seed_yz( nwy, nwz, bwy, bwz, fwo_yz, id_stg_left, id_stg_right )

          CALL stg_generate_seed_xz( nux, nuz, bux, buz, fuo_xz, id_stg_south, id_stg_north )
          CALL stg_generate_seed_xz( nvx, nvz, bvx, bvz, fvo_xz, id_stg_south, id_stg_north )
          CALL stg_generate_seed_xz( nwx, nwz, bwx, bwz, fwo_xz, id_stg_south, id_stg_north )
       ENDIF
    ENDIF

!
!-- Turbulence generation at left and/or right boundary
    IF ( myidx == id_stg_left  .OR.  myidx == id_stg_right )  THEN
!
!--    Calculate new set of perturbations. Do this only at last RK3-substep and when dt_stg_call is
!--    exceeded. Else the old set of perturbations is imposed
       IF ( stg_call )  THEN

          DO  j = nys, nyn
             DO  k = nzb, nzt + 1
!
!--             Update fu, fv, fw following Eq. 14 of Xie and Castro (2008)
                IF ( tu(k) == 0.0_wp  .OR.  adjustment_step )  THEN
                   fu_yz(k,j) = fuo_yz(k,j)
                ELSE
                   fu_yz(k,j) = fu_yz(k,j) * EXP( -pi * dt_stg * 0.5_wp / tu(k) ) +                &
                                fuo_yz(k,j) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tu(k) ) )
                ENDIF

                IF ( tv(k) == 0.0_wp  .OR.  adjustment_step )  THEN
                   fv_yz(k,j) = fvo_yz(k,j)
                ELSE
                   fv_yz(k,j) = fv_yz(k,j) * EXP( -pi * dt_stg * 0.5_wp / tv(k) ) +                &
                                fvo_yz(k,j) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tv(k) ) )
                ENDIF

                IF ( tw(k) == 0.0_wp  .OR.  adjustment_step )  THEN
                   fw_yz(k,j) = fwo_yz(k,j)
                ELSE
                   fw_yz(k,j) = fw_yz(k,j) * EXP( -pi * dt_stg * 0.5_wp / tw(k) ) +                &
                                fwo_yz(k,j) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tw(k) ) )
                ENDIF
             ENDDO
          ENDDO

          dist_yz(nzb,:,1) = 0.0_wp
          dist_yz(nzb,:,2) = 0.0_wp
          dist_yz(nzb,:,3) = 0.0_wp

          IF ( myidx == id_stg_left  )  i = nxl
          IF ( myidx == id_stg_right )  i = nxr+1
          DO  j = nys, nyn
             DO  k = nzb+1, nzt + 1
!
!--             Lund rotation following Eq. 17 in Xie and Castro (2008).
                dist_yz(k,j,1) = a11(k) * fu_yz(k,j) *                                             &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

             ENDDO
          ENDDO

          IF ( myidx == id_stg_left  )  i = nxl-1
          IF ( myidx == id_stg_right )  i = nxr+1
          DO  j = nys, nyn
             DO  k = nzb+1, nzt + 1
!
!--             Lund rotation following Eq. 17 in Xie and Castro (2008).
!--             Additional factors are added to improve the variance of v and w experimental test
!--             of 1.2
                dist_yz(k,j,2) = ( a21(k) * fu_yz(k,j) + a22(k) * fv_yz(k,j) ) *                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                dist_yz(k,j,3) = ( a31(k) * fu_yz(k,j) + a32(k) * fv_yz(k,j) +                     &
                                   a33(k) * fw_yz(k,j) ) *                                         &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )

             ENDDO
          ENDDO
       ENDIF
!
!--    Mass flux correction following Kim et al. (2013)
!--    This correction factor ensures that the mass flux is preserved at the inflow boundary. First,
!--    calculate mean value of the imposed perturbations at yz boundary. Note, this needs to be done
!--    only after the last RK3-substep, else the boundary values won't be accessed.
       IF ( intermediate_timestep_count == intermediate_timestep_count_max )  THEN
          mc_factor_l = 0.0_wp
          mc_factor   = 0.0_wp
!
!--       Sum up the original volume flows (with and without perturbations).
!--       Note, for non-normal components (here v and w) it is actually no volume flow.
          IF ( myidx == id_stg_left  )  i = nxl
          IF ( myidx == id_stg_right )  i = nxr+1

          mc_factor_l(1) = SUM( dist_yz(nzb:nzt,nys:nyn,1) *                                       &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzb:nzt,nys:nyn,i), 1 ) ) )

          IF ( myidx == id_stg_left  )  i = nxl-1
          IF ( myidx == id_stg_right )  i = nxr+1

          mc_factor_l(2) = SUM( dist_yz(nzb:nzt,nysv:nyn,2) *                                      &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzb:nzt,nysv:nyn,i), 2 ) ) )
          mc_factor_l(3) = SUM( dist_yz(nzb:nzt,nys:nyn,3) *                                       &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzb:nzt,nys:nyn,i), 3 ) ) )

#if defined( __parallel )
          CALL MPI_ALLREDUCE( mc_factor_l, mc_factor, 3, MPI_REAL, MPI_SUM, comm1dy, ierr )
#else
          mc_factor   = mc_factor_l
#endif
!
!--       Calculate correction factor and force zero mean perturbations.
          mc_factor = mc_factor / REAL( nr_non_topo_yz, KIND = wp )

          IF ( myidx == id_stg_left  )  i = nxl
          IF ( myidx == id_stg_right )  i = nxr+1

          dist_yz(:,nys:nyn,1) = ( dist_yz(:,nys:nyn,1) - mc_factor(1) ) *                         &
                                MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(:,nys:nyn,i), 1 ) )


          IF ( myidx == id_stg_left  )  i = nxl-1
          IF ( myidx == id_stg_right )  i = nxr+1

          dist_yz(:,nys:nyn,2) = ( dist_yz(:,nys:nyn,2) - mc_factor(2) ) *                         &
                                MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(:,nys:nyn,i), 2 ) )

          dist_yz(:,nys:nyn,3) = ( dist_yz(:,nys:nyn,3) - mc_factor(3) ) *                         &
                                MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(:,nys:nyn,i), 3 ) )
!
!--       Add disturbances
          IF ( myidx == id_stg_left  )  THEN
!
!--          For the left boundary distinguish between mesoscale offline / self nesting and
!--          turbulent inflow at the left boundary only. In the latter case turbulence is imposed on
!--          the mean inflow profiles.
             IF ( .NOT. nesting_offline  .AND.  .NOT. child_domain )  THEN
!
!--             Add disturbance at the inflow
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      u(k,j,-nbgp+1:0) = ( mean_inflow_profiles(k,1) + dist_yz(k,j,1) ) *          &
                                     MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,0), 1 ) )
                      v(k,j,-nbgp:-1)  = ( mean_inflow_profiles(k,2) + dist_yz(k,j,2) ) *          &
                                     MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,-1), 2 ) )
                      w(k,j,-nbgp:-1)  = dist_yz(k,j,3) *                                          &
                                     MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,-1), 3 ) )
                   ENDDO
                ENDDO
             ELSE

                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      u(k,j,0)   = ( u(k,j,0)  + dist_yz(k,j,1) ) *                                &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,0), 1 ) )
                      u(k,j,-1)  = u(k,j,0)
                      v(k,j,-1)  = ( v(k,j,-1) + dist_yz(k,j,2) ) *                                &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,-1), 2 ) )
                      w(k,j,-1)  = ( w(k,j,-1) + dist_yz(k,j,3) ) *                                &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,-1), 3 ) )
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
          IF ( myidx == id_stg_right  )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   u(k,j,nxr+1) = ( u(k,j,nxr+1) + dist_yz(k,j,1) ) *                              &
                                  MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,nxr+1), 1 ) )
                   v(k,j,nxr+1) = ( v(k,j,nxr+1) + dist_yz(k,j,2) ) *                              &
                                  MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,nxr+1), 2 ) )
                   w(k,j,nxr+1) = ( w(k,j,nxr+1) + dist_yz(k,j,3) ) *                              &
                                  MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,nxr+1), 3 ) )
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDIF
!
!-- Turbulence generation at north and south boundary
    IF ( myidy == id_stg_north  .OR.  myidy == id_stg_south )  THEN
!
!--    Calculate new set of perturbations. Do this only at last RK3-substep and when dt_stg_call is
!--    exceeded. Else the old set of perturbations is imposed
       IF ( stg_call )  THEN
          DO  i = nxl, nxr
             DO  k = nzb, nzt + 1
!
!--             Update fu, fv, fw following Eq. 14 of Xie and Castro (2008)
                IF ( tu(k) == 0.0_wp  .OR.  adjustment_step )  THEN
                   fu_xz(k,i) = fuo_xz(k,i)
                ELSE
                   fu_xz(k,i) = fu_xz(k,i) * EXP( -pi * dt_stg * 0.5_wp / tu(k) ) +                &
                                fuo_xz(k,i) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tu(k) ) )
                ENDIF

                IF ( tv(k) == 0.0_wp  .OR.  adjustment_step )  THEN
                   fv_xz(k,i) = fvo_xz(k,i)
                ELSE
                   fv_xz(k,i) = fv_xz(k,i) * EXP( -pi * dt_stg * 0.5_wp / tv(k) ) +                &
                                fvo_xz(k,i) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tv(k) ) )
                ENDIF

                IF ( tw(k) == 0.0_wp  .OR.  adjustment_step )  THEN
                   fw_xz(k,i) = fwo_xz(k,i)
                ELSE
                   fw_xz(k,i) = fw_xz(k,i) * EXP( -pi * dt_stg * 0.5_wp / tw(k) ) +                &
                                fwo_xz(k,i) * SQRT( 1.0_wp - EXP( -pi * dt_stg / tw(k) ) )
                ENDIF
             ENDDO
          ENDDO


          dist_xz(nzb,:,1) = 0.0_wp
          dist_xz(nzb,:,2) = 0.0_wp
          dist_xz(nzb,:,3) = 0.0_wp

          IF ( myidy == id_stg_south  ) j = nys
          IF ( myidy == id_stg_north )  j = nyn+1
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt + 1
!
!--             Lund rotation following Eq. 17 in Xie and Castro (2008).
!--             Additional factors are added to improve the variance of v and w
                !experimental test of 1.2
                dist_xz(k,i,2) = ( a21(k) * fu_xz(k,i) + a22(k) * fv_xz(k,i) ) *                   &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
             ENDDO
          ENDDO

          IF ( myidy == id_stg_south  ) j = nys-1
          IF ( myidy == id_stg_north )  j = nyn+1
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt + 1
!
!--             Lund rotation following Eq. 17 in Xie and Castro (2008).
!--             Additional factors are added to improve the variance of v and w
                dist_xz(k,i,1) = a11(k) * fu_xz(k,i) *                                             &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                dist_xz(k,i,3) = ( a31(k) * fu_xz(k,i) + a32(k) * fv_xz(k,i) +                     &
                                   a33(k) * fw_xz(k,i) ) *                                         &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
             ENDDO
          ENDDO
       ENDIF

!
!--    Mass flux correction following Kim et al. (2013). This correction factor ensures that the
!--    mass flux is preserved at the inflow boundary. First, calculate mean value of the imposed
!--    perturbations at yz boundary. Note, this needs to be done only after the last RK3-substep,
!--    else the boundary values won't be accessed.
       IF ( intermediate_timestep_count == intermediate_timestep_count_max )  THEN
          mc_factor_l = 0.0_wp
          mc_factor   = 0.0_wp

          IF ( myidy == id_stg_south  ) j = nys
          IF ( myidy == id_stg_north )  j = nyn+1

          mc_factor_l(2) = SUM( dist_xz(nzb:nzt,nxl:nxr,2) *                                       &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzb:nzt,j,nxl:nxr), 2 ) ) )

          IF ( myidy == id_stg_south ) j = nys-1
          IF ( myidy == id_stg_north ) j = nyn+1

          mc_factor_l(1) = SUM( dist_xz(nzb:nzt,nxlu:nxr,1) *                                      &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzb:nzt,j,nxlu:nxr), 1 ) ) )
          mc_factor_l(3) = SUM( dist_xz(nzb:nzt,nxl:nxr,3) *                                       &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(nzb:nzt,j,nxl:nxr), 3 ) ) )

#if defined( __parallel )
          CALL MPI_ALLREDUCE( mc_factor_l, mc_factor, 3, MPI_REAL, MPI_SUM, comm1dx, ierr )
#else
          mc_factor   = mc_factor_l
#endif

          mc_factor = mc_factor / REAL( nr_non_topo_xz, KIND = wp )

          IF ( myidy == id_stg_south ) j = nys
          IF ( myidy == id_stg_north ) j = nyn+1

          dist_xz(:,nxl:nxr,2) = ( dist_xz(:,nxl:nxr,2) - mc_factor(2) ) *                         &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(:,j,nxl:nxr), 2 ) )


          IF ( myidy == id_stg_south ) j = nys-1
          IF ( myidy == id_stg_north ) j = nyn+1

          dist_xz(:,nxl:nxr,1) = ( dist_xz(:,nxl:nxr,1) - mc_factor(1) ) *                         &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(:,j,nxl:nxr), 1 ) )

          dist_xz(:,nxl:nxr,3) = ( dist_xz(:,nxl:nxr,3) - mc_factor(3) ) *                         &
                                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(:,j,nxl:nxr), 3 ) )
!
!--       Add disturbances
          IF ( myidy == id_stg_south )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   u(k,-1,i) = ( u(k,-1,i) + dist_xz(k,i,1) )                                      &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,-1,i), 1 ) )
                   v(k,0,i)  = ( v(k,0,i)  + dist_xz(k,i,2) )                                      &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,0,i), 2 ) )
                   v(k,-1,i) = v(k,0,i)
                   w(k,-1,i) = ( w(k,-1,i) + dist_xz(k,i,3) )                                      &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,-1,i), 3 ) )
                ENDDO
             ENDDO
          ENDIF
          IF ( myidy == id_stg_north  )  THEN

             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   u(k,nyn+1,i) = ( u(k,nyn+1,i) + dist_xz(k,i,1) ) *                              &
                                  MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,nyn+1,i), 1 ) )
                   v(k,nyn+1,i) = ( v(k,nyn+1,i) + dist_xz(k,i,2) ) *                              &
                                  MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,nyn+1,i), 2 ) )
                   w(k,nyn+1,i) = ( w(k,nyn+1,i) + dist_xz(k,i,3) ) *                              &
                                  MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,nyn+1,i), 3 ) )
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDIF
!
!-- Exchange ghost points.
    CALL exchange_horiz( u, nbgp )
    CALL exchange_horiz( v, nbgp )
    CALL exchange_horiz( w, nbgp )
!
!-- Finally, set time counter for calling STG to zero
    IF ( stg_call )  time_stg_call = 0.0_wp
!
!-- Set adjustment step to False to indicate that time correlation can be
!-- switched-on again.
    adjustment_step = .FALSE.

    IF ( debug_output_timestep )  CALL debug_message( 'stg_main', 'end' )

 END SUBROUTINE stg_main

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Generate a set of random number rand_it wich is equal on each PE and calculate the velocity seed
!> f_n. f_n is splitted in vertical direction by the number of PEs in x-direction and each PE
!> calculates a vertical subsection of f_n. At the the end, all parts are collected to form the full
!> array.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_generate_seed_yz( n_y, n_z, b_y, b_z, f_n, id_left, id_right )

    INTEGER(iwp)           ::  id_left     !< core ids at respective boundaries
    INTEGER(iwp), OPTIONAL ::  id_right    !< core ids at respective boundaries
    INTEGER(iwp)           ::  j           !< loop index in y-direction
    INTEGER(iwp)           ::  jj          !< loop index in y-direction
    INTEGER(iwp)           ::  k           !< loop index in z-direction
    INTEGER(iwp)           ::  kk          !< loop index in z-direction
    INTEGER(iwp)           ::  send_count  !< send count for MPI_GATHERV

    INTEGER(iwp), DIMENSION(nzb:nzt+1) ::  n_y  !< length scale in y-direction
    INTEGER(iwp), DIMENSION(nzb:nzt+1) ::  n_z  !< length scale in z-direction

    REAL(wp) ::  nyz_inv         !< inverse of number of grid points in yz-slice
    REAL(wp) ::  rand_av         !< average of random number
    REAL(wp) ::  rand_sigma_inv  !< inverse of stdev of random number

    REAL(wp), DIMENSION(-mergp_y:mergp_y,nzb:nzt+1) :: b_y     !< filter function in y-direction
    REAL(wp), DIMENSION(-mergp_z:mergp_z,nzb:nzt+1) :: b_z     !< filter function in z-direction

    REAL(wp), DIMENSION(nzb_x_stg:nzt_x_stg+1,nys:nyn) ::  f_n_l  !<  local velocity seed
    REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn)             ::  f_n    !<  velocity seed

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rand_it  !< global array of random numbers
!
!-- Generate random numbers using the parallel random generator. The set of random numbers are
!-- modified to have an average of 0 and unit variance. Note, at the end the random number array
!-- must be defined globally in order to compute the correlation matrices. However, the calculation
!-- and normalization of random numbers is done locally, while the result is later distributed to
!-- all processes. Further, please note, a set of random numbers is only calculated for the west
!-- boundary, while the east boundary uses the same random numbers and thus also computes the same
!-- correlation matrix.
    ALLOCATE( rand_it(nzb-mergp_z:nzt+1+mergp_z,-mergp_y+nys:nyn+mergp_y) )
    rand_it = 0.0_wp

    rand_av        = 0.0_wp
    rand_sigma_inv = 0.0_wp
    nyz_inv        = 1.0_wp / REAL( ( nzt + 1 + mergp_z - ( nzb - mergp_z ) + 1 )                  &
                     * ( ny + mergp_y - ( 0 - mergp_y ) + 1 ), KIND = wp )
!
!-- Compute and normalize random numbers.
    DO  j = nys - mergp_y, nyn + mergp_y
!
!--    Put the random seeds at grid point j
       CALL random_seed_parallel( put=seq_rand_yz(:,j) )
       DO  k = nzb - mergp_z, nzt + 1 + mergp_z
          CALL random_number_parallel( random_dummy )
          rand_it(k,j) = random_dummy
       ENDDO
!
!--    Get the new random seeds from last call at grid point j
       CALL random_seed_parallel( get=seq_rand_yz(:,j) )
    ENDDO
!
!-- For normalization to zero mean, sum-up the global random numers. To normalize the global set of
!-- random numbers, the inner ghost layers mergp must not be summed-up, else the random numbers on
!-- the ghost layers will be stronger weighted as they also occur on the inner subdomains.
    DO  j = MERGE( nys, nys - mergp_y, nys /= 0 ), MERGE( nyn, nyn + mergp_y, nyn /= ny )
       DO  k = nzb - mergp_z, nzt + 1 + mergp_z
          rand_av = rand_av + rand_it(k,j)
       ENDDO
    ENDDO

#if defined( __parallel )
!
!-- Sum-up the local averages of the random numbers
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, rand_av, 1, MPI_REAL, MPI_SUM, comm1dy, ierr )
#endif
    rand_av = rand_av * nyz_inv
!
!-- Obtain zero mean
    rand_it= rand_it - rand_av
!
!-- Now, compute the variance
    DO  j = MERGE( nys, nys - mergp_y, nys /= 0 ), MERGE( nyn, nyn + mergp_y, nyn /= ny )
       DO  k = nzb - mergp_z, nzt + 1 + mergp_z
          rand_sigma_inv = rand_sigma_inv + rand_it(k,j)**2
       ENDDO
    ENDDO

#if defined( __parallel )
!
!-- Sum-up the local quadratic averages of the random numbers
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, rand_sigma_inv, 1, MPI_REAL, MPI_SUM, comm1dy, ierr )
#endif
!
!-- Compute standard deviation
    IF ( rand_sigma_inv /= 0.0_wp )  THEN
       rand_sigma_inv = 1.0_wp / SQRT( rand_sigma_inv * nyz_inv )
    ELSE
       rand_sigma_inv = 1.0_wp
    ENDIF
!
!-- Normalize with standard deviation to obtain unit variance
    rand_it = rand_it * rand_sigma_inv

    CALL cpu_log( log_point_s(31), 'STG f_n factors', 'start' )
!
!-- Generate velocity seed following Eq.6 of Xie and Castro (2008). There are two options. In the
!-- first one, the computation of the seeds is distributed to all processes along the communicator
!-- comm1dy and gathered on the leftmost and, if necessary, on the rightmost process. For huge
!-- length scales the computational effort can become quite huge (it scales with the turbulent
!-- length scales), so that gain by parallelization exceeds the costs by the subsequent
!-- communication. In the second option, which performs better when the turbulent length scales
!-- are parametrized and thus the loops are smaller, the seeds are computed locally and no
!-- communication is necessary.
    IF ( compute_velocity_seeds_local )  THEN

       f_n  = 0.0_wp
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
             DO  jj = -n_y(k), n_y(k)
                DO  kk = -n_z(k), n_z(k)
                   f_n(k,j) = f_n(k,j) + b_y(jj,k) * b_z(kk,k) * rand_it(k+kk,j+jj)
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSE

       f_n_l  = 0.0_wp
       DO  j = nys, nyn
          DO  k = nzb_x_stg, nzt_x_stg+1
             DO  jj = -n_y(k), n_y(k)
                DO  kk = -n_z(k), n_z(k)
                   f_n_l(k,j) = f_n_l(k,j) + b_y(jj,k) * b_z(kk,k) * rand_it(k+kk,j+jj)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!
!--    Gather velocity seeds of full subdomain
       send_count = nzt_x_stg - nzb_x_stg + 1
       IF ( nzt_x_stg == nzt )  send_count = send_count + 1

#if defined( __parallel )
!
!--    Gather the velocity seed matrix on left boundary mpi ranks.
       CALL MPI_GATHERV( f_n_l(nzb_x_stg,nys), send_count, stg_type_yz_small, f_n(nzb+1,nys),      &
                         recv_count_yz, displs_yz, stg_type_yz, id_left, comm1dx, ierr )
!
!--    If required, gather the same velocity seed matrix on right boundary mpi ranks (in offline
!--    nesting for example).
       IF ( PRESENT( id_right ) )  THEN
          CALL MPI_GATHERV( f_n_l(nzb_x_stg,nys), send_count, stg_type_yz_small, f_n(nzb+1,nys),   &
                            recv_count_yz, displs_yz, stg_type_yz, id_right, comm1dx, ierr )
       ENDIF
#else
       f_n(nzb+1:nzt+1,nys:nyn) = f_n_l(nzb_x_stg:nzt_x_stg+1,nys:nyn)
!
!--    Next line required to avoid compile errors because of unused dummy arguments
       IF ( id_left == 0 )  id_right = 0
#endif

    ENDIF

    DEALLOCATE( rand_it )

    CALL cpu_log( log_point_s(31), 'STG f_n factors', 'stop' )

 END SUBROUTINE stg_generate_seed_yz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Generate a set of random number rand_it wich is equal on each PE and calculate the velocity seed
!> f_n.
!> f_n is splitted in vertical direction by the number of PEs in y-direction and and each PE
!> calculates a vertical subsection of f_n. At the the end, all parts are collected to form the
!> full array.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_generate_seed_xz( n_x, n_z, b_x, b_z, f_n, id_south, id_north )

    INTEGER(iwp) :: i           !< loop index in x-direction
    INTEGER(iwp) :: id_north    !< core ids at respective boundaries
    INTEGER(iwp) :: id_south    !< core ids at respective boundaries
    INTEGER(iwp) :: ii          !< loop index in x-direction
    INTEGER(iwp) :: k           !< loop index in z-direction
    INTEGER(iwp) :: kk          !< loop index in z-direction
    INTEGER(iwp) :: send_count  !< send count for MPI_GATHERV

    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_x    !< length scale in x-direction
    INTEGER(iwp), DIMENSION(nzb:nzt+1) :: n_z    !< length scale in z-direction

    REAL(wp) :: nxz_inv         !< inverse of number of grid points in xz-slice
    REAL(wp) :: rand_av         !< average of random number
    REAL(wp) :: rand_sigma_inv  !< inverse of stdev of random number

    REAL(wp), DIMENSION(-mergp_x:mergp_x,nzb:nzt+1) :: b_x     !< filter function in x-direction
    REAL(wp), DIMENSION(-mergp_z:mergp_z,nzb:nzt+1) :: b_z     !< filter function in z-direction

    REAL(wp), DIMENSION(nzb_y_stg:nzt_y_stg+1,nxl:nxr) :: f_n_l   !<  local velocity seed
    REAL(wp), DIMENSION(nzb:nzt+1,nxl:nxr)             :: f_n     !<  velocity seed

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rand_it   !< global array of random numbers

!
!-- Generate random numbers using the parallel random generator. The set of random numbers are
!-- modified to have an average of 0 and unit variance. Note, at the end the random number array
!-- must be defined globally in order to compute the correlation matrices. However, the calculation
!-- and normalization of random numbers is done locally, while the result is later distributed to
!-- all processes. Further, please note, a set of random numbers is only calculated for the south
!-- boundary, while the north boundary uses the same random numbers and thus also computes the same
!-- correlation matrix.
    ALLOCATE( rand_it(nzb-mergp_z:nzt+1+mergp_z,-mergp_x+nxl:nxr+mergp_x) )
    rand_it = 0.0_wp

    rand_av        = 0.0_wp
    rand_sigma_inv = 0.0_wp
    nxz_inv        = 1.0_wp / REAL( ( nzt + 1 + mergp_z - ( nzb - mergp_z ) + 1 )                  &
                     * ( nx + mergp_x - ( 0 - mergp_x ) +1 ), KIND = wp )
!
!-- Compute and normalize random numbers.
    DO  i = nxl - mergp_x, nxr + mergp_x
!
!--    Put the random seeds at grid point ii
       CALL random_seed_parallel( put=seq_rand_xz(:,i) )
       DO  k = nzb - mergp_z, nzt + 1 + mergp_z
          CALL random_number_parallel( random_dummy )
          rand_it(k,i) = random_dummy
       ENDDO
!
!--    Get the new random seeds from last call at grid point ii
       CALL random_seed_parallel( get=seq_rand_xz(:,i) )
    ENDDO
!
!-- For normalization to zero mean, sum-up the global random numers.
!-- To normalize the global set of random numbers, the inner ghost layers mergp must not be
!-- summed-up, else the random numbers on the ghost layers will be stronger weighted as they
!-- also occur on the inner subdomains.
    DO  i = MERGE( nxl, nxl - mergp_x, nxl /= 0 ), MERGE( nxr, nxr + mergp_x, nxr /= nx )
       DO  k = nzb - mergp_z, nzt + 1 + mergp_z
          rand_av = rand_av + rand_it(k,i)
       ENDDO
    ENDDO

#if defined( __parallel )
!
!-- Sum-up the local averages of the random numbers
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, rand_av, 1, MPI_REAL, MPI_SUM, comm1dx, ierr )
#endif
    rand_av = rand_av * nxz_inv
!
!-- Obtain zero mean
    rand_it= rand_it - rand_av
!
!-- Now, compute the variance
    DO  i = MERGE( nxl, nxl - mergp_x, nxl /= 0 ), MERGE( nxr, nxr + mergp_x, nxr /= nx )
       DO  k = nzb - mergp_z, nzt + 1 + mergp_z
          rand_sigma_inv = rand_sigma_inv + rand_it(k,i)**2
       ENDDO
    ENDDO

#if defined( __parallel )
!
!-- Sum-up the local quadratic averages of the random numbers
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, rand_sigma_inv, 1, MPI_REAL, MPI_SUM, comm1dx, ierr )
#endif
!
!-- Compute standard deviation
    IF ( rand_sigma_inv /= 0.0_wp )  THEN
       rand_sigma_inv = 1.0_wp / SQRT( rand_sigma_inv * nxz_inv )
    ELSE
       rand_sigma_inv = 1.0_wp
    ENDIF
!
!-- Normalize with standard deviation to obtain unit variance
    rand_it = rand_it * rand_sigma_inv

    CALL cpu_log( log_point_s(31), 'STG f_n factors', 'start' )
!
!-- Generate velocity seed following Eq.6 of Xie and Castro (2008). There are two options. In the
!-- first one, the computation of the seeds is distributed to all processes along the communicator
!-- comm1dx and gathered on the southmost and, if necessary, on the northmost process. For huge
!-- length scales the computational effort can become quite huge (it scales with the turbulent
!-- length scales), so that gain by parallelization exceeds the costs by the subsequent communication.
!-- In the second option, which performs better when the turbulent length scales are parametrized
!-- and thus the loops are smaller, the seeds are computed locally and no communication is necessary.
    IF ( compute_velocity_seeds_local )  THEN

       f_n  = 0.0_wp
       DO  i = nxl, nxr
          DO  k = nzb, nzt+1
             DO  ii = -n_x(k), n_x(k)
                DO  kk = -n_z(k), n_z(k)
                   f_n(k,i) = f_n(k,i) + b_x(ii,k) * b_z(kk,k) * rand_it(k+kk,i+ii)
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSE

       f_n_l  = 0.0_wp
       DO  i = nxl, nxr
          DO  k = nzb_y_stg, nzt_y_stg+1
             DO  ii = -n_x(k), n_x(k)
                DO  kk = -n_z(k), n_z(k)
                   f_n_l(k,i) = f_n_l(k,i) + b_x(ii,k) * b_z(kk,k) * rand_it(k+kk,i+ii)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!
!--    Gather velocity seeds of full subdomain
       send_count = nzt_y_stg - nzb_y_stg + 1
       IF ( nzt_y_stg == nzt )  send_count = send_count + 1

#if defined( __parallel )
!
!--    Gather the processed velocity seed on south boundary mpi ranks.
       CALL MPI_GATHERV( f_n_l(nzb_y_stg,nxl), send_count, stg_type_xz_small, f_n(nzb+1,nxl),      &
                         recv_count_xz, displs_xz, stg_type_xz, id_south, comm1dy, ierr )
!
!--    Gather the processed velocity seed on north boundary mpi ranks.
       CALL MPI_GATHERV( f_n_l(nzb_y_stg,nxl), send_count, stg_type_xz_small, f_n(nzb+1,nxl),      &
                         recv_count_xz, displs_xz, stg_type_xz, id_north, comm1dy, ierr )
#else
       f_n(nzb+1:nzt+1,nxl:nxr) = f_n_l(nzb_y_stg:nzt_y_stg+1,nxl:nxr)
!
!--    Next line required to avoid compile errors because of unused dummy arguments
       IF ( id_north == 0 )  id_south = 0
#endif

    ENDIF

    DEALLOCATE( rand_it )

    CALL cpu_log( log_point_s(31), 'STG f_n factors', 'stop' )

 END SUBROUTINE stg_generate_seed_xz

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parametrization of the Reynolds-stress componentes, turbulent length- and time scales. The
!> parametrization follows Brost et al. (1982) with modifications described in Rotach et al. (1996),
!> which is applied in state-of-the-art dispserion modelling.
!--------------------------------------------------------------------------------------------------!
  SUBROUTINE parametrize_turbulence

    INTEGER(iwp) :: k                 !< loop index in z-direction

    REAL(wp) ::  corr_term_uh         !< correction term in parametrization of horizontal variances for unstable stratification
    REAL(wp) ::  d_zi                 !< inverse boundary-layer depth
    REAL(wp) ::  length_scale_lon     !< length scale in flow direction
    REAL(wp) ::  length_scale_lon_zi  !< length scale in flow direction at boundary-layer top
    REAL(wp) ::  length_scale_lat     !< length scale in crosswind direction
    REAL(wp) ::  length_scale_lat_zi  !< length scale in crosswind direction at boundary-layer top
    REAL(wp) ::  length_scale_vert    !< length scale in vertical direction
    REAL(wp) ::  length_scale_vert_zi !< length scale in vertical direction at boundary-layer top
    REAL(wp) ::  time_scale_zi        !< time scale at boundary-layer top
    REAL(wp) ::  r11_zi               !< longitudinal variance at boundary-layer top
    REAL(wp) ::  r22_zi               !< crosswind variance at boundary-layer top
    REAL(wp) ::  r33_zi               !< vertical variance at boundary-layer top
    REAL(wp) ::  r31_zi               !< r31 at boundary-layer top
    REAL(wp) ::  r32_zi               !< r32 at boundary-layer top
    REAL(wp) ::  rlat1                !< first dummy argument for crosswind compoment of reynolds stress
    REAL(wp) ::  rlat2                !< second dummy argument for crosswind compoment of reynolds stress
    REAL(wp) ::  rlon1                !< first dummy argument for longitudinal compoment of reynolds stress
    REAL(wp) ::  rlon2                !< second dummy argument for longitudinal compoment of reynolds stress
    REAL(wp) ::  zzi                  !< ratio of z/zi

    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  cos_phi    !< cosine of angle between mean flow direction and x-axis
    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  phi        !< angle between mean flow direction and x-axis
    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  sin_phi    !< sine of angle between mean flow direction and x-axis

!
!-- Calculate the boundary-layer top index. Boundary-layer top is calculated by Richardson-bulk
!-- criterion according to Heinze et al. (2017).
    k_zi = MAX( 1, MINLOC( ABS( zu - zi_ribulk ), DIM = 1 ) )
    IF ( zu(k_zi) > zi_ribulk )  k_zi = k_zi - 1
    k_zi = MIN( k_zi, nzt )
!
!-- Calculate angle between flow direction and x-axis.
    DO  k = nzb+1, nzt + 1
       IF ( u_init(k) /= 0.0_wp )  THEN
          phi(k) = ATAN( v_init(k) / u_init(k) )
       ELSE
          phi(k) = 0.5 * pi
       ENDIF
!
!--    Pre-calculate sine and cosine.
       cos_phi(k) = COS( phi(k) )
       sin_phi(k) = SIN( phi(k) )
    ENDDO
!
!-- Parametrize Reynolds-stress components. Please note, parametrization is formulated for the
!-- stream- and spanwise components, while stream- and spanwise direction do not necessarily
!-- coincide with the grid axis. Hence, components are rotated after computation.
    d_zi = 1.0_wp / zi_ribulk
    DO  k = nzb+1, k_zi

       corr_term_uh = MERGE( 0.35_wp * ABS(                                                        &
                        - zi_ribulk / ( kappa * scale_l - 10E-4_wp )                               &
                                       )**( 2.0_wp / 3.0_wp ),                                     &
                          0.0_wp,                                                                  &
                          scale_l < -5.0_wp )
!
!--    Determine normalized height coordinate
       zzi = zu(k) * d_zi
!
!--    Compute longitudinal and crosswind component of reynolds stress. Rotate these components by
!--    the wind direction to map onto the u- and v-component. Note, since reynolds stress for
!--    variances cannot become negative, take the absolute values.
       rlon1 = scale_us**2 * ( corr_term_uh + 5.0_wp - 4.0_wp * zzi )
       rlat1 = scale_us**2 * ( corr_term_uh + 2.0_wp - zzi )

       r11(k) = ABS(   cos_phi(k) * rlon1 + sin_phi(k) * rlat1 )
       r22(k) = ABS( - sin_phi(k) * rlon1 + cos_phi(k) * rlat1 )
!
!--    w'w'
       r33(k) = scale_wm**2 * (                                                                    &
                   1.5_wp * zzi**( 2.0_wp / 3.0_wp ) * EXP( -2.0_wp * zzi )                        &
                 + ( 1.7_wp - zzi ) * ( scale_us / scale_wm )**2                                   &
                              )
!
!--    u'w' and v'w'. After calculation of the longitudinal and crosswind component
!--    these are projected along the x- and y-direction. Note, it is assumed that
!--    the flux within the boundary points opposite to the vertical gradient.
       rlon2 = scale_us**2 * ( zzi - 1.0_wp )
       rlat2 = scale_us**2 * ( 0.4 * zzi * ( 1.0_wp - zzi ) )

       r31(k) = SIGN( ABS(   cos_phi(k) * rlon2 + sin_phi(k) * rlat2 ),                            &
                   - ( u_init(k+1) - u_init(k-1) ) )
       r32(k) = SIGN( ABS( - sin_phi(k) * rlon2 + cos_phi(k) * rlat2 ),                            &
                   - ( v_init(k+1) - v_init(k-1) ) )
!
!--    For u'v' no parametrization exist so far. For simplicity assume a similar profile as
!--    for the vertical transport.
       r21(k) = 0.5_wp * ( r31(k) + r32(k) )
!
!--    Compute turbulent time scales according to Brost et al. (1982). Note, time scales are
!--    limited to the adjustment time scales.
       tu(k)  = MIN( dt_stg_adjust,                                                                &
                     3.33_wp * zzi * ( 1.0 - 0.67_wp * zzi ) / scale_wm * zi_ribulk )

       tv(k)  = tu(k)
       tw(k)  = tu(k)

       length_scale_lon  = MIN( SQRT( r11(k) ) * tu(k), zi_ribulk )
       length_scale_lat  = MIN( SQRT( r22(k) ) * tv(k), zi_ribulk )
       length_scale_vert = MIN( SQRT( r33(k) ) * tw(k), zi_ribulk )
!
!--    Assume isotropic turbulence length scales
       nux(k) = MAX( INT( length_scale_lon  * ddx     ), 1 )
       nuy(k) = MAX( INT( length_scale_lat  * ddy     ), 1 )
       nvx(k) = MAX( INT( length_scale_lon  * ddx     ), 1 )
       nvy(k) = MAX( INT( length_scale_lat  * ddy     ), 1 )
       nwx(k) = MAX( INT( length_scale_lon  * ddx     ), 1 )
       nwy(k) = MAX( INT( length_scale_lat  * ddy     ), 1 )
       nuz(k) = MAX( INT( length_scale_vert * ddzw(k) ), 1 )
       nvz(k) = MAX( INT( length_scale_vert * ddzw(k) ), 1 )
       nwz(k) = MAX( INT( length_scale_vert * ddzw(k) ), 1 )

    ENDDO
!
!-- Above boundary-layer top length- and timescales as well as reynolds-stress components are
!-- reduced to zero by a blending function.
    length_scale_lon_zi  = SQRT( r11(k_zi) ) * tu(k_zi)
    length_scale_lat_zi  = SQRT( r22(k_zi) ) * tv(k_zi)
    length_scale_vert_zi = SQRT( r33(k_zi) ) * tu(k_zi)

    time_scale_zi        = tu(k_zi)
    r11_zi               = r11(k_zi)
    r22_zi               = r22(k_zi)
    r33_zi               = r33(k_zi)
    r31_zi               = r31(k_zi)
    r32_zi               = r32(k_zi)

    d_l = blend_coeff / MAX( length_scale_vert_zi, dx, dy, MINVAL( dzw ) )

    DO  k = k_zi+1, nzt+1
!
!--    Calculate function to gradually decrease Reynolds stress above ABL top.
       blend = MIN( 1.0_wp, EXP( d_l * zu(k) - d_l * zi_ribulk ) )
!
!--    u'u' and v'v'. Assume isotropy. Note, add a small negative number to the denominator, else
!--    the mergpe-function can crash if scale_l is zero.
       r11(k) = r11_zi * blend
       r22(k) = r22_zi * blend
       r33(k) = r33_zi * blend
       r31(k) = r31_zi * blend
       r32(k) = r32_zi * blend
       r21(k) = 0.5_wp * ( r31(k) + r32(k) )

!
!--    Compute turbulent time scales according to Brost et al. (1982).
!--    Note, time scales are limited to the adjustment time scales.
       tu(k)  = time_scale_zi * blend
       tv(k)  = time_scale_zi * blend
       tw(k)  = time_scale_zi * blend

       length_scale_lon  = length_scale_lon_zi  * blend
       length_scale_lat  = length_scale_lat_zi  * blend
       length_scale_vert = length_scale_vert_zi * blend
!
!--    Assume isotropic turbulence length scales
       nux(k) = MAX( INT( length_scale_lon  * ddx ),     1 )
       nuy(k) = MAX( INT( length_scale_lat  * ddy ),     1 )
       nvx(k) = MAX( INT( length_scale_lon  * ddx ),     1 )
       nvy(k) = MAX( INT( length_scale_lat  * ddy ),     1 )
       nwx(k) = MAX( INT( length_scale_lon  * ddx ),     1 )
       nwy(k) = MAX( INT( length_scale_lat  * ddy ),     1 )
       nuz(k) = MAX( INT( length_scale_vert * ddzw(k) ), 1 )
       nvz(k) = MAX( INT( length_scale_vert * ddzw(k) ), 1 )
       nwz(k) = MAX( INT( length_scale_vert * ddzw(k) ), 1 )
    ENDDO
!
!-- Set bottom boundary condition for reynolds stresses
    r11(nzb) = r11(nzb+1)
    r22(nzb) = r22(nzb+1)
    r33(nzb) = r33(nzb+1)

    r21(nzb) = r11(nzb+1)
    r31(nzb) = r31(nzb+1)
    r32(nzb) = r32(nzb+1)

!
!-- Set bottom boundary condition for the length and time scales
    nux(nzb) = nux(nzb+1)
    nuy(nzb) = nuy(nzb+1)
    nuz(nzb) = nuz(nzb+1)
    nvx(nzb) = nvx(nzb+1)
    nvy(nzb) = nvy(nzb+1)
    nvz(nzb) = nvz(nzb+1)
    nwx(nzb) = nwx(nzb+1)
    nwy(nzb) = nwy(nzb+1)
    nwz(nzb) = nwz(nzb+1)

    tu(nzb)  = tu(nzb+1)
    tv(nzb)  = tv(nzb+1)
    tw(nzb)  = tw(nzb+1)

    adjustment_step = .TRUE.

 END SUBROUTINE parametrize_turbulence

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the coefficient matrix from the Lund rotation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_coeff_matrix

    INTEGER(iwp) ::  k  !< loop index in z-direction

!
!-- Calculate coefficient matrix. Split loops to allow for loop vectorization.
    DO  k = nzb+1, nzt+1
       IF ( r11(k) > 10E-6_wp )  THEN
          a11(k) = SQRT( r11(k) )
          a21(k) = r21(k) / a11(k)
          a31(k) = r31(k) / a11(k)
       ELSE
          a11(k) = 10E-8_wp
          a21(k) = 10E-8_wp
          a31(k) = 10E-8_wp
       ENDIF
    ENDDO
    DO  k = nzb+1, nzt+1
       a22(k) = r22(k) - a21(k)**2
       IF ( a22(k) > 10E-6_wp )  THEN
          a22(k) = SQRT( a22(k) )
          a32(k) = ( r32(k) - a21(k) * a31(k) ) / a22(k)
       ELSE
          a22(k) = 10E-8_wp
          a32(k) = 10E-8_wp
       ENDIF
    ENDDO
    DO  k = nzb+1, nzt+1
       a33(k) = r33(k) - a31(k)**2 - a32(k)**2
       IF ( a33(k) > 10E-6_wp )  THEN
          a33(k) =  SQRT( a33(k) )
       ELSE
          a33(k) = 10E-8_wp
       ENDIF
    ENDDO
!
!-- Set bottom boundary condition
    a11(nzb) = a11(nzb+1)
    a22(nzb) = a22(nzb+1)
    a21(nzb) = a21(nzb+1)
    a33(nzb) = a33(nzb+1)
    a31(nzb) = a31(nzb+1)
    a32(nzb) = a32(nzb+1)

 END SUBROUTINE calc_coeff_matrix

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine controls the re-adjustment of the turbulence statistics used for generating
!> turbulence at the lateral boundaries.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stg_adjust

    IF ( debug_output_timestep )  CALL debug_message( 'stg_adjust', 'start' )
!
!-- In case of dirichlet inflow boundary conditions only at one lateral boundary, i.e. in the case
!-- no offline or self nesting is applied but synthetic turbulence shall be parametrized
!-- nevertheless, the boundary-layer depth need to determined first.
    IF ( .NOT. nesting_offline  .AND.  .NOT.  child_domain )                                       &
       CALL nesting_offl_calc_zi
!
!-- Compute scaling parameters (domain-averaged), such as friction velocity are calculated.
    CALL calc_scaling_variables
!
!-- Parametrize Reynolds-stress tensor. Parametrization follows Brost et al. (1982) with
!-- modifications described in Rotach et al. (1996) and is based on boundary-layer depth, friction
!-- velocity and velocity scale.
    CALL parametrize_turbulence
!
!-- Calculate coefficient matrix from Reynolds stress tensor
!-- (Lund rotation)
    CALL calc_coeff_matrix
!
!-- Setup filter functions according to updated length scales.
    CALL stg_setup_filter_function
!
!-- Reset time counter for controlling next adjustment to zero
    time_stg_adjust = 0.0_wp

    IF ( debug_output_timestep )  CALL debug_message( 'stg_adjust', 'end' )

 END SUBROUTINE stg_adjust


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate scaling variables which are used for turbulence parametrization according to Rotach
!> et al. (1996). Scaling variables are: friction velocity, boundary-layer depth, momentum velocity
!> scale, and Obukhov length.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_scaling_variables

    INTEGER(iwp) ::  i  !< loop index in x-direction
    INTEGER(iwp) ::  j  !< loop index in y-direction
    INTEGER(iwp) ::  k  !< loop index in z-direction
    INTEGER(iwp) ::  m  !< surface element index

    REAL(wp) ::  friction_vel_l  !< mean friction veloctiy on subdomain
    REAL(wp) ::  pt_surf_mean    !< mean near surface temperature (at 1st grid point)
    REAL(wp) ::  pt_surf_mean_l  !< mean near surface temperature (at 1st grid point) on subdomain
    REAL(wp) ::  scale_l_l       !< mean Obukhov lenght on subdomain
    REAL(wp) ::  shf_mean        !< mean surface sensible heat flux
    REAL(wp) ::  shf_mean_l      !< mean surface sensible heat flux on subdomain
    REAL(wp) ::  w_convective    !< convective velocity scale

!
!-- Calculate mean friction velocity, velocity scale, heat flux and near-surface temperature in the
!-- model domain.
!     pt_surf_mean_l = 0.0_wp
!     shf_mean_l     = 0.0_wp
!     scale_l_l      = 0.0_wp
!     friction_vel_l = 0.0_wp
!     DO  m = 1, surf_def_h(0)%ns
!        i = surf_def_h(0)%i(m)
!        j = surf_def_h(0)%j(m)
!        k = surf_def_h(0)%k(m)
!        friction_vel_l = friction_vel_l  + surf_def_h(0)%us(m)
!        shf_mean_l     = shf_mean_l      + surf_def_h(0)%shf(m) * drho_air(k)
!        scale_l_l      = scale_l_l       + surf_def_h(0)%ol(m)
!        pt_surf_mean_l = pt_surf_mean_l  + pt(k,j,i)
!     ENDDO
!     DO  m = 1, surf_lsm_h(0)%ns
!        i = surf_lsm_h(0)%i(m)
!        j = surf_lsm_h(0)%j(m)
!        k = surf_lsm_h(0)%k(m)
!        friction_vel_l = friction_vel_l  + surf_lsm_h(0)%us(m)
!        shf_mean_l     = shf_mean_l      + surf_lsm_h(0)%shf(m) * drho_air(k)
!        scale_l_l      = scale_l_l       + surf_lsm_h(0)%ol(m)
!        pt_surf_mean_l = pt_surf_mean_l  + pt(k,j,i)
!     ENDDO
!     DO  m = 1, surf_usm_h(0)%ns
!        i = surf_usm_h(0)%i(m)
!        j = surf_usm_h(0)%j(m)
!        k = surf_usm_h(0)%k(m)
!        friction_vel_l = friction_vel_l  + surf_usm_h(0)%us(m)
!        shf_mean_l     = shf_mean_l      + surf_usm_h(0)%shf(m) * drho_air(k)
!        scale_l_l      = scale_l_l       + surf_usm_h(0)%ol(m)
!        pt_surf_mean_l = pt_surf_mean_l  + pt(k,j,i)
!     ENDDO
!
! #if defined( __parallel )
!     CALL MPI_ALLREDUCE( friction_vel_l, scale_us,     1, MPI_REAL, MPI_SUM, comm2d, ierr )
!     CALL MPI_ALLREDUCE( shf_mean_l, shf_mean,         1, MPI_REAL, MPI_SUM, comm2d, ierr )
!     CALL MPI_ALLREDUCE( scale_l_l, scale_l,           1, MPI_REAL, MPI_SUM, comm2d, ierr )
!     CALL MPI_ALLREDUCE( pt_surf_mean_l, pt_surf_mean, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
! #else
!     scale_us     = friction_vel_l
!     shf_mean     = shf_mean_l
!     scale_l      = scale_l_l
!     pt_surf_mean = pt_surf_mean_l
! #endif

    scale_us     = 0.3_wp !( ( hom(topo_min_level+4,19,pr_palm,0) / momentumflux_output_conversion(topo_min_level+4) )**2 &
!                    + ( hom(topo_min_level+4,20,pr_palm,0) / momentumflux_output_conversion(topo_min_level+4) )**2 &
!                    )**0.25_wp ! scale_us     * d_nxy
    shf_mean     = hom(topo_min_level+4,1,18,0) / heatflux_output_conversion(topo_min_level+4) ! shf_mean     * d_nxy
!     scale_l      = scale_l      * d_nxy
    pt_surf_mean = hom(topo_min_level+4,1,4,0)  !pt_surf_mean * d_nxy
!
!-- Compute mean convective velocity scale. Note, in case the mean heat flux is negative, set
!-- convective velocity scale to zero.
    IF ( shf_mean > 0.0_wp )  THEN
       w_convective = ( g * shf_mean * zi_ribulk / pt_surf_mean )**( 1.0_wp / 3.0_wp )
       scale_l = - scale_us**3 * pt_surf_mean / (0.4_wp * 9.81_wp * shf_mean )
    ELSE
       w_convective = 0.0_wp
       scale_l = 0.1_wp
    ENDIF
!
!-- At the initial run the friction velocity is initialized close to zero, leading to almost zero
!-- disturbances at the boundaries. In order to obtain disturbances nevertheless, set a minium
!-- value of friction velocity for the reynolds-stress parametrization.
    IF ( time_since_reference_point <= 0.0_wp )  scale_us = MAX( scale_us, 0.05_wp )
!
!-- Finally, in order to consider also neutral or stable stratification, compute momentum velocity
!-- scale from u* and convective velocity scale, according to Rotach et al. (1996).
    scale_wm = ( scale_us**3 + 0.6_wp * w_convective**3 )**( 1.0_wp / 3.0_wp )

    write(9,*) "stg scales: ", scale_us, scale_wm, w_convective, scale_l, pt_surf_mean, shf_mean, &
                               zi_ribulk
    flush(9)

 END SUBROUTINE calc_scaling_variables

 END MODULE synthetic_turbulence_generator_mod
