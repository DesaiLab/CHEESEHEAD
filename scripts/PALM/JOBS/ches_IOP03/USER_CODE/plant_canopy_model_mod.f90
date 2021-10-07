!> @file plant_canopy_model_mod.f90
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
! Copyright 2017-2021 Institute of Computer Science of the
!                     Czech Academy of Sciences, Prague
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> 1) Initialization of the canopy model, e.g. construction of leaf area density profile
!>    (subroutine pcm_init).
!> 2) Calculation of sinks and sources of momentum, heat and scalar concentration due to canopy
!>    elements (subroutine pcm_tendency).
!
! @todo - precalculate constant terms in pcm_calc_transpiration_rate
! @todo - unify variable names (pcm_, pc_, ...)
! @todo - get rid-off dependency on radiation model
!--------------------------------------------------------------------------------------------------!
 MODULE plant_canopy_model_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzu, dzw, e, exner, hyp, pt, q, s, tend, u, v, w, zu, zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  c_p, degc_to_k, l_v, lv_d_cp, r_d, rd_d_rv

    USE bulk_cloud_model_mod,                                                                      &
        ONLY: bulk_cloud_model, microphysics_seifert

    USE control_parameters,                                                                        &
        ONLY: average_count_3d,                                                                    &
              coupling_char,                                                                       &
              debug_output,                                                                        &
              dt_3d,                                                                               &
              dz,                                                                                  &
              humidity,                                                                            &
              land_surface,                                                                        &
              length,                                                                              &
              message_string,                                                                      &
              ocean_mode,                                                                          &
              passive_scalar,                                                                      &
              plant_canopy,                                                                        &
              restart_data_format_output,                                                          &
              restart_string,                                                                      &
              urban_surface

    USE grid_variables,                                                                            &
        ONLY:  dx, dy

    USE indices,                                                                                   &
        ONLY:  nbgp, nxl, nxlg, nxlu, nxr, nxrg, nyn, nyng, nys, nysg, nysv, nz, nzb, nzt,         &
               topo_top_ind, topo_flags

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  char_fill,                                                                          &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               input_file_static,                                                                  &
               input_pids_static,                                                                  &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               pids_id,                                                                            &
               real_3d,                                                                            &
               vars_pids

    USE pegrid

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io,                                                                         &
               wrd_mpi_io

    USE surface_mod,                                                                               &
        ONLY: surf_def_h, surf_lsm_h, surf_usm_h


    IMPLICIT NONE

    CHARACTER (LEN=30) ::  canopy_mode = 'homogeneous'           !< canopy coverage
    INTEGER(iwp) ::  pch_index = 0                                 !< plant canopy height/top index

    INTEGER(iwp) ::  lad_vertical_gradient_level_ind(10) = -9999   !< lad-profile levels (index)

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  pch_index_ji     !< local plant canopy top

    LOGICAL ::  calc_beta_lad_profile = .FALSE.       !< switch for calc. of lad from beta func.
    LOGICAL ::  plant_canopy_transpiration = .FALSE.  !< flag to switch calculation of transpiration and corresponding latent heat
                                                      !< for resolved plant canopy inside radiation model
                                                      !< (calls subroutine pcm_calc_transpiration_rate from module plant_canopy_mod)

    REAL(wp) ::  alpha_lad = 9999999.9_wp                        !< coefficient for lad calculation
    REAL(wp) ::  beta_lad = 9999999.9_wp                         !< coefficient for lad calculation
    REAL(wp) ::  canopy_drag_coeff = 0.0_wp                      !< canopy drag coefficient (parameter)
    REAL(wp) ::  cthf = 0.0_wp                                   !< canopy top heat flux
    REAL(wp) ::  dt_plant_canopy = 0.0_wp                        !< timestep account. for canopy drag
    REAL(wp) ::  ext_coef = 0.6_wp                               !< extinction coefficient
    REAL(wp) ::  lad_surface = 0.0_wp                            !< lad surface value
    REAL(wp) ::  lad_type_coef(0:10) = 1.0_wp                    !< multiplicative coeficients for particular types
                                                                 !< of plant canopy (e.g. deciduous tree during winter)
    REAL(wp) ::  lad_vertical_gradient(10) = 0.0_wp              !< lad gradient
    REAL(wp) ::  lad_vertical_gradient_level(10) = -9999999.9_wp !< lad-prof. levels (in m)
    REAL(wp) ::  lai_beta = 0.0_wp                               !< leaf area index (lai) for lad calc.
    REAL(wp) ::  leaf_scalar_exch_coeff = 0.0_wp                 !< canopy scalar exchange coeff.
    REAL(wp) ::  leaf_surface_conc = 0.0_wp                      !< leaf surface concentration

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  lad            !< leaf area density
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pre_lad        !< preliminary lad

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  bad_s                    !< basal-area density on scalar-grid
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  cum_lai_hf               !< cumulative lai for heatflux calc.
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  lad_s                    !< lad on scalar-grid
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pcm_heating_rate         !< plant canopy heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pcm_heatrate_av          !< array for averaging plant canopy sensible heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pcm_latent_rate          !< plant canopy latent heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pcm_latentrate_av        !< array for averaging plant canopy latent heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pcm_transpiration_rate   !< plant canopy transpiration rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pcm_transpirationrate_av !< array for averaging plant canopy transpiration rate

    TYPE(real_3d) ::  basal_area_density_f    !< input variable for basal area density - resolved vegetation
    TYPE(real_3d) ::  leaf_area_density_f     !< input variable for leaf area density - resolved vegetation
    TYPE(real_3d) ::  root_area_density_lad_f !< input variable for root area density - resolved vegetation

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC pcm_calc_transpiration_rate,                                                            &
           pcm_check_data_output,                                                                  &
           pcm_check_parameters,                                                                   &
           pcm_3d_data_averaging,                                                                  &
           pcm_data_output_3d,                                                                     &
           pcm_define_netcdf_grid,                                                                 &
           pcm_header,                                                                             &
           pcm_init,                                                                               &
           pcm_parin,                                                                              &
           pcm_rrd_global,                                                                         &
           pcm_rrd_local,                                                                          &
           pcm_tendency,                                                                           &
           pcm_wrd_global,                                                                         &
           pcm_wrd_local

!
!-- Public variables and constants
    PUBLIC canopy_drag_coeff, pcm_heating_rate, pcm_transpiration_rate, pcm_latent_rate,           &
           canopy_mode, cthf, dt_plant_canopy, lad, lad_s, pch_index, plant_canopy_transpiration,  &
           pcm_heatrate_av, pcm_latentrate_av, pch_index_ji

    INTERFACE pcm_calc_transpiration_rate
       MODULE PROCEDURE pcm_calc_transpiration_rate
    END INTERFACE pcm_calc_transpiration_rate

    INTERFACE pcm_check_data_output
       MODULE PROCEDURE pcm_check_data_output
    END INTERFACE pcm_check_data_output

    INTERFACE pcm_check_parameters
       MODULE PROCEDURE pcm_check_parameters
    END INTERFACE pcm_check_parameters

    INTERFACE pcm_3d_data_averaging
       MODULE PROCEDURE pcm_3d_data_averaging
    END INTERFACE pcm_3d_data_averaging

    INTERFACE pcm_data_output_3d
       MODULE PROCEDURE pcm_data_output_3d
    END INTERFACE pcm_data_output_3d

    INTERFACE pcm_define_netcdf_grid
       MODULE PROCEDURE pcm_define_netcdf_grid
    END INTERFACE pcm_define_netcdf_grid

     INTERFACE pcm_header
       MODULE PROCEDURE pcm_header
    END INTERFACE pcm_header

    INTERFACE pcm_init
       MODULE PROCEDURE pcm_init
    END INTERFACE pcm_init

    INTERFACE pcm_parin
       MODULE PROCEDURE pcm_parin
    END INTERFACE pcm_parin

    INTERFACE pcm_read_plant_canopy_3d
       MODULE PROCEDURE pcm_read_plant_canopy_3d
    END INTERFACE pcm_read_plant_canopy_3d

    INTERFACE pcm_rrd_local
       MODULE PROCEDURE pcm_rrd_local_ftn
       MODULE PROCEDURE pcm_rrd_local_mpi
    END INTERFACE pcm_rrd_local

    INTERFACE pcm_rrd_global
       MODULE PROCEDURE pcm_rrd_global_ftn
       MODULE PROCEDURE pcm_rrd_global_mpi
    END INTERFACE pcm_rrd_global

    INTERFACE pcm_tendency
       MODULE PROCEDURE pcm_tendency
       MODULE PROCEDURE pcm_tendency_ij
    END INTERFACE pcm_tendency

    INTERFACE pcm_wrd_local
       MODULE PROCEDURE pcm_wrd_local
    END INTERFACE pcm_wrd_local

    INTERFACE pcm_wrd_global
       MODULE PROCEDURE pcm_wrd_global
    END INTERFACE pcm_wrd_global


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the plant canopy transpiration rate based on the Jarvis-Stewart with
!> parametrizations described in Daudet et al. (1999; Agricult. and Forest Meteorol. 97) and Ngao,
!> Adam and Saudreau (2017;  Agricult. and Forest Meteorol 237-238). Model functions f1-f4 were
!> adapted from Stewart (1998; Agric. and Forest. Meteorol. 43) instead, because they are valid for
!> broader intervals of values. Funcion f4 used in form present in van Wijk et al. (1998; Tree
!> Physiology 20).
!>
!> This subroutine is called from subroutine radiation_interaction after the calculation of
!> radiation in plant canopy boxes.
!> (arrays pcbinsw and pcbinlw).
!>
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_calc_transpiration_rate(i, j, k, kk, pcbsw, pcblw, pcbtr, pcblh)

!
!--  Input parameters
     INTEGER(iwp), INTENT(IN) ::  i, j, k, kk        !< indices of the pc gridbox
     REAL(wp), INTENT(IN)     ::  pcblw              !< lw radiation in gridbox (W)
     REAL(wp), INTENT(IN)     ::  pcbsw              !< sw radiation in gridbox (W)
     REAL(wp), INTENT(OUT)    ::  pcblh              !< latent heat from transpiration dT/dt (K/s)
     REAL(wp), INTENT(OUT)    ::  pcbtr              !< transpiration rate dq/dt (kg/kg/s)

!--  Variables and parameters for calculation of transpiration rate
     REAL(wp), PARAMETER ::  gama_psychr = 66.0_wp !< psychrometric constant (Pa/K)
     REAL(wp), PARAMETER ::  g_s_max = 0.01        !< maximum stomatal conductivity (m/s)
     REAL(wp), PARAMETER ::  m_soil = 0.4_wp       !< soil water content (needs to adjust or take from LSM)
     REAL(wp), PARAMETER ::  m_wilt = 0.01_wp      !< wilting point soil water content (needs to adjust or take from LSM)
     REAL(wp), PARAMETER ::  m_sat = 0.51_wp       !< saturation soil water content (needs to adjust or take from LSM)
     REAL(wp), PARAMETER ::  t2_min = 0.0_wp       !< minimal temperature for calculation of f2
     REAL(wp), PARAMETER ::  t2_max = 40.0_wp      !< maximal temperature for calculation of f2

     REAL(wp) ::  d_fact
     REAL(wp) ::  e_eq
     REAL(wp) ::  e_imp
     REAL(wp) ::  evapor_rate
     REAL(wp) ::  f1
     REAL(wp) ::  f2
     REAL(wp) ::  f3
     REAL(wp) ::  f4
     REAL(wp) ::  g_b
     REAL(wp) ::  g_s
     REAL(wp) ::  rad
     REAL(wp) ::  rswc
     REAL(wp) ::  sat_press
     REAL(wp) ::  sat_press_d
     REAL(wp) ::  temp
     REAL(wp) ::  v_lad
     REAL(wp) ::  vpd
     REAL(wp) ::  wind_speed

!
!--  Temperature (deg C)
     temp = pt(k,j,i) * exner(k) - degc_to_k
!
!--  Coefficient for conversion of radiation to grid to radiation to unit leaves surface
     v_lad = 1.0_wp / ( MAX( lad_s(kk,j,i), 1.0E-10_wp ) * dx * dy * dz(1) )
!
!--  Magnus formula for the saturation pressure (see Ngao, Adam and Saudreau (2017) eq. 1)
!--  There are updated formulas available, kept consistent with the rest of the parametrization
     sat_press = 610.8_wp * EXP( 17.27_wp * temp / ( temp + 237.3_wp ) )
!
!--  Saturation pressure derivative (derivative of the above)
     sat_press_d = sat_press * 17.27_wp * 237.3_wp / ( temp + 237.3_wp )**2
!
!--  Wind speed
     wind_speed = SQRT( ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) )**2 +                                &
                        ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) )**2 +                                &
                        ( 0.5_wp * ( w(k,j,i) + w(k-1,j,i) ) )**2 )
!
!--  Aerodynamic conductivity (Daudet et al. (1999) eq. 14
     g_b = 0.01_wp * wind_speed + 0.0071_wp
!
!--  Radiation flux per leaf surface unit
     rad = pcbsw * v_lad
!
!--  First function for calculation of stomatal conductivity (radiation dependency)
!--  Stewart (1988; Agric. and Forest. Meteorol. 43) eq. 17
     f1 = rad * ( 1000.0_wp + 42.1_wp ) / 1000.0_wp / ( rad + 42.1_wp )
!
!--  Second function for calculation of stomatal conductivity (temperature dependency)
!--  Stewart (1988; Agric. and Forest. Meteorol. 43) eq. 21
     f2 = MAX( t2_min, ( temp - t2_min ) * MAX( 0.0_wp, t2_max - temp )**( ( t2_max - 16.9_wp ) /  &
                                                                           ( 16.9_wp - t2_min ) )  &
               / ( ( 16.9_wp - t2_min ) * ( t2_max - 16.9_wp )**( ( t2_max - 16.9_wp ) /           &
                                                                  ( 16.9_wp - t2_min ) ) ) )
!
!--  Water pressure deficit
!--  Ngao, Adam and Saudreau (2017) eq. 6 but with water vapour partial pressure
     vpd = MAX( sat_press - q(k,j,i) * hyp(k) / rd_d_rv, 0._wp )
!
!--  Third function for calculation of stomatal conductivity (water pressure deficit dependency)
!--  Ngao, Adam and Saudreau (2017) Table 1, limited from below according to Stewart (1988)
!--  The coefficients of the linear dependence should better correspond to broad-leaved trees than
!--  the coefficients from Stewart (1988) which correspond to conifer trees.
     vpd = MIN( MAX( vpd, 770.0_wp ), 3820.0_wp )
     f3 = -2E-4_wp * vpd + 1.154_wp
!
!--  Fourth function for calculation of stomatal conductivity (soil moisture dependency)
!--  Residual soil water content
!--  van Wijk et al. (1998; Tree Physiology 20) eq. 7
!--  TODO - over LSM surface might be calculated from LSM parameters
     rswc = ( m_sat - m_soil ) / ( m_sat - m_wilt )
!
!--  van Wijk et al. (1998; Tree Physiology 20) eq. 5-6 (it is a reformulation of eq. 22-23 of
!--  Stewart(1988))
     f4 = MAX( 0.0_wp, MIN( 1.0_wp - 0.041_wp * EXP( 3.2_wp * rswc ), 1.0_wp - 0.041_wp ) )
!
!--  Stomatal conductivity
!--  Stewart (1988; Agric. and Forest. Meteorol. 43) eq. 12
!--  (notation according to Ngao, Adam and Saudreau (2017) and others)
     g_s = g_s_max * f1 * f2 * f3 * f4 + 1.0E-10_wp
!
!--  Decoupling factor
!--  Daudet et al. (1999) eq. 6
     d_fact = ( sat_press_d / gama_psychr + 2.0_wp ) /                                             &
              ( sat_press_d / gama_psychr + 2.0_wp + 2.0_wp * g_b / g_s )
!
!--  Equilibrium evaporation rate
!--  Daudet et al. (1999) eq. 4
     e_eq = ( pcbsw + pcblw ) * v_lad * sat_press_d /                                              &
            gama_psychr / ( sat_press_d / gama_psychr + 2.0_wp ) / l_v
!
!--  Imposed evaporation rate
!--  Daudet et al. (1999) eq. 5
     e_imp = r_d * pt(k,j,i) * exner(k) / hyp(k) * c_p * g_s * vpd / gama_psychr / l_v
!
!--  Evaporation rate
!--  Daudet et al. (1999) eq. 3
!--  (evaporation rate is limited to non-negative values)
     evapor_rate = MAX( d_fact * e_eq + ( 1.0_wp - d_fact ) * e_imp, 0.0_wp )
!
!--  Conversion of evaporation rate to q tendency in gridbox
!--  dq/dt = E * LAD * V_g / (rho_air * V_g)
     pcbtr = evapor_rate * r_d * pt(k,j,i) * exner(k) * lad_s(kk,j,i) / hyp(k)  !-- = dq/dt
!
!--  latent heat from evaporation
     pcblh = pcbtr * lv_d_cp  !-- = - dT/dt

 END SUBROUTINE pcm_calc_transpiration_rate


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for plant canopy model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_check_data_output( var, unit )

    CHARACTER (LEN=*) ::  unit  !<
    CHARACTER (LEN=*) ::  var   !<


    SELECT CASE ( TRIM( var ) )

       CASE ( 'pcm_heatrate' )
!
!--       Output of heatrate can be only done if it is explicitely set by cthf, or parametrized by
!--       absorption of radiation. The latter, however, is only available if radiation_interactions
!--       are on. Note, these are enabled if land-surface or urban-surface is switched-on. Using
!--       radiation_interactions_on directly is not possible since it belongs to the
!--       radition_model, which in turn depends on the plant-canopy model, creating circular
!--       dependencies.
          IF ( cthf == 0.0_wp  .AND.  ( .NOT. urban_surface .AND. .NOT. land_surface ) )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requi' //                         &
                              'res setting of parameter cthf /= 0.0'
             CALL message( 'pcm_check_data_output', 'PA0718', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'K s-1'

       CASE ( 'pcm_transpirationrate' )
          unit = 'kg kg-1 s-1'

       CASE ( 'pcm_latentrate' )
          unit = 'K s-1'

       CASE ( 'pcm_bad', 'pcm_lad' )
          unit = 'm2 m-3'


       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE pcm_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for plant canopy model
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_check_parameters

       IF ( ocean_mode )  THEN
          message_string = 'plant_canopy = .TRUE. is not allowed in the ocean'
          CALL message( 'pcm_check_parameters', 'PA0696', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( canopy_drag_coeff == 0.0_wp )  THEN
          message_string = 'plant_canopy = .TRUE. requires a non-zero drag ' //                    &
                           'coefficient & given value is canopy_drag_coeff = 0.0'
          CALL message( 'pcm_check_parameters', 'PA0041', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( ( alpha_lad /= 9999999.9_wp .AND. beta_lad == 9999999.9_wp )  .OR.                     &
            beta_lad /= 9999999.9_wp  .AND.  alpha_lad == 9999999.9_wp )  THEN
          message_string = 'using the beta function for the construction ' //                      &
                           'of the leaf area density profile requires '    //                      &
                           'both alpha_lad and beta_lad to be /= 9999999.9'
          CALL message( 'pcm_check_parameters', 'PA0118', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( calc_beta_lad_profile  .AND.  lai_beta == 0.0_wp )  THEN
          message_string = 'using the beta function for the construction ' //                      &
                           'of the leaf area density profile requires '    //                      &
                           'a non-zero lai_beta, but given value is '      //                      &
                           'lai_beta = 0.0'
          CALL message( 'pcm_check_parameters', 'PA0119', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( calc_beta_lad_profile  .AND.  lad_surface /= 0.0_wp )  THEN
          message_string = 'simultaneous setting of alpha_lad /= 9999999.9 '//                     &
                           'combined with beta_lad /= 9999999.9 '           //                     &
                           'and lad_surface /= 0.0 is not possible, '       //                     &
                           'use either vertical gradients or the beta '     //                     &
                           'function for the construction of the leaf area '//                     &
                           'density profile'
          CALL message( 'pcm_check_parameters', 'PA0120', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
          message_string = 'plant_canopy = .TRUE. requires cloud_scheme /= seifert_beheng'
          CALL message( 'pcm_check_parameters', 'PA0360', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE pcm_check_parameters


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_3d_data_averaging( mode, variable )

    CHARACTER (LEN=*) ::  mode    !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  k            !<


    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'pcm_heatrate' )
             IF ( .NOT. ALLOCATED( pcm_heatrate_av ) )  THEN
                ALLOCATE( pcm_heatrate_av(0:pch_index,nysg:nyng,nxlg:nxrg) )
             ENDIF
             pcm_heatrate_av = 0.0_wp


          CASE ( 'pcm_latentrate' )
             IF ( .NOT. ALLOCATED( pcm_latentrate_av ) )  THEN
                ALLOCATE( pcm_latentrate_av(0:pch_index,nysg:nyng,nxlg:nxrg) )
             ENDIF
             pcm_latentrate_av = 0.0_wp


          CASE ( 'pcm_transpirationrate' )
             IF ( .NOT. ALLOCATED( pcm_transpirationrate_av ) )  THEN
                ALLOCATE( pcm_transpirationrate_av(0:pch_index,nysg:nyng,nxlg:nxrg) )
             ENDIF
             pcm_transpirationrate_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSE IF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'pcm_heatrate' )
             IF ( ALLOCATED( pcm_heatrate_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      IF ( pch_index_ji(j,i) /= 0 )  THEN
                         DO  k = 0, pch_index_ji(j,i)
                            pcm_heatrate_av(k,j,i) = pcm_heatrate_av(k,j,i) +                      &
                                                     pcm_heating_rate(k,j,i)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'pcm_latentrate' )
             IF ( ALLOCATED( pcm_latentrate_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      IF ( pch_index_ji(j,i) /= 0 )  THEN
                         DO  k = 0, pch_index_ji(j,i)
                            pcm_latentrate_av(k,j,i) = pcm_latentrate_av(k,j,i) +                  &
                                                       pcm_latent_rate(k,j,i)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'pcm_transpirationrate' )
             IF ( ALLOCATED( pcm_transpirationrate_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      IF ( pch_index_ji(j,i) /= 0 )  THEN
                         DO  k = 0, pch_index_ji(j,i)
                            pcm_transpirationrate_av(k,j,i) = pcm_transpirationrate_av(k,j,i) +    &
                                                              pcm_transpiration_rate(k,j,i)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSE IF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'pcm_heatrate' )
             IF ( ALLOCATED( pcm_heatrate_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      IF ( pch_index_ji(j,i) /= 0 )  THEN
                         DO  k = 0, pch_index_ji(j,i)
                            pcm_heatrate_av(k,j,i) = pcm_heatrate_av(k,j,i)                        &
                                                     / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'pcm_latentrate' )
             IF ( ALLOCATED( pcm_latentrate_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      IF ( pch_index_ji(j,i) /= 0 )  THEN
                         DO  k = 0, pch_index_ji(j,i)
                            pcm_latentrate_av(k,j,i) = pcm_latentrate_av(k,j,i)                    &
                                                       / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'pcm_transpirationrate' )
             IF ( ALLOCATED( pcm_transpirationrate_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      IF ( pch_index_ji(j,i) /= 0 )  THEN
                         DO  k = 0, pch_index_ji(j,i)
                            pcm_transpirationrate_av(k,j,i) = pcm_transpirationrate_av(k,j,i)      &
                                                              / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE pcm_3d_data_averaging

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables.
!> Note, 3D plant-canopy output has it's own vertical output dimension, meaning that 3D output is
!> relative to the model surface now rather than at the actual grid point where the plant canopy is
!> located.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )

    CHARACTER (LEN=*) ::  variable !< treated variable

    INTEGER(iwp) ::  av     !< flag indicating instantaneous or averaged data output
    INTEGER(iwp) ::  i      !< grid index x-direction
    INTEGER(iwp) ::  j      !< grid index y-direction
    INTEGER(iwp) ::  k      !< grid index z-direction
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found  !< flag indicating if variable is found

    REAL(wp)     ::  fill_value !< fill value

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< data output array


    found = .TRUE.

    local_pf = REAL( fill_value, KIND = 4 )

    SELECT CASE ( TRIM( variable ) )
!
!--    Note, to save memory arrays for heating are allocated from 0:pch_index.
!--    Thus, output must be relative to these array indices. Further, check whether the output is
!--    within the vertical output range, i.e. nzb_do:nzt_do, which is necessary as local_pf is only
!--    allocated for this index space. Note, plant-canopy output has a separate vertical output
!--    coordinate zlad, so that output is mapped down to the surface.
       CASE ( 'pcm_heatrate' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = pcm_heating_rate(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = pcm_heatrate_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'pcm_latentrate' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = pcm_latent_rate(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = pcm_latentrate_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'pcm_transpirationrate' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = pcm_transpiration_rate(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = pcm_transpirationrate_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'pcm_lad' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = lad_s(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'pcm_bad' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = MAX( 1, nzb_do ), MIN( pch_index_ji(j,i), nzt_do )
                      local_pf(i,j,k) = bad_s(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE pcm_data_output_3d

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

     CHARACTER (LEN=*), INTENT(IN)  ::  var         !<
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
     CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

     LOGICAL, INTENT(OUT)           ::  found       !<


     found  = .TRUE.

!
!--  Check for the grid. zpc is zu(nzb:nzb+pch_index)
     SELECT CASE ( TRIM( var ) )

        CASE ( 'pcm_heatrate', 'pcm_bad', 'pcm_lad', 'pcm_transpirationrate', 'pcm_latentrate' )
           grid_x = 'x'
           grid_y = 'y'
           grid_z = 'zpc'

        CASE DEFAULT
           found  = .FALSE.
           grid_x = 'none'
           grid_y = 'none'
           grid_z = 'none'
     END SELECT

 END SUBROUTINE pcm_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for plant canopy model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_header ( io )

    CHARACTER (LEN=10) ::  coor_chr            !<
    CHARACTER (LEN=86) ::  coordinates         !<
    CHARACTER (LEN=86) ::  gradients           !<
    CHARACTER (LEN=86) ::  leaf_area_density   !<
    CHARACTER (LEN=86) ::  slices              !<

    INTEGER(iwp)              ::  i   !<
    INTEGER(iwp),  INTENT(IN) ::  io  !< Unit of the output file
    INTEGER(iwp)              ::  k   !<

    REAL(wp) ::  canopy_height       !< canopy height (in m)


    canopy_height = zw(pch_index)

    WRITE( io, 1 )  canopy_mode, canopy_height, pch_index, canopy_drag_coeff
    IF ( passive_scalar )  THEN
       WRITE( io, 2 )  leaf_scalar_exch_coeff, leaf_surface_conc
    ENDIF

!
!-- Heat flux at the top of vegetation
    WRITE( io, 3 )  cthf

!
!-- Leaf area density profile, calculated either from given vertical gradients or from beta
!-- probability density function.
    IF ( .NOT. calc_beta_lad_profile )  THEN

!      Building output strings, starting with surface value
       WRITE( leaf_area_density, '(F7.4)' )  lad_surface
       gradients = '------'
       slices = '     0'
       coordinates = '   0.0'
       DO  i = 1, UBOUND( lad_vertical_gradient_level_ind, DIM=1 )
          IF  ( lad_vertical_gradient_level_ind(i) /= -9999 )  THEN

             WRITE( coor_chr, '(F7.2)' ) lad(lad_vertical_gradient_level_ind(i))
             leaf_area_density = TRIM( leaf_area_density ) // ' ' // TRIM( coor_chr )

             WRITE( coor_chr, '(F7.2)' ) lad_vertical_gradient(i)
             gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

             WRITE( coor_chr, '(I7)' ) lad_vertical_gradient_level_ind(i)
             slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

             WRITE( coor_chr, '(F7.1)' ) lad_vertical_gradient_level(i)
             coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )
          ELSE
             EXIT
          ENDIF
       ENDDO

       WRITE( io, 4 )  TRIM( coordinates ), TRIM( leaf_area_density ), TRIM( gradients ),          &
                       TRIM( slices )

    ELSE

       WRITE( leaf_area_density, '(F7.4)' )  lad_surface
       coordinates = '   0.0'

       DO  k = 1, pch_index

          WRITE( coor_chr,'(F7.2)' )  lad(k)
          leaf_area_density = TRIM( leaf_area_density ) // ' ' // TRIM( coor_chr )

          WRITE(coor_chr,'(F7.1)')  zu(k)
          coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

       ENDDO

       WRITE( io, 5 ) TRIM( coordinates ), TRIM( leaf_area_density ), alpha_lad, beta_lad, lai_beta

    ENDIF

1   FORMAT (/ /' Vegetation canopy (drag) model:' / ' ------------------------------' //           &
              ' Canopy mode: ', A / ' Canopy height: ', F6.2, 'm (',I4,' grid points)' /           &
              ' Leaf drag coefficient: ', F6.2 /)
2   FORMAT (/ ' Scalar exchange coefficient: ',F6.2 /                                              &
              ' Scalar concentration at leaf surfaces in kg/m**3: ', F6.2 /)
3   FORMAT ( ' Predefined constant heatflux at the top of the vegetation: ', F6.2, ' K m/s')
4   FORMAT (/ ' Characteristic levels of the leaf area density:' //                                &
              ' Height:              ', A, '  m' /                                                 &
              ' Leaf area density:   ', A, '  m**2/m**3' /                                         &
              ' Gradient:            ', A, '  m**2/m**4' /                                         &
              ' Gridpoint:           ', A )
5   FORMAT (//' Characteristic levels of the leaf area density and coefficients:' //               &
              ' Height:              ', A, '  m' /                                                 &
              ' Leaf area density:   ', A, '  m**2/m**3' /                                         &
              ' Coefficient alpha: ',F6.2 /                                                        &
              ' Coefficient beta: ',F6.2 /                                                         &
              ' Leaf area index: ',F6.2,'  m**2/m**2' /)

    END SUBROUTINE pcm_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the plant canopy model
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_init

       USE exchange_horiz_mod,                                                                     &
           ONLY:  exchange_horiz

       INTEGER(iwp) ::  i   !< running index
       INTEGER(iwp) ::  j   !< running index
       INTEGER(iwp) ::  k   !< running index
       INTEGER(iwp) ::  m   !< running index

       LOGICAL ::  lad_on_top = .FALSE.  !< dummy flag to indicate that LAD is defined on a building roof
       LOGICAL ::  bad_on_top = .FALSE.  !< dummy flag to indicate that BAD is defined on a building roof

       REAL(wp) ::  canopy_height   !< canopy height for lad-profile construction
       REAL(wp) ::  gradient        !< gradient for lad-profile construction
       REAL(wp) ::  int_bpdf        !< vertical integral for lad-profile construction
       REAL(wp) ::  lad_max         !< maximum LAD value in the model domain, used to perform a check


       IF ( debug_output )  CALL debug_message( 'pcm_init', 'start' )
!
!--    Allocate one-dimensional arrays for the computation of the leaf area density (lad) profile
       ALLOCATE( lad(0:nz+1), pre_lad(0:nz+1) )
       lad = 0.0_wp
       pre_lad = 0.0_wp

!
!--    Set flag that indicates that the lad-profile shall be calculated by using a beta probability
!--    density function
       IF ( alpha_lad /= 9999999.9_wp  .AND.  beta_lad /= 9999999.9_wp )  THEN
          calc_beta_lad_profile = .TRUE.
       ENDIF


!
!--    Compute the profile of leaf area density used in the plant canopy model. The profile can
!--    either be constructed from prescribed vertical gradients of the leaf area density or by using
!--    a beta probability density function (see e.g. Markkanen et al., 2003: Boundary-Layer
!--    Meteorology, 106, 437-459)
       IF ( .NOT. calc_beta_lad_profile )  THEN

!
!--       Use vertical gradients for lad-profile construction
          i = 1
          gradient = 0.0_wp

          lad(0) = lad_surface
          lad_vertical_gradient_level_ind(1) = 0

          DO  k = 1, pch_index
             IF ( i < 11 )  THEN
                IF ( lad_vertical_gradient_level(i) < zu(k)  .AND.                                 &
                     lad_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = lad_vertical_gradient(i)
                   lad_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   lad(k) = lad(k-1) + dzu(k) * gradient
                ELSE
                   lad(k) = lad_surface + dzu(k) * gradient
                ENDIF
             ELSE
                lad(k) = lad(k-1)
             ENDIF
          ENDDO

!
!--       In case of no given leaf area density gradients, choose a vanishing gradient. This
!--       information is used for the HEADER and the RUN_CONTROL file.
          IF ( lad_vertical_gradient_level(1) == -9999999.9_wp )  THEN
             lad_vertical_gradient_level(1) = 0.0_wp
          ENDIF

       ELSE

!
!--       Use beta function for lad-profile construction
          int_bpdf = 0.0_wp
          canopy_height = zw(pch_index)

          DO  k = 0, pch_index
             int_bpdf = int_bpdf +                                                                 &
                                 ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) ) *           &
                                   ( ( 1.0_wp - ( zw(k) / canopy_height ) )**( beta_lad-1.0_wp ) ) &
                                   * ( ( zw(k+1)-zw(k) ) / canopy_height ) )
          ENDDO

!
!--       Preliminary lad profile (defined on w-grid)
          DO  k = 0, pch_index
             pre_lad(k) =  lai_beta *                                                              &
                              ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) )                &
                                * ( ( 1.0_wp - ( zw(k) / canopy_height ) )**( beta_lad-1.0_wp ) )  &
                                / int_bpdf                                                         &
                              ) / canopy_height
          ENDDO

!
!--       Final lad profile (defined on scalar-grid level, since most prognostic quantities are
!--       defined there, hence, less interpolation is required when calculating the canopy
!--       tendencies)
          lad(0) = pre_lad(0)
          DO  k = 1, pch_index
             lad(k) = 0.5 * ( pre_lad(k-1) + pre_lad(k) )
          ENDDO

       ENDIF

!
!--    Allocate 3D-array for the leaf-area density (lad_s) as well as for basal-area densitiy
!--    (bad_s). Note, by default bad_s is zero.
       ALLOCATE( lad_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( bad_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       bad_s = 0.0_wp
!
!--    Initialization of the canopy coverage in the model domain:
!--    Setting the parameter canopy_mode = 'homogeneous' initializes a canopy, which fully covers
!--    the domain surface
       SELECT CASE ( TRIM( canopy_mode ) )

          CASE ( 'homogeneous' )

             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   lad_s(:,j,i) = lad(:)
                ENDDO
             ENDDO

          CASE ( 'read_from_file' )
!
!--          Read plant canopy
             IF ( input_pids_static )  THEN
!
!--             Open the static input file
#if defined( __netcdf )
                CALL open_read_file( TRIM( input_file_static ) //                                  &
                                     TRIM( coupling_char ),                                        &
                                     pids_id )

                CALL inquire_num_variables( pids_id, num_var_pids )
!
!--             Allocate memory to store variable names and read them
                ALLOCATE( vars_pids(1:num_var_pids) )
                CALL inquire_variable_names( pids_id, vars_pids )
!
!--             Read leaf area density - resolved vegetation
                IF ( check_existence( vars_pids, 'lad' ) )  THEN
                   leaf_area_density_f%from_file = .TRUE.
                   CALL get_attribute( pids_id, char_fill,                                         &
                                       leaf_area_density_f%fill,                                   &
                                       .FALSE., 'lad' )
!
!--                Inquire number of vertical vegetation layer
                   CALL get_dimension_length( pids_id,                                             &
                                              leaf_area_density_f%nz,                              &
                                              'zlad' )
!
!--                Allocate variable for leaf-area density
                   ALLOCATE( leaf_area_density_f%var                                               &
                                                 (0:leaf_area_density_f%nz-1,nys:nyn,nxl:nxr) )

                   CALL get_variable( pids_id, 'lad', leaf_area_density_f%var, nxl, nxr, nys, nyn, &
                                      0, leaf_area_density_f%nz-1 )

                ELSE
                   leaf_area_density_f%from_file = .FALSE.
                ENDIF
!
!--             Read basal area density - resolved vegetation
                IF ( check_existence( vars_pids, 'bad' ) )  THEN
                   basal_area_density_f%from_file = .TRUE.
                   CALL get_attribute( pids_id, char_fill,                                         &
                                       basal_area_density_f%fill,                                  &
                                       .FALSE., 'bad' )
!
!--                Inquire number of vertical vegetation layer
                   CALL get_dimension_length( pids_id,                                             &
                                              basal_area_density_f%nz,                             &
                                              'zlad' )
!
!--                Allocate variable
                   ALLOCATE( basal_area_density_f%var                                              &
                                                  (0:basal_area_density_f%nz-1,nys:nyn,nxl:nxr) )

                   CALL get_variable( pids_id, 'bad', basal_area_density_f%var, nxl, nxr, nys, nyn,&
                                      0, basal_area_density_f%nz-1 )
                ELSE
                   basal_area_density_f%from_file = .FALSE.
                ENDIF
!
!--             Read root area density - resolved vegetation
                IF ( check_existence( vars_pids, 'root_area_dens_r' ) )  THEN
                   root_area_density_lad_f%from_file = .TRUE.
                   CALL get_attribute( pids_id, char_fill,                                         &
                                       root_area_density_lad_f%fill,                               &
                                       .FALSE., 'root_area_dens_r' )
!
!--                Inquire number of vertical soil layers
                   CALL get_dimension_length( pids_id,                                             &
                                              root_area_density_lad_f%nz,                          &
                                              'zsoil' )
!
!--                Allocate variable
                   ALLOCATE( root_area_density_lad_f%var                                           &
                                                 (0:root_area_density_lad_f%nz-1,nys:nyn,nxl:nxr) )

                   CALL get_variable( pids_id, 'root_area_dens_r', root_area_density_lad_f%var,    &
                                      nxl, nxr, nys, nyn, 0, root_area_density_lad_f%nz-1 )
                ELSE
                   root_area_density_lad_f%from_file = .FALSE.
                ENDIF

                DEALLOCATE( vars_pids )
!
!--             Finally, close the input file and deallocate temporary array
                CALL close_input_file( pids_id )
#endif
             ENDIF

!
!--          Initialize LAD with data from file. If LAD is given in NetCDF file, use these values,
!--          else take LAD profiles from ASCII file.
!--          Please note, in NetCDF file LAD is only given up to the maximum canopy top, indicated
!--          by leaf_area_density_f%nz.
             lad_s = 0.0_wp
             IF ( leaf_area_density_f%from_file )  THEN
!
!--             Set also pch_index, used to be the upper bound of the vertical loops. Therefore, use
!--             the global top of the canopy layer.
                pch_index = leaf_area_density_f%nz - 1

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = 0, leaf_area_density_f%nz - 1
                         IF ( leaf_area_density_f%var(k,j,i) /=  leaf_area_density_f%fill )        &
                            lad_s(k,j,i) = leaf_area_density_f%var(k,j,i)
                      ENDDO
!
!--                   Check if resolved vegetation is mapped onto buildings.
!--                   In general, this is allowed and also meaningful, e.g. when trees carry across
!--                   roofs. However, due to the topography filtering, new building grid points can
!--                   emerge at locations where also plant canopy is defined. As a result, plant
!--                   canopy is mapped on top of roofs, with siginficant impact on the downstream
!--                   flow field and the nearby surface radiation. In order to avoid that plant
!--                   canopy is mistakenly mapped onto building roofs, check for building grid
!--                   points (bit 6) that emerge from the filtering (bit 4) and set LAD to zero at
!--                   these artificially created building grid points. This case, an informative
!--                   message is given.
                      IF ( ANY( lad_s(:,j,i) /= 0.0_wp )         .AND.                             &
                           ANY( BTEST( topo_flags(:,j,i), 6 ) )  .AND.                             &
                           ANY( BTEST( topo_flags(:,j,i), 4 ) ) )  THEN
                         lad_s(:,j,i) = 0.0_wp
                         lad_on_top   = .TRUE.
                      ENDIF
                   ENDDO
                ENDDO
#if defined( __parallel )
               CALL MPI_ALLREDUCE( MPI_IN_PLACE, lad_on_top, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
                IF ( lad_on_top )  THEN
                   WRITE( message_string, * )                                                      &
                                        'Resolved plant-canopy is defined on top of an ' //        &
                                        'artificially created building grid point(s) '//           &
                                        'the filtering) - LAD/BAD profile is omitted at this / ' //&
                                        'these grid point(s).'
                   CALL message( 'pcm_init', 'PA0313', 0, 0, 0, 6, 0 )
                ENDIF
                CALL exchange_horiz( lad_s, nbgp )
!
!            ASCII file
!--          Initialize canopy parameters canopy_drag_coeff, leaf_scalar_exch_coeff,
!--          leaf_surface_conc from file which contains complete 3D data (separate vertical profiles
!--          for each location).
             ELSE
                CALL pcm_read_plant_canopy_3d
             ENDIF
!
!--          Initialize LAD with data from file. If LAD is given in NetCDF file, use these values,
!--          else take LAD profiles from ASCII file.
!--          Please note, in NetCDF file LAD is only given up to the maximum canopy top, indicated
!--          by basal_area_density_f%nz.
             bad_s = 0.0_wp
             IF ( basal_area_density_f%from_file )  THEN
!
!--             Set also pch_index, used to be the upper bound of the vertical loops. Therefore, use
!--             the global top of the canopy layer.
                pch_index = basal_area_density_f%nz - 1

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = 0, basal_area_density_f%nz - 1
                         IF ( basal_area_density_f%var(k,j,i) /= basal_area_density_f%fill )       &
                            bad_s(k,j,i) = basal_area_density_f%var(k,j,i)
                      ENDDO
!
!--                   Check if resolved vegetation is mapped onto buildings.
!--                   Please see comment for leaf_area density
                      IF ( ANY( bad_s(:,j,i) /= 0.0_wp )         .AND.                             &
                           ANY( BTEST( topo_flags(:,j,i), 6 ) )  .AND.                             &
                           ANY( BTEST( topo_flags(:,j,i), 4 ) ) )  THEN
                         bad_s(:,j,i) = 0.0_wp
                         bad_on_top   = .TRUE.
                      ENDIF
                   ENDDO
                ENDDO
#if defined( __parallel )
               CALL MPI_ALLREDUCE( MPI_IN_PLACE, bad_on_top, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
                IF ( bad_on_top )  THEN
                   WRITE( message_string, * )                                                      &
                                        'Resolved plant-canopy is defined on top of an ' //        &
                                        'artificially created building grid point(s) '//           &
                                        'the filtering) - LAD/BAD profile is omitted at this / ' //&
                                        'these grid point(s).'
                   CALL message( 'pcm_init', 'PA0313', 0, 0, 0, 6, 0 )
                ENDIF
                CALL exchange_horiz( bad_s, nbgp )
             ENDIF

          CASE DEFAULT
!
!--          The DEFAULT case is reached either if the parameter canopy mode contains a wrong
!--          character string or if the user has coded a special case in the user interface.
!--          There, the subroutine user_init_plant_canopy checks which of these two conditions
!--          applies.
             CALL user_init_plant_canopy

       END SELECT
!
!--    Check that at least one grid point has an LAD /= 0, else this may cause errors in the
!--    radiation model.
       lad_max = MAXVAL( lad_s )
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, lad_max, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
#endif
       IF ( lad_max <= 0.0_wp )  THEN
          message_string = 'Plant-canopy model is switched-on but no ' //                          &
                           'plant canopy is present in the model domain.'
          CALL message( 'pcm_init', 'PA0685', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Initialize 2D index array indicating canopy top index.
       ALLOCATE( pch_index_ji(nysg:nyng,nxlg:nxrg) )
       pch_index_ji = 0

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = 0, pch_index
                IF ( lad_s(k,j,i) /= 0.0_wp  .OR.  bad_s(k,j,i) /= 0.0_wp )  pch_index_ji(j,i) = k
             ENDDO
!
!--          Check whether topography and local vegetation on top exceed height of the model domain.
             IF ( topo_top_ind(j,i,0) + pch_index_ji(j,i) >= nzt + 1 )  THEN
                message_string =  'Local vegetation height on top of ' //                          &
                                  'topography exceeds height of model domain.'
                CALL message( 'pcm_init', 'PA0674', 2, 2, myid, 6, 0 )
             ENDIF

          ENDDO
       ENDDO
!
!--    Calculate global pch_index value (index of top of plant canopy from ground)
       pch_index = MAXVAL( pch_index_ji )
!
!--    Exchange pch_index from all processors
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, pch_index, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr)
#endif
!
!--    Allocation of arrays pcm_heating_rate, pcm_transpiration_rate and pcm_latent_rate
       ALLOCATE( pcm_heating_rate(0:pch_index,nysg:nyng,nxlg:nxrg) )
       pcm_heating_rate = 0.0_wp

       IF ( humidity )  THEN
          ALLOCATE( pcm_transpiration_rate(0:pch_index,nysg:nyng,nxlg:nxrg) )
          pcm_transpiration_rate = 0.0_wp
          ALLOCATE( pcm_latent_rate(0:pch_index,nysg:nyng,nxlg:nxrg) )
          pcm_latent_rate = 0.0_wp
       ENDIF
!
!--    Initialization of the canopy heat source distribution due to heating of the canopy layers by
!--    incoming solar radiation, in case that a non-zero
!--    value is set for the canopy top heat flux (cthf), which equals the available net radiation at
!--    canopy top.
!--    The heat source distribution is calculated by a decaying exponential function of the downward
!--    cumulative leaf area index (cum_lai_hf), assuming that the foliage inside the plant canopy is
!--    heated by solar radiation penetrating the canopy layers according to the distribution of net
!--    radiation as suggested by Brown & Covey (1966; Agric. Meteorol. 3, 73–96). This approach has
!--    been applied e.g. by Shaw & Schumann (1992; Bound.-Layer Meteorol. 61, 47–64).
!--    When using the radiation_interactions, canopy heating (pcm_heating_rate) and plant canopy
!--    transpiration (pcm_transpiration_rate, pcm_latent_rate) are calculated in the RTM after the
!--    calculation of radiation.
       IF ( cthf /= 0.0_wp )  THEN

          ALLOCATE( cum_lai_hf(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--       Piecewise calculation of the cumulative leaf area index by vertical integration of the
!--       leaf area density
          cum_lai_hf(:,:,:) = 0.0_wp
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = pch_index_ji(j,i)-1, 0, -1
                   IF ( k == pch_index_ji(j,i)-1 )  THEN
                      cum_lai_hf(k,j,i) = cum_lai_hf(k+1,j,i) +                                    &
                                          ( 0.5_wp * lad_s(k+1,j,i) *                              &
                                            ( zw(k+1) - zu(k+1) ) )  +                             &
                                          ( 0.5_wp * ( 0.5_wp * ( lad_s(k+1,j,i) +                 &
                                                                  lad_s(k,j,i) ) +                 &
                                                       lad_s(k+1,j,i) ) *                          &
                                            ( zu(k+1) - zw(k) ) )
                   ELSE
                      cum_lai_hf(k,j,i) = cum_lai_hf(k+1,j,i) +                                    &
                                          ( 0.5_wp * ( 0.5_wp * ( lad_s(k+2,j,i) +                 &
                                                                  lad_s(k+1,j,i) ) +               &
                                                       lad_s(k+1,j,i) ) *                          &
                                            ( zw(k+1) - zu(k+1) ) )  +                             &
                                          ( 0.5_wp * ( 0.5_wp * ( lad_s(k+1,j,i) +                 &
                                                                  lad_s(k,j,i) ) +                 &
                                                       lad_s(k+1,j,i) ) *                          &
                                            ( zu(k+1) - zw(k) ) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

!
!--       In areas with canopy the surface value of the canopy heat flux distribution overrides the
!--       surface heat flux (shf),
!--       Start with default surface type
          DO  m = 1, surf_def_h(0)%ns
             i = surf_def_h(0)%i(m)
             j = surf_def_h(0)%j(m)
             IF ( cum_lai_hf(0,j,i) /= 0.0_wp )                                                    &
                surf_def_h(0)%shf(m) = cthf * EXP( -ext_coef * cum_lai_hf(0,j,i) )
          ENDDO
!
!--       Natural surfaces
          DO  m = 1, surf_lsm_h(0)%ns
             i = surf_lsm_h(0)%i(m)
             j = surf_lsm_h(0)%j(m)
             IF ( cum_lai_hf(0,j,i) /= 0.0_wp )                                                    &
                surf_lsm_h(0)%shf(m) = cthf * EXP( -ext_coef * cum_lai_hf(0,j,i) )
          ENDDO
!
!--       Urban surfaces
          DO  m = 1, surf_usm_h(0)%ns
             i = surf_usm_h(0)%i(m)
             j = surf_usm_h(0)%j(m)
             IF ( cum_lai_hf(0,j,i) /= 0.0_wp )                                                    &
                surf_usm_h(0)%shf(m) = cthf * EXP( -ext_coef * cum_lai_hf(0,j,i) )
          ENDDO
!
!
!--       Calculation of the heating rate (K/s) within the different layers of the plant canopy.
!--       Calculation is only necessary in areas covered with canopy.
!--       Within the different canopy layers the plant-canopy heating rate (pcm_heating_rate) is
!--       calculated as the vertical divergence of the canopy heat fluxes at the top and bottom of
!--       the respective layer.
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = 1, pch_index_ji(j,i)
                   IF ( cum_lai_hf(0,j,i) /= 0.0_wp )  THEN
                      pcm_heating_rate(k,j,i) = cthf *                                             &
                                                ( EXP( -ext_coef * cum_lai_hf(k,j,i) ) -           &
                                                  EXP( -ext_coef * cum_lai_hf(k-1,j,i) ) ) / dzw(k)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ENDIF

       IF ( debug_output )  CALL debug_message( 'pcm_init', 'end' )

    END SUBROUTINE pcm_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &plant_canopy_parameters for plant canopy model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pcm_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /plant_canopy_parameters/  alpha_lad,                                                 &
                                        beta_lad,                                                  &
                                        canopy_drag_coeff,                                         &
                                        canopy_mode,                                               &
                                        cthf,                                                      &
                                        lad_surface,                                               &
                                        lad_type_coef,                                             &
                                        lad_vertical_gradient,                                     &
                                        lad_vertical_gradient_level,                               &
                                        lai_beta,                                                  &
                                        leaf_scalar_exch_coeff,                                    &
                                        leaf_surface_conc,                                         &
                                        pch_index,                                                 &
                                        plant_canopy_transpiration,                                &
                                        switch_off_module

!
!-- Move to the beginning of the namelist file and try to find and read the user-defined namelist
!-- plant_canopy_parameters.
    REWIND( 11 )
    READ( 11, plant_canopy_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    plant_canopy_parameters namelist was found and read correctly. Set flag that indicates that
!--    the plant-canopy model is switched on.
       IF ( .NOT. switch_off_module )  plant_canopy = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    plant_canopy_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'plant_canopy_parameters', line )

    ENDIF

 END SUBROUTINE pcm_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Loads 3D plant canopy data from file. File format is as follows:
!>
!> num_levels
!> dtype,x,y,pctype,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> dtype,x,y,pctype,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> dtype,x,y,pctype,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> ...
!>
!> i.e. first line determines number of levels and further lines represent plant canopy data, one
!> line per column and variable. In each data line, dtype represents variable to be set:
!>
!> dtype=1: leaf area density (lad_s)
!> dtype=2....n: some additional plant canopy input data quantity
!>
!> Zeros are added automatically above num_levels until top of domain. Any non-specified (x,y)
!> columns have zero values as default.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_read_plant_canopy_3d

       USE exchange_horiz_mod,                                                                     &
           ONLY:  exchange_horiz

       INTEGER(iwp) ::  dtype     !< type of input data (1=lad)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  kk        !<
       INTEGER(iwp) ::  nzp       !< number of vertical layers of plant canopy
       INTEGER(iwp) ::  nzpltop   !<
       INTEGER(iwp) ::  nzpl      !<
       INTEGER(iwp) ::  pctype    !< type of plant canopy (deciduous,non-deciduous,...)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  col   !< vertical column of input data

!
!--    Initialize lad_s array
       lad_s = 0.0_wp

!
!--    Open and read plant canopy input data
       OPEN( 152, FILE='PLANT_CANOPY_DATA_3D' // TRIM( coupling_char ), ACCESS='SEQUENTIAL',       &
                  ACTION='READ', STATUS='OLD', FORM='FORMATTED', ERR=515 )
       READ( 152, *, ERR=516, END=517 )  nzp   !< read first line = number of vertical layers
       nzpltop = MIN(nzt+1, nzb+nzp-1)
       nzpl = nzpltop - nzb + 1    !< no. of layers to assign
       ALLOCATE( col(0:nzp-1) )

       DO
          READ( 152, *, ERR=516, END=517 ) dtype, i, j, pctype, col(:)
          IF ( i < nxlg  .OR.  i > nxrg  .OR.  j < nysg  .OR.  j > nyng )  CYCLE

          SELECT CASE ( dtype )
             CASE( 1 )   !< leaf area density
!
!--             This is just the pure canopy layer assumed to be grounded to a flat domain surface.
!--             At locations where plant canopy sits on top of any kind of topography, the vertical
!--             plant column must be "lifted", which is done in SUBROUTINE pcm_tendency.
                IF ( pctype < 0  .OR.  pctype > 10 )  THEN   !< incorrect plant canopy type
                   WRITE( message_string, * ) 'Incorrect type of plant canopy. '   //              &
                                              'Allowed values 0 <= pctype <= 10, ' //              &
                                              'but pctype is ', pctype
                   CALL message( 'pcm_read_plant_canopy_3d', 'PA0349', 1, 2, 0, 6, 0 )
                ENDIF
                kk = topo_top_ind(j,i,0)
                lad_s(nzb:nzpltop-kk, j, i) = col(kk:nzpl-1)*lad_type_coef(pctype)
             CASE DEFAULT
                WRITE( message_string, '(a,i2,a)' )                                                &
                     'Unknown record type in file PLANT_CANOPY_DATA_3D: "', dtype, '"'
                CALL message( 'pcm_read_plant_canopy_3d', 'PA0530', 1, 2, 0, 6, 0 )
          END SELECT
       ENDDO

515    message_string = 'error opening file PLANT_CANOPY_DATA_3D'
       CALL message( 'pcm_read_plant_canopy_3d', 'PA0531', 1, 2, 0, 6, 0 )

516    message_string = 'error reading file PLANT_CANOPY_DATA_3D'
       CALL message( 'pcm_read_plant_canopy_3d', 'PA0532', 1, 2, 0, 6, 0 )

517    CLOSE( 152 )
       DEALLOCATE( col )

       CALL exchange_horiz( lad_s, nbgp )

    END SUBROUTINE pcm_read_plant_canopy_3d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_rrd_global_ftn( found )

       LOGICAL, INTENT(OUT) ::  found


       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'pch_index' )
             READ( 13 )  pch_index

          CASE DEFAULT

             found = .FALSE.

       END SELECT

    END SUBROUTINE pcm_rrd_global_ftn

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_rrd_global_mpi

       CALL rrd_mpi_io( 'pch_index', pch_index )

    END SUBROUTINE pcm_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync, &
                                  nyn_on_file, nysf, nysc, nys_on_file, found )

       INTEGER(iwp) ::  k               !<
       INTEGER(iwp) ::  nxl_on_file     !<
       INTEGER(iwp) ::  nxlc            !<
       INTEGER(iwp) ::  nxlf            !<
       INTEGER(iwp) ::  nxr_on_file     !<
       INTEGER(iwp) ::  nxrc            !<
       INTEGER(iwp) ::  nxrf            !<
       INTEGER(iwp) ::  nyn_on_file     !<
       INTEGER(iwp) ::  nync            !<
       INTEGER(iwp) ::  nynf            !<
       INTEGER(iwp) ::  nys_on_file     !<
       INTEGER(iwp) ::  nysc            !<
       INTEGER(iwp) ::  nysf            !<

       LOGICAL, INTENT(OUT) ::  found

       REAL(wp), DIMENSION( 0:pch_index,                                                           &
                            nys_on_file-nbgp:nyn_on_file+nbgp,                                     &
                            nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d2  !< temporary 3D array for entire vertical
                                                                           !< extension of canopy layer
       found = .TRUE.


       SELECT CASE ( restart_string(1:length) )

          CASE ( 'pcm_heatrate_av' )
             IF ( .NOT. ALLOCATED( pcm_heatrate_av ) )  THEN
                ALLOCATE( pcm_heatrate_av(0:pch_index,nysg:nyng,nxlg:nxrg) )
                pcm_heatrate_av = 0.0_wp
             ENDIF
             IF ( k == 1 )  READ( 13 )  tmp_3d2
             pcm_heatrate_av(0:pch_index,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                &
                        tmp_3d2(0:pch_index,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'pcm_latentrate_av' )
             IF ( .NOT. ALLOCATED( pcm_latentrate_av ) )  THEN
                ALLOCATE( pcm_latentrate_av(0:pch_index,nysg:nyng,nxlg:nxrg) )
                pcm_latentrate_av = 0.0_wp
             ENDIF
             IF ( k == 1 )  READ( 13 )  tmp_3d2
             pcm_latentrate_av(0:pch_index,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                        tmp_3d2(0:pch_index,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'pcm_transpirationrate_av' )
             IF ( .NOT. ALLOCATED( pcm_transpirationrate_av ) )  THEN
                ALLOCATE( pcm_transpirationrate_av(0:pch_index,nysg:nyng,nxlg:nxrg) )
                pcm_transpirationrate_av = 0.0_wp
             ENDIF
             IF ( k == 1 )  READ( 13 )  tmp_3d2
             pcm_transpirationrate_av(0:pch_index,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =       &
                        tmp_3d2(0:pch_index,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE DEFAULT

             found = .FALSE.

       END SELECT

    END SUBROUTINE pcm_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_rrd_local_mpi

       IMPLICIT NONE

       LOGICAL ::  array_found  !<

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp_3d  !< temporary array to store pcm data with
                                                           !< non-standard vertical index bounds

!
!--    Plant canopy arrays have non standard reduced vertical index bounds. They are stored with
!--    full vertical bounds (bzb:nzt+1) in the restart file and must be re-stored after reading.
       CALL rd_mpi_io_check_array( 'pcm_heatrate_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( pcm_heatrate_av ) )  THEN
             ALLOCATE( pcm_heatrate_av(nzb:pch_index,nysg:nyng,nxlg:nxrg) )
          ENDIF
          ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'pcm_heatrate_av', tmp_3d )
          pcm_heatrate_av = tmp_3d(nzb:pch_index,:,:)
          DEALLOCATE( tmp_3d )
       ENDIF

       CALL rd_mpi_io_check_array( 'pcm_latentrate_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( pcm_latentrate_av ) )  THEN
             ALLOCATE( pcm_latentrate_av(nzb:pch_index,nysg:nyng,nxlg:nxrg) )
          ENDIF
          ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'pcm_latentrate_av', tmp_3d )
          pcm_latentrate_av = tmp_3d(nzb:pch_index,:,:)
          DEALLOCATE( tmp_3d )
       ENDIF

       CALL rd_mpi_io_check_array( 'pcm_transpirationrate_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( pcm_transpirationrate_av ) )  THEN
             ALLOCATE( pcm_transpirationrate_av(nzb:pch_index,nysg:nyng,nxlg:nxrg) )
          ENDIF
          ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'pcm_transpirationrate_av', tmp_3d )
          pcm_transpirationrate_av = tmp_3d(nzb:pch_index,:,:)
          DEALLOCATE( tmp_3d )
       ENDIF

    END SUBROUTINE pcm_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the tendency terms, accounting for the effect of the plant canopy on momentum and
!> scalar quantities.
!>
!> The canopy is located where the leaf area density lad_s(k,j,i) > 0.0 (defined on scalar grid), as
!> initialized in subroutine pcm_init.
!> The lad on the w-grid is vertically interpolated from the surrounding lad_s. The upper boundary
!> of the canopy is defined on the w-grid at k = pch_index. Here, the lad is zero.
!>
!> The canopy drag must be limited (previously accounted for by calculation of a limiting canopy
!> timestep for the determination of the maximum LES timestep in subroutine timestep), since it is
!> physically impossible that the canopy drag alone can locally change the sign of a velocity
!> component. This limitation is realized by calculating preliminary tendencies and velocities. It
!> is subsequently checked if the preliminary new velocity has a different sign than the current
!> velocity. If so, the tendency is limited in a way that the velocity can at maximum be reduced to
!> zero by the canopy drag.
!>
!>
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_tendency( component )

       INTEGER(iwp) ::  component !< prognostic variable (u,v,w,pt,q,e)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  k         !< running index
       INTEGER(iwp) ::  kk        !< running index for flat lad arrays

       LOGICAL ::  building_edge_e !< control flag indicating an eastward-facing building edge
       LOGICAL ::  building_edge_n !< control flag indicating a north-facing building edge
       LOGICAL ::  building_edge_s !< control flag indicating a south-facing building edge
       LOGICAL ::  building_edge_w !< control flag indicating a westward-facing building edge

       REAL(wp) ::  bad_local !< local bad value
       REAL(wp) ::  ddt_3d    !< inverse of the LES timestep (dt_3d)
       REAL(wp) ::  lad_local !< local lad value
       REAL(wp) ::  pre_tend  !< preliminary tendency
       REAL(wp) ::  pre_u     !< preliminary u-value
       REAL(wp) ::  pre_v     !< preliminary v-value
       REAL(wp) ::  pre_w     !< preliminary w-value


       ddt_3d = 1.0_wp / dt_3d

!
!--    Compute drag for the three velocity components and the SGS-TKE:
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  i = nxlu, nxr
                DO  j = nys, nyn
!
!--                Set control flags indicating east- and westward-orientated building edges. Note,
!--                building_egde_w is set from the perspective of the potential rooftop grid point,
!--                while building_edge_e is set from the perspective of the non-building grid point.
                   building_edge_w =              ANY( BTEST( topo_flags(:,j,i),   6 ) )           &
                                     .AND.  .NOT. ANY( BTEST( topo_flags(:,j,i-1), 6 ) )          
                   building_edge_e =              ANY( BTEST( topo_flags(:,j,i-1), 6 ) )           &
                                     .AND.  .NOT. ANY( BTEST( topo_flags(:,j,i),   6 ) )
!
!--                Determine topography-top index on u-grid
                   DO  k = topo_top_ind(j,i,1)+1, topo_top_ind(j,i,1) + pch_index_ji(j,i)

                      kk = k - topo_top_ind(j,i,1)   !- lad arrays are defined flat
!
!--                   In order to create sharp boundaries of the plant canopy, the lad on the u-grid
!--                   at index (k,j,i) is equal to lad_s(k,j,i), rather than being interpolated from
!--                   the surrounding lad_s, because this would yield smaller lad at the canopy
!--                   boundaries than inside of the canopy.
!--                   For the same reason, the lad at the rightmost(i+1)canopy boundary on the
!--                   u-grid equals lad_s(k,j,i), which is considered in the next if-statement.
!--                   Note, at left-sided building edges this is not applied, here the LAD equals
!--                   the LAD at grid point (k,j,i), in order to avoid that LAD is mistakenly mapped
!--                   on top of a roof where (usually) no LAD is defined. The same is also valid for
!--                   bad_s.
                      lad_local = lad_s(kk,j,i)
                      IF ( lad_local == 0.0_wp  .AND.  lad_s(kk,j,i-1) > 0.0_wp                    &
                           .AND.  .NOT. building_edge_w )  lad_local = lad_s(kk,j,i-1)

                      bad_local = bad_s(kk,j,i)
                      IF ( bad_local == 0.0_wp  .AND.  bad_s(kk,j,i-1) > 0.0_wp                    &
                           .AND.  .NOT. building_edge_w )  bad_local = bad_s(kk,j,i-1)
!
!--                   In order to avoid that LAD is mistakenly considered at right-sided building
!--                   edges (here the topography-top index for the u-component at index j,i is still
!--                   on the building while the topography top for the scalar isn't), LAD is taken
!--                   from grid point (j,i-1). The same is also valid for bad_s.
                      IF ( lad_local > 0.0_wp  .AND.  lad_s(kk,j,i-1) == 0.0_wp                    &
                           .AND.  building_edge_e )  lad_local = lad_s(kk,j,i-1)
                      IF ( bad_local > 0.0_wp  .AND.  bad_s(kk,j,i-1) == 0.0_wp                    &
                           .AND.  building_edge_e )  bad_local = bad_s(kk,j,i-1)

                      pre_tend = 0.0_wp
                      pre_u = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - canopy_drag_coeff *                                             &
                                   ( lad_local + bad_local ) *                                     &
                                   SQRT( u(k,j,i)**2 +                                             &
                                         ( 0.25_wp * ( v(k,j,i-1) +                                &
                                                       v(k,j,i)   +                                &
                                                       v(k,j+1,i) +                                &
                                                       v(k,j+1,i-1) )                              &
                                         )**2 +                                                    &
                                         ( 0.25_wp * ( w(k-1,j,i-1) +                              &
                                                       w(k-1,j,i)   +                              &
                                                       w(k,j,i-1)   +                              &
                                                       w(k,j,i) )                                  &
                                         )**2                                                      &
                                       ) *                                                         &
                                   u(k,j,i)

!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_u = u(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN( pre_u,u(k,j,i) ) /= pre_u )  THEN
                         pre_tend = - u(k,j,i) * ddt_3d
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  i = nxl, nxr
                DO  j = nysv, nyn
!
!--                Set control flags indicating north- and southward-orientated building edges.
!--                Note, building_egde_s is set from the perspective of the potential rooftop grid
!--                point, while building_edge_n is set from the perspective of the non-building grid
!--                point.
                   building_edge_s =              ANY( BTEST( topo_flags(:,j,i),   6 ) )           &
                                     .AND.  .NOT. ANY( BTEST( topo_flags(:,j-1,i), 6 ) )           
                   building_edge_n =              ANY( BTEST( topo_flags(:,j-1,i), 6 ) )           &
                                     .AND.  .NOT. ANY( BTEST( topo_flags(:,j,i),   6 ) )
!
!--                Determine topography-top index on v-grid
                   DO  k = topo_top_ind(j,i,2)+1, topo_top_ind(j,i,2) + pch_index_ji(j,i)

                      kk = k - topo_top_ind(j,i,2)   !- lad arrays are defined flat
!
!--                   In order to create sharp boundaries of the plant canopy, the lad on the v-grid
!--                   at index (k,j,i) is equal to lad_s(k,j,i), rather than being interpolated from
!--                   the surrounding lad_s, because this would yield smaller lad at the canopy
!--                   boundaries  than inside of the canopy.
!--                   For the same reason, the lad at the northmost (j+1) canopy boundary on the
!--                   v-grid equals lad_s(k,j,i), which is considered in the next if-statement.
!--                   Note, at left-sided building edges this is not applied, here the LAD equals
!--                   the LAD at grid point (k,j,i), in order to avoid that LAD is mistakenly mapped
!--                   on top of a roof where (usually) no LAD is defined.
!--                   The same is also valid for bad_s.
                      lad_local = lad_s(kk,j,i)
                      IF ( lad_local == 0.0_wp  .AND.  lad_s(kk,j-1,i) > 0.0_wp                    &
                           .AND.  .NOT. building_edge_s )  lad_local = lad_s(kk,j-1,i)

                      bad_local = bad_s(kk,j,i)
                      IF ( bad_local == 0.0_wp  .AND.  bad_s(kk,j-1,i) > 0.0_wp                    &
                           .AND.  .NOT. building_edge_s )  bad_local = bad_s(kk,j-1,i)
!
!--                   In order to avoid that LAD is mistakenly considered at right-sided building
!--                   edges (here the topography-top index for the u-component at index j,i is still
!--                   on the building while the topography top for the scalar isn't), LAD is taken
!--                   from grid point (j,i-1). The same is also valid for bad_s.
                      IF ( lad_local > 0.0_wp  .AND.  lad_s(kk,j-1,i) == 0.0_wp                    &
                           .AND.  building_edge_n )  lad_local = lad_s(kk,j-1,i)

                      IF ( bad_local > 0.0_wp  .AND.  bad_s(kk,j-1,i) == 0.0_wp                    &
                           .AND.  building_edge_n )  bad_local = bad_s(kk,j-1,i)

                      pre_tend = 0.0_wp
                      pre_v = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - canopy_drag_coeff *                                             &
                                   ( lad_local + bad_local ) *                                     &
                                   SQRT( ( 0.25_wp * ( u(k,j-1,i)   +                              &
                                                       u(k,j-1,i+1) +                              &
                                                       u(k,j,i)     +                              &
                                                       u(k,j,i+1) )                                &
                                         )**2 +                                                    &
                                         v(k,j,i)**2 +                                             &
                                         ( 0.25_wp * ( w(k-1,j-1,i) +                              &
                                                       w(k-1,j,i)   +                              &
                                                       w(k,j-1,i)   +                              &
                                                       w(k,j,i) )                                  &
                                         )**2                                                      &
                                       ) *                                                         &
                                   v(k,j,i)

!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_v = v(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity, and in case the
!--                   signs are different, limit the tendency.
                      IF ( SIGN( pre_v,v(k,j,i) ) /= pre_v )  THEN
                         pre_tend = - v(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on w-grid
                   DO  k = topo_top_ind(j,i,3)+1, topo_top_ind(j,i,3) + pch_index_ji(j,i) - 1

                      kk = k - topo_top_ind(j,i,3)   !- lad arrays are defined flat

                      pre_tend = 0.0_wp
                      pre_w = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - canopy_drag_coeff *                                             &
                                   ( 0.5_wp * ( lad_s(kk+1,j,i) + lad_s(kk,j,i) ) +                &
                                     0.5_wp * ( bad_s(kk+1,j,i) + bad_s(kk,j,i) ) ) *              &
                                   SQRT( ( 0.25_wp * ( u(k,j,i)   +                                &
                                                       u(k,j,i+1) +                                &
                                                       u(k+1,j,i) +                                &
                                                       u(k+1,j,i+1) )                              &
                                         )**2 +                                                    &
                                         ( 0.25_wp * ( v(k,j,i)   +                                &
                                                       v(k,j+1,i) +                                &
                                                       v(k+1,j,i) +                                &
                                                       v(k+1,j+1,i) )                              &
                                         )**2 +                                                    &
                                         w(k,j,i)**2                                               &
                                       ) *                                                         &
                                   w(k,j,i)
!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_w = w(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity, and in case the
!--                   signs are different, limit the tendency
                      IF ( SIGN( pre_w,w(k,j,i) ) /= pre_w )  THEN
                         pre_tend = - w(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       potential temperature
          CASE ( 4 )
             IF ( humidity )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
!--                   Determine topography-top index on scalar-grid
                      DO  k = topo_top_ind(j,i,0)+1, topo_top_ind(j,i,0) + pch_index_ji(j,i)
                         kk = k - topo_top_ind(j,i,0)   !- lad arrays are defined flat
                         tend(k,j,i) = tend(k,j,i) + pcm_heating_rate(kk,j,i) - pcm_latent_rate(kk,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
!--                   Determine topography-top index on scalar-grid
                      DO  k = topo_top_ind(j,i,0)+1, topo_top_ind(j,i,0) + pch_index_ji(j,i)
                         kk = k - topo_top_ind(j,i,0)   !- lad arrays are defined flat
                         tend(k,j,i) = tend(k,j,i) + pcm_heating_rate(kk,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

!
!--       humidity
          CASE ( 5 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   DO  k = topo_top_ind(j,i,0)+1, topo_top_ind(j,i,0) + pch_index_ji(j,i)

                      kk = k - topo_top_ind(j,i,0)   !- lad arrays are defined flat

                      IF ( .NOT. plant_canopy_transpiration )  THEN
                         ! pcm_transpiration_rate is calculated in radiation model
                         ! in case of plant_canopy_transpiration = .T.
                         ! to include also the dependecy to the radiation
                         ! in the plant canopy box
                         pcm_transpiration_rate(kk,j,i) = - leaf_scalar_exch_coeff                 &
                                                            * lad_s(kk,j,i) *                      &
                                                            SQRT( ( 0.5_wp * ( u(k,j,i) +          &
                                                                               u(k,j,i+1) )        &
                                                                  )**2 +                           &
                                                                  ( 0.5_wp * ( v(k,j,i) +          &
                                                                               v(k,j+1,i) )        &
                                                                  )**2 +                           &
                                                                  ( 0.5_wp * ( w(k-1,j,i) +        &
                                                                               w(k,j,i) )          &
                                                                  )**2                             &
                                                                ) *                                &
                                                            ( q(k,j,i) - leaf_surface_conc )
                      ENDIF

                      tend(k,j,i) = tend(k,j,i) + pcm_transpiration_rate(kk,j,i)
                   ENDDO
                ENDDO
             ENDDO

!
!--       sgs-tke
          CASE ( 6 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   DO  k = topo_top_ind(j,i,0)+1, topo_top_ind(j,i,0) + pch_index_ji(j,i)

                      kk = k - topo_top_ind(j,i,0)   !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -  2.0_wp * canopy_drag_coeff *                    &
                                                   ( lad_s(kk,j,i) + bad_s(kk,j,i) ) *             &
                                                   SQRT( ( 0.5_wp * ( u(k,j,i) +                   &
                                                                      u(k,j,i+1) )                 &
                                                         )**2 +                                    &
                                                         ( 0.5_wp * ( v(k,j,i) +                   &
                                                                      v(k,j+1,i) )                 &
                                                         )**2 +                                    &
                                                         ( 0.5_wp * ( w(k,j,i) +                   &
                                                                      w(k+1,j,i) )                 &
                                                         )**2                                      &
                                                       ) *                                         &
                                                   e(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
!
!--       scalar concentration
          CASE ( 7 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Determine topography-top index on scalar-grid
                   DO  k = topo_top_ind(j,i,0)+1, topo_top_ind(j,i,0) + pch_index_ji(j,i)

                      kk = k - topo_top_ind(j,i,0)   !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -  leaf_scalar_exch_coeff *                        &
                                                   lad_s(kk,j,i) *                                 &
                                                   SQRT( ( 0.5_wp * ( u(k,j,i) +                   &
                                                                      u(k,j,i+1) )                 &
                                                         )**2 +                                    &
                                                         ( 0.5_wp * ( v(k,j,i) +                   &
                                                                      v(k,j+1,i) )                 &
                                                         )**2 +                                    &
                                                         ( 0.5_wp * ( w(k-1,j,i) +                 &
                                                                      w(k,j,i) )                   &
                                                         )**2                                      &
                                                       ) *                                         &
                                                   ( s(k,j,i) - leaf_surface_conc )
                   ENDDO
                ENDDO
             ENDDO



          CASE DEFAULT

             WRITE( message_string, * ) 'wrong component: ', component
             CALL message( 'pcm_tendency', 'PA0279', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE pcm_tendency


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the tendency terms, accounting for the effect of the plant canopy on momentum and
!> scalar quantities.
!>
!> The canopy is located where the leaf area density lad_s(k,j,i) > 0.0 (defined on scalar grid), as
!> initialized in subroutine pcm_init.
!> The lad on the w-grid is vertically interpolated from the surrounding lad_s. The upper boundary
!> of the canopy is defined on the w-grid at k = pch_index. Here, the lad is zero.
!>
!> The canopy drag must be limited (previously accounted for by calculation of a limiting canopy
!> timestep for the determination of the maximum LES timestep in subroutine timestep), since it is
!> physically impossible that the canopy drag alone can locally change the sign of a velocity
!> component. This limitation is realized by calculating preliminary tendencies and velocities. It
!> is subsequently checked if the preliminary new velocity has a different sign than the current
!> velocity. If so, the tendency is limited in a way that the velocity can at maximum be reduced to
!> zero by the canopy drag.
!>
!>
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_tendency_ij( i, j, component )

       INTEGER(iwp) ::  component !< prognostic variable (u,v,w,pt,q,e)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  k         !< running index
       INTEGER(iwp) ::  kk        !< running index for flat lad arrays

       LOGICAL ::  building_edge_e !< control flag indicating an eastward-facing building edge
       LOGICAL ::  building_edge_n !< control flag indicating a north-facing building edge
       LOGICAL ::  building_edge_s !< control flag indicating a south-facing building edge
       LOGICAL ::  building_edge_w !< control flag indicating a westward-facing building edge

       REAL(wp) ::  bad_local !< local lad value
       REAL(wp) ::  ddt_3d    !< inverse of the LES timestep (dt_3d)
       REAL(wp) ::  lad_local !< local lad value
       REAL(wp) ::  pre_tend  !< preliminary tendency
       REAL(wp) ::  pre_u     !< preliminary u-value
       REAL(wp) ::  pre_v     !< preliminary v-value
       REAL(wp) ::  pre_w     !< preliminary w-value


       ddt_3d = 1.0_wp / dt_3d
!
!--    Compute drag for the three velocity components and the SGS-TKE
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
!
!--          Set control flags indicating east- and westward-orientated building edges. Note,
!--          building_egde_w is set from the perspective of the potential rooftop grid point, while
!--          building_edge_e is set from the perspective of the non-building grid point.
             building_edge_w =       ANY( BTEST( topo_flags(:,j,i),   6 ) )  .AND.                 &
                               .NOT. ANY( BTEST( topo_flags(:,j,i-1), 6 ) )                       
             building_edge_e =       ANY( BTEST( topo_flags(:,j,i-1), 6 ) )  .AND.                 &
                               .NOT. ANY( BTEST( topo_flags(:,j,i),   6 ) )
!
!--          Determine topography-top index on u-grid
             DO  k = topo_top_ind(j,i,1) + 1, topo_top_ind(j,i,1) + pch_index_ji(j,i)

                kk = k - topo_top_ind(j,i,1)  !- lad arrays are defined flat

!
!--             In order to create sharp boundaries of the plant canopy, the lad on the u-grid at
!--             index (k,j,i) is equal to lad_s(k,j,i), rather than being interpolated from the
!--             surrounding lad_s, because this would yield smaller lad at the canopy boundaries
!--             than inside of the canopy.
!--             For the same reason, the lad at the rightmost(i+1)canopy boundary on the u-grid
!--             equals lad_s(k,j,i), which is considered in the next if-statement. Note, at
!--             left-sided building edges this is not applied, here the LAD is equals the LAD at
!--             grid point (k,j,i), in order to avoid that LAD is mistakenly mapped on top of a roof
!--             where (usually) is no LAD is defined.
!--             The same is also valid for bad_s.
                lad_local = lad_s(kk,j,i)
                IF ( lad_local == 0.0_wp  .AND.  lad_s(kk,j,i-1) > 0.0_wp  .AND.                   &
                     .NOT. building_edge_w )  lad_local = lad_s(kk,j,i-1)

                bad_local = bad_s(kk,j,i)
                IF ( bad_local == 0.0_wp  .AND.  bad_s(kk,j,i-1) > 0.0_wp  .AND.                   &
                     .NOT. building_edge_w )  bad_local = bad_s(kk,j,i-1)
!
!--             In order to avoid that LAD is mistakenly considered at right-sided building edges
!--             (here the topography-top index for the u-component at index j,i is still on the
!--             building while the topography top for the scalar isn't), LAD is taken from grid
!--             point (j,i-1). The same is also valid for bad_s.
                IF ( lad_local > 0.0_wp  .AND.  lad_s(kk,j,i-1) == 0.0_wp  .AND.                   &
                     building_edge_e )  lad_local = lad_s(kk,j,i-1)
                IF ( bad_local > 0.0_wp  .AND.  bad_s(kk,j,i-1) == 0.0_wp  .AND.                   &
                     building_edge_e )  bad_local = bad_s(kk,j,i-1)

                pre_tend = 0.0_wp
                pre_u = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - canopy_drag_coeff *                                                   &
                             ( lad_local + bad_local ) *                                           &
                             SQRT( u(k,j,i)**2 +                                                   &
                                   ( 0.25_wp * ( v(k,j,i-1)  +                                     &
                                                 v(k,j,i)    +                                     &
                                                 v(k,j+1,i)  +                                     &
                                                 v(k,j+1,i-1) )                                    &
                                   )**2 +                                                          &
                                   ( 0.25_wp * ( w(k-1,j,i-1) +                                    &
                                                 w(k-1,j,i)   +                                    &
                                                 w(k,j,i-1)   +                                    &
                                                 w(k,j,i) )                                        &
                                   )**2                                                            &
                                 ) *                                                               &
                             u(k,j,i)

!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_u = u(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity, and in case the signs are
!--             different, limit the tendency.
                IF ( SIGN( pre_u,u(k,j,i) ) /= pre_u )  THEN
                   pre_tend = - u(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO


!
!--       v-component
          CASE ( 2 )
!
!--          Set control flags indicating north- and southward-orientated building edges. Note,
!--          building_egde_s is set from the perspective of the potential rooftop grid point, while
!--          building_edge_n is set from the perspective of the non-building grid point.
             building_edge_s =       ANY( BTEST( topo_flags(:,j,i),   6 ) )  .AND.                 &
                               .NOT. ANY( BTEST( topo_flags(:,j-1,i), 6 ) )                       
             building_edge_n =       ANY( BTEST( topo_flags(:,j-1,i), 6 ) )  .AND.                 &
                               .NOT. ANY( BTEST( topo_flags(:,j,i),   6 ) )
!
!--          Determine topography-top index on v-grid
             DO  k = topo_top_ind(j,i,2) + 1, topo_top_ind(j,i,2) + pch_index_ji(j,i)

                kk = k - topo_top_ind(j,i,2)  !- lad arrays are defined flat
!
!--             In order to create sharp boundaries of the plant canopy, the lad on the v-grid at
!--             index (k,j,i) is equal to lad_s(k,j,i), rather than being interpolated from the
!--             surrounding lad_s, because this would yield smaller lad at the canopy boundaries
!--             than inside of the canopy.
!--             For the same reason, the lad at the northmost (j+1) canopy boundary on the v-grid
!--             equals lad_s(k,j,i), which is considered in the next if-statement. Note, at
!--             left-sided building edges this is not applied, here the LAD is equals the LAD at
!--             grid point (k,j,i), in order to avoid that LAD is mistakenly mapped on top of a roof
!--             where (usually) is no LAD is defined.
!--             The same is also valid for bad_s.
                lad_local = lad_s(kk,j,i)
                IF ( lad_local == 0.0_wp  .AND.  lad_s(kk,j-1,i) > 0.0_wp  .AND.                   &
                     .NOT. building_edge_s )  lad_local = lad_s(kk,j-1,i)

                bad_local = bad_s(kk,j,i)
                IF ( bad_local == 0.0_wp  .AND.  bad_s(kk,j-1,i) > 0.0_wp  .AND.                   &
                     .NOT. building_edge_s )  bad_local = bad_s(kk,j-1,i)
!
!--             In order to avoid that LAD is mistakenly considered at right-sided building edges
!--             (here the topography-top index for the u-component at index j,i is still on the
!--             building while the topography top for the scalar isn't), LAD is taken from grid
!--             point (j,i-1). The same is also valid for bad_s.
                IF ( lad_local > 0.0_wp  .AND.  lad_s(kk,j-1,i) == 0.0_wp  .AND.                   &
                     building_edge_n )  lad_local = lad_s(kk,j-1,i)
                IF ( bad_local > 0.0_wp  .AND.  bad_s(kk,j-1,i) == 0.0_wp  .AND.                   &
                     building_edge_n )  bad_local = bad_s(kk,j-1,i)

                pre_tend = 0.0_wp
                pre_v = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - canopy_drag_coeff *                                                   &
                             ( lad_local + bad_local ) *                                           &
                             SQRT( ( 0.25_wp * ( u(k,j-1,i)   +                                    &
                                                 u(k,j-1,i+1) +                                    &
                                                 u(k,j,i)     +                                    &
                                                 u(k,j,i+1) )                                      &
                                   )**2 +                                                          &
                                   v(k,j,i)**2 +                                                   &
                                   ( 0.25_wp * ( w(k-1,j-1,i) +                                    &
                                                 w(k-1,j,i)   +                                    &
                                                 w(k,j-1,i)   +                                    &
                                                 w(k,j,i) )                                        &
                                   )**2                                                            &
                                 ) *                                                               &
                             v(k,j,i)

!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_v = v(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN( pre_v,v(k,j,i) ) /= pre_v )  THEN
                   pre_tend = - v(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO


!
!--       w-component
          CASE ( 3 )
!
!--          Determine topography-top index on w-grid
             DO  k = topo_top_ind(j,i,3) + 1, topo_top_ind(j,i,3) + pch_index_ji(j,i) - 1

                kk = k - topo_top_ind(j,i,3)  !- lad arrays are defined flat

                pre_tend = 0.0_wp
                pre_w = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - canopy_drag_coeff *                                                   &
                             ( 0.5_wp * ( lad_s(kk+1,j,i) + lad_s(kk,j,i) ) +                      &
                               0.5_wp * ( bad_s(kk+1,j,i) + bad_s(kk,j,i) ) ) *                    &
                             SQRT( ( 0.25_wp * ( u(k,j,i)    +                                     &
                                                 u(k,j,i+1)  +                                     &
                                                 u(k+1,j,i)  +                                     &
                                                 u(k+1,j,i+1) )                                    &
                                   )**2 +                                                          &
                                   ( 0.25_wp * ( v(k,j,i)    +                                     &
                                                 v(k,j+1,i)  +                                     &
                                                 v(k+1,j,i)  +                                     &
                                                 v(k+1,j+1,i) )                                    &
                                   )**2 +                                                          &
                                   w(k,j,i)**2                                                     &
                                 ) *                                                               &
                             w(k,j,i)
!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_w = w(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity, and in case the signs are
!--             different, limit the tendency.
                IF ( SIGN( pre_w,w(k,j,i) ) /= pre_w )  THEN
                   pre_tend = - w(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO

!
!--       potential temperature
          CASE ( 4 )
!
!--          Determine topography-top index on scalar grid
             IF ( humidity )  THEN
                DO  k = topo_top_ind(j,i,0) + 1, topo_top_ind(j,i,0) + pch_index_ji(j,i)
                   kk = k - topo_top_ind(j,i,0)  !- lad arrays are defined flat
                   tend(k,j,i) = tend(k,j,i) + pcm_heating_rate(kk,j,i) - pcm_latent_rate(kk,j,i)
                ENDDO
             ELSE
                DO  k = topo_top_ind(j,i,0) + 1, topo_top_ind(j,i,0) + pch_index_ji(j,i)
                   kk = k - topo_top_ind(j,i,0)  !- lad arrays are defined flat
                   tend(k,j,i) = tend(k,j,i) + pcm_heating_rate(kk,j,i)
                ENDDO
             ENDIF

!
!--       humidity
          CASE ( 5 )
!
!--          Determine topography-top index on scalar grid
             DO  k = topo_top_ind(j,i,0) + 1, topo_top_ind(j,i,0) + pch_index_ji(j,i)
                kk = k - topo_top_ind(j,i,0)  !- lad arrays are defined flat
                IF ( .NOT. plant_canopy_transpiration )  THEN
                   ! pcm_transpiration_rate is calculated in radiation model
                   ! in case of plant_canopy_transpiration = .T.
                   ! to include also the dependecy to the radiation
                   ! in the plant canopy box
                   pcm_transpiration_rate(kk,j,i) = - leaf_scalar_exch_coeff                       &
                                                    * lad_s(kk,j,i) *                              &
                                                    SQRT( ( 0.5_wp * ( u(k,j,i) +                  &
                                                                       u(k,j,i+1) )                &
                                                          )**2  +                                  &
                                                          ( 0.5_wp * ( v(k,j,i) +                  &
                                                                       v(k,j+1,i) )                &
                                                          )**2 +                                   &
                                                          ( 0.5_wp * ( w(k-1,j,i) +                &
                                                                       w(k,j,i) )                  &
                                                          )**2                                     &
                                                        ) *                                        &
                                                    ( q(k,j,i) - leaf_surface_conc )
                ENDIF

                tend(k,j,i) = tend(k,j,i) + pcm_transpiration_rate(kk,j,i)

             ENDDO

!
!--       sgs-tke
          CASE ( 6 )
!
!--          Determine topography-top index on scalar grid
             DO  k = topo_top_ind(j,i,0) + 1, topo_top_ind(j,i,0) + pch_index_ji(j,i)

                kk = k - topo_top_ind(j,i,0)
                tend(k,j,i) = tend(k,j,i) -  2.0_wp * canopy_drag_coeff *                          &
                                             ( lad_s(kk,j,i) + bad_s(kk,j,i) ) *                   &
                                             SQRT( ( 0.5_wp * ( u(k,j,i) +                         &
                                                                u(k,j,i+1) )                       &
                                                   )**2 +                                          &
                                                   ( 0.5_wp * ( v(k,j,i) +                         &
                                                                v(k,j+1,i) )                       &
                                                   )**2 +                                          &
                                                   ( 0.5_wp * ( w(k,j,i) +                         &
                                                                w(k+1,j,i) )                       &
                                                   )**2                                            &
                                                 ) *                                               &
                                             e(k,j,i)
             ENDDO
!
!--       scalar concentration
          CASE ( 7 )
!
!--          Determine topography-top index on scalar grid
             DO  k = topo_top_ind(j,i,0) + 1, topo_top_ind(j,i,0) + pch_index_ji(j,i)

                kk = k - topo_top_ind(j,i,0)
                tend(k,j,i) = tend(k,j,i) -  leaf_scalar_exch_coeff *                              &
                                             lad_s(kk,j,i) *                                       &
                                             SQRT( ( 0.5_wp * ( u(k,j,i) +                         &
                                                                u(k,j,i+1) )                       &
                                                   )**2  +                                         &
                                                   ( 0.5_wp * ( v(k,j,i) +                         &
                                                                v(k,j+1,i) )                       &
                                                   )**2 +                                          &
                                                   ( 0.5_wp * ( w(k-1,j,i) +                       &
                                                                w(k,j,i) )                         &
                                                   )**2                                            &
                                                 ) *                                               &
                                             ( s(k,j,i) - leaf_surface_conc )
             ENDDO

       CASE DEFAULT

          WRITE( message_string, * ) 'wrong component: ', component
          CALL message( 'pcm_tendency', 'PA0279', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE pcm_tendency_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes global restart data
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_wrd_global

       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          CALL wrd_write_string( 'pch_index' )
          WRITE( 14 )  pch_index

       ELSE IF ( restart_data_format_output(1:3) == 'mpi' )  THEN

          CALL wrd_mpi_io( 'pch_index', pch_index )

       ENDIF

    END SUBROUTINE pcm_wrd_global

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes local (subdomain) restart data
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pcm_wrd_local

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp_3d  !< temporary array to store pcm data with
                                                           !< non-standard vertical index bounds


       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          IF ( ALLOCATED( pcm_heatrate_av ) )  THEN
             CALL wrd_write_string( 'pcm_heatrate_av' )
             WRITE( 14 )  pcm_heatrate_av
          ENDIF

          IF ( ALLOCATED( pcm_latentrate_av ) )  THEN
             CALL wrd_write_string( 'pcm_latentrate_av' )
             WRITE( 14 )  pcm_latentrate_av
          ENDIF

          IF ( ALLOCATED( pcm_transpirationrate_av ) )  THEN
             CALL wrd_write_string( 'pcm_transpirationrate_av' )
             WRITE( 14 )  pcm_transpirationrate_av
          ENDIF

       ELSE IF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!
!--       Plant canopy arrays have non standard reduced vertical index bounds. They are stored with
!--       full vertical bounds (bzb:nzt+1) in the restart file and must be re-stored before writing.
          IF ( ALLOCATED( pcm_heatrate_av ) )  THEN
             ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             tmp_3d(nzb:pch_index,:,:)     = pcm_heatrate_av
             tmp_3d(pch_index+1:nzt+1,:,:) = 0.0_wp
             CALL wrd_mpi_io( 'pcm_heatrate_av', tmp_3d )
             DEALLOCATE( tmp_3d )
          ENDIF
          IF ( ALLOCATED( pcm_latentrate_av ) )  THEN
             ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             tmp_3d(nzb:pch_index,:,:)     = pcm_latentrate_av
             tmp_3d(pch_index+1:nzt+1,:,:) = 0.0_wp
             CALL wrd_mpi_io( 'pcm_latentrate_av', tmp_3d )
             DEALLOCATE( tmp_3d )
          ENDIF
          IF ( ALLOCATED( pcm_transpirationrate_av ) )  THEN
             ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             tmp_3d(nzb:pch_index,:,:)     = pcm_transpirationrate_av
             tmp_3d(pch_index+1:nzt+1,:,:) = 0.0_wp
             CALL wrd_mpi_io( 'pcm_transpirationrate_av', tmp_3d )
             DEALLOCATE( tmp_3d )
          ENDIF

       ENDIF

    END SUBROUTINE pcm_wrd_local


 END MODULE plant_canopy_model_mod
