#--------------------------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along with PALM. If not, see
# <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2021  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------------------------#
#
# WARNING: don't write filenames with extension .f90 in this header!!!!
#
# Description:
# ------------
# Makefile for the parallelized LES model (PALM) to be used by the
# shellscript "mrun"
#------------------------------------------------------------------------------!

PROG  =  palm
PROG1 =  combine_plot_fields.x

SOURCES = \
	advec_s_bc.f90 \
	advec_s_pw.f90 \
	advec_s_up.f90 \
	advec_u_pw.f90 \
	advec_u_up.f90 \
	advec_v_pw.f90 \
	advec_v_up.f90 \
	advec_w_pw.f90 \
	advec_ws.f90 \
	advec_w_up.f90 \
	average_3d_data.f90 \
	basic_constants_and_equations_mod.f90 \
	biometeorology_mod.f90 \
	buoyancy.f90 \
	calc_mean_profile.f90 \
	check_for_restart.f90 \
	check_open.f90 \
	check_parameters.f90 \
	chem_emis_domestic_mod.f90 \
	chem_emis_generic_mod.f90 \
	chem_emis_traffic_mod.f90 \
	chem_emis_vsrc_mod.f90 \
	chem_emissions_mod.f90 \
	chem_gasphase_mod.f90 \
	chemistry_model_mod.f90 \
	chem_modules.f90 \
	chem_photolysis_mod.f90 \
	close_file.f90 \
	compute_vpt.f90 \
	coriolis.f90 \
	cpulog_mod.f90 \
	cuda_fft_interfaces.f90 \
	data_log.f90 \
	data_output_binary_module.f90 \
	data_output_module.f90 \
	data_output_netcdf4_module.f90 \
	data_output_2d.f90 \
	data_output_3d.f90 \
	data_output_flight.f90\
	data_output_mask.f90 \
	data_output_particle_mod.f90 \
	data_output_profiles.f90 \
	data_output_spectra.f90 \
	data_output_topo_and_surface_setup_mod.f90 \
	data_output_tseries.f90 \
	diagnostic_output_quantities_mod.f90 \
	diffusion_s.f90 \
	diffusion_u.f90 \
	diffusion_v.f90 \
	diffusion_w.f90 \
	disturb_field.f90 \
	disturb_heatflux.f90 \
	dynamics_mod.f90 \
	exchange_horiz_mod.f90 \
	fft_xy_mod.f90 \
	flow_statistics.f90 \
	global_min_max.f90 \
	gust_mod.f90 \
	header.f90 \
	indoor_model_mod.f90 \
	inflow_turbulence.f90 \
	init_3d_model.f90 \
	init_advec.f90 \
	init_coupling.f90 \
	init_grid.f90 \
	init_masks.f90 \
	init_pegrid.f90 \
	init_pt_anomaly.f90 \
	init_rankine.f90 \
	init_slope.f90 \
	init_vertical_profiles.f90 \
	lagrangian_particle_model_mod.f90 \
	land_surface_model_mod.f90 \
	large_scale_forcing_nudging_mod.f90 \
	local_stop.f90 \
	local_system.f90 \
	local_tremain.f90 \
	local_tremain_ini.f90 \
	message.f90 \
	bulk_cloud_model_mod.f90 \
	model_1d_mod.f90 \
	module_interface.f90 \
	mod_kinds.f90 \
	mod_particle_attributes.f90 \
	modules.f90 \
	multi_agent_system_mod.f90 \
	netcdf_data_input_mod.f90 \
	netcdf_interface_mod.f90 \
	nesting_offl_mod.f90 \
	ocean_mod.f90 \
	outflow_turbulence.f90 \
	palm.f90 \
	palm_date_time_mod.f90 \
	parin.f90 \
	plant_canopy_model_mod.f90 \
	pmc_child_mod.f90 \
	pmc_general_mod.f90 \
	pmc_handle_communicator_mod.f90 \
	pmc_interface_mod.f90 \
	pmc_mpi_wrapper_mod.f90 \
	pmc_parent_mod.f90 \
	pmc_particle_interface.f90 \
	poisfft_mod.f90 \
	poismg_mod.f90 \
	poismg_noopt_mod.f90 \
	posix_interface_mod.f90 \
	pres.f90 \
	print_1d.f90 \
	prognostic_equations.f90 \
	progress_bar_mod.f90 \
	radiation_model_mod.f90 \
	random_function_mod.f90 \
	random_gauss.f90 \
	random_generator_parallel_mod.f90 \
	read_restart_data_mod.f90 \
	restart_data_mpi_io_mod.f90 \
	run_control.f90 \
	salsa_mod.f90 \
	shared_memory_io_mod.f90 \
	singleton_mod.f90 \
	sor.f90 \
	spectra_mod.f90 \
	subsidence_mod.f90 \
	sum_up_3d_data.f90 \
	surface_coupler.f90 \
	surface_layer_fluxes_mod.f90 \
	surface_mod.f90 \
	surface_data_output_mod.f90 \
	synthetic_turbulence_generator_mod.f90 \
	temperton_fft_mod.f90 \
	time_integration.f90 \
	time_integration_spinup.f90 \
	timestep.f90 \
	timestep_scheme_steering.f90 \
	time_to_string.f90 \
	transpose_mod.f90 \
	tridia_solver_mod.f90 \
	turbulence_closure_mod.f90 \
	urban_surface_mod.f90 \
	user_data_output_mask.f90 \
	user_flight.f90\
	user_init_3d_model.f90 \
	user_init_flight_mod.f90 \
	user_init_grid.f90 \
	user_init_land_surface.f90 \
	user_init_plant_canopy.f90 \
	user_init_radiation.f90 \
	user_init_urban_surface.f90 \
	user_lpm_advec.f90 \
	user_lpm_init.f90 \
	user_module.f90 \
	user_spectra.f90 \
	vdi_internal_controls.f90 \
	virtual_flight_mod.f90 \
	virtual_measurement_mod.f90 \
	wind_turbine_model_mod.f90 \
	wrd_write_string.f90 \
	write_restart_data_mod.f90


OBJS  = $(SOURCES:.f90=.o)
OBJS1 = combine_plot_fields.o

CC = cc
CFLAGS = -O

F90 =
F90_SER =
COPT =
F90FLAGS =
LDFLAGS =


# The following line is needed for palm_simple_install, don't remove it!
#to_be_replaced_by_include

.SUFFIXES:
.SUFFIXES: .o .f90


all: $(PROG) $(PROG1)

$(PROG): $(OBJS)
	$(F90) -o $(PROG) $(OBJS) $(LDFLAGS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.smod *.i *.lst

.f90.o:
	$(F90) $(F90FLAGS) $(COPT) -c $<


advec_s_bc.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
advec_s_pw.o: \
	mod_kinds.o \
	modules.o
advec_s_up.o: \
	mod_kinds.o \
	modules.o
advec_u_pw.o: \
	mod_kinds.o \
	modules.o
advec_u_up.o: \
	mod_kinds.o \
	modules.o
advec_v_pw.o: \
	mod_kinds.o \
	modules.o
advec_v_up.o: \
	mod_kinds.o \
	modules.o
advec_ws.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o
advec_w_pw.o: \
	mod_kinds.o \
	modules.o
advec_w_up.o: \
	mod_kinds.o \
	modules.o
average_3d_data.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	module_interface.o \
	modules.o
basic_constants_and_equations_mod.o: \
	mod_kinds.o
biometeorology_mod.o: \
	basic_constants_and_equations_mod.o \
	modules.o \
	mod_kinds.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	radiation_model_mod.o \
	restart_data_mpi_io_mod.o
bulk_cloud_model_mod.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	advec_ws.o \
	advec_s_bc.o \
	advec_s_pw.o \
	advec_s_up.o \
	diffusion_s.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
buoyancy.o: \
	basic_constants_and_equations_mod.o \
	mod_kinds.o \
	modules.o
calc_mean_profile.o: \
	mod_kinds.o \
	modules.o
check_for_restart.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
check_open.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	posix_interface_mod.o
check_parameters.o: \
	basic_constants_and_equations_mod.o \
	init_vertical_profiles.o \
	model_1d_mod.o \
	module_interface.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	pmc_interface_mod.o \
	subsidence_mod.o
chem_emis_domestic_mod.o: \
	chem_emis_generic_mod.o \
	chem_emis_vsrc_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	netcdf_data_input_mod.o \
	mod_kinds.o \
	modules.o \
	palm_date_time_mod.o
chem_emis_generic_mod.o: \
	chem_emis_vsrc_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	netcdf_data_input_mod.o \
	mod_kinds.o \
	modules.o \
	palm_date_time_mod.o
chem_emis_traffic_mod.o: \
	chem_emis_generic_mod.o \
	chem_emis_vsrc_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	netcdf_data_input_mod.o \
	mod_kinds.o \
	modules.o \
	palm_date_time_mod.o
chem_emis_vsrc_mod.o: \
	chem_gasphase_mod.o \
	chem_modules.o \
	mod_kinds.o \
	modules.o
chem_emissions_mod.o: \
	chem_gasphase_mod.o \
	chem_modules.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	surface_mod.o
chemistry_model_mod.o: \
	advec_ws.o \
	chem_emis_domestic_mod.o \
	chem_emis_generic_mod.o \
	chem_emis_traffic_mod.o \
	chem_emis_vsrc_mod.o \
	chem_emissions_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	chem_photolysis_mod.o \
	diffusion_s.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	radiation_model_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
chem_gasphase_mod.o: \
	mod_kinds.o \
	modules.o
chem_modules.o: \
	mod_kinds.o \
	modules.o
chem_photolysis_mod.o: \
	basic_constants_and_equations_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	mod_kinds.o \
	modules.o \
	palm_date_time_mod.o \
	radiation_model_mod.o
close_file.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
compute_vpt.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	mod_kinds.o \
	modules.o
coriolis.o: \
	mod_kinds.o \
	modules.o
cpulog_mod.o: \
	mod_kinds.o \
	modules.o
cpu_statistics.o: \
	mod_kinds.o \
	modules.o
cuda_fft_interfaces.o: \
	modules.o \
	mod_kinds.o
data_log.o: \
	mod_kinds.o \
	modules.o
data_output_binary_module.o: \
	mod_kinds.o
data_output_mask.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	chemistry_model_mod.o \
	cpulog_mod.o \
	diagnostic_output_quantities_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_interface_mod.o \
	salsa_mod.o
data_output_module.o: \
	data_output_binary_module.o \
	data_output_netcdf4_module.o \
	mod_kinds.o
data_output_netcdf4_module.o: \
	mod_kinds.o
data_output_particle_mod.o: \
	data_output_netcdf4_module.o \
	modules.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	shared_memory_io_mod.o
data_output_profiles.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
data_output_spectra.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	spectra_mod.o
data_output_topo_and_surface_setup_mod.o: \
	data_output_module.o \
	netcdf_data_input_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
data_output_tseries.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o
data_output_flight.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	virtual_flight_mod.o
data_output_2d.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	module_interface.o \
	modules.o \
	netcdf_interface_mod.o \
	surface_mod.o \
	urban_surface_mod.o
data_output_3d.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	module_interface.o \
	modules.o \
	netcdf_interface_mod.o \
	plant_canopy_model_mod.o \
	radiation_model_mod.o \
	urban_surface_mod.o
diagnostic_output_quantities_mod.o: \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	restart_data_mpi_io_mod.o \
	surface_layer_fluxes_mod.o
diffusion_s.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
diffusion_u.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
diffusion_v.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
diffusion_w.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o
disturb_field.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	random_function_mod.o \
	random_generator_parallel_mod.o
disturb_heatflux.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	random_generator_parallel_mod.o \
	surface_mod.o
dynamics_mod.o: \
	basic_constants_and_equations_mod.o \
	mod_kinds.o \
	surface_mod.o \
	pmc_interface_mod.o \
	modules.o
exchange_horiz_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
fft_xy_mod.o: \
	cuda_fft_interfaces.o \
	mod_kinds.o \
	modules.o \
	singleton_mod.o \
	temperton_fft_mod.o \
	transpose_mod.o
flow_statistics.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	chem_modules.o \
	cpulog_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	mod_kinds.o \
	module_interface.o \
	modules.o \
	netcdf_interface_mod.o \
	radiation_model_mod.o \
	surface_mod.o
global_min_max.o: \
	mod_kinds.o \
	modules.o
gust_mod.o: \
	mod_kinds.o \
	modules.o
header.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	model_1d_mod.o \
	module_interface.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	ocean_mod.o \
	palm_date_time_mod.o \
	plant_canopy_model_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_interface_mod.o \
	subsidence_mod.o \
	surface_mod.o
indoor_model_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	surface_mod.o \
	urban_surface_mod.o
inflow_turbulence.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
init_3d_model.o: \
	advec_ws.o \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	chem_modules.o \
	cpulog_mod.o \
	disturb_heatflux.o \
	exchange_horiz_mod.o \
	large_scale_forcing_nudging_mod.o \
	model_1d_mod.o \
	mod_kinds.o \
	modules.o \
	multi_agent_system_mod.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	nesting_offl_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	read_restart_data_mod.o \
	palm_date_time_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o \
	surface_data_output_mod.o \
	synthetic_turbulence_generator_mod.o
init_advec.o: \
	mod_kinds.o \
	modules.o
init_coupling.o: \
	mod_kinds.o \
	modules.o
init_grid.o: \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o
init_masks.o: \
	bulk_cloud_model_mod.o \
	mod_kinds.o \
	module_interface.o \
	modules.o \
	netcdf_interface_mod.o
init_pegrid.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	spectra_mod.o \
	synthetic_turbulence_generator_mod.o \
	transpose_mod.o
init_pt_anomaly.o: \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o
init_rankine.o: \
	basic_constants_and_equations_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o
init_slope.o: \
	basic_constants_and_equations_mod.o \
	mod_kinds.o \
	modules.o
init_vertical_profiles.o: \
	modules.o \
	mod_kinds.o
lagrangian_particle_model_mod.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	data_output_particle_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	pmc_interface_mod.o \
	pmc_particle_interface.o \
	netcdf_interface_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
land_surface_model_mod.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	calc_mean_profile.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	surface_mod.o
large_scale_forcing_nudging_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
local_stop.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
local_tremain.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
local_tremain_ini.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
message.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
model_1d_mod.o: \
	basic_constants_and_equations_mod.o \
	mod_kinds.o \
	modules.o \
	time_to_string.o
module_interface.o: \
	mod_kinds.o \
	modules.o \
	biometeorology_mod.o \
	bulk_cloud_model_mod.o \
	chemistry_model_mod.o \
	data_output_module.o \
	data_output_topo_and_surface_setup_mod.o \
	diagnostic_output_quantities_mod.o \
	dynamics_mod.o \
	gust_mod.o \
	indoor_model_mod.o \
	lagrangian_particle_model_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	mod_particle_attributes.o \
	multi_agent_system_mod.o \
	nesting_offl_mod.o \
	ocean_mod.o \
	plant_canopy_model_mod.o \
	poismg_noopt_mod.o \
	radiation_model_mod.o \
	salsa_mod.o \
	spectra_mod.o \
	surface_data_output_mod.o \
	surface_mod.o \
	synthetic_turbulence_generator_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o \
	user_module.o \
	vdi_internal_controls.o \
	virtual_flight_mod.o \
	virtual_measurement_mod.o \
	wind_turbine_model_mod.o
modules.o: \
	mod_kinds.o
multi_agent_system_mod.o: \
	basic_constants_and_equations_mod.o \
	biometeorology_mod.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	random_function_mod.o
mod_kinds.o: \

mod_particle_attributes.o: \
	modules.o \
	mod_kinds.o
nesting_offl_mod.o: \
	chem_modules.o \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	salsa_mod.o
netcdf_data_input_mod.o: \
	chem_modules.o \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o   \
	surface_mod.o
netcdf_interface_mod.o: \
	basic_constants_and_equations_mod.o \
	biometeorology_mod.o \
	chemistry_model_mod.o \
	diagnostic_output_quantities_mod.o \
	gust_mod.o \
	indoor_model_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	ocean_mod.o \
	palm_date_time_mod.o \
	radiation_model_mod.o \
	salsa_mod.o \
	spectra_mod.o \
	turbulence_closure_mod.o \
	urban_surface_mod.o \
	user_module.o
ocean_mod.o: \
	advec_s_pw.o \
	advec_s_up.o \
	advec_ws.o \
	basic_constants_and_equations_mod.o \
	diffusion_s.o \
	init_vertical_profiles.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
outflow_turbulence.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
palm.o: \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	diagnostic_output_quantities_mod.o \
	exchange_horiz_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	module_interface.o \
	modules.o \
	multi_agent_system_mod.o \
	netcdf_data_input_mod.o \
	pmc_interface_mod.o \
	pmc_particle_interface.o \
	restart_data_mpi_io_mod.o \
	surface_layer_fluxes_mod.o \
	surface_data_output_mod.o \
	time_to_string.o \
	write_restart_data_mod.o
palm_date_time_mod.o: \
	mod_kinds.o \
	modules.o
parin.o: \
	chem_modules.o \
	cpulog_mod.o \
	model_1d_mod.o \
	mod_kinds.o \
	module_interface.o \
	modules.o \
	netcdf_interface_mod.o \
	pmc_interface_mod.o \
	progress_bar_mod.o \
	read_restart_data_mod.o \
	turbulence_closure_mod.o
plant_canopy_model_mod.o: \
	bulk_cloud_model_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
pmc_particle_interface.o: \
	mod_particle_attributes.o \
	modules.o \
	pmc_child_mod.o \
	pmc_interface_mod.o \
	pmc_parent_mod.o
pmc_interface_mod.o: \
	bulk_cloud_model_mod.o \
	chem_modules.o \
	chemistry_model_mod.o \
	chem_gasphase_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	pmc_child_mod.o \
	pmc_general_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_mpi_wrapper_mod.o \
	pmc_parent_mod.o \
	salsa_mod.o \
	surface_mod.o
pmc_child_mod.o: \
	mod_kinds.o \
	pmc_general_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_mpi_wrapper_mod.o
pmc_general_mod.o: \
	mod_kinds.o
pmc_handle_communicator_mod.o: \
	mod_kinds.o \
	modules.o \
	pmc_general_mod.o
pmc_mpi_wrapper_mod.o: \
	pmc_handle_communicator_mod.o
pmc_parent_mod.o: \
	pmc_general_mod.o \
	pmc_handle_communicator_mod.o \
	pmc_mpi_wrapper_mod.o
poisfft_mod.o: \
	cpulog_mod.o \
	fft_xy_mod.o \
	mod_kinds.o \
	modules.o \
	transpose_mod.o \
	tridia_solver_mod.o
poismg_mod.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o
poismg_noopt_mod.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o
posix_interface_mod.o: \
	mod_kinds.o
pres.o: \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o \
	poisfft_mod.o \
	poismg_mod.o \
	poismg_noopt_mod.o \
	surface_mod.o
print_1d.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
prognostic_equations.o: \
	advec_s_bc.o \
	advec_s_pw.o \
	advec_s_up.o \
	advec_u_pw.o \
	advec_u_up.o \
	advec_v_pw.o \
	advec_v_up.o \
	advec_w_pw.o \
	advec_ws.o \
	advec_w_up.o \
	buoyancy.o \
	coriolis.o \
	cpulog_mod.o \
	cpulog_mod.o \
	diffusion_s.o \
	diffusion_u.o \
	diffusion_v.o \
	diffusion_w.o \
	large_scale_forcing_nudging_mod.o \
	mod_kinds.o \
	mod_kinds.o \
	module_interface.o \
	modules.o \
	ocean_mod.o \
	plant_canopy_model_mod.o \
	radiation_model_mod.o \
	subsidence_mod.o \
	surface_mod.o \
	surface_mod.o \
	wind_turbine_model_mod.o
progress_bar_mod.o: \
	mod_kinds.o \
	modules.o
radiation_model_mod.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	calc_mean_profile.o \
	exchange_horiz_mod.o \
	mod_particle_attributes.o \
	modules.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	plant_canopy_model_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
random_function_mod.o: \
	mod_kinds.o
random_gauss.o: \
	mod_kinds.o \
	random_function_mod.o \
	random_generator_parallel_mod.o
random_generator_parallel_mod.o: \
	mod_kinds.o \
	modules.o \
	posix_interface_mod.o
read_restart_data_mod.o: \
	bulk_cloud_model_mod.o \
	chem_modules.o \
	cpulog_mod.o \
	gust_mod.o \
	indoor_model_mod.o \
	mod_kinds.o \
	model_1d_mod.o \
	module_interface.o \
	modules.o \
	mod_particle_attributes.o \
	netcdf_interface_mod.o \
	palm_date_time_mod.o \
	radiation_model_mod.o \
	random_function_mod.o\
	random_generator_parallel_mod.o \
	restart_data_mpi_io_mod.o \
	spectra_mod.o \
	surface_mod.o \
	user_module.o \
	virtual_measurement_mod.o
restart_data_mpi_io_mod.o: \
	exchange_horiz_mod.o \
	modules.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	posix_interface_mod.o \
	shared_memory_io_mod.o
run_control.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
salsa_mod.o: \
	modules.o \
	mod_kinds.o \
	advec_ws.o \
	advec_s_pw.o \
	advec_s_up.o \
	basic_constants_and_equations_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	diffusion_s.o \
	exchange_horiz_mod.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	plant_canopy_model_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o \
	netcdf_data_input_mod.o
singleton_mod.o: \
	mod_kinds.o
shared_memory_io_mod.o: \
	modules.o \
	transpose_mod.o
sor.o: \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o
spectra_mod.o: \
	cpulog_mod.o \
	fft_xy_mod.o \
	mod_kinds.o \
	modules.o \
	transpose_mod.o
subsidence_mod.o: \
	mod_kinds.o \
	modules.o
sum_up_3d_data.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	mod_kinds.o \
	module_interface.o \
	modules.o \
	surface_mod.o \
	urban_surface_mod.o
surface_coupler.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
surface_layer_fluxes_mod.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	exchange_horiz_mod.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	urban_surface_mod.o
surface_mod.o: \
	chem_gasphase_mod.o \
	chem_modules.o \
	model_1d_mod.o \
	mod_kinds.o \
	modules.o \
	restart_data_mpi_io_mod.o
surface_data_output_mod.o: \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	netcdf_interface_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
synthetic_turbulence_generator_mod.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	modules.o \
	nesting_offl_mod.o \
	pmc_interface_mod.o \
	random_generator_parallel_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o \
	transpose_mod.o
temperton_fft_mod.o: \
	mod_kinds.o \
	modules.o
time_integration.o: \
	advec_ws.o \
	biometeorology_mod.o \
	bulk_cloud_model_mod.o \
	buoyancy.o \
	calc_mean_profile.o \
	chem_emissions_mod.o \
	chem_gasphase_mod.o \
	chem_modules.o \
	chemistry_model_mod.o \
	cpulog_mod.o \
	data_output_flight.o \
	diagnostic_output_quantities_mod.o \
	disturb_heatflux.o \
	exchange_horiz_mod.o \
	indoor_model_mod.o \
	lagrangian_particle_model_mod.o \
	land_surface_model_mod.o \
	large_scale_forcing_nudging_mod.o \
	mod_kinds.o \
	modules.o \
	multi_agent_system_mod.o \
	nesting_offl_mod.o \
	ocean_mod.o \
	palm_date_time_mod.o \
	pmc_interface_mod.o \
	prognostic_equations.o \
	progress_bar_mod.o \
	radiation_model_mod.o \
	salsa_mod.o \
	spectra_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o \
	surface_data_output_mod.o \
	synthetic_turbulence_generator_mod.o \
	time_to_string.o \
	turbulence_closure_mod.o\
	urban_surface_mod.o \
	virtual_flight_mod.o \
	virtual_measurement_mod.o \
	wind_turbine_model_mod.o
time_integration_spinup.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	diagnostic_output_quantities_mod.o \
	exchange_horiz_mod.o \
	land_surface_model_mod.o \
	indoor_model_mod.o \
	mod_kinds.o \
	modules.o \
	palm_date_time_mod.o \
	pmc_interface_mod.o \
	radiation_model_mod.o \
	surface_data_output_mod.o \
	surface_layer_fluxes_mod.o \
	time_to_string.o \
	urban_surface_mod.o
time_to_string.o: \
	mod_kinds.o
timestep.o: \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
timestep_scheme_steering.o: \
	mod_kinds.o \
	modules.o
transpose_mod.o: \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o
tridia_solver_mod.o: \
	basic_constants_and_equations_mod.o \
	mod_kinds.o \
	modules.o \
	transpose_mod.o
turbulence_closure_mod.o: \
	advec_s_bc.o \
	advec_s_pw.o \
	advec_s_up.o \
	advec_ws.o \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	exchange_horiz_mod.o \
	mod_kinds.o \
	mod_particle_attributes.o \
	modules.o \
	ocean_mod.o \
	plant_canopy_model_mod.o \
	pmc_interface_mod.o \
	surface_layer_fluxes_mod.o \
	surface_mod.o
urban_surface_mod.o: \
	basic_constants_and_equations_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	palm_date_time_mod.o \
	radiation_model_mod.o \
	restart_data_mpi_io_mod.o \
	surface_mod.o
user_data_output_mask.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_flight.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_init_3d_model.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	user_module.o
user_init_flight_mod.o: \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	user_module.o
user_init_grid.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_init_land_surface.o: \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	surface_mod.o
user_init_plant_canopy.o: \
	mod_kinds.o \
	modules.o \
	plant_canopy_model_mod.o \
	user_module.o
user_init_radiation.o: \
	mod_kinds.o \
	modules.o \
	radiation_model_mod.o \
	user_module.o
user_init_urban_surface.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	urban_surface_mod.o \
	user_module.o
user_lpm_advec.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_lpm_init.o: \
	mod_kinds.o \
	modules.o \
	user_module.o
user_module.o: \
	mod_kinds.o \
	modules.o \
	surface_mod.o \
	nesting_offl_mod.o \
	netcdf_data_input_mod.o \
	plant_canopy_model_mod.o
user_spectra.o: \
	mod_kinds.o \
	modules.o \
	spectra_mod.o \
	user_module.o
vdi_internal_controls.o: \
	mod_kinds.o \
	modules.o \
	pmc_interface_mod.o
virtual_flight_mod.o: \
	basic_constants_and_equations_mod.o \
	bulk_cloud_model_mod.o \
	cpulog_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_interface_mod.o \
	restart_data_mpi_io_mod.o \
	user_flight.o \
	user_init_flight_mod.o
virtual_measurement_mod.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	chem_modules.o \
	chem_gasphase_mod.o \
	data_output_module.o \
	land_surface_model_mod.o \
	mod_kinds.o \
	modules.o \
	netcdf_data_input_mod.o \
	radiation_model_mod.o \
	surface_mod.o \
	urban_surface_mod.o
wind_turbine_model_mod.o: \
	basic_constants_and_equations_mod.o \
	cpulog_mod.o \
	data_output_module.o \
	exchange_horiz_mod.o \
	netcdf_data_input_mod.o \
	mod_kinds.o \
	modules.o \
	restart_data_mpi_io_mod.o
wrd_write_string.o: \
	mod_kinds.o
write_restart_data_mod.o: \
	bulk_cloud_model_mod.o \
	gust_mod.o \
	indoor_model_mod.o \
	mod_kinds.o \
	model_1d_mod.o \
	module_interface.o \
	modules.o \
	mod_particle_attributes.o \
	netcdf_interface_mod.o \
	palm_date_time_mod.o \
	radiation_model_mod.o \
	random_function_mod.o \
	random_generator_parallel_mod.o \
	restart_data_mpi_io_mod.o \
	spectra_mod.o \
	surface_mod.o \
	user_module.o \
	virtual_measurement_mod.o


$(OBJS1): combine_plot_fields.f90
	$(F90_SER) $(F90FLAGS) $(COPT) -c combine_plot_fields.f90

$(PROG1): $(OBJS1)
	$(F90_SER) -o $(PROG1) $(OBJS1) $(LDFLAGS)
