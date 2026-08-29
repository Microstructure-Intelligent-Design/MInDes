#include "LatticeBoltzmann.h"
namespace pf {
	namespace lattice_boltzmann {
		namespace default_functions {
			// fluid flow
			// fi = fi + mi
			static void collision_SRT_d2q9(long long x, long long y, long long z) {
				LBMPoint& point = external_physical_field::lbm_field(x, y, z);
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double tau = lbm_boundary_condition::tau(x, y, z);
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++) {
						double f_eq = f_eq_i(index, point.pressure, point.F_MACRO, point.velocity);
						point.M[index] = (f_eq - point.F[index]) / tau + time_parameters::delt_t * fluid_source_i(x, y, z, tau, index);
					}
				}
				else {
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++) {
						point.M[index] = 0.0;
					}
				}
			}
			static void collision_SRT_d3q19(long long x, long long y, long long z) {
				LBMPoint& point = external_physical_field::lbm_field(x, y, z);
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double tau = lbm_boundary_condition::tau(x, y, z);
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++) {
						double f_eq = f_eq_i(index, point.pressure, point.F_MACRO, point.velocity);
						point.M[index] = (f_eq - point.F[index]) / tau + time_parameters::delt_t * fluid_source_i(x, y, z, tau, index);
					}
				}
				else {
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++) {
						point.M[index] = 0.0;
					}
				}
			}
		}
		void init_distribution_functions_d2q9() {
#pragma omp parallel for
			for (long long x = external_physical_field::lbm_field.COMP_X_BGN(); x <= external_physical_field::lbm_field.COMP_X_END(); x++)
				for (long long y = external_physical_field::lbm_field.COMP_Y_BGN(); y <= external_physical_field::lbm_field.COMP_Y_END(); y++)
					for (long long z = external_physical_field::lbm_field.COMP_Z_BGN(); z <= external_physical_field::lbm_field.COMP_Z_END(); z++) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double density = lbm_boundary_condition::density(x, y, z);
						point.mass = density;
						point.F_MACRO = density;
						for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++)
							point.F[index] = LBM_d2q9_w[index] * density;
					}
		}
		void init_distribution_functions_d3q19() {
#pragma omp parallel for
			for (long long x = external_physical_field::lbm_field.COMP_X_BGN(); x <= external_physical_field::lbm_field.COMP_X_END(); x++)
				for (long long y = external_physical_field::lbm_field.COMP_Y_BGN(); y <= external_physical_field::lbm_field.COMP_Y_END(); y++)
					for (long long z = external_physical_field::lbm_field.COMP_Z_BGN(); z <= external_physical_field::lbm_field.COMP_Z_END(); z++) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double density = lbm_boundary_condition::density(x, y, z);
						point.mass = density;
						point.F_MACRO = density;
						for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++)
							point.F[index] = LBM_d3q19_w[index] * density;
					}
		}
		void lbm_properties_automatically_change() {
			lbm_boundary_condition::lbm_properties_automatically_change();
			lbm_equilibrium_distribution_function::lbm_properties_automatically_change();
			lbm_source::lbm_properties_automatically_change();
			lbm_macro_variable::lbm_properties_automatically_change();
		}
		void init() {
			fluid_lbm_solver.init(external_physical_field::lbm_field, mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ,
				mesh_parameters::delt_r, mesh_parameters::x_down, mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up,
				mesh_parameters::z_down, mesh_parameters::z_up);
			infile_reader::read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.max_iterate_steps", max_iterate_steps, true);
			if (infile_reader::read_bool_value("Postprocess.FluidDynamics.LatticeBoltzmann.debug_solver", debug_solver, true))
				infile_reader::read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.debug_output_step", debug_output_step, true);
			infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.momentum_accuracy", momentum_accuracy, true);
			lbm_boundary_condition::init(fluid_lbm_solver);
			lbm_source::init(fluid_lbm_solver);
			lbm_equilibrium_distribution_function::init(fluid_lbm_solver);
			lbm_macro_variable::init(fluid_lbm_solver);
			InputFileReader::get_instance()->read_bool_value("Postprocess.FluidDynamics.LatticeBoltzmann.two_phase_flow", is_two_phase_flow, true);
			if (is_two_phase_flow) {
				// 
			}
			fluid_source_i = lbm_source::fluid_source_i;
			f_eq_i = lbm_equilibrium_distribution_function::f_eq_i;
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				fluid_lbm_solver._collision = default_functions::collision_SRT_d2q9;
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				fluid_lbm_solver._collision = default_functions::collision_SRT_d3q19;
			}
		}
		void exec_pre() {
			stringstream output;
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				init_distribution_functions_d2q9();
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				init_distribution_functions_d3q19();
			}
			lbm_boundary_condition::cal_fluid_domain();
			fluid_lbm_solver.do_boundary_condition();
			fluid_lbm_solver.cal_macro_variables();
			if (is_two_phase_flow) {
				//
			}
			if (debug_solver) {
				output << "> Fluid field solver debug:" << std::endl;
				WriteLog(output.str());
			}
			MACRO_MAX_VARIATION MAX_CHANGE;
			int istep = 0;
			for (istep = 1; istep <= max_iterate_steps; istep++) {
				fluid_lbm_solver.do_collision();
				fluid_lbm_solver.do_streaming();
				fluid_lbm_solver.do_boundary_condition();
				MAX_CHANGE = fluid_lbm_solver.cal_macro_variables();
				if (is_two_phase_flow) {
					// field_lbm_two_phase_solver.collision();
					// field_lbm_two_phase_solver.streaming();
					// field_lbm_two_phase_solver.boundary_condition();
					// field_lbm_two_phase_solver.cal_macro_variables(MAX_F_MACRO_CHANGE);
				}
				if (debug_solver && (istep % debug_output_step == 0 || istep == 1)) {
					output.str("");
					output << "> Fluid field iterate step:" << istep << std::endl;
					output << "		MAX_MOMENTUM_CHANGE_X = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[0] << std::endl;
					output << "		MAX_MOMENTUM_CHANGE_Y = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[1] << std::endl;
					output << "		MAX_MOMENTUM_CHANGE_Z = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[2] << std::endl;
					// if (is_two_phase_flow)
					// 	output << "		TWO_PHASE_VARIATION   = " << MAX_F_MACRO_CHANGE << endl;
					WriteLog(output.str());
				}
				if (MAX_CHANGE.FV_MACRO_MAX_VARIATION[0] < momentum_accuracy
					&& MAX_CHANGE.FV_MACRO_MAX_VARIATION[1] < momentum_accuracy
					&& MAX_CHANGE.FV_MACRO_MAX_VARIATION[2] < momentum_accuracy)
					break;
			}
			istep--;
			output.str("");
			output << "> Fluid field solver:" << std::endl;
			output << "		MAX_ITERATE_TIMES = " << istep << std::endl;
			output << "		MAX_MOMENTUM_CHANGE_X = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[0] << std::endl;
			output << "		MAX_MOMENTUM_CHANGE_Y = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[1] << std::endl;
			output << "		MAX_MOMENTUM_CHANGE_Z = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[2] << std::endl;
			// if (is_two_phase_flow)
			// 	output << "		TWO_PHASE_VARIATION   = " << MAX_F_MACRO_CHANGE << endl;
			WriteLog(output.str());
		}
		std::string exec_loop() {
			lbm_properties_automatically_change();
			stringstream report;
			MACRO_MAX_VARIATION MAX_CHANGE;
			lbm_boundary_condition::cal_fluid_domain();
			int istep = 0;
			for (istep = 1; istep <= max_iterate_steps; istep++) {
				fluid_lbm_solver.do_collision();
				fluid_lbm_solver.do_streaming();
				fluid_lbm_solver.do_boundary_condition();
				MAX_CHANGE = fluid_lbm_solver.cal_macro_variables();
				if (is_two_phase_flow) {
					// field_lbm_two_phase_solver.collision();
					// field_lbm_two_phase_solver.streaming();
					// field_lbm_two_phase_solver.boundary_condition();
					// field_lbm_two_phase_solver.cal_macro_variables(MAX_F_MACRO_CHANGE);
				}
				if (MAX_CHANGE.FV_MACRO_MAX_VARIATION[0] < momentum_accuracy
					&& MAX_CHANGE.FV_MACRO_MAX_VARIATION[1] < momentum_accuracy
					&& MAX_CHANGE.FV_MACRO_MAX_VARIATION[2] < momentum_accuracy)
					break;
			}
			istep--;
			report.str("");
			report << "> Fluid field solver:" << std::endl;
			report << "		MAX_ITERATE_TIMES = " << istep << std::endl;
			report << "		MAX_MOMENTUM_CHANGE_X = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[0] << std::endl;
			report << "		MAX_MOMENTUM_CHANGE_Y = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[1] << std::endl;
			report << "		MAX_MOMENTUM_CHANGE_Z = " << MAX_CHANGE.FV_MACRO_MAX_VARIATION[2] << std::endl;
			// if (is_two_phase_flow)
			// 	report << "		TWO_PHASE_VARIATION   = " << MAX_F_MACRO_CHANGE << endl;
			return report.str();
		}
		void deinit() {
			fluid_lbm_solver.free();
			lbm_boundary_condition::deinit();
		}
		void write_velocity(std::ofstream& fout) {
			fout << "<DataArray type = \"Float64\" Name = \"fluid_velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
			for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						LBMPoint& point = external_physical_field::lbm_field(i, j, k);
						fout << point.velocity[0] << " "
							<< point.velocity[1] << " "
							<< point.velocity[2] << std::endl;
					}
			fout << "</DataArray>" << std::endl;
		}
		void write_abs_velocity(std::ofstream& fout) {
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_abs_velocity" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
			for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						fout << external_physical_field::lbm_field(i, j, k).velocity.abs() << std::endl;
					}
			fout << "</DataArray>" << std::endl;
		}
		void write_pressure(std::ofstream& fout) {
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_pressure" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
			for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						fout << external_physical_field::lbm_field(i, j, k).pressure << std::endl;
					}
			fout << "</DataArray>" << std::endl;
		}
		void write_density(std::ofstream& fout) {
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_density" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
			for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						fout << external_physical_field::lbm_field(i, j, k).mass << std::endl;
					}
			fout << "</DataArray>" << std::endl;
		}
	}
}