#include "MacroVar.h"
namespace pf {
	namespace lbm_macro_variable {
		namespace macro_variable_funcs {
			void load_forces() {
				fluid_force_list.clear();
				std::string force_key = "Postprocess.FluidDynamics.LatticeBoltzmann.force", force_input = "()";
				infile_reader::read_string_value(force_key, force_input, false);
				std::vector<input_value> force_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, force_key, force_input, false);
				for (auto force = force_value.begin(); force < force_value.end(); force++) {
					bool is_already_load = false;
					for (auto old_force = force_value.begin(); old_force < force; old_force++)
						if (force->int_value == old_force->int_value)
							is_already_load = true;
					if (!is_already_load) {
						switch (LBM_Force_Type(force->int_value))
						{
						case LBM_FM_ThermalExpansion:
							fluid_force_list.push_back(lbm_source::force_funcs::Fluid_Force_Thermal_Expansion);
							break;
						case LBM_FM_Gravity:
							fluid_force_list.push_back(lbm_source::force_funcs::Fluid_Force_Gravity);
							break;
						case LBM_FM_H_Liang_SurfaceTension:
							fluid_force_list.push_back(lbm_source::force_funcs::Fluid_Force_H_Liang_Surface_Tension);
							break;
						default:
							break;
						}
					}
				}
			}

			static void cal_macro_variables_d2q9_standard(long long x, long long y, long long z) {
				LBMPoint& point = main_field::lbm_field(x, y, z);
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					point.F_MACRO = 0.0;
					point.FV_MACRO.set_to_zero();
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++) {
						point.F_MACRO += point.F[index];
						point.FV_MACRO += LBM_d2q9_w_vec[index] * point.F[index];
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(x, y, z);
					point.FV_MACRO += force * time_parameters::delt_t * 0.5; // force
					point.velocity = point.FV_MACRO / point.F_MACRO;
					point.pressure = point.F_MACRO * Cs2;
					point.mass = point.F_MACRO;
				}
				else {
					point.F_MACRO = 0.0;
					point.FV_MACRO.set_to_zero();
					point.velocity.set_to_zero();
					point.pressure = 0.0;
					point.mass = 0.0;
				}
			}
			static void cal_macro_variables_d3q19_standard(long long x, long long y, long long z) {
				LBMPoint& point = main_field::lbm_field(x, y, z);
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					point.F_MACRO = 0.0;
					point.FV_MACRO.set_to_zero();
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++) {
						point.F_MACRO += point.F[index];
						point.FV_MACRO += LBM_d3q19_w_vec[index] * point.F[index];
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(x, y, z);
					point.FV_MACRO += force * time_parameters::delt_t * 0.5; // force
					point.velocity = point.FV_MACRO / point.F_MACRO;
					point.pressure = point.F_MACRO * Cs2;
					point.mass = point.F_MACRO;
				}
				else {
					point.F_MACRO = 0.0;
					point.FV_MACRO.set_to_zero();
					point.velocity.set_to_zero();
					point.pressure = 0.0;
					point.mass = 0.0;
				}
			}
			static double s_0_d2q9_two_phase_flow(Vector3& U) {
				return -w0_d2q9 * (U * U) / Cs2 / 2.0;
			}
			static void cal_macro_variables_d2q9_H_Liang(long long x, long long y, long long z) {
				LBMPoint& point = main_field::lbm_field(x, y, z);
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double buff = 0.0, density = lbm_boundary_condition::density(x, y, z);
					point.F_MACRO = density;
					point.FV_MACRO.set_to_zero();
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++) {
						buff += point.F[index];
						point.FV_MACRO += LBM_d2q9_w_vec[index] * point.F[index];
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(x, y, z);
					point.FV_MACRO += force * time_parameters::delt_t * 0.5; // force
					Vector3 Velocity = point.FV_MACRO / density;
					point.velocity = Velocity;
					point.mass = density;
					Vector3 grad_density = Vector3((main_field::lbm_field(x + 1, y, z).mass - main_field::lbm_field(x - 1, y, z).mass) / 2 / main_field::lbm_field.DeltR(),
						(main_field::lbm_field(x, y + 1, z).mass - main_field::lbm_field(x, y - 1, z).mass) / 2 / main_field::lbm_field.DeltR(),
						(main_field::lbm_field(x, y, z + 1).mass - main_field::lbm_field(x, y, z - 1).mass) / 2 / main_field::lbm_field.DeltR());
					buff -= point.F[LBM_0];
					point.pressure = Cs2 / (1.0 - w0_d2q9)
						* (buff + time_parameters::delt_t * 0.5 * (Velocity * grad_density) + density * s_0_d2q9_two_phase_flow(Velocity));
				}
				else {
					point.F_MACRO = 0.0;
					point.FV_MACRO.set_to_zero();
					point.velocity.set_to_zero();
					point.pressure = 0.0;
					point.mass = 0.0;
				}
			}
			static double s_0_d3q19_two_phase_flow(Vector3& U) {
				return -w0_d3q19 * (U * U) / Cs2 / 2.0;
			}
			static void cal_macro_variables_d3q19_H_Liang(long long x, long long y, long long z) {
				LBMPoint& point = main_field::lbm_field(x, y, z);
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double buff = 0.0, density = lbm_boundary_condition::density(x, y, z);
					point.F_MACRO = density;
					point.FV_MACRO.set_to_zero();
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++) {
						buff += point.F[index];
						point.FV_MACRO += LBM_d3q19_w_vec[index] * point.F[index];
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(x, y, z);
					point.FV_MACRO += force * time_parameters::delt_t * 0.5; // force
					Vector3 Velocity = point.FV_MACRO / density;
					point.velocity = Velocity;
					point.mass = density;
					Vector3 grad_density = Vector3((main_field::lbm_field(x + 1, y, z).mass - main_field::lbm_field(x - 1, y, z).mass) / 2 / main_field::lbm_field.DeltR(),
						(main_field::lbm_field(x, y + 1, z).mass - main_field::lbm_field(x, y - 1, z).mass) / 2 / main_field::lbm_field.DeltR(),
						(main_field::lbm_field(x, y, z + 1).mass - main_field::lbm_field(x, y, z - 1).mass) / 2 / main_field::lbm_field.DeltR());
					buff -= point.F[LBM_0];
					point.pressure = Cs2 / (1.0 - w0_d3q19)
						* (buff + time_parameters::delt_t * 0.5 * (Velocity * grad_density) + density * s_0_d3q19_two_phase_flow(Velocity));
				}
				else {
					point.F_MACRO = 0.0;
					point.FV_MACRO.set_to_zero();
					point.velocity.set_to_zero();
					point.pressure = 0.0;
					point.mass = 0.0;
				}
			}
			static void cal_macro_variables_d2q9_two_phase(long long x, long long y, long long z) {
				LBMPoint& point = main_field::lbm_field(x, y, z);
				point.F_MACRO = 0.0;
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++)
						point.F_MACRO += point.F[index];
				}
			}
			static void cal_macro_variables_d3q19_two_phase(long long x, long long y, long long z) {
				LBMPoint& point = main_field::lbm_field(x, y, z);
				point.F_MACRO = 0.0;
				if (point.fluid_region >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++)
						point.F_MACRO += point.F[index];
				}
			}
		}

		void lbm_properties_automatically_change() {
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
		}

		void init(LBM& fluid_lbm_solver) {
			// init parameters
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d2q9_standard;
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d3q19_standard;
			}
			// init cal_macro_variables function
			macro_variable_funcs::load_forces();
			// init solver
			fluid_lbm_solver._cal_macro_variables = cal_macro_variables;
		}

		void init_two_phase(LBM& fluid_lbm_solver, LBM& field_lbm_two_phase_solver) {
			// init parameters
			if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				cal_macro_variables_two_phase = macro_variable_funcs::cal_macro_variables_d2q9_two_phase;
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d2q9_H_Liang;
			}
			else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				cal_macro_variables_two_phase = macro_variable_funcs::cal_macro_variables_d3q19_two_phase;
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d3q19_H_Liang;
			}
			// init solver
			fluid_lbm_solver._cal_macro_variables = cal_macro_variables;
			field_lbm_two_phase_solver._cal_macro_variables = cal_macro_variables_two_phase;
		}
	}
}