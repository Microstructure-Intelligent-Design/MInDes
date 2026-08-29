#include "SourceF.h"
namespace pf {
	namespace lbm_source {
		namespace force_funcs {
			// fluid force
			Vector3 Fluid_Force_Gravity(long long x, long long y, long long z) {
				return gravitational_acceleration * (external_physical_field::lbm_field(x, y, z).mass - ref_density);
			}
			Vector3 Fluid_Force_H_Liang_Surface_Tension(long long x, long long y, long long z) {
				double phi = external_physical_field::lbm_field(x, y, z).F_MACRO,
					lap_phi = (external_physical_field::lbm_field(x + 1, y, z).F_MACRO + external_physical_field::lbm_field(x - 1, y, z).F_MACRO
						+ external_physical_field::lbm_field(x, y + 1, z).F_MACRO + external_physical_field::lbm_field(x, y - 1, z).F_MACRO
						+ external_physical_field::lbm_field(x, y, z + 1).F_MACRO + external_physical_field::lbm_field(x, y, z - 1).F_MACRO - 6.0 * phi)
					/ external_physical_field::lbm_field.DeltR() / external_physical_field::lbm_field.DeltR();
				Vector3 grad_phi = Vector3((external_physical_field::lbm_field(x + 1, y, z).F_MACRO - external_physical_field::lbm_field(x - 1, y, z).F_MACRO) / 2 / external_physical_field::lbm_field.DeltR(),
					(external_physical_field::lbm_field(x, y + 1, z).F_MACRO - external_physical_field::lbm_field(x, y - 1, z).F_MACRO) / 2 / external_physical_field::lbm_field.DeltR(),
					(external_physical_field::lbm_field(x, y, z + 1).F_MACRO - external_physical_field::lbm_field(x, y, z - 1).F_MACRO) / 2 / external_physical_field::lbm_field.DeltR());
				return grad_phi * (4.0 * beta * phi * (phi - 1) * (phi - 0.5) - kappa * lap_phi);
			}
			Vector3 Fluid_Force_Thermal_Expansion(long long x, long long y, long long z) {
				return gravitational_acceleration * thermal_expansion_parameter
					* (ref_temp - main_field::temperature_field(x, y, z)) * external_physical_field::lbm_field(x, y, z).mass;
			}

			// fluid force model
			static double Fluid_Modified_LGA_d2q9_Force_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
				Vector3 force(0.0, 0.0, 0.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(x, y, z);
				return LBM_d2q9_w[LBM_i] * (LBM_d2q9_w_vec[LBM_i] * force) / Cs2;
			}
			static double Fluid_Modified_LGA_d3q19_Force_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
				Vector3 force(0.0, 0.0, 0.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(x, y, z);
				return LBM_d3q19_w[LBM_i] * (LBM_d3q19_w_vec[LBM_i] * force) / Cs2;
			}
			static double Fluid_Modified_GZS_d2q9_Force_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = LBM_d2q9_w_vec[LBM_i], u = external_physical_field::lbm_field(x, y, z).velocity;
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(x, y, z);
				return (1.0 - 1.0 / 2.0 / tau) * LBM_d2q9_w[LBM_i] * (((ci - u) / Cs2 + ci * (ci * u) / Cs4 / 2.0) * force);
			}
			static double Fluid_Modified_GZS_d3q19_Force_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = LBM_d3q19_w_vec[LBM_i], u = external_physical_field::lbm_field(x, y, z).velocity;
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(x, y, z);
				return (1.0 - 1.0 / 2.0 / tau) * LBM_d3q19_w[LBM_i] * (((ci - u) / Cs2 + ci * (ci * u) / Cs4 / 2.0) * force);
			}
			static double Fluid_Modified_HL_d2q9_Force_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = LBM_d2q9_w_vec[LBM_i], u = external_physical_field::lbm_field(x, y, z).velocity,
					grad_density = Vector3((external_physical_field::lbm_field(x + 1, y, z).mass - external_physical_field::lbm_field(x - 1, y, z).mass) / 2 / external_physical_field::lbm_field.DeltR(),
						(external_physical_field::lbm_field(x, y + 1, z).mass - external_physical_field::lbm_field(x, y - 1, z).mass) / 2 / external_physical_field::lbm_field.DeltR(),
						(external_physical_field::lbm_field(x, y, z + 1).mass - external_physical_field::lbm_field(x, y, z - 1).mass) / 2 / external_physical_field::lbm_field.DeltR()) * (-1.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(x, y, z);
				return (1.0 - 1.0 / 2.0 / tau) * LBM_d2q9_w[LBM_i] * ((ci * force) / Cs2 + (u * ci) * (grad_density * ci) / Cs2);
			}
			static double Fluid_Modified_HL_d3q19_Force_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = LBM_d3q19_w_vec[LBM_i], u = external_physical_field::lbm_field(x, y, z).velocity,
					grad_density = Vector3((external_physical_field::lbm_field(x + 1, y, z).mass - external_physical_field::lbm_field(x - 1, y, z).mass) / 2 / external_physical_field::lbm_field.DeltR(),
						(external_physical_field::lbm_field(x, y + 1, z).mass - external_physical_field::lbm_field(x, y - 1, z).mass) / 2 / external_physical_field::lbm_field.DeltR(),
						(external_physical_field::lbm_field(x, y, z + 1).mass - external_physical_field::lbm_field(x, y, z - 1).mass) / 2 / external_physical_field::lbm_field.DeltR()) * (-1.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(x, y, z);
				return (1.0 - 1.0 / 2.0 / tau) * LBM_d3q19_w[LBM_i] * ((ci * force) / Cs2 + (u * ci) * (grad_density * ci) / Cs2);
			}
			static double Fluid_Modified_two_phase_d2q9_Force_i(long long x, long long y, long long z, double tau, Vector3 prefactor, size_t LBM_i) {
				Vector3 ci = LBM_d2q9_w_vec[LBM_i];
				return (1.0 - 1.0 / 2.0 / tau) * LBM_d2q9_w[LBM_i] * (ci * prefactor) / Cs2;
			}
			static double Fluid_Modified_two_phase_d3q19_Force_i(long long x, long long y, long long z, double tau, Vector3 prefactor, size_t LBM_i) {
				Vector3 ci = LBM_d3q19_w_vec[LBM_i];
				return (1.0 - 1.0 / 2.0 / tau) * LBM_d3q19_w[LBM_i] * (ci * prefactor) / Cs2;
			}
			void load_forces() {
				fluid_force_list.clear();
				WriteDebugFile("# Postprocess.FluidDynamics.LatticeBoltzmann.force = () \n");
				WriteDebugFile("#             0 - ThermalExpansion, 1 - Gravity, 2 - H_Liang_SurfaceTension \n");
				std::string force_key = "Postprocess.FluidDynamics.LatticeBoltzmann.force", force_input = "()";
				infile_reader::read_string_value(force_key, force_input, true);
				std::vector<input_value> force_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, force_key, force_input, true);
				for (auto force = force_value.begin(); force < force_value.end(); force++) {
					bool is_already_load = false;
					for (auto old_force = force_value.begin(); old_force < force; old_force++)
						if (force->int_value == old_force->int_value)
							is_already_load = true;
					if (!is_already_load) {
						switch (LBM_Force_Type(force->int_value))
						{
						case LBM_FM_ThermalExpansion:
							fluid_force_list.push_back(Fluid_Force_Thermal_Expansion);
							infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.ThermalExpansionForce.reference_temperature", ref_temp, true);
							infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.ThermalExpansionForce.thermal_expansion_parameter", thermal_expansion_parameter, true);
							break;
						case LBM_FM_Gravity:
							infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.Gravity.reference_density", ref_density, true);
							fluid_force_list.push_back(Fluid_Force_Gravity);
							break;
						case LBM_FM_H_Liang_SurfaceTension:
							fluid_force_list.push_back(Fluid_Force_H_Liang_Surface_Tension);
							break;
						default:
							break;
						}
					}
				}
			}
		}

		// source model
		double fluid_source_i(long long x, long long y, long long z, double tau, size_t LBM_i) {
			double buff = 0.0;
			for (auto func = fluid_source_list.begin(); func < fluid_source_list.end(); func++)
				buff += (*func)(x, y, z, tau, LBM_i);
			return buff;
		}
		double fluid_two_phase_source_i(long long x, long long y, long long z, double tau, Vector3 prefactor, size_t LBM_i) {
			return fluid_two_phase_source_list(x, y, z, tau, prefactor, LBM_i);
		}
		void lbm_properties_automatically_change() {
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
		}
		void init(LBM& fluid_lbm_solver) {
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			// load source
			WriteDebugFile("# Postprocess.FluidDynamics.LatticeBoltzmann.source = () \n");
			WriteDebugFile("#             0 - Forces \n");
			std::string source_key = "Postprocess.FluidDynamics.LatticeBoltzmann.source", source_input = "()";
			infile_reader::read_string_value(source_key, source_input, true);
			std::vector<input_value> source_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, source_key, source_input, true);
			int force_model_type = LBM_Force_Term_Model::LBM_FTM_NONE;
			for (auto source = source_value.begin(); source < source_value.end(); source++) {
				bool is_already_load = false;
				for (auto old_source = source_value.begin(); old_source < source; old_source++)
					if (source->int_value == old_source->int_value)
						is_already_load = true;
				if (!is_already_load) {
					int gravity_direction = 0;
					switch (LBM_Source_Type(source->int_value))
					{
					case pf::LBM_ST_Forces:
						WriteDebugFile("# .gravity_direction : 0 - x_down, 1 - x_up, 2 - y_down, 3 - y_up, 4 - z_dowm, 5 - z_up \n");
						infile_reader::read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.Force.gravity_direction", gravity_direction, true);
						switch (gravity_direction)
						{
						case 0:  // down_x
							force_funcs::gravitational_acceleration = Vector3(-9.8, 0.0, 0.0);
							break;
						case 1:  // up_x
							force_funcs::gravitational_acceleration = Vector3(9.8, 0.0, 0.0);
							break;
						case 2:  // down_y
							force_funcs::gravitational_acceleration = Vector3(0.0, -9.8, 0.0);
							break;
						case 3:  // up_y
							force_funcs::gravitational_acceleration = Vector3(0.0, 9.8, 0.0);
							break;
						case 4:  // down_z
							force_funcs::gravitational_acceleration = Vector3(0.0, 0.0, -9.8);
							break;
						case 5:  // up_z
							force_funcs::gravitational_acceleration = Vector3(0.0, 0.0, 9.8);
							break;
						default:
							// down_x
							force_funcs::gravitational_acceleration = Vector3(-9.8, 0.0, 0.0);
							break;
						}
						WriteDebugFile("# .LatticeBoltzmann.force_model : 0 - NONE, 1 - LGA, 2 - ZS_Guo, 3 - H_Liang \n");
						infile_reader::read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.force_model", force_model_type, true);
						if (force_model_type == LBM_Force_Term_Model::LBM_FTM_LGA) {
							if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_LGA_d2q9_Force_i);
							else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_LGA_d3q19_Force_i);
						}
						else if (force_model_type == LBM_Force_Term_Model::LBM_FTM_ZS_Guo) {
							if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_GZS_d2q9_Force_i);
							else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_GZS_d3q19_Force_i);
						}
						else if (force_model_type == LBM_Force_Term_Model::LBM_FTM_H_Liang) {
							if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_HL_d2q9_Force_i);
							else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_HL_d3q19_Force_i);
						}
						force_funcs::load_forces();
						break;
					default:
						break;
					}
				}
			}
		}
		void init_two_phase_solver(LBM& field_lbm_two_phase_solver) {
			WriteDebugFile("# kappa = 3.0 / 2.0 * surface_tension * interface_thickness \n");
			WriteDebugFile("# beta  = 12.0 * surface_tension / interface_thickness \n");

			infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.surface_tension", surface_tension, true);
			infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.interface_thickness", interface_thickness, true);
			beta = 12.0 * surface_tension / interface_thickness;
			kappa = 3.0 / 2.0 * surface_tension * interface_thickness;

			if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
				fluid_two_phase_source_list = force_funcs::Fluid_Modified_two_phase_d2q9_Force_i;
			else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
				fluid_two_phase_source_list = force_funcs::Fluid_Modified_two_phase_d3q19_Force_i;
		}
	}
}