#include "PlasticSolver.h"
namespace pf {
	namespace plastic_solver {
		namespace default_functions {
			static void cal_plastic_parameters_norm(long long x, long long y, long long z, 
				double& yield_stress, double& hardening_modulus, double& shear_modulus) {
				Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
				for (size_t index = 0; index < phi_number; index++) {
					yield_stress += phi[index] * phi_yield_stress[index];
					hardening_modulus += phi[index] * phi_hardening_modulus[index];
					shear_modulus += phi[index] * phi_shear_modulus[index];
				}
			}
			static void cal_plasticstrain(long long x, long long y, long long z, vStrain& eigenstrain) {
				eigenstrain += external_physical_field::plastic_field(x, y, z).PlasticStrain;
			}
		}

		double Prandtl_Reuss_formula(double deviatoric_part_of_mises_stress_norm, double yield_stress
			, double hardening_modulus, double ave_plastic_strain) {
			// sqrt(2.0 / 3.0) = 0.8164965809
			return deviatoric_part_of_mises_stress_norm - 0.8164965809 * (yield_stress + hardening_modulus * ave_plastic_strain);
		}

		void solve_plastic_flow(double& max_dplastic_strain, double& max_dave_plastic_strain) {
			// sqrt(2.0 / 3.0) = 0.8164965809
			max_dplastic_strain = 0.0;
			max_dave_plastic_strain = 0.0;
#pragma omp parallel for
			for (long long x = external_physical_field::plastic_field.COMP_X_BGN(); x <= external_physical_field::plastic_field.COMP_X_END(); x++)
				for (long long y = external_physical_field::plastic_field.COMP_Y_BGN(); y <= external_physical_field::plastic_field.COMP_Y_END(); y++)
					for (long long z = external_physical_field::plastic_field.COMP_Z_BGN(); z <= external_physical_field::plastic_field.COMP_Z_END(); z++) {
						ElasticPoint& elastic_point = external_physical_field::elastic_field(x, y, z);
						PlasticPoint& plastic_point = external_physical_field::plastic_field(x, y, z);
						Vector6 deviatoric_stress = elastic_point.Stress;
						double mean_stress = (deviatoric_stress[0] + deviatoric_stress[1] + deviatoric_stress[2]) / 3.0, 
							yield_stress = 0.0, hardening_modulus = 0.0, shear_modulus = 0.0;
						deviatoric_stress[0] -= mean_stress;
						deviatoric_stress[1] -= mean_stress;
						deviatoric_stress[2] -= mean_stress;
						double deviatoric_stress_norm = deviatoric_stress.norm();
						// - calculate plastic parameters
						cal_plastic_parameters(x, y, z, yield_stress, hardening_modulus, shear_modulus);
						// - 
						double plastic_judgement = Prandtl_Reuss_formula(deviatoric_stress_norm, yield_stress, hardening_modulus, plastic_point.AvePlasticStrain);
						// - 
						Vector6 old_plastic_strain(plastic_point.PlasticStrain);
						double old_ave_plastic_strain = plastic_point.AvePlasticStrain;
						if (plastic_judgement > 0.0) {
							double magnitude_plastic_flow = plastic_judgement / (2.0 * shear_modulus + 2.0 / 3.0 * hardening_modulus);
							plastic_point.PlasticStrain += deviatoric_stress / deviatoric_stress_norm * magnitude_plastic_flow;
							plastic_point.AvePlasticStrain += 0.8164965809 * magnitude_plastic_flow;
						}
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							for (int i = 0; i < 6; i++) {
								double delt_strain = abs(old_plastic_strain[i] - plastic_point.PlasticStrain[i]);
								if (delt_strain > max_dplastic_strain)
									max_dplastic_strain = delt_strain;
							}
							double delt_ave_strain = abs(old_ave_plastic_strain - plastic_point.AvePlasticStrain);
							if (delt_ave_strain > max_dave_plastic_strain)
								max_dave_plastic_strain = delt_ave_strain;
						}
					}
		}

		void init() {
			phi_number = main_field::phi_number;
			phi_yield_stress.resize(phi_number, 0);
			phi_hardening_modulus.resize(phi_number, 0);
			phi_shear_modulus.resize(phi_number, 0);

			if (PhiProperties::instance().is_init()) {
				phi_property_number = PhiProperties::instance().phi_property_number();
				phase_yield_stress.resize(phi_property_number, 0);
				phase_hardening_modulus.resize(phi_property_number, 0);
				phase_shear_modulus.resize(phi_property_number, 0);
				WriteDebugFile("# Postprocess.SolidMechanics.Plasticity.yield_stress = ( phi_property_0_stress, ... ) \n");
				WriteDebugFile("#                                      .hardening_modulus = ( phi_property_0_modulus, ...) \n");
				WriteDebugFile("#                                      .shear_modulus = ( phi_property_0_modulus, ...) \n");
				std::string yield_stress_key = "Postprocess.SolidMechanics.Plasticity.yield_stress", yield_stress_input = "()";
				infile_reader::read_string_value(yield_stress_key, yield_stress_input, true);
				std::vector<input_value> yield_stress_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, yield_stress_key, yield_stress_input, true);
				if (yield_stress_value.size() != PhiProperties::instance().phi_property_number()) {
					WriteDebugFile("# ERROR : Phi Property Number is " + std::to_string(PhiProperties::instance().phi_property_number()) + ", size of yield_stress mismatch \n");
					SYS_PROGRAM_STOP;
				}
				std::string hardening_modulus_key = "Postprocess.SolidMechanics.Plasticity.hardening_modulus", hardening_modulus_input = "()";
				infile_reader::read_string_value(hardening_modulus_key, hardening_modulus_input, true);
				std::vector<input_value> hardening_modulus_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, hardening_modulus_key, hardening_modulus_input, true);
				if (hardening_modulus_value.size() != PhiProperties::instance().phi_property_number()) {
					WriteDebugFile("# ERROR : Phi Property Number is " + std::to_string(PhiProperties::instance().phi_property_number()) + ", size of hardening_modulus_value mismatch \n");
					SYS_PROGRAM_STOP;
				}
				std::string shear_modulus_key = "Postprocess.SolidMechanics.Plasticity.shear_modulus", shear_modulus_input = "()";
				infile_reader::read_string_value(shear_modulus_key, shear_modulus_input, true);
				std::vector<input_value> shear_modulus_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, shear_modulus_key, shear_modulus_input, true);
				if (shear_modulus_value.size() != PhiProperties::instance().phi_property_number()) {
					WriteDebugFile("# ERROR : Phi Property Number is " + std::to_string(PhiProperties::instance().phi_property_number()) + ", size of shear_modulus_value mismatch \n");
					SYS_PROGRAM_STOP;
				}
				for (size_t index = 0; index < phi_property_number; index++) {
					phase_yield_stress[index] = yield_stress_value[index].REAL_value;
					phase_hardening_modulus[index] = hardening_modulus_value[index].REAL_value;
					phase_shear_modulus[index] = shear_modulus_value[index].REAL_value;
				}
				for (size_t index = 0; index < phi_number; index++) {
					phi_yield_stress[index] = phase_yield_stress[PhiProperties::instance().phi_property(index)];
					phi_hardening_modulus[index] = phase_hardening_modulus[PhiProperties::instance().phi_property(index)];
					phi_shear_modulus[index] = phase_shear_modulus[PhiProperties::instance().phi_property(index)];
				}
			}
			else {
				WriteDebugFile("# Postprocess.SolidMechanics.Plasticity.yield_stress = ( phi_index_0_stress, ... ) \n");
				WriteDebugFile("#                                      .hardening_modulus = ( phi_index_0_modulus, ...) \n");
				WriteDebugFile("#                                      .shear_modulus = ( phi_index_0_modulus, ...) \n");
				std::string yield_stress_key = "Postprocess.SolidMechanics.Plasticity.yield_stress", yield_stress_input = "()";
				infile_reader::read_string_value(yield_stress_key, yield_stress_input, true);
				std::vector<input_value> yield_stress_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, yield_stress_key, yield_stress_input, true);
				if (yield_stress_value.size() != main_field::phi_number) {
					WriteDebugFile("# ERROR : Phi Number is " + std::to_string(main_field::phi_number) + ", size of yield_stress_value mismatch \n");
					SYS_PROGRAM_STOP;
				}
				std::string hardening_modulus_key = "Postprocess.SolidMechanics.Plasticity.hardening_modulus", hardening_modulus_input = "()";
				infile_reader::read_string_value(hardening_modulus_key, hardening_modulus_input, true);
				std::vector<input_value> hardening_modulus_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, hardening_modulus_key, hardening_modulus_input, true);
				if (hardening_modulus_value.size() != main_field::phi_number) {
					WriteDebugFile("# ERROR : Phi Number is " + std::to_string(main_field::phi_number) + ", size of hardening_modulus_value mismatch \n");
					SYS_PROGRAM_STOP;
				}
				std::string shear_modulus_key = "Postprocess.SolidMechanics.Plasticity.shear_modulus", shear_modulus_input = "()";
				infile_reader::read_string_value(shear_modulus_key, shear_modulus_input, true);
				std::vector<input_value> shear_modulus_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, shear_modulus_key, shear_modulus_input, true);
				if (shear_modulus_value.size() != main_field::phi_number) {
					WriteDebugFile("# ERROR : Phi Number is " + std::to_string(main_field::phi_number) + ", size of shear_modulus_value mismatch \n");
					SYS_PROGRAM_STOP;
				}
				for (size_t index = 0; index < phi_number; index++) {
					phi_yield_stress[index] = yield_stress_value[index].REAL_value;
					phi_hardening_modulus[index] = hardening_modulus_value[index].REAL_value;
					phi_shear_modulus[index] = shear_modulus_value[index].REAL_value;
				}
			}

			external_physical_field::init_plastic_field();
			stiffness_eigenstrain::eigenstrain_list.push_back(default_functions::cal_plasticstrain);
			WriteLog("> MODULE INIT : Plastic Explicit Solver On ! \n");
		}

		void deinit() {
			external_physical_field::plastic_field.clear();
		}
	}
}