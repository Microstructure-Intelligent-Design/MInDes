#include "StiffnessEigenStrain.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../model_modules/GrainsOrientations.h"
#include "../../model_modules/PhiProperties.h"
namespace pf {
	namespace stiffness_eigenstrain {
		void do_phi_index_orientation(size_t phi_index, size_t phi_property) {
			phi_index_stiffness[phi_index] = phi_property_stiffness[phi_property].
				get_rotated_matrix(GrainsOrientations::instance().get_phi_rotationMatrix(phi_index));
		}
		namespace default_functions {
			inline void cal_eigenstrain_phi_dependent_norm(long long x, long long y, long long z, vStrain& eigenstrain) {
				std::vector<REAL>& phi = phase_field->at(x, y, z);
				for (size_t index = 0; index < phi_number; index++) {
					if (phi[index] < SYS_EPSILON)
						continue;
					eigenstrain += phi_index_eigen_strain[index] * phi[index];
				}
			}
			inline void cal_stiffness_phi_dependent_norm(long long x, long long y, long long z, Matrix6x6& stiffness) {
				std::vector<REAL>& phi = phase_field->at(x, y, z);
				for (size_t index = 0; index < phi_number; index++) {
					if (phi[index] < SYS_EPSILON)
						continue;
					stiffness += phi_index_stiffness[index] * phi[index];
				}
			}
			inline void cal_stiffness_phi_dependent_with_crack_norm(long long x, long long y, long long z, Matrix6x6& stiffness) {
				std::vector<REAL>& phi = phase_field->at(x, y, z);
				REAL sum_phi = 0, crack_phi = 0;
				for (size_t index = 0; index < phi_number; index++)
					if (PhiProperties::instance().phi_property(index) != crack_phi_property) {
						if (phi[index] > SYS_EPSILON) {
							sum_phi += phi[index];
							stiffness += phi_index_stiffness[index] * phi[index];
						}
					}
					else {
						crack_phi = phi[index];
					}
				if (sum_phi > SYS_EPSILON)
					stiffness = stiffness / sum_phi * (1 - crack_phi) * (1 - crack_phi);
			}
			inline void cal_eigenstrain_normal(long long x, long long y, long long z, vStrain& eigenstrain) {
				for (auto func = eigenstrain_list.begin(); func < eigenstrain_list.end(); func++)
					(*func)(x, y, z, eigenstrain);
			}
		}
		void init(size_t _phi_number, size_t _phi_property_number) {
			phi_number = _phi_number;
			phi_property_number = _phi_property_number;
			Matrix6x6 Cij; Cij.set_to_zero();
			phi_property_stiffness.resize(phi_property_number, Cij);
			phi_index_stiffness.resize(phi_number, Cij);
			vStrain Eij; Eij.set_to_zero();
			phi_property_eigen_strain.resize(phi_property_number, Eij);
			phi_index_eigen_strain.resize(phi_number, Eij);

			stiffness = default_functions::cal_stiffness_phi_dependent_norm;
			eigen_strain = default_functions::cal_eigenstrain_normal;
			eigenstrain_list.push_back(default_functions::cal_eigenstrain_phi_dependent_norm);
			std::string solid_phi_key = "Postprocess.SolidMechanics.solid_phases", solid_phi_input = "()";
			infile_reader::read_string_value(solid_phi_key, solid_phi_input, true);
			std::vector<input_value> solid_phi_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, solid_phi_key, solid_phi_input, true);

			bool is_es_phi_dependent = false;
			for (auto solid_phi = solid_phi_value.begin(); solid_phi < solid_phi_value.end(); solid_phi++)
				if (PhiProperties::instance().is_phi_property(solid_phi->string_value)) {
					WriteDebugFile("# Postprocess.SolidMechanics.Stiffness.PhaseName = [(C11,C12,C13,C14,C15,C16),\\\\ \n");
					WriteDebugFile("                                                    (C21,C22,C23,C24,C25,C26),\\\\ \n");
					WriteDebugFile("                                                    (C31,C32,C33,C34,C35,C36),\\\\ \n");
					WriteDebugFile("                                                    (C41,C42,C43,C44,C45,C46),\\\\ \n");
					WriteDebugFile("                                                    (C51,C52,C53,C54,C55,C56),\\\\ \n");
					WriteDebugFile("                                                    (C61,C62,C63,C64,C65,C66)] \n");
					string stiffness_key = "Postprocess.SolidMechanics.Stiffness." + solid_phi->string_value, stiffness_input = "[(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0)]";
					infile_reader::read_string_value(stiffness_key, stiffness_input, true);
					vector<vector<input_value>> stiffness_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, stiffness_key, stiffness_input, true);
					Matrix6x6 Cij;
					for (size_t xx = 0; xx < 6; xx++)
						for (size_t yy = 0; yy < 6; yy++)
							Cij(xx, yy) = stiffness_value[xx][yy].REAL_value;
					phi_property_stiffness[PhiProperties::instance().phi_property(solid_phi->string_value)] = Cij;
					WriteDebugFile("# Postprocess.SolidMechanics.EigenStrain.PhaseName = [(E11,E22,E33,2E23,2E13,2E12)] \n");
					string eigenstrain_key = "Postprocess.SolidMechanics.EigenStrain." + solid_phi->string_value, eigenstrain_input = "(0,0,0,0,0,0)";
					if (infile_reader::read_string_value(eigenstrain_key, eigenstrain_input, true))
						is_es_phi_dependent = true;
					vector<input_value> eigenstrain_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, eigenstrain_key, eigenstrain_input, true);
					Vector6 Eij;
					for (size_t xx = 0; xx < 6; xx++)
						Eij[xx] = eigenstrain_value[xx].REAL_value;
					phi_property_eigen_strain[PhiProperties::instance().phi_property(solid_phi->string_value)] = Eij;
				}
			for (size_t index = 0; index < phi_number; index++)
				do_phi_index_orientation(index, PhiProperties::instance().phi_property(index));
			WriteLog("> MODULE INIT : Stiffness & EigenStrain \n");
		}
	}
}