#pragma once
#include "DDC_Functions.h"
namespace pf {
	namespace data_driven_complex_model {
		namespace phase_field_functions {
			bool is_interphase(size_t x, size_t y, size_t z, std::vector<REAL>& phi, size_t phi_index) {
				if (phi[phi_index] >= Phi_Num_Cut_Off && phi[phi_index] <= (1.0 - Phi_Num_Cut_Off))
					return true;
				if (phi[phi_index] < Phi_Num_Cut_Off) {
					if (parameters::diff_method == parameters::DifferenceMethod::FIVE_POINT) {
						if (parameters::phase_field->at(x + 1, y, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x - 1, y, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y + 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y - 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y, z + 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y, z - 1)[phi_index] >= Phi_Num_Cut_Off)
							return true;
					}
					else {
						if (parameters::phase_field->at(x + 1, y, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x - 1, y, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y + 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y - 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y, z + 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y, z - 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x + 1, y + 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x + 1, y - 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x + 1, y, z + 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x + 1, y, z - 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x - 1, y + 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x - 1, y - 1, z)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x - 1, y, z + 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x - 1, y, z - 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y + 1, z + 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y + 1, z - 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y - 1, z + 1)[phi_index] >= Phi_Num_Cut_Off
							|| parameters::phase_field->at(x, y - 1, z - 1)[phi_index] >= Phi_Num_Cut_Off)
							return true;
					}
				}
				else if (phi[phi_index] > Phi_Num_Cut_Off_R) {
					if (parameters::diff_method == parameters::DifferenceMethod::FIVE_POINT) {
						if (parameters::phase_field->at(x + 1, y, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x - 1, y, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y + 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y - 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y, z + 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y, z - 1)[phi_index] <= Phi_Num_Cut_Off_R)
							return true;
					}
					else {
						if (parameters::phase_field->at(x + 1, y, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x - 1, y, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y + 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y - 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y, z + 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y, z - 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x + 1, y + 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x + 1, y - 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x + 1, y, z + 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x + 1, y, z - 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x - 1, y + 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x - 1, y - 1, z)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x - 1, y, z + 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x - 1, y, z - 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y + 1, z + 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y + 1, z - 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y - 1, z + 1)[phi_index] <= Phi_Num_Cut_Off_R
							|| parameters::phase_field->at(x, y - 1, z - 1)[phi_index] <= Phi_Num_Cut_Off_R)
							return true;
					}
				}
				return false;
			}
			static void normalize_phi(FIELD_VAR& var) {
				REAL all_phi = 0;
				for (size_t index = 0; index < parameters::phi_number; index++)
					all_phi += var.new_phi[index];
				if (all_phi < SYS_EPSILON) {
					return;
				}
				else {
					for (size_t index = 0; index < parameters::phi_number; index++)
						var.new_phi[index] /= all_phi;
				}
			}
			// - 
			Vector3 normals(REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				Vector3 normal;
				normal[0] = alpha_grad[0] * beta_phi - alpha_phi * beta_grad[0];
				normal[1] = alpha_grad[1] * beta_phi - alpha_phi * beta_grad[1];
				normal[2] = alpha_grad[2] * beta_phi - alpha_phi * beta_grad[2];
				REAL length = sqrt(normal * normal);
				if (length < SYS_EPSILON) {
					normal[0] = 0.0;
					normal[1] = 0.0;
					normal[2] = 0.0;
				}
				else {
					normal[0] = normal[0] / length;
					normal[1] = normal[1] / length;
					normal[2] = normal[2] / length;
				};
				return normal;
			};
			// - interface mobility
			REAL Lij(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				return parameters::Lij[alpha_index][beta_index];
			};
			REAL Lij_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index];
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 - parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property]
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return Lij;
			};
			REAL Lij_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index];
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 + parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
				return Lij;
			};
			REAL Lij_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index];
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 - parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- parameters::intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- parameters::intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- parameters::intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
				return Lij;
			};
			REAL Lij_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index];
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 - parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- parameters::intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + parameters::intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
				return Lij;
			};
			REAL Lij_temp(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				return parameters::Lij[alpha_index][beta_index] * std::exp(-parameters::Qij[alpha_index][beta_index] / parameters::R / temperature);
			};
			REAL Lij_temp_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index] * std::exp(-parameters::Qij[alpha_index][beta_index] / parameters::R / temperature);
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 - parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property]
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return Lij;
			};
			REAL Lij_temp_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index] * std::exp(-parameters::Qij[alpha_index][beta_index] / parameters::R / temperature);
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 + parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
				return Lij;
			};
			REAL Lij_temp_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index] * std::exp(-parameters::Qij[alpha_index][beta_index] / parameters::R / temperature);
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 - parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- parameters::intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- parameters::intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- parameters::intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
				return Lij;
			};
			REAL Lij_temp_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij[alpha_index][beta_index] * std::exp(-parameters::Qij[alpha_index][beta_index] / parameters::R / temperature);
				// - anisotropic
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				Lij = Lij * REAL(1.0 - parameters::intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- parameters::intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + parameters::intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
				return Lij;
			};
			// - interface energy
			REAL xi_ab(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				return parameters::xi_ab[parameters::phi_property[alpha_index]][parameters::phi_property[beta_index]];
			};
			REAL xi_ab_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				REAL loc_xi_ab = parameters::xi_ab[phi_alpha_property][phi_beta_property];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + parameters::intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return loc_xi_ab;
			};
			REAL xi_ab_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				REAL loc_xi_ab = parameters::xi_ab[phi_alpha_property][phi_beta_property];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				loc_xi_ab = loc_xi_ab * REAL(1.0 - parameters::intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] +
					15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
				return loc_xi_ab;
			};
			REAL xi_ab_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				REAL loc_xi_ab = parameters::xi_ab[phi_alpha_property][phi_beta_property];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + parameters::intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					+ parameters::intEnAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0) - 30.0 * norm[2] * norm[2] + 3.0)
					+ parameters::intEnAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0) - 315.0 * std::pow(norm[2], 4.0)
						+ 105.0 * norm[2] * norm[2] - 5.0)
					+ parameters::intEnAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0) - std::pow(norm[1], 6.0)));
				return loc_xi_ab;
			};
			REAL xi_ab_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = parameters::phi_property[alpha_index], phi_beta_property = parameters::phi_property[beta_index];
				REAL loc_xi_ab = parameters::xi_ab[phi_alpha_property][phi_beta_property];
				Vector3 norm = parameters::grains_orientation.get_phi_rotationMatrix(beta_index) * (parameters::grains_orientation.get_phi_rotationMatrix(alpha_index)
					* normals(alpha_phi, beta_phi, alpha_grad, beta_grad));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + parameters::intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3.0 * norm[2] * norm[2] - 1.0), 2.0)
					+ parameters::intEnAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + parameters::intEnAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
				return loc_xi_ab;
			};
			REAL xi_abc(size_t alpha_index, size_t beta_index, size_t gamma_index) {
				return parameters::xi_abc[parameters::phi_property[alpha_index]][parameters::phi_property[beta_index]][parameters::phi_property[gamma_index]];
			};
			// - interface energy model
			REAL dfint_dphi_grad_S1996_acc(FIELD_VAR& point, size_t phi_index) {
				REAL grad = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (phi_bIndex != parameters::PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
						grad += 2 * parameters::interface_width 
							* _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index], point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex])
							* (point.grad_phi[phi_bIndex] * point.grad_phi[phi_bIndex] * point.old_phi[phi_index]
							- point.grad_phi[phi_index] * point.grad_phi[phi_bIndex] * point.old_phi[phi_bIndex] + point.old_phi[phi_index] * point.old_phi[phi_bIndex] * point.lap_phi[phi_bIndex]
							- point.old_phi[phi_bIndex] * point.old_phi[phi_bIndex] * point.lap_phi[phi_index]);
					}
				}
				return grad;
			};
			REAL dfint_dphi_grad_S1999_acc(FIELD_VAR& point, size_t phi_index) {
				REAL grad = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (phi_bIndex != parameters::PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
						grad += parameters::interface_width 
							* _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index], point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex]) 
							* point.lap_phi[phi_bIndex];
					}
				}
				return grad;
			};
			REAL dfint_dphi_pot_Nobstacle_acc(FIELD_VAR& point, size_t phi_index) {
				REAL pot = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (phi_bIndex != parameters::PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
						pot += 16 / parameters::interface_width / PI / PI * _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index], point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex]) * point.old_phi[phi_bIndex];
						for (size_t index2 = index + 1; index2 < parameters::PHI_ACC_NUMBER; index2++) {
							size_t phi_gIndex = point.active_index[index2];
							if (phi_gIndex != parameters::PAIRWISE_ACC_STOP && phi_gIndex != phi_index)
								pot += _xi_abc(phi_index, phi_bIndex, phi_gIndex) * point.old_phi[phi_bIndex] * point.old_phi[phi_gIndex] / parameters::interface_width;
						}
					}
				}
				return pot;
			};
			REAL dfint_dphi_pot_Nwell_acc(FIELD_VAR& point, size_t phi_index) {
				REAL pot = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (phi_bIndex != parameters::PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
						pot += 18 / parameters::interface_width * _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index], point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex]) * point.old_phi[phi_index] * point.old_phi[phi_bIndex] * point.old_phi[phi_bIndex];
						for (size_t index2 = index + 1; index2 < parameters::PHI_ACC_NUMBER; index2++) {
							size_t phi_gIndex = point.active_index[index2];
							if (phi_gIndex != parameters::PAIRWISE_ACC_STOP && phi_gIndex != phi_index)
								pot += 2 / parameters::interface_width * _xi_abc(phi_index, phi_bIndex, phi_gIndex) * point.old_phi[phi_index]
								* point.old_phi[phi_bIndex] * point.old_phi[phi_gIndex] * point.old_phi[phi_bIndex] * point.old_phi[phi_gIndex];
						}
					}
				}
				return pot;
			};
			// - 
			REAL dfint_dphi_pairwise_S2009(FIELD_VAR& point, size_t alpha_index, size_t beta_index) {
				return _xi_ab(alpha_index, beta_index, point.old_phi[alpha_index], point.old_phi[beta_index], point.grad_phi[alpha_index], point.grad_phi[beta_index]) 
					* (PI * PI / 2 / parameters::interface_width * (point.old_phi[alpha_index] - point.old_phi[beta_index])
					+ parameters::interface_width * (point.old_phi[beta_index] * point.lap_phi[alpha_index] - point.old_phi[alpha_index] * point.lap_phi[beta_index]));
			};
			REAL dfint_dphi_pairwise_S1996_NO(FIELD_VAR& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1996_acc(point, beta_index) + dfint_dphi_pot_Nobstacle_acc(point, beta_index) -
					dfint_dphi_grad_S1996_acc(point, alpha_index) - dfint_dphi_pot_Nobstacle_acc(point, alpha_index);
			};
			REAL dfint_dphi_pairwise_S1996_NW(FIELD_VAR& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1996_acc(point, beta_index) + dfint_dphi_pot_Nwell_acc(point, beta_index) -
					dfint_dphi_grad_S1996_acc(point, alpha_index) - dfint_dphi_pot_Nwell_acc(point, alpha_index);
			};
			REAL dfint_dphi_pairwise_S1999_NO(FIELD_VAR& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1999_acc(point, beta_index) + dfint_dphi_pot_Nobstacle_acc(point, beta_index) -
					dfint_dphi_grad_S1999_acc(point, alpha_index) - dfint_dphi_pot_Nobstacle_acc(point, alpha_index);
			};
			REAL dfint_dphi_pairwise_S1999_NW(FIELD_VAR& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1999_acc(point, beta_index) + dfint_dphi_pot_Nwell_acc(point, beta_index) -
					dfint_dphi_grad_S1999_acc(point, alpha_index) - dfint_dphi_pot_Nwell_acc(point, alpha_index);
			};
			REAL dfbulk_dphi_pairwise_acc(long long x, long long y, long long z, size_t phi_index) {
				REAL dfbulk_sum = 0.0;
				for (auto dfbulk = parameters::delt_Fbulk_delt_phi.begin(); dfbulk < parameters::delt_Fbulk_delt_phi.end(); dfbulk++)
					dfbulk_sum += (*dfbulk)(x, y, z, phi_index);
				return dfbulk_sum;
			}
			REAL source_pairwise_acc(long long x, long long y, long long z, size_t alpha_index, size_t beta_index) {
				REAL source_sum = 0.0;
				for (auto source = parameters::source_alpha_beta.begin(); source < parameters::source_alpha_beta.end(); source++)
					source_sum += (*source)(x, y, z, alpha_index, beta_index);
				return source_sum;
			}
			void pairwise_normalize_acc(std::vector<size_t>& active_index, std::vector<bool>& interphase, 
				std::vector<REAL>& old_phi, std::vector<REAL>& phi_increment) {
				REAL scale = 1.0, increment = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
					size_t acc_index = active_index[index];
					if (acc_index != parameters::PAIRWISE_ACC_STOP && interphase[acc_index]) {
						increment = *parameters::delt_t * phi_increment[acc_index];
						if (isTwoREALEquality(increment, 0))
							continue;
						else if ((old_phi[acc_index] + increment) > 1) {
							REAL p_scale = (1 - old_phi[acc_index]) / increment;
							if (p_scale < scale)
								scale = p_scale;
						}
						else if ((old_phi[acc_index] + increment) < 0) {
							REAL p_scale = abs(old_phi[acc_index] / increment);
							if (p_scale < scale)
								scale = p_scale;
						}
					}
				}
				for (size_t index = 0; index < parameters::PAIRWISE_ACC_STOP; index++) {
					size_t acc_index = active_index[index];
					if (acc_index != parameters::PAIRWISE_ACC_STOP && interphase[acc_index])
						phi_increment[acc_index] *= scale;
				}
			}
			//=========================================================================================================================================
			void init_phi_pair_wise() {
#pragma omp parallel for collapse(3)
				for (long long x = 0; x < parameters::phase_field->Nx(); x++)
					for (long long y = 0; y < parameters::phase_field->Ny(); y++)
						for (long long z = 0; z < parameters::phase_field->Nz(); z++) {
							std::vector<REAL>& phi = parameters::phase_field->at(x, y, z);
							FIELD_VAR& field_var = parameters::field_variable(x, y, z);
							if (parameters::is_phi_normalized) {
								REAL sum_phi = 0.0;
								for (size_t index = 0; index < parameters::phi_number; index++)
									sum_phi += phi[index];
								if (sum_phi > SYS_EPSILON) {
									for (size_t index = 0; index < parameters::phi_number; index++)
										phi[index] /= sum_phi;
								}
							}
							for (int index = 0; index < parameters::phi_number; index++) {
								field_var.new_phi[index] = phi[index];
								field_var.interphase[index] = is_interphase(x, y, z, field_var.new_phi, index);
							}
						}
			}
			void pre_calculation_phi_pair_wise() {
#pragma omp parallel for collapse(3)
				for (long long x = parameters::phase_field->COMP_X_BGN(); x <= parameters::phase_field->COMP_X_END(); x++)
					for (long long y = parameters::phase_field->COMP_Y_BGN(); y <= parameters::phase_field->COMP_Y_END(); y++)
						for (long long z = parameters::phase_field->COMP_Z_BGN(); z <= parameters::phase_field->COMP_Z_END(); z++) {
							FIELD_VAR& field_var = parameters::field_variable(x, y, z);
							for (size_t index = 0; index < parameters::phi_number; index++) {
								field_var.old_phi[index] = field_var.new_phi[index]; // save phi 
								field_var.new_phi[index] = 0;
							}
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++)
								field_var.active_index[index] = parameters::PAIRWISE_ACC_STOP;
							size_t actIndex = 0;
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
								if (field_var.interphase[index] || field_var.old_phi[index] > Phi_Num_Cut_Off) {
									field_var.active_index[actIndex] = index; // save the active phi index
									actIndex++;
								}
								if (field_var.interphase[index]) {
									field_var.grad_phi[index][0] = (parameters::field_variable(x + 1, y, z).old_phi[index] - parameters::field_variable(x - 1, y, z).old_phi[index]) / 2 / parameters::delt_r;
									field_var.grad_phi[index][1] = (parameters::field_variable(x, y + 1, z).old_phi[index] - parameters::field_variable(x, y - 1, z).old_phi[index]) / 2 / parameters::delt_r;
									field_var.grad_phi[index][2] = (parameters::field_variable(x, y, z + 1).old_phi[index] - parameters::field_variable(x, y, z - 1).old_phi[index]) / 2 / parameters::delt_r;
									if (parameters::diff_method == parameters::DifferenceMethod::FIVE_POINT) {
										field_var.lap_phi[index] =
											(parameters::field_variable(x + 1, y, z).old_phi[index] + parameters::field_variable(x - 1, y, z).old_phi[index]
												+ parameters::field_variable(x, y + 1, z).old_phi[index] + parameters::field_variable(x, y - 1, z).old_phi[index]
												+ parameters::field_variable(x, y, z + 1).old_phi[index] + parameters::field_variable(x, y, z - 1).old_phi[index]
												- 6 * field_var.old_phi[index]) / parameters::delt_r / parameters::delt_r;
									}
									else if (parameters::diff_method == parameters::DifferenceMethod::NINE_POINT) {
										field_var.lap_phi[index] =
											((parameters::field_variable(x + 1, y, z).old_phi[index] + parameters::field_variable(x - 1, y, z).old_phi[index]
												+ parameters::field_variable(x, y + 1, z).old_phi[index] + parameters::field_variable(x, y - 1, z).old_phi[index]
												+ parameters::field_variable(x, y, z + 1).old_phi[index] + parameters::field_variable(x, y, z - 1).old_phi[index]) * 4
												+ parameters::field_variable(x - 1, y - 1, z).old_phi[index] + parameters::field_variable(x - 1, y + 1, z).old_phi[index]
												+ parameters::field_variable(x + 1, y - 1, z).old_phi[index] + parameters::field_variable(x + 1, y + 1, z).old_phi[index]
												+ parameters::field_variable(x - 1, y, z - 1).old_phi[index] + parameters::field_variable(x - 1, y, z + 1).old_phi[index]
												+ parameters::field_variable(x + 1, y, z - 1).old_phi[index] + parameters::field_variable(x + 1, y, z + 1).old_phi[index]
												+ parameters::field_variable(x, y - 1, z - 1).old_phi[index] + parameters::field_variable(x, y - 1, z + 1).old_phi[index]
												+ parameters::field_variable(x, y + 1, z - 1).old_phi[index] + parameters::field_variable(x, y + 1, z + 1).old_phi[index]
												- 36 * field_var.old_phi[index]) / 6 / parameters::delt_r / parameters::delt_r;
									}
								}
								else {
									field_var.grad_phi[index][0] = 0.0;
									field_var.grad_phi[index][1] = 0.0;
									field_var.grad_phi[index][2] = 0.0;
									field_var.lap_phi[index] = 0.0;
								}
							}
						}
#pragma omp parallel for collapse(3)
				for (long long x = parameters::phase_field->COMP_X_BGN(); x <= parameters::phase_field->COMP_X_END(); x++)
					for (long long y = parameters::phase_field->COMP_Y_BGN(); y <= parameters::phase_field->COMP_Y_END(); y++)
						for (long long z = parameters::phase_field->COMP_Z_BGN(); z <= parameters::phase_field->COMP_Z_END(); z++) {
							FIELD_VAR& field_var = parameters::field_variable(x, y, z);
							std::vector<REAL>& phi_increment = field_var.new_phi;
							REAL N = 0;
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
								size_t phi_index = field_var.active_index[index];
								if (phi_index != parameters::PAIRWISE_ACC_STOP && field_var.interphase[phi_index])
									N = N + 1;
								else
									continue;
							}
							// interface energy
							for (size_t aIndex = 0; aIndex < parameters::PHI_ACC_NUMBER; aIndex++) {
								size_t alpha_index = field_var.active_index[aIndex];
								if (alpha_index == parameters::PAIRWISE_ACC_STOP || !field_var.interphase[alpha_index])
									continue;
								for (size_t bIndex = aIndex + 1; bIndex < parameters::PHI_ACC_NUMBER; bIndex++) {
									size_t beta_index = field_var.active_index[bIndex];
									if (beta_index == parameters::PAIRWISE_ACC_STOP || !field_var.interphase[beta_index])
										continue;
									REAL int_incre_b_a = _Lij(alpha_index, beta_index, field_var.old_phi[alpha_index], field_var.old_phi[beta_index], 
										field_var.grad_phi[alpha_index], field_var.grad_phi[beta_index], field_var.old_temp) / parameters::interface_width / N
										* dfint_dphi_pairwise_acc(field_var, alpha_index, beta_index);
									if (parameters::is_phi_normalized) {
										if ((int_incre_b_a > SYS_EPSILON && (field_var.old_phi[alpha_index] > SYS_EPSILON_R || field_var.old_phi[beta_index] < SYS_EPSILON))
											|| (int_incre_b_a < -SYS_EPSILON && (field_var.old_phi[alpha_index] < SYS_EPSILON || field_var.old_phi[beta_index] > SYS_EPSILON_R)))
											int_incre_b_a = 0.0;
									}
									phi_increment[alpha_index] += int_incre_b_a;
									phi_increment[beta_index] -= int_incre_b_a;
#ifdef _DEBUG
									if (_isnan(int_incre_b_a)) {
										std::cout << "DEBUG: interface energy (pair-wise functions) error !" << std::endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
							// driving force and source term
							for (size_t aIndex = 0; aIndex < parameters::PHI_ACC_NUMBER; aIndex++) {
								size_t alpha_index = field_var.active_index[aIndex];
								if (alpha_index == parameters::PAIRWISE_ACC_STOP || !field_var.interphase[alpha_index])
									continue;
								for (size_t bIndex = aIndex + 1; bIndex < parameters::PHI_ACC_NUMBER; bIndex++) {
									size_t beta_index = field_var.active_index[bIndex];
									if (beta_index == parameters::PAIRWISE_ACC_STOP || !field_var.interphase[beta_index])
										continue;
									if (field_var.old_phi[alpha_index] > PhiCon_Num_Cut_Off && field_var.old_phi[beta_index] > PhiCon_Num_Cut_Off) {
										REAL bulk_increment_b_a = _Lij(alpha_index, beta_index, field_var.old_phi[alpha_index], field_var.old_phi[beta_index],
											field_var.grad_phi[alpha_index], field_var.grad_phi[beta_index], field_var.old_temp) / parameters::interface_width / N
											* (dfbulk_dphi_pairwise_acc(x, y, z, beta_index) - dfbulk_dphi_pairwise_acc(x, y, z, alpha_index))
											+ source_pairwise_acc(x, y, z, alpha_index, beta_index);
										if (parameters::is_phi_normalized) {
											if ((bulk_increment_b_a > SYS_EPSILON && (field_var.old_phi[alpha_index] > SYS_EPSILON_R || field_var.old_phi[beta_index] < SYS_EPSILON))
												|| (bulk_increment_b_a < -SYS_EPSILON && (field_var.old_phi[alpha_index] < SYS_EPSILON || field_var.old_phi[beta_index] > SYS_EPSILON_R)))
												bulk_increment_b_a = 0.0;
										}
										phi_increment[alpha_index] += bulk_increment_b_a;
										phi_increment[beta_index] -= bulk_increment_b_a;
									}
								}
							}
							// numerical treatment to nomalize the phi
							if (parameters::is_phi_normalized)
								pairwise_normalize_acc(field_var.active_index, field_var.interphase, field_var.old_phi, phi_increment);
						}
			}
			REAL solve_phi_pair_wise() {
				REAL MAX_PHI_INCREMENT = 0.0;
				bool phi_change = false;
#pragma omp parallel for collapse(3)
				for (long long x = parameters::phase_field->COMP_X_BGN(); x <= parameters::phase_field->COMP_X_END(); x++)
					for (long long y = parameters::phase_field->COMP_Y_BGN(); y <= parameters::phase_field->COMP_Y_END(); y++)
						for (long long z = parameters::phase_field->COMP_Z_BGN(); z <= parameters::phase_field->COMP_Z_END(); z++) {
							FIELD_VAR& field_var = parameters::field_variable(x, y, z);
							std::vector<REAL>& phi_increment = field_var.new_phi;
							phi_change = false;
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++) {
								size_t phi_index = field_var.active_index[index];
								if (phi_index != parameters::PAIRWISE_ACC_STOP && field_var.interphase[phi_index]) {
									phi_change = true;
									field_var.new_phi[phi_index] = field_var.old_phi[phi_index] + (*parameters::delt_t) * phi_increment[phi_index];
#ifdef _OPENMP
#pragma omp critical
#endif
									{
										if (abs(field_var.old_phi[phi_index] - field_var.new_phi[phi_index]) > MAX_PHI_INCREMENT)
											MAX_PHI_INCREMENT = abs(field_var.old_phi[phi_index] - field_var.new_phi[phi_index]);
									}
#ifdef _DEBUG
									if (_isnan(field_var.new_phi[phi_index])) {
										std::string error_report = "ERROR : Phase fraction is NaN in x = " + std::to_string(x) + " , y = " + std::to_string(y)
											+ " , z = " + std::to_string(z) + " position, the phase index is " + std::to_string(phi_index) + "\n";
										std::cout << error_report << std::endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
							if (phi_change) {
								// normalize the phi
								if (parameters::is_phi_normalized)
									normalize_phi(field_var);
								for (size_t index = 0; index < parameters::phi_number; index++) {
									// modify phi
									parameters::phase_field->at(x, y, z)[index] = field_var.new_phi[index];
									// change _flag
									if (field_var.new_phi[index] >= Phi_Num_Cut_Off && field_var.new_phi[index] <= Phi_Num_Cut_Off_R) {
										if (!field_var.interphase[index]) {
											field_var.interphase[index] = true;
											parameters::field_variable(x + 1, y, z).interphase[index] = true;
											parameters::field_variable(x - 1, y, z).interphase[index] = true;
											parameters::field_variable(x, y + 1, z).interphase[index] = true;
											parameters::field_variable(x, y - 1, z).interphase[index] = true;
											parameters::field_variable(x, y, z + 1).interphase[index] = true;
											parameters::field_variable(x, y, z - 1).interphase[index] = true;
											if (parameters::diff_method == parameters::DifferenceMethod::NINE_POINT) {
												parameters::field_variable(x + 1, y + 1, z).interphase[index] = true;
												parameters::field_variable(x + 1, y - 1, z).interphase[index] = true;
												parameters::field_variable(x + 1, y, z + 1).interphase[index] = true;
												parameters::field_variable(x + 1, y, z - 1).interphase[index] = true;
												parameters::field_variable(x - 1, y + 1, z).interphase[index] = true;
												parameters::field_variable(x - 1, y - 1, z).interphase[index] = true;
												parameters::field_variable(x - 1, y, z + 1).interphase[index] = true;
												parameters::field_variable(x - 1, y, z - 1).interphase[index] = true;
												parameters::field_variable(x, y + 1, z + 1).interphase[index] = true;
												parameters::field_variable(x, y + 1, z - 1).interphase[index] = true;
												parameters::field_variable(x, y - 1, z + 1).interphase[index] = true;
												parameters::field_variable(x, y - 1, z - 1).interphase[index] = true;
											}
										}
									}
									else if (field_var.new_phi[index] < Phi_Num_Cut_Off) {
										if (parameters::diff_method == parameters::DifferenceMethod::FIVE_POINT) {
											if (parameters::field_variable(x + 1, y, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x - 1, y, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y + 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y - 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y, z + 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y, z - 1).new_phi[index] >= Phi_Num_Cut_Off)
												field_var.interphase[index] = true;
											else
												field_var.interphase[index] = false;
										}
										else {
											if (parameters::field_variable(x + 1, y, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x - 1, y, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y + 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y - 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y, z + 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y, z - 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x + 1, y + 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x + 1, y - 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x + 1, y, z + 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x + 1, y, z - 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x - 1, y + 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x - 1, y - 1, z).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x - 1, y, z + 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x - 1, y, z - 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y + 1, z + 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y + 1, z - 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y - 1, z + 1).new_phi[index] >= Phi_Num_Cut_Off
												|| parameters::field_variable(x, y - 1, z - 1).new_phi[index] >= Phi_Num_Cut_Off)
												field_var.interphase[index] = true;
											else
												field_var.interphase[index] = false;
										}
									}
									else if (field_var.new_phi[index] > Phi_Num_Cut_Off_R) {
										if (parameters::diff_method == parameters::DifferenceMethod::FIVE_POINT) {
											if (parameters::field_variable(x + 1, y, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x - 1, y, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y + 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y - 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y, z + 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y, z - 1).new_phi[index] <= Phi_Num_Cut_Off_R)
												field_var.interphase[index] = true;
											else
												field_var.interphase[index] = false;
										}
										else {
											if (parameters::field_variable(x + 1, y, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x - 1, y, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y + 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y - 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y, z + 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y, z - 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x + 1, y + 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x + 1, y - 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x + 1, y, z + 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x + 1, y, z - 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x - 1, y + 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x - 1, y - 1, z).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x - 1, y, z + 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x - 1, y, z - 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y + 1, z + 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y + 1, z - 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y - 1, z + 1).new_phi[index] <= Phi_Num_Cut_Off_R
												|| parameters::field_variable(x, y - 1, z - 1).new_phi[index] <= Phi_Num_Cut_Off_R)
												field_var.interphase[index] = true;
											else
												field_var.interphase[index] = false;
										}
									}
								}
							}
						}
				return MAX_PHI_INCREMENT;
			}
		}
		void exec_pre_iii() {
			void init_phi_pair_wise();
		}
		void exec_i() {
			void pre_calculation_phi_pair_wise();
			REAL solve_phi_pair_wise();
		}
		void deinit() {

		}
	}
}