/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include "../Base.h"

namespace pf {
	enum Int_Gradient { Steinbach_1996, Steinbach_1999, Steinbach_G2009, Int_GStandard };
	enum Int_Potential { Nestler_Well, Nestler_Obstacle, Steinbach_P2009 };
	enum Int_AnisoModel { No_Aniso, Aniso_Cos_E1, Aniso_Cos_E2 };

	namespace interface_energy_funcs {
		static double xi_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return 0.0;
		};
		static double xi_abc(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma) {
			return 0.0;
		};
		static double dfint_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, double int_width) {
			return 0.0;
		}
		static Vector3 normals(pf::PhaseEntry& phase) {
			Vector3 normal;
			normal = phase.phi_grad;
			double length = sqrt(normal * normal);
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
		static Vector3 normals(pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			Vector3 normal;
			normal[0] = alpha.phi_grad[0] * beta.phi - alpha.phi * beta.phi_grad[0];
			normal[1] = alpha.phi_grad[1] * beta.phi - alpha.phi * beta.phi_grad[1];
			normal[2] = alpha.phi_grad[2] * beta.phi - alpha.phi * beta.phi_grad[2];
			double length = sqrt(normal * normal);
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
		static double abs_grad_phi_S1996(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			Vector3 grad;
			grad[0] = alpha.phi_grad[0] * beta.phi - alpha.phi * beta.phi_grad[0];
			grad[1] = alpha.phi_grad[1] * beta.phi - alpha.phi * beta.phi_grad[1];
			grad[2] = alpha.phi_grad[2] * beta.phi - alpha.phi * beta.phi_grad[2];
			return grad.abs();
		}
		static double abs_grad_phi_S1999(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return fabs(alpha.phi_grad * beta.phi_grad);
		}
		static double abs_grad_phi_standard(pf::PhaseNode& node, pf::PhaseEntry& alpha) {
			return alpha.phi_grad.abs();
		}
	}
	namespace interface_energy {
		// models
		static int interface_gradient = Int_Gradient::Steinbach_G2009;
		static int interface_potential = Int_Potential::Steinbach_P2009;
		static double(*_xi_a)(pf::PhaseNode& node, pf::PhaseEntry& alpha);
		static double(*_xi_ab)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta);
		static double(*_xi_abc)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma);
		static double (*_abs_grad_phi_standard)(pf::PhaseNode& node, pf::PhaseEntry& alpha);
		static double (*_abs_grad_phi_pairwise)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta);
		// parameters
		static double interface_width = 4.0;
		static double const_xi_a = 0.0;
		static double const_xi_ab = 0.0;
		static double const_xi_abc = 0.0;
		static double_box matirx_xi_a;
		static tensor2_double matirx_xi_ab;
		static tensor3_double matirx_xi_abc;
		// anisotropy
		static bool is_anisotropy_on = false;
		static int anisotropy_model = Int_AnisoModel::No_Aniso;
		static double aniso_strength = 0.0;
		static double aniso_module_num = 1.0;
		static double aniso_angle = 0.0;
		//static;
		// output
		int_box interface_energy_standard_output;
		PairFlag interface_energy_pairwise_output;

		static Int_Gradient get_interface_gradient() {
			return Int_Gradient(interface_gradient);
		}
		static Int_Potential get_interface_potential() {
			return Int_Potential(interface_potential);
		}
		static double get_interface_width() {
			return interface_width;
		}
		static double xi_a_const(pf::PhaseNode& node, pf::PhaseEntry& alpha) {
			return const_xi_a;
		};
		static double xi_ab_const(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return const_xi_ab;
		}
		static double xi_ab_const_aniso_cos(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return const_xi_ab;
		};
		static double xi_abc_const(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma) {
			return const_xi_abc;
		};

		static double xi_a_matrix(pf::PhaseNode& node, pf::PhaseEntry& alpha) {
			// to be defined
			double _xi = 0.0;
			for (auto xi_a = matirx_xi_a.begin(); xi_a < matirx_xi_a.end(); xi_a++)
				if (xi_a->index == alpha.index)
					_xi = xi_a->value;
			return _xi;
		};

		static double xi_ab_matrix(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			// to be defined
			double _xi = 0.0;
			for (auto xi_a = matirx_xi_ab.begin(); xi_a < matirx_xi_ab.end(); xi_a++)
				for (auto xi_b = xi_a->begin(); xi_b < xi_a->end(); xi_b++)
					if ((xi_a->index == alpha.index && xi_b->index == beta.index)
						|| (xi_a->index == beta.index && xi_b->index == alpha.index))
						_xi = xi_b->val;
			return _xi;
		};

		static double xi_abc_matrix(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma) {
			// to be defined
			double _xi = 0.0;
			for (auto xi_a = matirx_xi_abc.begin(); xi_a < matirx_xi_abc.end(); xi_a++)
				for (auto xi_b = xi_a->begin(); xi_b < xi_a->end(); xi_b++)
					for (auto xi_c = xi_b->begin(); xi_c < xi_b->end(); xi_c++)
						if ((xi_a->index == alpha.index && xi_b->index == beta.index && xi_c->index == gamma.index)
							|| (xi_a->index == beta.index && xi_b->index == alpha.index && xi_c->index == gamma.index)
							|| (xi_a->index == alpha.index && xi_b->index == gamma.index && xi_c->index == beta.index)
							|| (xi_a->index == gamma.index && xi_b->index == beta.index && xi_c->index == alpha.index))
							_xi = xi_c->val;
			return _xi;
		};

		// store epsilon and epsilon*Depsilon_Dtheta inside nodes
		static void _xi_a_aniso_cos_e1(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double theta = std::atan2(phase.phi_grad.getY(), phase.phi_grad.getX());
			double angle = aniso_module_num * (theta - aniso_angle);
			double sigma = 1 + aniso_strength * std::cos(angle);
			double xi = _xi_a(node, phase);
			double epsilon = xi * sigma;
			double De_Dtheta = -1.0 * xi * aniso_strength * aniso_module_num * std::sin(angle);
			node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff, epsilon);
			node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff_Derived, De_Dtheta);
		}

		// store epsilon and epsilon*Depsilon_Dtheta inside nodes
		static void _xi_a_aniso_cos_e2(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double theta = std::atan2(phase.phi_grad.getY(), phase.phi_grad.getX());
			double angle = aniso_module_num * (theta - aniso_angle);
			double sigma = 1 + aniso_strength * std::cos(angle);
			double xi = _xi_a(node, phase);
			double epsilon = xi * sigma;
			double e_de = -1.0 * epsilon * xi * aniso_strength * aniso_module_num * std::sin(angle);
			node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff, epsilon);
			node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff_Derived, e_de);
		}

		// -----------------------------------------------------------
		static double dfint_dphi_S2009(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double intface_width) {
			return _xi_ab(node, alpha, beta) * (intface_width * (beta.phi * alpha.laplacian - alpha.phi * beta.laplacian)
				+ PI * PI / 2.0 / intface_width * (alpha.phi - beta.phi));
		};
		static double dfint_dphi_grad_S1996_acc(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double grad = 0.0;
			for (auto beta = node.Goast_Phase.begin(); beta < node.Goast_Phase.end(); beta++)
				if ((*beta)->index != phase.index) {
					grad += 2.0 * intface_width * _xi_ab(node, phase, **beta) * ((*beta)->phi_grad * (*beta)->phi_grad * phase.phi
						- phase.phi_grad * (*beta)->phi_grad * (*beta)->phi + phase.phi * (*beta)->phi * (*beta)->laplacian - (*beta)->phi * (*beta)->phi * phase.laplacian);
				}
			return grad;
		};
		static double dfint_dphi_grad_S1996_norm(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double grad = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++)
				if (beta->index != phase.index && beta->_flag) {
					grad += 2.0 * intface_width * _xi_ab(node, phase, *beta) * (beta->phi_grad * beta->phi_grad * phase.phi
						- phase.phi_grad * beta->phi_grad * beta->phi + phase.phi * beta->phi * beta->laplacian - beta->phi * beta->phi * phase.laplacian);
				}
			return grad;
		};
		static double dfint_dphi_grad_S1999_acc(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double grad = 0.0;
			for (auto beta = node.Goast_Phase.begin(); beta < node.Goast_Phase.end(); beta++)
				if ((*beta)->index != phase.index) {
					grad += intface_width * _xi_ab(node, phase, **beta) * (*beta)->laplacian;
				}
			return grad;
		};
		static double dfint_dphi_grad_S1999_norm(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double grad = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++)
				if (beta->index != phase.index && beta->_flag) {
					grad += intface_width * _xi_ab(node, phase, *beta) * beta->laplacian;
				}
			return grad;
		};
		static double dfint_dphi_pot_Nobstacle_acc(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double pot = 0.0;
			for (auto beta = node.Goast_Phase.begin(); beta < node.Goast_Phase.end(); beta++) {
				if ((*beta)->index == phase.index)
					continue;
				pot += 16.0 / intface_width / PI / PI * _xi_ab(node, phase, **beta) * (*beta)->phi;
				for (auto gamma = beta + 1; gamma < node.Goast_Phase.end(); gamma++) {
					if ((*gamma)->index == phase.index)
						continue;
					pot += _xi_abc(node, phase, **beta, **gamma) * (*beta)->phi * (*gamma)->phi / intface_width;
				}
			}
			return pot;
		};
		static double dfint_dphi_pot_Nobstacle_norm(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double pot = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++) {
				if (beta->index == phase.index || beta->_flag == pf_BULK)
					continue;
				pot += 16.0 / intface_width / PI / PI * _xi_ab(node, phase, *beta) * beta->phi;
				for (auto gamma = beta + 1; gamma < node.end(); gamma++) {
					if (gamma->index == phase.index || gamma->_flag == pf_BULK)
						continue;
					pot += _xi_abc(node, phase, *beta, *gamma) * beta->phi * gamma->phi / intface_width;
				}
			}
			return pot;
		};
		static double dfint_dphi_pot_Nwell_acc(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double pot = 0.0;
			for (auto beta = node.Goast_Phase.begin(); beta < node.Goast_Phase.end(); beta++) {
				if ((*beta)->index == phase.index)
					continue;
				pot += 18.0 / intface_width * _xi_ab(node, phase, **beta) * phase.phi * (*beta)->phi * (*beta)->phi;
				for (auto gamma = beta + 1; gamma < node.Goast_Phase.end(); gamma++) {
					if ((*gamma)->index == phase.index)
						continue;
					pot += 2.0 / intface_width * _xi_abc(node, phase, **beta, **gamma) * phase.phi
						* (*beta)->phi * (*gamma)->phi * (*beta)->phi * (*gamma)->phi;
				}
			}
			return pot;
		};
		static double dfint_dphi_pot_Nwell_norm(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			double pot = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++) {
				if (beta->index == phase.index || beta->_flag == pf_BULK)
					continue;
				pot += 18.0 / intface_width * _xi_ab(node, phase, *beta) * phase.phi * beta->phi * beta->phi;
				for (auto gamma = beta + 1; gamma < node.end(); gamma++) {
					if (gamma->index == phase.index || gamma->_flag == pf_BULK)
						continue;
					pot += 2.0 / intface_width * _xi_abc(node, phase, *beta, *gamma) * phase.phi
						* beta->phi * gamma->phi * beta->phi * gamma->phi;
				}
			}
			return pot;
		};
		static double dfint_dphi_grad_standard(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return -_xi_a(node, phase) * phase.laplacian;
		};
		static double dfint_dphi_pairwise_acc(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1996_acc(node, phase, intface_width) + dfint_dphi_pot_Nobstacle_acc(node, phase, intface_width);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1996_acc(node, phase, intface_width) + dfint_dphi_pot_Nwell_acc(node, phase, intface_width);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1999_acc(node, phase, intface_width) + dfint_dphi_pot_Nobstacle_acc(node, phase, intface_width);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1999_acc(node, phase, intface_width) + dfint_dphi_pot_Nwell_acc(node, phase, intface_width);
			}
			else {
				Solvers::get_instance()->writer.add_string_to_txt_and_screen("> ERROR: Phi Pair-Wise models hasn't been defined !\n", LOG_FILE_NAME);
				exit(0);
			}
		};
		static double dfint_dphi_pairwise_norm(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width) {
			if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1996_norm(node, phase, intface_width) + dfint_dphi_pot_Nobstacle_norm(node, phase, intface_width);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1996_norm(node, phase, intface_width) + dfint_dphi_pot_Nwell_norm(node, phase, intface_width);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1999_norm(node, phase, intface_width) + dfint_dphi_pot_Nobstacle_norm(node, phase, intface_width);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1999_norm(node, phase, intface_width) + dfint_dphi_pot_Nwell_norm(node, phase, intface_width);
			}
			else {
				Solvers::get_instance()->writer.add_string_to_txt_and_screen("> ERROR: Phi Pair-Wise models hasn't been defined !\n", LOG_FILE_NAME);
				exit(0);
			}
		}

		// From DOI 10.1007/978-3-319-41196-5 CHAP 4 SEC 7
		// calculate gradient energy terms (three terms)
		static double dfint_dphi_aniso_cos_e1(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			pf::Vector3 e_de_grad = node.cal_customValues_gradient(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff_Derived);
			double epsilon = node.customValues[pf::ExternalFieldsPlus::EFP_AnisoGradCoeff];
			return 0.5 * (phase.phi_grad.getY() * e_de_grad.getX() - phase.phi_grad.getX() * e_de_grad.getY()) - epsilon * phase.laplacian;
		}
		static double dfint_dphi_aniso_cos_e2(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			pf::Vector3 e_de_grad = node.cal_customValues_gradient(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff_Derived);
			double epsilon = node.customValues[pf::ExternalFieldsPlus::EFP_AnisoGradCoeff];
			return phase.phi_grad.getY() * e_de_grad.getX() - phase.phi_grad.getX() * e_de_grad.getY() - epsilon * epsilon * phase.laplacian;
		}

		static void cal_interface_increment_ac_pair_wise_normal(pf::PhaseNode& node, bool adjust_phi_0_1) {
			//double first_term, second_term;
			auto mobility = Solvers::get_instance()->Phi_Solver_AC.Lij;
			///< simplified version, to avoid issues
			if (interface_gradient == Int_Gradient::Steinbach_G2009 || interface_potential == Int_Potential::Steinbach_P2009) {
				for (auto alpha = node.begin(); alpha < node.end() - 1; alpha++)
					for (auto beta = alpha + 1; beta < node.end(); beta++)
						if (alpha->_flag && beta->_flag) {
							double int_incre_b_a = 0.0;
							int_incre_b_a = mobility(node, *alpha, *beta) / interface_width * dfint_dphi_S2009(node, *alpha, *beta, interface_width);
							if (adjust_phi_0_1) {
								if ((int_incre_b_a > SYS_EPSILON && (alpha->phi > (1.0 - SYS_EPSILON) || beta->phi < SYS_EPSILON))
									|| (int_incre_b_a < -SYS_EPSILON && (alpha->phi < SYS_EPSILON || beta->phi >(1.0 - SYS_EPSILON))))
									int_incre_b_a = 0.0;
							}
							alpha->int_increment += int_incre_b_a;
							beta->int_increment -= int_incre_b_a;
#ifdef _DEBUG
							if (_isnan(int_incre_b_a)) {
								cout << "DEBUG: int_incre_b_a interface error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
			}
			else {
				for (auto alpha = node.begin(); alpha < node.end() - 1; alpha++)
					for (auto beta = alpha + 1; beta < node.end(); beta++)
						if (alpha->_flag && beta->_flag) {
							double int_incre_b_a = mobility(node, *alpha, *beta) / interface_width
								* (dfint_dphi_pairwise_norm(node, *beta, interface_width) - dfint_dphi_pairwise_norm(node, *alpha, interface_width));
							if (adjust_phi_0_1) {
								bool check1 = int_incre_b_a > SYS_EPSILON && (alpha->phi > (1.0 - SYS_EPSILON) || beta->phi < SYS_EPSILON),
									check2 = int_incre_b_a < (-SYS_EPSILON) && (alpha->phi < SYS_EPSILON || beta->phi >(1.0 - SYS_EPSILON));
								if (check1 || check2) {
									int_incre_b_a = 0.0;
								}
							}
							alpha->int_increment += int_incre_b_a;
							beta->int_increment -= int_incre_b_a;
#ifdef _DEBUG
							if (_isnan(int_incre_b_a)) {
								cout << "DEBUG: int_incre_b_a interface error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
			}
		}
		static void cal_interface_increment_ac_pair_wise_accelarate(pf::PhaseNode& node, bool adjust_phi_0_1) {
			//double first_term, second_term;
			auto mobility = Solvers::get_instance()->Phi_Solver_AC.Lij;
			///< simplified version, to avoid issues
			if (interface_gradient == Int_Gradient::Steinbach_G2009 || interface_potential == Int_Potential::Steinbach_P2009) {
				for (auto alpha = node.Goast_Phase.begin(); alpha < node.Goast_Phase.end() - 1; alpha++)
					for (auto beta = alpha + 1; beta < node.Goast_Phase.end(); beta++) {
						double int_incre_b_a = 0.0;
						int_incre_b_a = mobility(node, **alpha, **beta) / interface_width * dfint_dphi_S2009(node, **alpha, **beta, interface_width);
						if (adjust_phi_0_1) {
							if ((int_incre_b_a > SYS_EPSILON && ((*alpha)->phi > (1.0 - SYS_EPSILON) || (*beta)->phi < SYS_EPSILON))
								|| (int_incre_b_a < -SYS_EPSILON && ((*alpha)->phi < SYS_EPSILON || (*beta)->phi >(1.0 - SYS_EPSILON))))
								int_incre_b_a = 0.0;
						}
						(*alpha)->int_increment += int_incre_b_a;
						(*beta)->int_increment -= int_incre_b_a;
#ifdef _DEBUG
						if (_isnan(int_incre_b_a)) {
							cout << "DEBUG: int_incre_b_a interface error !" << endl;
							SYS_PROGRAM_STOP;
						}
#endif
					}
			}
			else {
				for (auto alpha = node.Goast_Phase.begin(); alpha < node.Goast_Phase.end() - 1; alpha++)
					for (auto beta = alpha + 1; beta < node.Goast_Phase.end(); beta++) {
						double int_incre_b_a = mobility(node, **alpha, **beta) / interface_width
							* (dfint_dphi_pairwise_acc(node, **beta, interface_width) - dfint_dphi_pairwise_acc(node, **alpha, interface_width));
						if (adjust_phi_0_1) {
							bool check1 = int_incre_b_a > SYS_EPSILON && ((*alpha)->phi > (1.0 - SYS_EPSILON) || (*beta)->phi < SYS_EPSILON),
								check2 = int_incre_b_a < (-SYS_EPSILON) && ((*alpha)->phi < SYS_EPSILON || (*beta)->phi >(1.0 - SYS_EPSILON));
							if (check1 || check2) {
								int_incre_b_a = 0.0;
							}
						}
						(*alpha)->int_increment += int_incre_b_a;
						(*beta)->int_increment -= int_incre_b_a;
#ifdef _DEBUG
						if (_isnan(int_incre_b_a)) {
							cout << "DEBUG: int_incre_b_a interface error !" << endl;
							SYS_PROGRAM_STOP;
						}
#endif
					}
			}

		}
		static void cal_interface_increment_ac_standard(pf::PhaseNode& node, bool adjust_phi_0_1) {
			auto mobility = Solvers::get_instance()->Phi_Solver_AC.Lij;
			for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
				alpha->int_increment = 0.0;
				if (alpha->laplacian > SYS_EPSILON || alpha->laplacian < -SYS_EPSILON) {
					for (auto beta = node.begin(); beta < node.end(); beta++)
						if (beta->laplacian > SYS_EPSILON || beta->laplacian < -SYS_EPSILON) {
							alpha->int_increment += -mobility(node, *alpha, *beta) * dfint_dphi_grad_standard(node, *beta);
						}
				}
			}
		}
		static void cal_interface_aniso_cos_e1(pf::PhaseNode& node, bool adjust_phi_0_1) {
			auto mobility = Solvers::get_instance()->Phi_Solver_AC.Lij;
			// calculate dfint_dphi for every nodes
			for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
				alpha->int_increment = 0.0;
				if (alpha->laplacian > SYS_EPSILON || alpha->laplacian < -SYS_EPSILON) {
					for (auto beta = node.begin(); beta < node.end(); beta++)
						if (beta->laplacian > SYS_EPSILON || beta->laplacian < -SYS_EPSILON) {
							alpha->int_increment += -mobility(node, *alpha, *beta) * dfint_dphi_aniso_cos_e1(node, *beta);
						}
				}
			}
		}
		static void cal_interface_aniso_cos_e2(pf::PhaseNode& node, bool adjust_phi_0_1) {
			auto mobility = Solvers::get_instance()->Phi_Solver_AC.Lij;
			// calculate dfint_dphi for every nodes
			for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
				alpha->int_increment = 0.0;
				if (alpha->laplacian > SYS_EPSILON || alpha->laplacian < -SYS_EPSILON) {
					for (auto beta = node.begin(); beta < node.end(); beta++)
						if (beta->laplacian > SYS_EPSILON || beta->laplacian < -SYS_EPSILON) {
							alpha->int_increment += -mobility(node, *alpha, *beta) * dfint_dphi_aniso_cos_e2(node, *beta);
						}
				}
			}
		}

		static double cal_interface_increment_ch_standard(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return dfint_dphi_grad_standard(node, phase);
		}

		static void load_interface_anisotropic_model(bool infile_debug) {
			InputFileReader::get_instance()->read_bool_value("ModelsManager.Phi.InterfaceEnergy.is_anisotropy_on", is_anisotropy_on, infile_debug);
			if (is_anisotropy_on) {
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.anisotropy_model = 0: no anisotropic; 1: 1+\\delta\\cos(n\\theta), e^1; 2: 1+\\delta\\cos(n\\theta), e^2\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.Phi.InterfaceEnergy.anisotropy_model", anisotropy_model, infile_debug);
				std::string aniso_model_key{ "ModelsManager.Phi.InterfaceEnergy.cos_model_parameters" };
				std::vector<input_value> aniso_model_para;
				std::string aniso_model_value{ "()" };
				switch (anisotropy_model)
				{
				case(Int_AnisoModel::Aniso_Cos_E1): {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.cos_model_paras = (aniso_strength, aniso_module_num, aniso_angle) \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->read_string_value(aniso_model_key, aniso_model_value, infile_debug);
					aniso_model_para = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, aniso_model_key, aniso_model_value, infile_debug);
					aniso_strength = aniso_model_para[0].double_value;
					aniso_module_num = aniso_model_para[1].double_value;
					aniso_angle = AngleToRadians(aniso_model_para[2].double_value);
					Solvers::get_instance()->Phi_Solver_AC.dfint_dphi = cal_interface_aniso_cos_e1;
					break;
				}
				case(Int_AnisoModel::Aniso_Cos_E2): {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.cos_model_paras = (aniso_strength, aniso_module_num, aniso_angle) \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->read_string_value(aniso_model_key, aniso_model_value, infile_debug);
					aniso_model_para = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, aniso_model_key, aniso_model_value, infile_debug);
					aniso_strength = aniso_model_para[0].double_value;
					aniso_module_num = aniso_model_para[1].double_value;
					aniso_angle = AngleToRadians(aniso_model_para[2].double_value);
					Solvers::get_instance()->Phi_Solver_AC.dfint_dphi = cal_interface_aniso_cos_e2;
					break;
				}
				default:
					InputFileReader::get_instance()->debug_writer->add_string_to_txt_and_screen("Error, have you indicated the correct anisotropy model?", LOG_FILE_NAME);
					exit(0);
					break;
				}
			}
		}

		static void load_interface_energy_model(bool infile_debug) {
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.int_gradient : 0 - Steinbach_1996 , 1 - Steinbach_1999 , 2 - Steinbach_G2009\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.Phi.InterfaceEnergy.int_gradient", interface_gradient, infile_debug);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.int_potential : 0 - Nestler_Well , 1 - Nestler_Obstacle , 2 - Steinbach_P2009\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.Phi.InterfaceEnergy.int_potential", interface_potential, infile_debug);

				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.InterfaceEnergy.int_width", interface_width, infile_debug);

				Solvers::get_instance()->Phi_Solver_AC.dfint_dphi = cal_interface_increment_ac_pair_wise_normal;

				if (Solvers::get_instance()->parameters.PhiEType == ConEquationType::CEType_PhaseX) {
					Solvers::get_instance()->C_Solver.abs_grad_phi_AB = _abs_grad_phi_pairwise;
					if (interface_gradient == Int_Gradient::Steinbach_1996)
						_abs_grad_phi_pairwise = interface_energy_funcs::abs_grad_phi_S1996;
					else if (interface_gradient == Int_Gradient::Steinbach_1999 || interface_gradient == Int_Gradient::Steinbach_G2009)
						_abs_grad_phi_pairwise = interface_energy_funcs::abs_grad_phi_S1999;
				}

				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.xi_ab.const  = xi_ab \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                        .matrix = [(phi_a, phi_b, xi_ab_value), ...] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key1 = "ModelsManager.Phi.xi_ab.matrix", matrix_input1 = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.xi_ab.const", const_xi_ab, infile_debug)) {
					_xi_ab = xi_ab_const;
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key1, matrix_input1, infile_debug)) {
					_xi_ab = xi_ab_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key1, matrix_input1, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_xi_ab.add_double(matrix_value[index][0].int_value, matrix_value[index][1].int_value, matrix_value[index][2].double_value);
				}
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.xi_abc.const  = xi_ab \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                         .matrix = [(phi_a, phi_b, phi_c, xi_abc_value), ...] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key2 = "ModelsManager.Phi.xi_abc.matrix", matrix_input2 = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.xi_abc.const", const_xi_abc, infile_debug)) {
					_xi_abc = xi_abc_const;
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key2, matrix_input2, infile_debug)) {
					_xi_abc = xi_abc_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT)
						; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key2, matrix_input2, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_xi_abc.add_double(matrix_value[index][0].int_value, matrix_value[index][1].int_value, matrix_value[index][2].int_value, matrix_value[index][3].double_value);
				}
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.int_gradient : 3 - Int_GStandard\n", InputFileReader::get_instance()->debug_file);
				interface_gradient = 3;
				InputFileReader::get_instance()->read_int_value("ModelsManager.Phi.InterfaceEnergy.int_gradient", interface_gradient, infile_debug);

				Solvers::get_instance()->Phi_Solver_AC.dfint_dphi = cal_interface_increment_ac_standard;

				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.xi_a.const  = xi_a \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                       .matrix = [(phi_a, xi_a_value), ...] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key1 = "ModelsManager.Phi.xi_a.matrix", matrix_input1 = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.xi_a.const", const_xi_a, infile_debug)) {
					_xi_a = xi_a_const;
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key1, matrix_input1, infile_debug)) {
					_xi_a = xi_a_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key1, matrix_input1, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_xi_a.add_double(matrix_value[index][0].int_value, matrix_value[index][1].double_value);
				}
				// borrow xi_a from before
				load_interface_anisotropic_model(infile_debug);
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.int_gradient : 3 - Int_GStandard\n", InputFileReader::get_instance()->debug_file);
				interface_gradient = 3;
				InputFileReader::get_instance()->read_int_value("ModelsManager.Phi.InterfaceEnergy.int_gradient", interface_gradient, infile_debug);

				Solvers::get_instance()->Phi_Solver_CH.dfint_dphi = cal_interface_increment_ch_standard;

				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.xi_a.const  = xi_a \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                       .matrix = [(phi_a, xi_a_value), ...] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key1 = "ModelsManager.Phi.xi_a.matrix", matrix_input1 = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.xi_a.const", const_xi_a, infile_debug)) {
					_xi_a = xi_a_const;
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key1, matrix_input1, infile_debug)) {
					_xi_a = xi_a_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key1, matrix_input1, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_xi_a.add_double(matrix_value[index][0].int_value, matrix_value[index][1].double_value);
				}
			}

		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			_xi_a = xi_a_const;
			_xi_ab = xi_ab_const;
			_xi_abc = xi_abc_const;
			_abs_grad_phi_pairwise = interface_energy_funcs::abs_grad_phi_S1996;
			_abs_grad_phi_standard = interface_energy_funcs::abs_grad_phi_standard;

			load_interface_energy_model(infile_debug);

			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				string df_dphi_key = "ModelsManager.Phi.InterfaceEnergy.vts_output", df_dphi_input = "[()]";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.vts_output = [(phi_index_0, phi_index_1), ... ]\n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value(df_dphi_key, df_dphi_input, infile_debug)) {
					vector<vector<input_value>> df_dphi_value;
					df_dphi_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, df_dphi_key, df_dphi_input, infile_debug);
					for (int i = 0; i < df_dphi_value.size(); i++)
						interface_energy_pairwise_output.set(df_dphi_value[i][0].int_value, df_dphi_value[i][1].int_value, 1);
				}
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard ||
				Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				string df_dphi_key = "ModelsManager.Phi.InterfaceEnergy.vts_output", df_dphi_input = "()";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.InterfaceEnergy.vts_output = (phi_index_0, phi_index_1, ... )\n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value(df_dphi_key, df_dphi_input, infile_debug)) {
					vector<input_value> df_dphi_value;
					df_dphi_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, df_dphi_key, df_dphi_input, infile_debug);
					for (int i = 0; i < df_dphi_value.size(); i++)
						interface_energy_standard_output.add_int(df_dphi_value[i].int_value, 1);
				}
			}

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : InterfaceEnergy !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_anisotropy_on) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = (phaseMesh)(x, y, z);
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								double xi = _xi_a(node, *phase);
								switch (anisotropy_model) {
								case(Int_AnisoModel::Aniso_Cos_E1): {
									node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff, xi);
									break;
								}
								case(Int_AnisoModel::Aniso_Cos_E2): {
									node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff, xi * xi);
									break;
								}
								default: {
									break;
								}
								}
								node.customValues.add_double(pf::ExternalFieldsPlus::EFP_AnisoGradCoeff_Derived, 0.0);
							}
						}
			}
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_anisotropy_on) {
				switch (anisotropy_model)
				{
				case(Int_AnisoModel::Aniso_Cos_E1): {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = (phaseMesh)(x, y, z);
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									_xi_a_aniso_cos_e1(node, *phase);
								}
							}
					break;
				}
				case(Int_AnisoModel::Aniso_Cos_E2): {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = (phaseMesh)(x, y, z);
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									_xi_a_aniso_cos_e2(node, *phase);
								}
							}
					break;
				}
				default:
					break;
				}

			}
			string report = "";
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (interface_energy_pairwise_output.pairValue.size() != 0
				&& Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "dfintdphi_merge" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double val = 0.0, merge = 0.0;
							for (auto f = interface_energy_pairwise_output.begin(); f < interface_energy_pairwise_output.end(); f++) {
								if (node[f->pairIndex_1]._flag && node[f->pairIndex_2]._flag) {
									if (interface_gradient == Int_Gradient::Steinbach_G2009 || interface_potential == Int_Potential::Steinbach_P2009) {
										val += dfint_dphi_S2009(node, node[f->pairIndex_1], node[f->pairIndex_2], interface_width)
											* (node[f->pairIndex_1].phi + node[f->pairIndex_2].phi);
										merge += node[f->pairIndex_1].phi + node[f->pairIndex_2].phi;
									}
									else {
										val += (dfint_dphi_pairwise_norm(node, node[f->pairIndex_2], interface_width) - dfint_dphi_pairwise_norm(node, node[f->pairIndex_1], interface_width))
											* (node[f->pairIndex_1].phi + node[f->pairIndex_2].phi);
										merge += node[f->pairIndex_1].phi + node[f->pairIndex_2].phi;
									}
								}
							}
							if (merge > SYS_EPSILON)
								fout << val / merge << endl;
							else
								fout << 0.0 << endl;
						}
				fout << "</DataArray>" << endl;
			}
			else if (interface_energy_standard_output.size() != 0
				&& (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard
					|| Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard)) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "dfintdphi_merge" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double val = 0.0, merge = 0.0;
							for (auto f = interface_energy_standard_output.begin(); f < interface_energy_standard_output.end(); f++) {
								if (abs(node[f->index].laplacian) > SYS_EPSILON) {
									val += dfint_dphi_grad_standard(node, node[f->index]) * node[f->index].phi;
									merge += node[f->index].phi;
								}
							}
							if (merge > SYS_EPSILON)
								fout << val / merge << endl;
							else
								fout << 0.0 << endl;
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}

	};
}