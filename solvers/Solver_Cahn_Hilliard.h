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
#include"base.h"
namespace pf {
	enum CH_RETURN { CH_MAX_DIFFUSION_FLUX, CH_MAX_REACTION_FLUX, CH_MAX_PHI_TRANS_FLUX, CH_MAX_X_INCREMENT_FLUX};
	namespace cahnHilliardSolver {
		static void dfbulk_dx(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static void Mbulk_ij(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static double abs_grad_phi_AB(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			Vector3 vec;
			vec[0] = alpha.phi_grad[0] * beta.phi - alpha.phi * beta.phi_grad[0];
			vec[1] = alpha.phi_grad[1] * beta.phi - alpha.phi * beta.phi_grad[1];
			vec[2] = alpha.phi_grad[2] * beta.phi - alpha.phi * beta.phi_grad[2];
			return vec.abs();
		}
		static void Source_A(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static void Source_AB(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double AbsGradPhi_Phi) {
			// AbsGradPhi_Phi * ( standard_reaction_flux_on_interface )
			return;
		}
		static void boundary_condition_phi(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static double dfint_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		static double dfbulk_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		static double M_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return 0.0;
		}
		static double Source_a(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}

		static double df_dx(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		static double M_ij(pf::PhaseNode& node, int con_i, int con_j) {
			return 0.0;
		}
		static double int_flux(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		static double Source(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		static double dphase_x_du(pf::PhaseNode& node, pf::PhaseEntry& phase, int) {
			return 0.0;
		}
		static void phase_x(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static void boundary_condition_total_x(pf::PhaseNode& node, int comp_index) {
			return;
		}
		static void boundary_condition_phase_x(pf::PhaseNode& node, pf::PhaseEntry& phase, int comp_index) {
			return;
		}
		static void init_phase_con_on_moving_interface(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			if (phase.phi > Phi_Num_Cut_Off) {
				// con diffuse/generate to new generate node
				if (phase.old_phi < Phi_Num_Cut_Off) {
					/*bool is_con_diffuse = false;
					double sum_eff_phis = 0.0;
					for(int r_x = -phase_con_phase_transition_optimize_range; r_x <= phase_con_phase_transition_optimize_range; r_x++)
						for (int r_y = -phase_con_phase_transition_optimize_range; r_y <= phase_con_phase_transition_optimize_range; r_y++)
							for (int r_z = -phase_con_phase_transition_optimize_range; r_z <= phase_con_phase_transition_optimize_range; r_z++) {
								if ((r_x == 0 && r_y == 0 && r_z == 0) ||
									(node._x + r_x) < 0 || (node._x + r_x) > phaseMesh->limit_x - 1 ||
									(node._y + r_y) < 0 || (node._y + r_y) > phaseMesh->limit_y - 1 ||
									(node._z + r_z) < 0 || (node._z + r_z) > phaseMesh->limit_z - 1)
									continue;
								PhaseEntry& r_phase = node.get_long_range_node(r_x, r_y, r_z)[phase->index];
								if (r_phase.phi > Phi_Num_Cut_Off && r_phase._cflag == CFlag::pf_Comp_Cal) {
									is_con_diffuse = true;
									sum_eff_phis += r_phase.phi;
									for (auto c = phase->x.begin(); c < phase->x.end(); c++)
										if(c->index != solvent)
											c->value = r_phase.x[c->index].value * r_phase.phi;
								}
							}
					if (is_con_diffuse) {
						double solute = 0.0;
						for (auto c = phase->x.begin(); c < phase->x.end(); c++)
							if (c->index != solvent) {
								c->value /= (sum_eff_phis + phase->phi);
								solute += c->value;
							}
						for (auto c = phase->x.begin(); c < phase->x.end(); c++)
							if (c->index == solvent)
								c->value = 1.0 - solute;
						for (int r_x = -phase_con_phase_transition_optimize_range; r_x <= phase_con_phase_transition_optimize_range; r_x++)
							for (int r_y = -phase_con_phase_transition_optimize_range; r_y <= phase_con_phase_transition_optimize_range; r_y++)
								for (int r_z = -phase_con_phase_transition_optimize_range; r_z <= phase_con_phase_transition_optimize_range; r_z++) {
									if ((r_x == 0 && r_y == 0 && r_z == 0) ||
										(node._x + r_x) < 0 || (node._x + r_x) > phaseMesh->limit_x - 1 ||
										(node._y + r_y) < 0 || (node._y + r_y) > phaseMesh->limit_y - 1 ||
										(node._z + r_z) < 0 || (node._z + r_z) > phaseMesh->limit_z - 1)
										continue;
									PhaseEntry& r_phase = node.get_long_range_node(r_x, r_y, r_z)[phase->index];
									if (r_phase.phi > Phi_Num_Cut_Off && r_phase._cflag == CFlag::pf_Comp_Cal) {
										double r_solute = 0.0;
										for (auto c = r_phase.x.begin(); c < r_phase.x.end(); c++)
											if (c->index != solvent) {
												c->value -= phase->x[c->index].value * phase->phi / sum_eff_phis;
												r_solute += c->value;
											}
										for (auto c = r_phase.x.begin(); c < r_phase.x.end(); c++)
											if (c->index == solvent)
												c->value = 1.0 - r_solute;
									}
								}
					}
					else {
						double_box sum_eff_phis;
						for (auto c = phase->x.begin(); c < phase->x.end(); c++) {
							c->value = 0.0;
							sum_eff_phis.add_double(c->index, 0.0);
						}
						for (auto r_phase = node.begin(); r_phase < node.end(); r_phase++) {
							if (r_phase->index == phase->index)
								continue;
							if (r_phase->phi > Phi_Num_Cut_Off && r_phase->_cflag == CFlag::pf_Comp_Cal) {
								for (auto r_c = r_phase->x.begin(); r_c < r_phase->x.end(); r_c++)
									for (auto c = phase->x.begin(); c < phase->x.end(); c++)
										if (r_c->index == c->index && r_c->index != solvent) {
											sum_eff_phis[c->index] += r_phase->phi;
											c->value += r_c->value * r_phase->phi;
										}
							}
						}
						double solute = 0.0;
						for (auto c = phase->x.begin(); c < phase->x.end(); c++)
							if (c->index != solvent && sum_eff_phis[c->index] > Phi_Num_Cut_Off) {
								c->value /= (sum_eff_phis[c->index] + phase->phi);
								solute += c->value;
							}
						for (auto c = phase->x.begin(); c < phase->x.end(); c++)
							if (c->index == solvent)
								c->value = 1.0 - solute;
						for (auto r_phase = node.begin(); r_phase < node.end(); r_phase++) {
							if (r_phase->index == phase->index)
								continue;
							double r_solute = 0.0;
							if (r_phase->phi > Phi_Num_Cut_Off && r_phase->_cflag == CFlag::pf_Comp_Cal) {
								for (auto r_c = r_phase->x.begin(); r_c < r_phase->x.end(); r_c++) {
									for (auto c = phase->x.begin(); c < phase->x.end(); c++)
										if (r_c->index == c->index && r_c->index != solvent) {
											r_c->value -= c->value * phase->phi / sum_eff_phis[c->index];
										}
									r_solute += r_c->value;
								}
								for (auto r_c = r_phase->x.begin(); r_c < r_phase->x.end(); r_c++)
									if (r_c->index == solvent)
										r_c->value = 1.0 - r_solute;
							}
						}
					}*/
					PhaseEntry& upx_phase = node.get_neighbor_node(Direction::x_up)[phase.index];
					PhaseEntry& downx_phase = node.get_neighbor_node(Direction::x_down)[phase.index];
					PhaseEntry& upy_phase = node.get_neighbor_node(Direction::y_up)[phase.index];
					PhaseEntry& downy_phase = node.get_neighbor_node(Direction::y_down)[phase.index];
					PhaseEntry& upz_phase = node.get_neighbor_node(Direction::z_up)[phase.index];
					PhaseEntry& downz_phase = node.get_neighbor_node(Direction::z_down)[phase.index];
					if (upx_phase.phi > Phi_Num_Cut_Off && upx_phase.old_phi > Phi_Num_Cut_Off) {
						for (auto c = phase.x.begin(); c < phase.x.end(); c++)
							c->value = upx_phase.x[c->index].value;
					}
					else if (downx_phase.phi > Phi_Num_Cut_Off && downx_phase.old_phi > Phi_Num_Cut_Off) {
						for (auto c = phase.x.begin(); c < phase.x.end(); c++)
							c->value = downx_phase.x[c->index].value;
					}
					else if (upy_phase.phi > Phi_Num_Cut_Off && upy_phase.old_phi > Phi_Num_Cut_Off) {
						for (auto c = phase.x.begin(); c < phase.x.end(); c++)
							c->value = upy_phase.x[c->index].value;
					}
					else if (downy_phase.phi > Phi_Num_Cut_Off && downy_phase.old_phi > Phi_Num_Cut_Off) {
						for (auto c = phase.x.begin(); c < phase.x.end(); c++)
							c->value = downy_phase.x[c->index].value;
					}
					else if (upz_phase.phi > Phi_Num_Cut_Off && upz_phase.old_phi > Phi_Num_Cut_Off) {
						for (auto c = phase.x.begin(); c < phase.x.end(); c++)
							c->value = upz_phase.x[c->index].value;
					}
					else if (downz_phase.phi > Phi_Num_Cut_Off && downz_phase.old_phi > Phi_Num_Cut_Off) {
						for (auto c = phase.x.begin(); c < phase.x.end(); c++)
							c->value = downz_phase.x[c->index].value;
					}
				}
			}
			else {
				// con collect from delete node
				if (phase.old_phi > Phi_Num_Cut_Off) {
					/*bool is_con_diffuse = false;
					double sum_eff_phis = 0.0;
					for (int r_x = -phase_con_phase_transition_optimize_range; r_x <= -phase_con_phase_transition_optimize_range; r_x++)
						for (int r_y = -phase_con_phase_transition_optimize_range; r_y <= -phase_con_phase_transition_optimize_range; r_y++)
							for (int r_z = -phase_con_phase_transition_optimize_range; r_z <= -phase_con_phase_transition_optimize_range; r_z++) {
								if (r_x == 0 && r_y == 0 && r_z == 0)
									continue;
								PhaseEntry& r_phase = node.get_long_range_node(r_x, r_y, r_z)[phase->index];
								if (r_phase.phi > Phi_Num_Cut_Off && r_phase._cflag == CFlag::pf_Comp_Cal) {
									is_con_diffuse = true;
									sum_eff_phis += r_phase.phi;
								}
							}
					if (is_con_diffuse) {
						for (int r_x = -phase_con_phase_transition_optimize_range; r_x <= -phase_con_phase_transition_optimize_range; r_x++)
							for (int r_y = -phase_con_phase_transition_optimize_range; r_y <= -phase_con_phase_transition_optimize_range; r_y++)
								for (int r_z = -phase_con_phase_transition_optimize_range; r_z <= -phase_con_phase_transition_optimize_range; r_z++) {
									if (r_x == 0 && r_y == 0 && r_z == 0)
										continue;
									PhaseEntry& r_phase = node.get_long_range_node(r_x, r_y, r_z)[phase->index];
									if (r_phase.phi > Phi_Num_Cut_Off && r_phase._cflag == CFlag::pf_Comp_Cal) {
										double r_solute = 0.0;
										for (auto c = r_phase.x.begin(); c < r_phase.x.end(); c++)
											if (c->index != solvent) {
												c->value += phase->x[c->index].value * phase->phi / sum_eff_phis;
												r_solute += c->value;
											}
										for (auto c = r_phase.x.begin(); c < r_phase.x.end(); c++)
											if (c->index == solvent)
												c->value = 1.0 - r_solute;
									}
								}
					}
					else {
						double_box sum_eff_phis;
						for (auto c = phase->x.begin(); c < phase->x.end(); c++)
							sum_eff_phis.add_double(c->index, 0.0);
						for (auto r_phase = node.begin(); r_phase < node.end(); r_phase++) {
							if (r_phase->index == phase->index)
								continue;
							if (r_phase->phi > Phi_Num_Cut_Off && r_phase->_cflag == CFlag::pf_Comp_Cal) {
								for (auto r_c = r_phase->x.begin(); r_c < r_phase->x.end(); r_c++)
									for (auto c = phase->x.begin(); c < phase->x.end(); c++)
										if (r_c->index == c->index && r_c->index != solvent)
											sum_eff_phis[c->index] += r_phase->phi;
							}
						}
						for (auto r_phase = node.begin(); r_phase < node.end(); r_phase++) {
							if (r_phase->index == phase->index)
								continue;
							double r_solute = 0.0;
							if (r_phase->phi > Phi_Num_Cut_Off && r_phase->_cflag == CFlag::pf_Comp_Cal) {
								for (auto r_c = r_phase->x.begin(); r_c < r_phase->x.end(); r_c++) {
									for (auto c = phase->x.begin(); c < phase->x.end(); c++)
										if (r_c->index == c->index && r_c->index != solvent) {
											r_c->value += c->value * phase->phi / sum_eff_phis[c->index];
										}
									r_solute += r_c->value;
								}
								for (auto r_c = r_phase->x.begin(); r_c < r_phase->x.end(); r_c++)
									if (r_c->index == solvent)
										r_c->value = 1.0 - r_solute;
							}
						}
					}*/
					for (auto c = phase.x.begin(); c < phase.x.end(); c++)
						c->value = 0.0;
				}
			}
		}
		static void init_con_on_moving_interface(pf::PhaseNode& node, ConEquationDomain _domain, double threshold) {
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
					// con diffuse/generate to new generate node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
						PhaseNode& upx_node = node.get_neighbor_node(Direction::x_up);
						PhaseNode& downx_node = node.get_neighbor_node(Direction::x_down);
						PhaseNode& upy_node = node.get_neighbor_node(Direction::y_up);
						PhaseNode& downy_node = node.get_neighbor_node(Direction::y_down);
						PhaseNode& upz_node = node.get_neighbor_node(Direction::z_up);
						PhaseNode& downz_node = node.get_neighbor_node(Direction::z_down);
						if (upx_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && upx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = upx_node.x[c->index].value;
						}
						else if (downx_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && downx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = downx_node.x[c->index].value;
						}
						else if (upy_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && upy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = upy_node.x[c->index].value;
						}
						else if (downy_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && downy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = downy_node.x[c->index].value;
						}
						else if (upz_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && upz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = upz_node.x[c->index].value;
						}
						else if (downz_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && downz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = downz_node.x[c->index].value;
						}
					}
				}
				else {
					// con collect from delete node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
						for (auto c = node.x.begin(); c < node.x.end(); c++)
							c->value = 0.0;
					}
				}
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
					// con diffuse/generate to new generate node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
						PhaseNode& upx_node = node.get_neighbor_node(Direction::x_up);
						PhaseNode& downx_node = node.get_neighbor_node(Direction::x_down);
						PhaseNode& upy_node = node.get_neighbor_node(Direction::y_up);
						PhaseNode& downy_node = node.get_neighbor_node(Direction::y_down);
						PhaseNode& upz_node = node.get_neighbor_node(Direction::z_up);
						PhaseNode& downz_node = node.get_neighbor_node(Direction::z_down);
						if (upx_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = upx_node.x[c->index].value;
						}
						else if (downx_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = downx_node.x[c->index].value;
						}
						else if (upy_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = upy_node.x[c->index].value;
						}
						else if (downy_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = downy_node.x[c->index].value;
						}
						else if (upz_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = upz_node.x[c->index].value;
						}
						else if (downz_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto c = node.x.begin(); c < node.x.end(); c++)
								c->value = downz_node.x[c->index].value;
						}
					}
				}
				else {
					// con collect from delete node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
						for (auto c = node.x.begin(); c < node.x.end(); c++)
							c->value = 0.0;
					}
				}
			}
		}
		static void init_grand_potential_on_moving_interface(pf::PhaseNode& node, ConEquationDomain _domain, double threshold) {
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
					// con diffuse/generate to new generate node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
						PhaseNode& upx_node = node.get_neighbor_node(Direction::x_up);
						PhaseNode& downx_node = node.get_neighbor_node(Direction::x_down);
						PhaseNode& upy_node = node.get_neighbor_node(Direction::y_up);
						PhaseNode& downy_node = node.get_neighbor_node(Direction::y_down);
						PhaseNode& upz_node = node.get_neighbor_node(Direction::z_up);
						PhaseNode& downz_node = node.get_neighbor_node(Direction::z_down);
						if (upx_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && upx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = upx_node.potential[p->index].value;
						}
						else if (downx_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && downx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = downx_node.potential[p->index].value;
						}
						else if (upy_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && upy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = upy_node.potential[p->index].value;
						}
						else if (downy_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && downy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = downy_node.potential[p->index].value;
						}
						else if (upz_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && upz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = upz_node.potential[p->index].value;
						}
						else if (downz_node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && downz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = downz_node.potential[p->index].value;
						}
					}
				}
				else {
					// con collect from delete node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
						for (auto p = node.potential.begin(); p < node.potential.end(); p++)
							p->value = 0.0;
					}
				}
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
					// con diffuse/generate to new generate node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
						PhaseNode& upx_node = node.get_neighbor_node(Direction::x_up);
						PhaseNode& downx_node = node.get_neighbor_node(Direction::x_down);
						PhaseNode& upy_node = node.get_neighbor_node(Direction::y_up);
						PhaseNode& downy_node = node.get_neighbor_node(Direction::y_down);
						PhaseNode& upz_node = node.get_neighbor_node(Direction::z_up);
						PhaseNode& downz_node = node.get_neighbor_node(Direction::z_down);
						if (upx_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = upx_node.potential[p->index].value;
						}
						else if (downx_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = downx_node.potential[p->index].value;
						}
						else if (upy_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = upy_node.potential[p->index].value;
						}
						else if (downy_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = downy_node.potential[p->index].value;
						}
						else if (upz_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = upz_node.potential[p->index].value;
						}
						else if (downz_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++)
								p->value = downz_node.potential[p->index].value;
						}
					}
				}
				else {
					// con collect from delete node
					if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
						for (auto p = node.potential.begin(); p < node.potential.end(); p++)
							p->value = 0.0;
					}
				}
			}
		}
	}

	class CahnHilliardSolver
	{
	public:
		CahnHilliardSolver() {
			grand_potential_range.resize(2);
			grand_potential_range[0] = -SYS_EPSILON;
			grand_potential_range[1] = SYS_EPSILON;
		};
		CahnHilliardSolver(FieldStorage_forPhaseNode& _phaseMesh) {
			init(_phaseMesh);
		}
		~CahnHilliardSolver() {
			free();
		}
		void init(FieldStorage_forPhaseNode& _phaseMesh) {
			phaseMesh = &_phaseMesh;

			Boundary_Condition_TotalX = cahnHilliardSolver::boundary_condition_total_x;
			Boundary_Condition_PhaseX = cahnHilliardSolver::boundary_condition_phase_x;
			Boundary_Condition_Phi = cahnHilliardSolver::boundary_condition_phi;

			dfint_dphi = cahnHilliardSolver::dfint_dphi;
			dfbulk_dphi = cahnHilliardSolver::dfbulk_dphi;
			M_ab = cahnHilliardSolver::M_ab;
			Source_a = cahnHilliardSolver::Source_a;

			dfbulk_dx = cahnHilliardSolver::dfbulk_dx;
			Mbulk_ij = cahnHilliardSolver::Mbulk_ij;
			abs_grad_phi_AB = cahnHilliardSolver::abs_grad_phi_AB;
			Source_A = cahnHilliardSolver::Source_A;
			Source_AB = cahnHilliardSolver::Source_AB;

			df_dx = cahnHilliardSolver::df_dx;
			M_ij = cahnHilliardSolver::M_ij;
			int_flux = cahnHilliardSolver::int_flux;
			Source = cahnHilliardSolver::Source;
			dphase_x_du = cahnHilliardSolver::dphase_x_du;
			phase_x = cahnHilliardSolver::phase_x;

			init_phase_con_on_moving_interface = cahnHilliardSolver::init_phase_con_on_moving_interface;
			init_con_on_moving_interface = cahnHilliardSolver::init_con_on_moving_interface;
			init_grand_potential_on_moving_interface = cahnHilliardSolver::init_grand_potential_on_moving_interface;

			phase_con_phase_transition_optimize_range = 1;
			solvent = SOLVENT_NONE;
			threshold = 0.1;
			diff_method = DifferenceMethod::FIVE_POINT;
		}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		void define_funcs_for_phase_concentration(int _solvent = SOLVENT_NONE, DifferenceMethod _diff_method = DifferenceMethod::FIVE_POINT,
			int _phase_con_phase_transition_optimize_range = 1,
			void(*_dfbulk_dx)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::dfbulk_dx,
			void(*_M_ij)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::Mbulk_ij,
			double (*_abs_grad_phi_AB)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&) = cahnHilliardSolver::abs_grad_phi_AB,
			void (*_Source_A)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::Source_A,
			void (*_Source_AB)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, double) = cahnHilliardSolver::Source_AB,
			void(*_Boundary_Condition_PhaseX)(pf::PhaseNode&, pf::PhaseEntry&, int) = cahnHilliardSolver::boundary_condition_phase_x) {
			solvent = _solvent;
			diff_method = _diff_method;
			phase_con_phase_transition_optimize_range = _phase_con_phase_transition_optimize_range;

			dfbulk_dx = _dfbulk_dx;
			Mbulk_ij = _M_ij;
			abs_grad_phi_AB = _abs_grad_phi_AB;
			Source_A = _Source_A;
			Source_AB = _Source_AB;
			Boundary_Condition_PhaseX = _Boundary_Condition_PhaseX;
		}

		void pre_calculation_phase_concentration(double dt);

		// [CH_MAX_DIFFUSION_FLUX, CH_MAX_REACTION_FLUX, CH_MAX_PHI_TRANS_FLUX, CH_MAX_X_INCREMENT_FLUX]
		vector<double> solve_phase_concentration(double dt, bool is_normalized);

		string print_phase_con_model() {
			stringstream rep;
			rep << "Phi Phase-Concentration Model : Vm^-1 * dCai_dt = 1 / Pa * [ //labla * ( Pa * FLUXai ) - //labla(Pa) * FLUXai + Pa * Sai - Vm^-1 * Cai * dPa_dt ]" << endl;
			rep << "                              : FLUXai = \\sum_j[ Maij * \\labla( dfbulk_dCaj )]" << endl;
			return rep.str();
		}

		//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		void define_funcs_for_phi(DifferenceMethod _diff_method = DifferenceMethod::FIVE_POINT,
			double (*_dfint_dphi)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::dfint_dphi,
			double (*_dfbulk_dphi)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::dfbulk_dphi,
			double (*_M_ab)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&) = cahnHilliardSolver::M_ab,
			double (*_Source)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::Source_a,
			void(*_Boundary_Condition_Phi)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::boundary_condition_phi) {
			diff_method = _diff_method;
			dfint_dphi = _dfint_dphi;
			dfbulk_dphi = _dfbulk_dphi;
			M_ab = _M_ab;
			Source_a = _Source;
			Boundary_Condition_Phi = _Boundary_Condition_Phi;
		}

		// int_increment in phase will be used as buff to save df_dx(df_dphi)
		void pre_calculation_phis();

		// MAX_PHI_INCREMENT
		double solve_phis(double dt, bool is_normalized = false);

		string print_phi_model_cahnhilliard() {
			stringstream rep;
			rep << "Phi Cahn-Hilliard Model : dPa_dt = \\labla{ \\sum_b[ Mab * \\labla( dfint_dPb + dfbulk_dPb )]} + Sa" << endl;
			return rep.str();
		}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		void define_funcs_for_total_concentration(vector<int> phiIndexes, int _solvent = SOLVENT_NONE, DifferenceMethod _diff_method = DifferenceMethod::FIVE_POINT,  double _threshold = 0.1,
			double (*_df_dx)(pf::PhaseNode&, int) = cahnHilliardSolver::df_dx,
			double (*_M_ij)(pf::PhaseNode&, int, int) = cahnHilliardSolver::M_ij,
			double (*_int_flux)(pf::PhaseNode&, int) = cahnHilliardSolver::int_flux,
			double (*_Source)(pf::PhaseNode&, int) = cahnHilliardSolver::Source,
			void(*_Boundary_Condition_TotalX)(pf::PhaseNode&, int) = cahnHilliardSolver::boundary_condition_total_x) {
			solvent = _solvent;
			diff_method = _diff_method;
			phase_indexes = phiIndexes;
			threshold = _threshold;
			df_dx = _df_dx;
			M_ij = _M_ij;
			int_flux = _int_flux;
			Source = _Source;
			Boundary_Condition_TotalX = _Boundary_Condition_TotalX;
		}

		void pre_calculation_total_concentration_inside_phis(double dt);

		// [CH_MAX_DIFFUSION_FLUX, CH_MAX_REACTION_FLUX, CH_MAX_PHI_TRANS_FLUX, CH_MAX_X_INCREMENT_FLUX]
		vector<double> solve_total_concentration_inside_phis(double dt, bool is_normalized);

		void pre_calculation_total_concentration_outside_phis(double dt);

		// [CH_MAX_DIFFUSION_FLUX, CH_MAX_REACTION_FLUX, CH_MAX_PHI_TRANS_FLUX, CH_MAX_X_INCREMENT_FLUX]
		vector<double> solve_total_concentration_outside_phis(double dt, bool is_normalized);

		string print_total_con_model() {
			stringstream rep;
			rep << "Phi Total-Concentration Model : Vm^-1 * dCi_dt = 1 / Peff * [ //labla * ( Peff * FLUXi ) - //labla(Peff) * FLUXi + Peff * Si ]" << endl;
			rep << "                              : FLUXi = \\sum_j[ Mij * \\labla( dfbulk_dCj )] , Peff : effective phis " << endl;
			return rep.str();
		}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		void define_funcs_for_grand_potential_functional(vector<int> phiIndexes, DifferenceMethod _diff_method = DifferenceMethod::FIVE_POINT, double _threshold = 0.1,
			double (*_M_ij)(pf::PhaseNode&, int, int) = cahnHilliardSolver::M_ij,
			double (*_int_flux)(pf::PhaseNode&, int) = cahnHilliardSolver::int_flux,
			double (*_Source)(pf::PhaseNode&, int) = cahnHilliardSolver::Source,
			void (*_df_dx)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::dfbulk_dx,
			double (*_dphase_x_du)(pf::PhaseNode&, pf::PhaseEntry&, int) = cahnHilliardSolver::dphase_x_du,
			void (*_phase_x)(pf::PhaseNode&, pf::PhaseEntry&) = cahnHilliardSolver::phase_x,
			void(*_Boundary_Condition_TotalX)(pf::PhaseNode&, int) = cahnHilliardSolver::boundary_condition_total_x) {
			diff_method = _diff_method;
			phase_indexes = phiIndexes;
			threshold = _threshold;
			dfbulk_dx = _df_dx;
			M_ij = _M_ij;
			int_flux = _int_flux;
			Source = _Source;
			dphase_x_du = _dphase_x_du;
			phase_x = _phase_x;
			Boundary_Condition_TotalX = _Boundary_Condition_TotalX;
		}

		// [CH_MAX_DIFFUSION_FLUX, CH_MAX_REACTION_FLUX, CH_MAX_PHI_TRANS_FLUX]
		vector<double> pre_calculation_grand_potential_functional_inside_phis(double dt);

		// CH_MAX_X_INCREMENT_FLUX
		double solve_grand_potential_functional_inside_phis(double dt);

		// [CH_MAX_DIFFUSION_FLUX, CH_MAX_REACTION_FLUX, CH_MAX_PHI_TRANS_FLUX]
		vector<double> pre_calculation_grand_potential_functional_outside_phis(double dt);

		// CH_MAX_X_INCREMENT_FLUX
		double solve_grand_potential_functional_outside_phis(double dt);

		string print_grand_potential_model() {
			stringstream rep;
			rep << "Phi Grand-Potential Model : dUi_dt = [ //sum_a( Pa * Xaij ) ]^-1 * { //labla * \\sum_j[ Mij * \\labla( Uj )] - //abs[ //labla(Peff) ] * REACi + Peff * Si }" << endl;
			rep << "                          : REACi : reaction on Peff interface , Peff : effective phis , Xaij = dCai_dUj " << endl;
			return rep.str();
		}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		void free() {
			phaseMesh = nullptr;
			dfbulk_dx = nullptr;
			Mbulk_ij = nullptr;
			abs_grad_phi_AB = nullptr;
			Source_A = nullptr;
			Source_AB = nullptr;
			phase_indexes.clear();
			df_dx = nullptr;
			M_ij = nullptr;
			int_flux = nullptr;
			Source = nullptr;
			dphase_x_du = nullptr;
			phase_x = nullptr;
		}

		void summary_phix_to_x();
		void summary_phip_to_p();

		int solvent;
		vector<int> phase_indexes;
		double threshold;
		vector<double> grand_potential_range;
		DifferenceMethod diff_method;
		int phase_con_phase_transition_optimize_range;
		void(*Boundary_Condition_Phi)(pf::PhaseNode&, pf::PhaseEntry&);
		void(*Boundary_Condition_TotalX)(pf::PhaseNode&, int);
		void(*Boundary_Condition_PhaseX)(pf::PhaseNode&, pf::PhaseEntry&, int);
		// phi
		double (*dfint_dphi)(pf::PhaseNode&, pf::PhaseEntry&);
		double (*dfbulk_dphi)(pf::PhaseNode&, pf::PhaseEntry&);
		double (*M_ab)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&);
		double (*Source_a)(pf::PhaseNode&, pf::PhaseEntry&);
		// total concentration + grand potential functional
		double (*df_dx)(pf::PhaseNode&, int);
		double (*M_ij)(pf::PhaseNode&, int, int);
		double (*int_flux)(pf::PhaseNode&, int);
		double (*Source)(pf::PhaseNode&, int);
		// grand potential functional
		double (*dphase_x_du)(pf::PhaseNode&, pf::PhaseEntry&, int);
		void (*phase_x)(pf::PhaseNode&, pf::PhaseEntry&);
		// phase concentration
		void(*dfbulk_dx)(pf::PhaseNode&, pf::PhaseEntry&);
		void(*Mbulk_ij)(pf::PhaseNode&, pf::PhaseEntry&);
		double (*abs_grad_phi_AB)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&);
		void (*Source_A)(pf::PhaseNode&, pf::PhaseEntry&);
		void (*Source_AB)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, double);
		// init con & phase con & grand potential on interface
		void (*init_phase_con_on_moving_interface)(pf::PhaseNode&, pf::PhaseEntry&);
		void (*init_con_on_moving_interface)(pf::PhaseNode&, ConEquationDomain, double);
		void (*init_grand_potential_on_moving_interface)(pf::PhaseNode&, ConEquationDomain, double);
	private:
		FieldStorage_forPhaseNode* phaseMesh;
	};
}