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


#include"Solver_Cahn_Hilliard.h"
namespace pf {

	void CahnHilliardSolver::pre_calculation_phis() {
		double MAX_PHI_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						//intphase and near intphase
						phase->phi_grad[0] = (node.get_neighbor_node(Direction::x_up)[phase->index].phi -
							node.get_neighbor_node(Direction::x_down)[phase->index].phi) / 2.0 / phaseMesh->dr;
						phase->phi_grad[1] = (node.get_neighbor_node(Direction::y_up)[phase->index].phi -
							node.get_neighbor_node(Direction::y_down)[phase->index].phi) / 2.0 / phaseMesh->dr;
						phase->phi_grad[2] = (node.get_neighbor_node(Direction::z_up)[phase->index].phi -
							node.get_neighbor_node(Direction::z_down)[phase->index].phi) / 2.0 / phaseMesh->dr;
						if (diff_method == DifferenceMethod::FIVE_POINT) {
							phase->laplacian = (node.get_neighbor_node(Direction::x_down)[phase->index].phi
								+ node.get_neighbor_node(Direction::x_up)[phase->index].phi
								+ node.get_neighbor_node(Direction::y_down)[phase->index].phi
								+ node.get_neighbor_node(Direction::y_up)[phase->index].phi
								+ node.get_neighbor_node(Direction::z_down)[phase->index].phi
								+ node.get_neighbor_node(Direction::z_up)[phase->index].phi - 6 * phase->phi) / phaseMesh->dr / phaseMesh->dr;
						}
						else if (diff_method == DifferenceMethod::NINE_POINT) {
							phase->laplacian = (4.0 * node.get_neighbor_node(Direction::x_down)[phase->index].phi
								+ 4.0 * node.get_neighbor_node(Direction::x_up)[phase->index].phi
								+ 4.0 * node.get_neighbor_node(Direction::y_down)[phase->index].phi
								+ 4.0 * node.get_neighbor_node(Direction::y_up)[phase->index].phi
								+ 4.0 * node.get_neighbor_node(Direction::z_down)[phase->index].phi
								+ 4.0 * node.get_neighbor_node(Direction::z_up)[phase->index].phi
								+ node.get_long_range_node(-1, -1, 0)[phase->index].phi
								+ node.get_long_range_node(-1, 1, 0)[phase->index].phi
								+ node.get_long_range_node(1, -1, 0)[phase->index].phi
								+ node.get_long_range_node(1, 1, 0)[phase->index].phi
								+ node.get_long_range_node(-1, 0, -1)[phase->index].phi
								+ node.get_long_range_node(-1, 0, 1)[phase->index].phi
								+ node.get_long_range_node(1, 0, -1)[phase->index].phi
								+ node.get_long_range_node(1, 0, 1)[phase->index].phi
								+ node.get_long_range_node(0, -1, -1)[phase->index].phi
								+ node.get_long_range_node(0, -1, 1)[phase->index].phi
								+ node.get_long_range_node(0, 1, -1)[phase->index].phi
								+ node.get_long_range_node(0, 1, 1)[phase->index].phi - 36 * phase->phi) / 6.0 / phaseMesh->dr / phaseMesh->dr;
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
							// buff
							phase->int_increment = dfint_dphi(node, *phase);
							phase->bulk_increment = dfbulk_dphi(node, *phase);
						}
					}
				}
	}

	double CahnHilliardSolver::solve_phis(double dt, bool is_normalized) {
		double MAX_PHI_INCREMENT = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					vec3_box grad_dfdphis;
					double_box lap_dfdphis;
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
							//intphase and near intphase
							PhaseEntry& phase_xdown = node.get_neighbor_node(Direction::x_down)[phase->index], phase_xup = node.get_neighbor_node(Direction::x_up)[phase->index],
								phase_ydown = node.get_neighbor_node(Direction::y_down)[phase->index], phase_yup = node.get_neighbor_node(Direction::y_up)[phase->index],
								phase_zdown = node.get_neighbor_node(Direction::z_down)[phase->index], phase_zup = node.get_neighbor_node(Direction::z_up)[phase->index];
							Vector3 vec((phase_xup.int_increment + phase_xup.bulk_increment - phase_xdown.int_increment - phase_xdown.bulk_increment) / 2.0 / phaseMesh->dr,
								(phase_yup.int_increment + phase_yup.bulk_increment - phase_ydown.int_increment - phase_ydown.bulk_increment) / 2.0 / phaseMesh->dr,
								(phase_zup.int_increment + phase_zup.bulk_increment - phase_zdown.int_increment - phase_zdown.bulk_increment) / 2.0 / phaseMesh->dr);
							grad_dfdphis.add_vec(phase->index, vec);
							if (diff_method == DifferenceMethod::FIVE_POINT) {
								lap_dfdphis.add_double(phase->index,
									(phase_xdown.int_increment + phase_xdown.bulk_increment
										+ phase_xup.int_increment + phase_xup.bulk_increment
										+ phase_ydown.int_increment + phase_ydown.bulk_increment
										+ phase_yup.int_increment + phase_yup.bulk_increment
										+ phase_zdown.int_increment + phase_zdown.bulk_increment
										+ phase_zup.int_increment + phase_zup.bulk_increment
										- 6 * (phase->int_increment + phase->bulk_increment)) / phaseMesh->dr / phaseMesh->dr);
							}
							else if (diff_method == DifferenceMethod::NINE_POINT) {
								PhaseEntry& phase_xdown_ydown = node.get_long_range_node(-1, -1, 0)[phase->index],
									phase_xdown_yup = node.get_long_range_node(-1, 1, 0)[phase->index],
									phase_xup_ydown = node.get_long_range_node(1, -1, 0)[phase->index],
									phase_xup_yup = node.get_long_range_node(1, 1, 0)[phase->index],
									phase_xdown_zdown = node.get_long_range_node(-1, 0, -1)[phase->index],
									phase_xdown_zup = node.get_long_range_node(-1, 0, 1)[phase->index],
									phase_xup_zdown = node.get_long_range_node(1, 0, -1)[phase->index],
									phase_xup_zup = node.get_long_range_node(1, 0, 1)[phase->index],
									phase_ydown_zdown = node.get_long_range_node(0, -1, -1)[phase->index],
									phase_ydown_zup = node.get_long_range_node(0, -1, 1)[phase->index],
									phase_yup_zdown = node.get_long_range_node(0, 1, -1)[phase->index],
									phase_yup_zup = node.get_long_range_node(0, 1, 1)[phase->index];
								lap_dfdphis.add_double(phase->index,
									(4.0 * (phase_xdown.int_increment + phase_xdown.bulk_increment)
										+ 4.0 * (phase_xup.int_increment + phase_xup.bulk_increment)
										+ 4.0 * (phase_ydown.int_increment + phase_ydown.bulk_increment)
										+ 4.0 * (phase_yup.int_increment + phase_yup.bulk_increment)
										+ 4.0 * (phase_zdown.int_increment + phase_zdown.bulk_increment)
										+ 4.0 * (phase_zup.int_increment + phase_zup.bulk_increment)
										+ phase_xdown_ydown.int_increment + phase_xdown_ydown.bulk_increment
										+ phase_xdown_yup.int_increment + phase_xdown_yup.bulk_increment
										+ phase_xup_ydown.int_increment + phase_xup_ydown.bulk_increment
										+ phase_xup_yup.int_increment + phase_xup_yup.bulk_increment
										+ phase_xdown_zdown.int_increment + phase_xdown_zdown.bulk_increment
										+ phase_xdown_zup.int_increment + phase_xdown_zup.bulk_increment
										+ phase_xup_zdown.int_increment + phase_xup_zdown.bulk_increment
										+ phase_xup_zup.int_increment + phase_xup_zup.bulk_increment
										+ phase_ydown_zdown.int_increment + phase_ydown_zdown.bulk_increment
										+ phase_ydown_zup.int_increment + phase_ydown_zup.bulk_increment
										+ phase_yup_zdown.int_increment + phase_yup_zdown.bulk_increment
										+ phase_yup_zup.int_increment + phase_yup_zup.bulk_increment
										- 36 * (phase->int_increment + phase->bulk_increment)) / 6.0 / phaseMesh->dr / phaseMesh->dr);
							}
						}
					}
					// assignment
					for (auto phi_a = node.begin(); phi_a < node.end(); phi_a++) {
						phi_a->old_phi = phi_a->phi;
						if (phi_a->laplacian > SYS_EPSILON || phi_a->laplacian < -SYS_EPSILON) {
							double phi_a_increment = 0.0;
							for (auto phi_b = node.begin(); phi_b < node.end(); phi_b++) {
								if (phi_b->laplacian > SYS_EPSILON || phi_b->laplacian < -SYS_EPSILON) {
									Vector3 vec_Mab((M_ab(node.get_neighbor_node(Direction::x_up), *phi_a, *phi_b) - M_ab(node.get_neighbor_node(Direction::x_down), *phi_a, *phi_b)) / 2.0 / phaseMesh->dr,
										            (M_ab(node.get_neighbor_node(Direction::y_up), *phi_a, *phi_b) - M_ab(node.get_neighbor_node(Direction::y_down), *phi_a, *phi_b)) / 2.0 / phaseMesh->dr,
										            (M_ab(node.get_neighbor_node(Direction::z_up), *phi_a, *phi_b) - M_ab(node.get_neighbor_node(Direction::z_down), *phi_a, *phi_b)) / 2.0 / phaseMesh->dr);
									double Mab = M_ab(node, *phi_a, *phi_b);
									phi_a_increment += vec_Mab * grad_dfdphis[phi_b->index] + Mab * lap_dfdphis[phi_b->index];
#ifdef _DEBUG
									if (_isnan(phi_a_increment)) {
										cout << "DEBUG: phi_a_increment interaction term error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
							// source term
							phi_a_increment += Source_a(node, *phi_a);
#ifdef _DEBUG
							if (_isnan(phi_a_increment)) {
								cout << "DEBUG: phi_a_increment source term error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif

							phi_a->phi += phi_a_increment * dt;
							Boundary_Condition_Phi(node, *phi_a);
							// normalize the phi
							if (is_normalized) {
								if (phi_a->phi > (1.0 - SYS_EPSILON))
									phi_a->phi = 1.0;
								else if (phi_a->phi < SYS_EPSILON)
									phi_a->phi = 0.0;
							}
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(phi_a->old_phi - phi_a->phi) > MAX_PHI_INCREMENT) // MAX_phi_increment
									MAX_PHI_INCREMENT = abs(phi_a->old_phi - phi_a->phi);
							}
#ifdef _DEBUG
							if (_isnan(phi_a->phi)) {
								cout << "DEBUG: phi_a->phi error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
				}
		return MAX_PHI_INCREMENT;
	}

	void CahnHilliardSolver::pre_calculation_phase_concentration(double dt) {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						for (auto con = phase->x.begin(); con < phase->x.end(); con++) {
							con->increment = 0.0;
							con->DiffusionFlux = 0.0;
							con->ChemicalReactionFlux = 0.0;
							con->PhaseTransitionFlux = 0.0;
						}
						for (auto pot = phase->potential.begin(); pot < phase->potential.end(); pot++)
							pot->value = 0.0;
						for (auto m = phase->kinetics_coeff.begin(); m < phase->kinetics_coeff.end(); m++)
							m->value = 0.0;
						init_phase_con_on_moving_interface(node, *phase);
						if (phase->phi > Phi_Num_Cut_Off) {
							dfbulk_dx(node, *phase);
							Mbulk_ij(node, *phase);
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if (phase->phi > Phi_Num_Cut_Off) {
							pf::PhaseEntry& phase_mid = node[phase->index];
							if (phase->laplacian < SYS_EPSILON) {
								if (diff_method == DifferenceMethod::FIVE_POINT) {
									pf::PhaseEntry& phase_upx = node.get_neighbor_node(Direction::x_up)[phase->index];
									pf::PhaseEntry& phase_downx = node.get_neighbor_node(Direction::x_down)[phase->index];
									pf::PhaseEntry& phase_upy = node.get_neighbor_node(Direction::y_up)[phase->index];
									pf::PhaseEntry& phase_downy = node.get_neighbor_node(Direction::y_down)[phase->index];
									pf::PhaseEntry& phase_upz = node.get_neighbor_node(Direction::z_up)[phase->index];
									pf::PhaseEntry& phase_downz = node.get_neighbor_node(Direction::z_down)[phase->index];
									for (auto chem = phase_mid.potential.begin(); chem < phase_mid.potential.end(); chem++) {
										chem->gradient[0] = (phase_upx.potential[chem->index].value - phase_downx.potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[1] = (phase_upy.potential[chem->index].value - phase_downy.potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[2] = (phase_upz.potential[chem->index].value - phase_downz.potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->laplacian = (phase_downx.potential[chem->index].value
											+ phase_upx.potential[chem->index].value
											+ phase_downy.potential[chem->index].value
											+ phase_upy.potential[chem->index].value
											+ phase_downz.potential[chem->index].value
											+ phase_upz.potential[chem->index].value - 6 * chem->value) / phaseMesh->dr / phaseMesh->dr;
									}
									for (auto coef = phase_mid.kinetics_coeff.begin(); coef < phase_mid.kinetics_coeff.end(); coef++) {
										coef->gradient[grad_x] = (phase_upx.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value - phase_downx.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_y] = (phase_upy.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value - phase_downy.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_z] = (phase_upz.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value - phase_downz.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
									}
								}
								else if (diff_method == DifferenceMethod::NINE_POINT) {
									pf::PhaseEntry& phase_upx = node.get_neighbor_node(Direction::x_up)[phase->index];
									pf::PhaseEntry& phase_downx = node.get_neighbor_node(Direction::x_down)[phase->index];
									pf::PhaseEntry& phase_upy = node.get_neighbor_node(Direction::y_up)[phase->index];
									pf::PhaseEntry& phase_downy = node.get_neighbor_node(Direction::y_down)[phase->index];
									pf::PhaseEntry& phase_upz = node.get_neighbor_node(Direction::z_up)[phase->index];
									pf::PhaseEntry& phase_downz = node.get_neighbor_node(Direction::z_down)[phase->index];
									pf::PhaseEntry& phase_upxupy = node.get_long_range_node(1, 1, 0)[phase->index];
									pf::PhaseEntry& phase_downxdowny = node.get_long_range_node(-1, -1, 0)[phase->index];
									pf::PhaseEntry& phase_upydownx = node.get_long_range_node(-1, 1, 0)[phase->index];
									pf::PhaseEntry& phase_downyupx = node.get_long_range_node(1, -1, 0)[phase->index];
									pf::PhaseEntry& phase_upxupz = node.get_long_range_node(1, 0, 1)[phase->index];
									pf::PhaseEntry& phase_downxdownz = node.get_long_range_node(-1, 0, -1)[phase->index];
									pf::PhaseEntry& phase_upzdownx = node.get_long_range_node(-1, 0, 1)[phase->index];
									pf::PhaseEntry& phase_downzupx = node.get_long_range_node(1, 0, -1)[phase->index];
									pf::PhaseEntry& phase_upzupy = node.get_long_range_node(0, 1, 1)[phase->index];
									pf::PhaseEntry& phase_downzdowny = node.get_long_range_node(0, -1, -1)[phase->index];
									pf::PhaseEntry& phase_upydownz = node.get_long_range_node(0, 1, -1)[phase->index];
									pf::PhaseEntry& phase_downyupz = node.get_long_range_node(0, -1, 1)[phase->index];
									for (auto chem = phase_mid.potential.begin(); chem < phase_mid.potential.end(); chem++) {
										chem->gradient[0] = (phase_upx.potential[chem->index].value - phase_downx.potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[1] = (phase_upy.potential[chem->index].value - phase_downy.potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[2] = (phase_upz.potential[chem->index].value - phase_downz.potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->laplacian = (4 * phase_downx.potential[chem->index].value + 4 * phase_upx.potential[chem->index].value
											+ 4 * phase_downy.potential[chem->index].value + 4 * phase_upy.potential[chem->index].value
											+ 4 * phase_downz.potential[chem->index].value + 4 * phase_upz.potential[chem->index].value
											+ phase_upxupy.potential[chem->index].value + phase_downxdowny.potential[chem->index].value
											+ phase_upydownx.potential[chem->index].value + phase_downyupx.potential[chem->index].value
											+ phase_upxupz.potential[chem->index].value + phase_downxdownz.potential[chem->index].value
											+ phase_upzdownx.potential[chem->index].value + phase_downzupx.potential[chem->index].value
											+ phase_upzupy.potential[chem->index].value + phase_downzdowny.potential[chem->index].value
											+ phase_upydownz.potential[chem->index].value + phase_downyupz.potential[chem->index].value - 36 * chem->value) / 6.0 / phaseMesh->dr / phaseMesh->dr;
									}
									for (auto coef = phase_mid.kinetics_coeff.begin(); coef < phase_mid.kinetics_coeff.end(); coef++) {
										coef->gradient[grad_x] = (phase_upx.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value - phase_downx.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_y] = (phase_upy.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value - phase_downy.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_z] = (phase_upz.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value - phase_downz.kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
									}
								}
							}
							else { //the Node in interface
								if (diff_method == DifferenceMethod::FIVE_POINT) {
									pf::PhaseEntry* phase_upx = &node.get_neighbor_node(Direction::x_up)[phase->index];
									pf::PhaseEntry* phase_downx = &node.get_neighbor_node(Direction::x_down)[phase->index];
									pf::PhaseEntry* phase_upy = &node.get_neighbor_node(Direction::y_up)[phase->index];
									pf::PhaseEntry* phase_downy = &node.get_neighbor_node(Direction::y_down)[phase->index];
									pf::PhaseEntry* phase_upz = &node.get_neighbor_node(Direction::z_up)[phase->index];
									pf::PhaseEntry* phase_downz = &node.get_neighbor_node(Direction::z_down)[phase->index];
									if (phase_upx->phi < Phi_Num_Cut_Off)
										phase_upx = &node[phase->index];
									if (phase_downx->phi < Phi_Num_Cut_Off)
										phase_downx = &node[phase->index];
									if (phase_upy->phi < Phi_Num_Cut_Off)
										phase_upy = &node[phase->index];
									if (phase_downy->phi < Phi_Num_Cut_Off)
										phase_downy = &node[phase->index];
									if (phase_upz->phi < Phi_Num_Cut_Off)
										phase_upz = &node[phase->index];
									if (phase_downz->phi < Phi_Num_Cut_Off)
										phase_downz = &node[phase->index];
									for (auto chem = phase_mid.potential.begin(); chem < phase_mid.potential.end(); chem++) {
										chem->gradient[0] = (phase_upx->potential[chem->index].value -
											phase_downx->potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[1] = (phase_upy->potential[chem->index].value -
											phase_downy->potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[2] = (phase_upz->potential[chem->index].value -
											phase_downz->potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->laplacian = (phase_downx->potential[chem->index].value
											+ phase_upx->potential[chem->index].value
											+ phase_downy->potential[chem->index].value
											+ phase_upy->potential[chem->index].value
											+ phase_downz->potential[chem->index].value
											+ phase_upz->potential[chem->index].value - 6 * chem->value) / phaseMesh->dr / phaseMesh->dr;
									}
									for (auto coef = phase_mid.kinetics_coeff.begin(); coef < phase_mid.kinetics_coeff.end(); coef++) {
										coef->gradient[grad_x] = (phase_upx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
											phase_downx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_y] = (phase_upy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
											phase_downy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_z] = (phase_upz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
											phase_downz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
									}
								}
								else if (diff_method == DifferenceMethod::NINE_POINT) {
									pf::PhaseEntry* phase_upx = &node.get_neighbor_node(Direction::x_up)[phase->index];
									pf::PhaseEntry* phase_downx = &node.get_neighbor_node(Direction::x_down)[phase->index];
									pf::PhaseEntry* phase_upy = &node.get_neighbor_node(Direction::y_up)[phase->index];
									pf::PhaseEntry* phase_downy = &node.get_neighbor_node(Direction::y_down)[phase->index];
									pf::PhaseEntry* phase_upz = &node.get_neighbor_node(Direction::z_up)[phase->index];
									pf::PhaseEntry* phase_downz = &node.get_neighbor_node(Direction::z_down)[phase->index];
									if (phase_upx->phi < Phi_Num_Cut_Off)
										phase_upx = &node[phase->index];
									if (phase_downx->phi < Phi_Num_Cut_Off)
										phase_downx = &node[phase->index];
									if (phase_upy->phi < Phi_Num_Cut_Off)
										phase_upy = &node[phase->index];
									if (phase_downy->phi < Phi_Num_Cut_Off)
										phase_downy = &node[phase->index];
									if (phase_upz->phi < Phi_Num_Cut_Off)
										phase_upz = &node[phase->index];
									if (phase_downz->phi < Phi_Num_Cut_Off)
										phase_downz = &node[phase->index];
									pf::PhaseEntry* phase_upxupy = &node.get_long_range_node(1, 1, 0)[phase->index];
									pf::PhaseEntry* phase_downxdowny = &node.get_long_range_node(-1, -1, 0)[phase->index];
									pf::PhaseEntry* phase_upydownx = &node.get_long_range_node(-1, 1, 0)[phase->index];
									pf::PhaseEntry* phase_downyupx = &node.get_long_range_node(1, -1, 0)[phase->index];
									pf::PhaseEntry* phase_upxupz = &node.get_long_range_node(1, 0, 1)[phase->index];
									pf::PhaseEntry* phase_downxdownz = &node.get_long_range_node(-1, 0, -1)[phase->index];
									pf::PhaseEntry* phase_upzdownx = &node.get_long_range_node(-1, 0, 1)[phase->index];
									pf::PhaseEntry* phase_downzupx = &node.get_long_range_node(1, 0, -1)[phase->index];
									pf::PhaseEntry* phase_upzupy = &node.get_long_range_node(0, 1, 1)[phase->index];
									pf::PhaseEntry* phase_downzdowny = &node.get_long_range_node(0, -1, -1)[phase->index];
									pf::PhaseEntry* phase_upydownz = &node.get_long_range_node(0, 1, -1)[phase->index];
									pf::PhaseEntry* phase_downyupz = &node.get_long_range_node(0, -1, 1)[phase->index];
									if (phase_upxupy->phi < Phi_Num_Cut_Off)
										phase_upxupy = &node[phase->index];
									if (phase_downxdowny->phi < Phi_Num_Cut_Off)
										phase_downxdowny = &node[phase->index];
									if (phase_upydownx->phi < Phi_Num_Cut_Off)
										phase_upydownx = &node[phase->index];
									if (phase_downyupx->phi < Phi_Num_Cut_Off)
										phase_downyupx = &node[phase->index];
									if (phase_upxupz->phi < Phi_Num_Cut_Off)
										phase_upxupz = &node[phase->index];
									if (phase_downxdownz->phi < Phi_Num_Cut_Off)
										phase_downxdownz = &node[phase->index];
									if (phase_upzdownx->phi < Phi_Num_Cut_Off)
										phase_upzdownx = &node[phase->index];
									if (phase_downzupx->phi < Phi_Num_Cut_Off)
										phase_downzupx = &node[phase->index];
									if (phase_upzupy->phi < Phi_Num_Cut_Off)
										phase_upzupy = &node[phase->index];
									if (phase_downzdowny->phi < Phi_Num_Cut_Off)
										phase_downzdowny = &node[phase->index];
									if (phase_upydownz->phi < Phi_Num_Cut_Off)
										phase_upydownz = &node[phase->index];
									if (phase_downyupz->phi < Phi_Num_Cut_Off)
										phase_downyupz = &node[phase->index];
									for (auto chem = phase_mid.potential.begin(); chem < phase_mid.potential.end(); chem++) {
										chem->gradient[0] = (phase_upx->potential[chem->index].value -
											phase_downx->potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[1] = (phase_upy->potential[chem->index].value -
											phase_downy->potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->gradient[2] = (phase_upz->potential[chem->index].value -
											phase_downz->potential[chem->index].value) / 2.0 / phaseMesh->dr;
										chem->laplacian = (4 * phase_downx->potential[chem->index].value + 4 * phase_upx->potential[chem->index].value
											+ 4 * phase_downy->potential[chem->index].value + 4 * phase_upy->potential[chem->index].value
											+ 4 * phase_downz->potential[chem->index].value + 4 * phase_upz->potential[chem->index].value
											+ phase_upxupy->potential[chem->index].value + phase_downxdowny->potential[chem->index].value
											+ phase_upydownx->potential[chem->index].value + phase_downyupx->potential[chem->index].value
											+ phase_upxupz->potential[chem->index].value + phase_downxdownz->potential[chem->index].value
											+ phase_upzdownx->potential[chem->index].value + phase_downzupx->potential[chem->index].value
											+ phase_upzupy->potential[chem->index].value + phase_downzdowny->potential[chem->index].value
											+ phase_upydownz->potential[chem->index].value + phase_downyupz->potential[chem->index].value - 36 * chem->value) / 6.0 / phaseMesh->dr / phaseMesh->dr;
									}
									for (auto coef = phase_mid.kinetics_coeff.begin(); coef < phase_mid.kinetics_coeff.end(); coef++) {
										coef->gradient[grad_x] = (phase_upx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
											phase_downx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_y] = (phase_upy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
											phase_downy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
										coef->gradient[grad_z] = (phase_upz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
											phase_downz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
									}
								}
							}

						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					// diffusion term
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if (phase->phi > Phi_Num_Cut_Off) {
							for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
								for (auto comp2 = phase->x.begin(); comp2 < phase->x.end(); comp2++) {
									if (comp->index != solvent && comp2->index != solvent) {
										Vector3 Dij_grad = phase->kinetics_coeff.get_gradientVec3(comp->index, comp2->index);
										double Dij = phase->kinetics_coeff(comp->index, comp2->index).value;
										ChemEntry& chem = phase->potential[comp2->index];
										comp->DiffusionFlux += Dij * (phase->phi_grad * chem.gradient) / phase->phi +
											(Dij_grad * chem.gradient + Dij * chem.laplacian);
#ifdef _DEBUG
										if (_isnan(comp->DiffusionFlux)) {
											cout << "DEBUG: comp->DiffusionFlux error !" << endl;
											SYS_PROGRAM_STOP;
										}
#endif
									}
								}
							}
						}
					}
					// reaction term
					for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
						if (alpha->phi > Phi_Num_Cut_Off)
							for (auto beta = node.begin(); beta < node.end(); beta++) {
								if (beta->phi > Phi_Num_Cut_Off && alpha->index != beta->index) {
									Source_AB(node, *alpha, *beta, abs_grad_phi_AB(node, *alpha, *beta) / alpha->phi);
								}
							}
					}
					for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
						if (alpha->phi > Phi_Num_Cut_Off) {
							Source_A(node, *alpha);
						}
					}
					// phase transtion term
					for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
						if (!alpha->_flag || alpha->phi < Phi_Num_Cut_Off)
							continue;
						for (auto x = alpha->x.begin(); x < alpha->x.end(); x++) {
							x->PhaseTransitionFlux -= x->value / alpha->phi * (alpha->phi - alpha->old_phi) / dt;
#ifdef _DEBUG
							if (_isnan(x->PhaseTransitionFlux)) {
								cout << "DEBUG: x->PhaseTransitionFlux error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
					// summary
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if (phase->phi > Phi_Num_Cut_Off || phase->_flag)
							for (auto x = phase->x.begin(); x < phase->x.end(); x++) {
								x->increment = x->DiffusionFlux + x->ChemicalReactionFlux + x->PhaseTransitionFlux;
#ifdef _DEBUG
								if (_isnan(x->increment)) {
									cout << "DEBUG: x->increment error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
							}
					}
				}
	}

	vector<double> CahnHilliardSolver::solve_phase_concentration(double dt, bool is_normalized) {
		vector<double> MAX_VARIATION; MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0);
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					// assignment
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if (phase->phi > Phi_Num_Cut_Off) {  //this phaseCon need to be calculated
							double sum_x = 0.0;
							for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
								if (comp->index != solvent) {
									double incre_DF = comp->DiffusionFlux * dt;
									double incre_CRF = comp->ChemicalReactionFlux * dt;
									double incre_PTF = comp->PhaseTransitionFlux * dt;
									double old_x = comp->value;
									comp->value += incre_CRF + incre_DF + incre_PTF;
									Boundary_Condition_PhaseX(node, *phase, comp->index);
#ifdef _OPENMP
#pragma omp critical
#endif
									{
										if (abs(incre_DF) > MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX]) // MAX_diffusionFlux
											MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX] = abs(incre_DF);
										if (abs(incre_CRF) > MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX]) // MAX_reactionFlux
											MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX] = abs(incre_CRF);
										if (abs(incre_PTF) > MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX]) // MAX_phaseTransFlux
											MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] = abs(incre_PTF);
										if (abs(old_x - comp->value) > MAX_VARIATION[CH_RETURN::CH_MAX_X_INCREMENT_FLUX]) // MAX_x_increment
											MAX_VARIATION[CH_RETURN::CH_MAX_X_INCREMENT_FLUX] = abs(old_x - comp->value);
									}

									
#ifdef _DEBUG
									if (_isnan(comp->value)) {
										cout << "DEBUG: comp->value error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
									if (is_normalized) {
										if (comp->value < SYS_EPSILON)
											comp->value = SYS_EPSILON;
										else if (comp->value > (1.0 - SYS_EPSILON))
											comp->value = 1.0 - SYS_EPSILON;
									}
									sum_x += comp->value;
								}
							}
							if (solvent != SOLVENT_NONE) {
								if (sum_x > (1.0 - SYS_EPSILON) && is_normalized) {
									phase->x[solvent].value = SYS_EPSILON;
									double scale = (1.0 - SYS_EPSILON) / sum_x;
									for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
										if (comp->index != solvent)
											comp->value *= scale;
								}
								else {
									phase->x[solvent].value = 1.0 - sum_x;
								}
							}
						}
					}
				}
				return MAX_VARIATION;
	}

	void CahnHilliardSolver::pre_calculation_total_concentration_inside_phis(double dt) {
		double MAX_COMP_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.customValues[ExternalFields::CON_Smooth_Phi] = node.cal_phases_fraction_by_index(phase_indexes);
					for (auto con = node.x.begin(); con < node.x.end(); con++) {
						con->increment = 0.0;
						con->DiffusionFlux = 0.0;
						con->ChemicalReactionFlux = 0.0;
						con->PhaseTransitionFlux = 0.0;
					}
					for (auto pot = node.potential.begin(); pot < node.potential.end(); pot++)
						pot->value = 0.0;
					for (auto m = node.kinetics_coeff.begin(); m < node.kinetics_coeff.end(); m++)
						m->value = 0.0;
					init_con_on_moving_interface(node, ConEquationDomain::CEDomain_Standard, threshold);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
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

						for (auto pot = node.potential.begin(); pot < node.potential.end(); pot++)
							pot->value = df_dx(node, pot->index);
						for (auto m = node.kinetics_coeff.begin(); m < node.kinetics_coeff.end(); m++)
							m->value = M_ij(node, m->pairIndex_1, m->pairIndex_2);

					}
					else {
						if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
							for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
								comp->value = 0.0;
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
						if (diff_method == DifferenceMethod::FIVE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downz = &node;
							for (auto chem = node.potential.begin(); chem < node.potential.end(); chem++) {
								chem->gradient[0] = (node_upx->potential[chem->index].value -
									node_downx->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[1] = (node_upy->potential[chem->index].value -
									node_downy->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[2] = (node_upz->potential[chem->index].value -
									node_downz->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->laplacian = (node_downx->potential[chem->index].value
									+ node_upx->potential[chem->index].value
									+ node_downy->potential[chem->index].value
									+ node_upy->potential[chem->index].value
									+ node_downz->potential[chem->index].value
									+ node_upz->potential[chem->index].value - 6 * chem->value) / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto coef = node.kinetics_coeff.begin(); coef < node.kinetics_coeff.end(); coef++) {
								coef->gradient[grad_x] = (node_upx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_y] = (node_upy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_z] = (node_upz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
							}
						}
						else if (diff_method == DifferenceMethod::NINE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downz = &node;
							pf::PhaseNode* node_upxupy = &node.get_long_range_node(1, 1, 0);
							pf::PhaseNode* node_downxdowny = &node.get_long_range_node(-1, -1, 0);
							pf::PhaseNode* node_upydownx = &node.get_long_range_node(-1, 1, 0);
							pf::PhaseNode* node_downyupx = &node.get_long_range_node(1, -1, 0);
							pf::PhaseNode* node_upxupz = &node.get_long_range_node(1, 0, 1);
							pf::PhaseNode* node_downxdownz = &node.get_long_range_node(-1, 0, -1);
							pf::PhaseNode* node_upzdownx = &node.get_long_range_node(-1, 0, 1);
							pf::PhaseNode* node_downzupx = &node.get_long_range_node(1, 0, -1);
							pf::PhaseNode* node_upzupy = &node.get_long_range_node(0, 1, 1);
							pf::PhaseNode* node_downzdowny = &node.get_long_range_node(0, -1, -1);
							pf::PhaseNode* node_upydownz = &node.get_long_range_node(0, 1, -1);
							pf::PhaseNode* node_downyupz = &node.get_long_range_node(0, -1, 1);
							if (node_upxupy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upxupy = &node;
							if (node_downxdowny->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downxdowny = &node;
							if (node_upydownx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upydownx = &node;
							if (node_downyupx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downyupx = &node;
							if (node_upxupz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upxupz = &node;
							if (node_downxdownz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downxdownz = &node;
							if (node_upzdownx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upzdownx = &node;
							if (node_downzupx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downzupx = &node;
							if (node_upzupy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upzupy = &node;
							if (node_downzdowny->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downzdowny = &node;
							if (node_upydownz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upydownz = &node;
							if (node_downyupz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downyupz = &node;
							for (auto chem = node.potential.begin(); chem < node.potential.end(); chem++) {
								chem->gradient[0] = (node_upx->potential[chem->index].value -
									node_downx->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[1] = (node_upy->potential[chem->index].value -
									node_downy->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[2] = (node_upz->potential[chem->index].value -
									node_downz->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->laplacian = (4 * node_downx->potential[chem->index].value + 4 * node_upx->potential[chem->index].value
									+ 4 * node_downy->potential[chem->index].value + 4 * node_upy->potential[chem->index].value
									+ 4 * node_downz->potential[chem->index].value + 4 * node_upz->potential[chem->index].value
									+ node_upxupy->potential[chem->index].value + node_downxdowny->potential[chem->index].value
									+ node_upydownx->potential[chem->index].value + node_downyupx->potential[chem->index].value
									+ node_upxupz->potential[chem->index].value + node_downxdownz->potential[chem->index].value
									+ node_upzdownx->potential[chem->index].value + node_downzupx->potential[chem->index].value
									+ node_upzupy->potential[chem->index].value + node_downzdowny->potential[chem->index].value
									+ node_upydownz->potential[chem->index].value + node_downyupz->potential[chem->index].value - 36 * chem->value) / 6.0 / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto coef = node.kinetics_coeff.begin(); coef < node.kinetics_coeff.end(); coef++) {
								coef->gradient[grad_x] = (node_upx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_y] = (node_upy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_z] = (node_upz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
							}
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
						// diffusion term
						Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh->dr);
						for (auto comp = node.x.begin(); comp < node.x.end(); comp++) {
							for (auto comp2 = node.x.begin(); comp2 < node.x.end(); comp2++) {
								if (comp->index != solvent && comp2->index != solvent) {
									Vector3 D_KineticCoef = node.kinetics_coeff.get_gradientVec3(comp->index, comp2->index);
									double kinetic_coef = node.kinetics_coeff(comp->index, comp2->index).value;
									ChemEntry& diff_c = node.potential[comp2->index];
									comp->DiffusionFlux += kinetic_coef * (s_phi_grad * diff_c.gradient) / node.customValues[ExternalFields::CON_Smooth_Phi] +
										(D_KineticCoef * diff_c.gradient + kinetic_coef * diff_c.laplacian);
#ifdef _DEBUG
									if (_isnan(comp->DiffusionFlux)) {
										cout << "DEBUG: comp->DiffusionFlux error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
						}
						// reaction term
						double s_phi_grad_norm = s_phi_grad.abs();
						if (s_phi_grad_norm > SYS_EPSILON)
							for (auto x = node.x.begin(); x < node.x.end(); x++)
								x->ChemicalReactionFlux = s_phi_grad_norm / node.customValues[ExternalFields::CON_Smooth_Phi] * int_flux(node, x->index);

						for (auto x = node.x.begin(); x < node.x.end(); x++) {
							x->ChemicalReactionFlux += Source(node, x->index);
						}
						// phase transtion term
						for (auto x = node.x.begin(); x < node.x.end(); x++) {
							x->PhaseTransitionFlux -= x->value / node.customValues[ExternalFields::CON_Smooth_Phi] 
								* (node.customValues[ExternalFields::CON_Smooth_Phi] - node.customValues[ExternalFields::CON_Smooth_Old_Phi]) / dt;
						}
						// summary
						for (auto x = node.x.begin(); x < node.x.end(); x++) {
							x->increment = x->DiffusionFlux + x->ChemicalReactionFlux + x->PhaseTransitionFlux;
#ifdef _DEBUG
							if (_isnan(x->increment)) {
								cout << "DEBUG: x->increment error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
				}
	}

	vector<double> CahnHilliardSolver::solve_total_concentration_inside_phis(double dt, bool is_normalized) {
		vector<double> MAX_VARIATION; MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0);
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					// assignment
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
						double sum_x = 0.0;
						for (auto comp = node.x.begin(); comp < node.x.end(); comp++) {
							if (comp->index != solvent) {
								double incre_DF = comp->DiffusionFlux * dt;
								double incre_CRF = comp->ChemicalReactionFlux * dt;
								double incre_PTF = comp->PhaseTransitionFlux * dt;
								double old_x = comp->value;
								comp->value += incre_CRF + incre_DF + incre_PTF;
								Boundary_Condition_TotalX(node, comp->index);
#ifdef _OPENMP
#pragma omp critical
#endif
								{
									if (abs(incre_DF) > MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX]) // MAX_diffusionFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX] = abs(incre_DF);
									if (abs(incre_CRF) > MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX]) // MAX_reactionFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX] = abs(incre_CRF);
									if (abs(incre_PTF) > MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX]) // MAX_phaseTransFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] = abs(incre_PTF);
									if (abs(old_x - comp->value) > MAX_VARIATION[CH_RETURN::CH_MAX_X_INCREMENT_FLUX]) // MAX_x_increment
										MAX_VARIATION[CH_RETURN::CH_MAX_X_INCREMENT_FLUX] = abs(old_x - comp->value);
								}
#ifdef _DEBUG
								if (_isnan(comp->value)) {
									cout << "DEBUG: comp->value error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
								if (is_normalized) {
									if (comp->value < SYS_EPSILON)
										comp->value = SYS_EPSILON;
									else if (comp->value > (1.0 - SYS_EPSILON))
										comp->value = 1.0 - SYS_EPSILON;
								}
								sum_x += comp->value;
							}
						}
						if (solvent != SOLVENT_NONE) {
							if (sum_x > (1.0 - SYS_EPSILON) && is_normalized) {
								node.x[solvent].value = SYS_EPSILON;
								double scale = (1.0 - SYS_EPSILON) / sum_x;
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (comp->index != solvent)
										comp->value *= scale;
							}
							else {
								node.x[solvent].value = 1.0 - sum_x;
							}
						}
					}
					node.customValues[ExternalFields::CON_Smooth_Old_Phi] = node.customValues[ExternalFields::CON_Smooth_Phi];
				}
		return MAX_VARIATION;
	}

	void CahnHilliardSolver::pre_calculation_total_concentration_outside_phis(double dt) {
		double MAX_COMP_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.customValues[ExternalFields::CON_Smooth_Phi] = node.cal_phases_fraction_by_index(phase_indexes);
					for (auto con = node.x.begin(); con < node.x.end(); con++) {
						con->increment = 0.0;
						con->DiffusionFlux = 0.0;
						con->ChemicalReactionFlux = 0.0;
						con->PhaseTransitionFlux = 0.0;
					}
					for (auto pot = node.potential.begin(); pot < node.potential.end(); pot++)
						pot->value = 0.0;
					for (auto m = node.kinetics_coeff.begin(); m < node.kinetics_coeff.end(); m++)
						m->value = 0.0;
					init_con_on_moving_interface(node, ConEquationDomain::CEDomain_Reverse, threshold);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
							PhaseNode& upx_node = node.get_neighbor_node(Direction::x_up);
							PhaseNode& downx_node = node.get_neighbor_node(Direction::x_down);
							PhaseNode& upy_node = node.get_neighbor_node(Direction::y_up);
							PhaseNode& downy_node = node.get_neighbor_node(Direction::y_down);
							PhaseNode& upz_node = node.get_neighbor_node(Direction::z_up);
							PhaseNode& downz_node = node.get_neighbor_node(Direction::z_down);
							if (upx_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
								for (auto c = node.x.begin(); c < node.x.end(); c++)
									c->value = upx_node.x[c->index].value;
							}
							else if (downx_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downx_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
								for (auto c = node.x.begin(); c < node.x.end(); c++)
									c->value = downx_node.x[c->index].value;
							}
							else if (upy_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
								for (auto c = node.x.begin(); c < node.x.end(); c++)
									c->value = upy_node.x[c->index].value;
							}
							else if (downy_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downy_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
								for (auto c = node.x.begin(); c < node.x.end(); c++)
									c->value = downy_node.x[c->index].value;
							}
							else if (upz_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && upz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
								for (auto c = node.x.begin(); c < node.x.end(); c++)
									c->value = upz_node.x[c->index].value;
							}
							else if (downz_node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && downz_node.customValues[ExternalFields::CON_Smooth_Old_Phi] > threshold) {
								for (auto c = node.x.begin(); c < node.x.end(); c++)
									c->value = downz_node.x[c->index].value;
							}
						}

						for (auto pot = node.potential.begin(); pot < node.potential.end(); pot++)
							pot->value = df_dx(node, pot->index);
						for (auto m = node.kinetics_coeff.begin(); m < node.kinetics_coeff.end(); m++)
							m->value = M_ij(node, m->pairIndex_1, m->pairIndex_2);

					}
					else {
						if (node.customValues[ExternalFields::CON_Smooth_Old_Phi] < threshold) {
							for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
								comp->value = 0.0;
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						if (diff_method == DifferenceMethod::FIVE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downz = &node;
							for (auto chem = node.potential.begin(); chem < node.potential.end(); chem++) {
								chem->gradient[0] = (node_upx->potential[chem->index].value -
									node_downx->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[1] = (node_upy->potential[chem->index].value -
									node_downy->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[2] = (node_upz->potential[chem->index].value -
									node_downz->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->laplacian = (node_downx->potential[chem->index].value
									+ node_upx->potential[chem->index].value
									+ node_downy->potential[chem->index].value
									+ node_upy->potential[chem->index].value
									+ node_downz->potential[chem->index].value
									+ node_upz->potential[chem->index].value - 6 * chem->value) / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto coef = node.kinetics_coeff.begin(); coef < node.kinetics_coeff.end(); coef++) {
								coef->gradient[grad_x] = (node_upx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_y] = (node_upy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_z] = (node_upz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
							}
						}
						else if (diff_method == DifferenceMethod::NINE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downz = &node;
							pf::PhaseNode* node_upxupy = &node.get_long_range_node(1, 1, 0);
							pf::PhaseNode* node_downxdowny = &node.get_long_range_node(-1, -1, 0);
							pf::PhaseNode* node_upydownx = &node.get_long_range_node(-1, 1, 0);
							pf::PhaseNode* node_downyupx = &node.get_long_range_node(1, -1, 0);
							pf::PhaseNode* node_upxupz = &node.get_long_range_node(1, 0, 1);
							pf::PhaseNode* node_downxdownz = &node.get_long_range_node(-1, 0, -1);
							pf::PhaseNode* node_upzdownx = &node.get_long_range_node(-1, 0, 1);
							pf::PhaseNode* node_downzupx = &node.get_long_range_node(1, 0, -1);
							pf::PhaseNode* node_upzupy = &node.get_long_range_node(0, 1, 1);
							pf::PhaseNode* node_downzdowny = &node.get_long_range_node(0, -1, -1);
							pf::PhaseNode* node_upydownz = &node.get_long_range_node(0, 1, -1);
							pf::PhaseNode* node_downyupz = &node.get_long_range_node(0, -1, 1);
							if (node_upxupy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upxupy = &node;
							if (node_downxdowny->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downxdowny = &node;
							if (node_upydownx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upydownx = &node;
							if (node_downyupx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downyupx = &node;
							if (node_upxupz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upxupz = &node;
							if (node_downxdownz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downxdownz = &node;
							if (node_upzdownx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upzdownx = &node;
							if (node_downzupx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downzupx = &node;
							if (node_upzupy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upzupy = &node;
							if (node_downzdowny->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downzdowny = &node;
							if (node_upydownz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upydownz = &node;
							if (node_downyupz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downyupz = &node;
							for (auto chem = node.potential.begin(); chem < node.potential.end(); chem++) {
								chem->gradient[0] = (node_upx->potential[chem->index].value -
									node_downx->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[1] = (node_upy->potential[chem->index].value -
									node_downy->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->gradient[2] = (node_upz->potential[chem->index].value -
									node_downz->potential[chem->index].value) / 2.0 / phaseMesh->dr;
								chem->laplacian = (4 * node_downx->potential[chem->index].value + 4 * node_upx->potential[chem->index].value
									+ 4 * node_downy->potential[chem->index].value + 4 * node_upy->potential[chem->index].value
									+ 4 * node_downz->potential[chem->index].value + 4 * node_upz->potential[chem->index].value
									+ node_upxupy->potential[chem->index].value + node_downxdowny->potential[chem->index].value
									+ node_upydownx->potential[chem->index].value + node_downyupx->potential[chem->index].value
									+ node_upxupz->potential[chem->index].value + node_downxdownz->potential[chem->index].value
									+ node_upzdownx->potential[chem->index].value + node_downzupx->potential[chem->index].value
									+ node_upzupy->potential[chem->index].value + node_downzdowny->potential[chem->index].value
									+ node_upydownz->potential[chem->index].value + node_downyupz->potential[chem->index].value - 36 * chem->value) / 6.0 / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto coef = node.kinetics_coeff.begin(); coef < node.kinetics_coeff.end(); coef++) {
								coef->gradient[grad_x] = (node_upx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downx->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_y] = (node_upy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downy->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
								coef->gradient[grad_z] = (node_upz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value -
									node_downz->kinetics_coeff(coef->pairIndex_1, coef->pairIndex_2).value) / 2.0 / phaseMesh->dr;
							}
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					//double s_phi = node.cal_phases_fraction(phase_indexes);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						// diffusion term
						Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh->dr) * (-1.0);
						for (auto comp = node.x.begin(); comp < node.x.end(); comp++) {
							for (auto comp2 = node.x.begin(); comp2 < node.x.end(); comp2++) {
								if (comp->index != solvent && comp2->index != solvent) {
									Vector3 D_KineticCoef = node.kinetics_coeff.get_gradientVec3(comp->index, comp2->index);
									double kinetic_coef = node.kinetics_coeff(comp->index, comp2->index).value;
									ChemEntry& diff_c = node.potential[comp2->index];
									comp->DiffusionFlux += kinetic_coef * (s_phi_grad * diff_c.gradient) / (1.0 - node.customValues[ExternalFields::CON_Smooth_Phi]) +
										(D_KineticCoef * diff_c.gradient + kinetic_coef * diff_c.laplacian);
#ifdef _DEBUG
									if (_isnan(comp->DiffusionFlux)) {
										cout << "DEBUG: comp->DiffusionFlux error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
						}
						// reaction term
						double smooth_grad_norm = s_phi_grad.abs();
						if (smooth_grad_norm > SYS_EPSILON)
							for (auto x = node.x.begin(); x < node.x.end(); x++)
								x->ChemicalReactionFlux = smooth_grad_norm / node.customValues[ExternalFields::CON_Smooth_Phi] * int_flux(node, x->index);
						for (auto x = node.x.begin(); x < node.x.end(); x++) {
							x->ChemicalReactionFlux += Source(node, x->index);
						}
						// phase transtion term
						for (auto x = node.x.begin(); x < node.x.end(); x++) {
							x->PhaseTransitionFlux -= x->value / node.customValues[ExternalFields::CON_Smooth_Phi]
								* (node.customValues[ExternalFields::CON_Smooth_Phi] - node.customValues[ExternalFields::CON_Smooth_Old_Phi]) / dt;
						}
						// summary
						for (auto x = node.x.begin(); x < node.x.end(); x++) {
							x->increment = x->DiffusionFlux + x->ChemicalReactionFlux + x->PhaseTransitionFlux;
#ifdef _DEBUG
							if (_isnan(x->increment)) {
								cout << "DEBUG: x->increment error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
				}
	}

	vector<double> CahnHilliardSolver::solve_total_concentration_outside_phis(double dt, bool is_normalized) {
		vector<double> MAX_VARIATION; MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0);
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						// assignment
						double sum_x = 0.0;
						for (auto comp = node.x.begin(); comp < node.x.end(); comp++) {
							if (comp->index != solvent) {
								double incre_DF = comp->DiffusionFlux * dt;
								double incre_CRF = comp->ChemicalReactionFlux * dt;
								double incre_PTF = comp->PhaseTransitionFlux * dt;
								double old_x = comp->value;
								comp->value += incre_CRF + incre_DF + incre_PTF;
								Boundary_Condition_TotalX(node, comp->index);
#ifdef _OPENMP
#pragma omp critical
#endif
								{
									if (abs(incre_DF) > MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX]) // MAX_diffusionFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX] = abs(incre_DF);
									if (abs(incre_CRF) > MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX]) // MAX_reactionFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX] = abs(incre_CRF);
									if (abs(incre_PTF) > MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX]) // MAX_phaseTransFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] = abs(incre_PTF);
									if (abs(old_x - comp->value) > MAX_VARIATION[CH_RETURN::CH_MAX_X_INCREMENT_FLUX]) // MAX_x_increment
										MAX_VARIATION[CH_RETURN::CH_MAX_X_INCREMENT_FLUX] = abs(old_x - comp->value);
								}

#ifdef _DEBUG
								if (_isnan(comp->value)) {
									cout << "DEBUG: comp->value error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
								if (is_normalized) {
									if (comp->value < SYS_EPSILON)
										comp->value = SYS_EPSILON;
									else if (comp->value > (1.0 - SYS_EPSILON))
										comp->value = 1.0 - SYS_EPSILON;
								}
								sum_x += comp->value;
							}
						}
						if (solvent != SOLVENT_NONE) {
							if (sum_x > (1.0 - SYS_EPSILON) && is_normalized) {
								node.x[solvent].value = SYS_EPSILON;
								double scale = (1.0 - SYS_EPSILON) / sum_x;
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (comp->index != solvent)
										comp->value *= scale;
							}
							else {
								node.x[solvent].value = 1.0 - sum_x;
							}
						}
					}
					node.customValues[ExternalFields::CON_Smooth_Old_Phi] = node.customValues[ExternalFields::CON_Smooth_Phi];
				}
		return MAX_VARIATION;
	}

	vector<double> CahnHilliardSolver::pre_calculation_grand_potential_functional_inside_phis(double dt) {
		vector<double> MAX_VARIATION; MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0);
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.customValues[ExternalFields::CON_Smooth_Phi] = node.cal_phases_fraction_by_index(phase_indexes);
					for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
						comp->value = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
							comp->value = 0.0;
						phase_x(node, *phase);
						for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
							node.x[comp->index].value += phase->phi * comp->value;
							phase->potential.add_con(comp->index, node.potential[comp->index].value);
						}
					}
					for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
						p->gradient.set_to_zero();
						p->increment = 0.0;
						p->laplacian = 0.0;
						for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++)
							node.kinetics_coeff.set(p->index, p2->index, M_ij(node, p->index, p2->index));
					}
					init_grand_potential_on_moving_interface(node, ConEquationDomain::CEDomain_Standard, threshold);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
						if (diff_method == DifferenceMethod::FIVE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downz = &node;
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								p->gradient[0] = (node_upx->potential[p->index].value -
									node_downx->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[1] = (node_upy->potential[p->index].value -
									node_downy->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[2] = (node_upz->potential[p->index].value -
									node_downz->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->laplacian = (node_downx->potential[p->index].value
									+ node_upx->potential[p->index].value
									+ node_downy->potential[p->index].value
									+ node_upy->potential[p->index].value
									+ node_downz->potential[p->index].value
									+ node_upz->potential[p->index].value - 6 * p->value) / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
									node.kinetics_coeff(p->index, p2->index).gradient[grad_x] = (node_upx->kinetics_coeff(p->index, p2->index).value -
										node_downx->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_y] = (node_upy->kinetics_coeff(p->index, p2->index).value -
										node_downy->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_z] = (node_upz->kinetics_coeff(p->index, p2->index).value -
										node_downz->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
								}
							}
						}
						else if (diff_method == DifferenceMethod::NINE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downz = &node;
							pf::PhaseNode* node_upxupy = &node.get_long_range_node(1, 1, 0);
							pf::PhaseNode* node_downxdowny = &node.get_long_range_node(-1, -1, 0);
							pf::PhaseNode* node_upydownx = &node.get_long_range_node(-1, 1, 0);
							pf::PhaseNode* node_downyupx = &node.get_long_range_node(1, -1, 0);
							pf::PhaseNode* node_upxupz = &node.get_long_range_node(1, 0, 1);
							pf::PhaseNode* node_downxdownz = &node.get_long_range_node(-1, 0, -1);
							pf::PhaseNode* node_upzdownx = &node.get_long_range_node(-1, 0, 1);
							pf::PhaseNode* node_downzupx = &node.get_long_range_node(1, 0, -1);
							pf::PhaseNode* node_upzupy = &node.get_long_range_node(0, 1, 1);
							pf::PhaseNode* node_downzdowny = &node.get_long_range_node(0, -1, -1);
							pf::PhaseNode* node_upydownz = &node.get_long_range_node(0, 1, -1);
							pf::PhaseNode* node_downyupz = &node.get_long_range_node(0, -1, 1);
							if (node_upxupy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upxupy = &node;
							if (node_downxdowny->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downxdowny = &node;
							if (node_upydownx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upydownx = &node;
							if (node_downyupx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downyupx = &node;
							if (node_upxupz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upxupz = &node;
							if (node_downxdownz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downxdownz = &node;
							if (node_upzdownx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upzdownx = &node;
							if (node_downzupx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downzupx = &node;
							if (node_upzupy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upzupy = &node;
							if (node_downzdowny->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downzdowny = &node;
							if (node_upydownz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_upydownz = &node;
							if (node_downyupz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
								node_downyupz = &node;

							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								p->gradient[0] = (node_upx->potential[p->index].value -
									node_downx->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[1] = (node_upy->potential[p->index].value -
									node_downy->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[2] = (node_upz->potential[p->index].value -
									node_downz->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->laplacian = (node_downx->potential[p->index].value
									+ node_upx->potential[p->index].value
									+ node_downy->potential[p->index].value
									+ node_upy->potential[p->index].value
									+ node_downz->potential[p->index].value
									+ node_upz->potential[p->index].value - 6 * p->value) / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
									node.kinetics_coeff(p->index, p2->index).gradient[grad_x] = (node_upx->kinetics_coeff(p->index, p2->index).value -
										node_downx->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_y] = (node_upy->kinetics_coeff(p->index, p2->index).value -
										node_downy->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_z] = (node_upz->kinetics_coeff(p->index, p2->index).value -
										node_downz->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
								}
							}
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
						Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh->dr);
						double s_phi_grad_norm = s_phi_grad.abs();
						for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
							double DIFF_FLUX = 0.0, REAC_FLUX = 0.0, PHASETRANS_FLUX = 0.0;
							// diffusion term
							for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
								Vector3 D_KineticCoef = node.kinetics_coeff.get_gradientVec3(p->index, p2->index);
								double kinetic_coef = node.kinetics_coeff(p->index, p2->index).value;

								DIFF_FLUX += kinetic_coef * (s_phi_grad * p2->gradient)
									+ (D_KineticCoef * p2->gradient + kinetic_coef * p2->laplacian) * node.customValues[ExternalFields::CON_Smooth_Phi];

							}
#ifdef _DEBUG
							if (_isnan(DIFF_FLUX)) {
								cout << "DEBUG: DIFF_FLUX error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							// reaction term
							if (s_phi_grad_norm > SYS_EPSILON)
								REAC_FLUX = s_phi_grad_norm * int_flux(node, p->index);
							REAC_FLUX += Source(node, p->index) * node.customValues[ExternalFields::CON_Smooth_Phi];

							// phase transtion term
							for (auto phase = node.begin(); phase < node.end(); phase++)
								PHASETRANS_FLUX -= (phase->phi - phase->old_phi) / dt * phase->x[p->index].value;
#ifdef _DEBUG
							if (_isnan(PHASETRANS_FLUX)) {
								cout << "DEBUG: PHASETRANS_FLUX error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							// summary
							p->increment = DIFF_FLUX + REAC_FLUX + PHASETRANS_FLUX;
#ifdef _DEBUG
							if (_isnan(p->increment)) {
								cout << "DEBUG: node.potential[x].increment error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							double pre_factor_con_i = 0.0;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								pre_factor_con_i += phase->phi * dphase_x_du(node, *phase, p->index);
#ifdef _DEBUG
							if (Is_Equality(pre_factor_con_i, 0.0)) {
								cout << "DEBUG: pre_factor_con_i error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							
							p->increment /= pre_factor_con_i;

#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(DIFF_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX]) // MAX_diffusionFlux
									MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX] = abs(DIFF_FLUX * dt);
								if (abs(REAC_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX]) // MAX_reactionFlux
									MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX] = abs(REAC_FLUX * dt);
								if (abs(PHASETRANS_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX]) // MAX_phaseTransFlux
									MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] = abs(PHASETRANS_FLUX * dt);
							}

						}

					}
				}
		return MAX_VARIATION;
	}
	
	double CahnHilliardSolver::solve_grand_potential_functional_inside_phis(double dt) {
		double MAX_VARIATION = 0.0; 
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
						for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
							double old_p = p->value;
							p->value += p->increment * dt;
							Boundary_Condition_TotalX(node, p->index);
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(old_p - p->value) > MAX_VARIATION) // MAX_x_increment
									MAX_VARIATION = abs(old_p - p->value);
							}
#ifdef _DEBUG
							if (_isnan(p->value)) {
								cout << "DEBUG: node.potential[x].value error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
					node.customValues[ExternalFields::CON_Smooth_Old_Phi] = node.customValues[ExternalFields::CON_Smooth_Phi];
				}
		return MAX_VARIATION;
	}

	vector<double> CahnHilliardSolver::pre_calculation_grand_potential_functional_outside_phis(double dt) {
		vector<double> MAX_VARIATION; MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0);
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.customValues[ExternalFields::CON_Smooth_Phi] = node.cal_phases_fraction_by_index(phase_indexes);
					for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
						comp->value = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
							comp->value = 0.0;
						phase_x(node, *phase);
						for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
							node.x[comp->index].value += phase->phi * comp->value;
							phase->potential.add_con(comp->index, node.potential[comp->index].value);
						}
					}
					for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
						p->gradient.set_to_zero();
						p->increment = 0.0;
						p->laplacian = 0.0;
						for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++)
							node.kinetics_coeff.set(p->index, p2->index, M_ij(node, p->index, p2->index));
					}
					init_grand_potential_on_moving_interface(node, ConEquationDomain::CEDomain_Reverse, threshold);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						if (diff_method == DifferenceMethod::FIVE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downz = &node;
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								p->gradient[0] = (node_upx->potential[p->index].value -
									node_downx->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[1] = (node_upy->potential[p->index].value -
									node_downy->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[2] = (node_upz->potential[p->index].value -
									node_downz->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->laplacian = (node_downx->potential[p->index].value
									+ node_upx->potential[p->index].value
									+ node_downy->potential[p->index].value
									+ node_upy->potential[p->index].value
									+ node_downz->potential[p->index].value
									+ node_upz->potential[p->index].value - 6 * p->value) / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
									node.kinetics_coeff(p->index, p2->index).gradient[grad_x] = (node_upx->kinetics_coeff(p->index, p2->index).value -
										node_downx->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_y] = (node_upy->kinetics_coeff(p->index, p2->index).value -
										node_downy->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_z] = (node_upz->kinetics_coeff(p->index, p2->index).value -
										node_downz->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
								}
							}
						}
						else if (diff_method == DifferenceMethod::NINE_POINT) {
							pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
							pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
							pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
							pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
							pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
							pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
							if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upx = &node;
							if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downx = &node;
							if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upy = &node;
							if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downy = &node;
							if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upz = &node;
							if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downz = &node;
							pf::PhaseNode* node_upxupy = &node.get_long_range_node(1, 1, 0);
							pf::PhaseNode* node_downxdowny = &node.get_long_range_node(-1, -1, 0);
							pf::PhaseNode* node_upydownx = &node.get_long_range_node(-1, 1, 0);
							pf::PhaseNode* node_downyupx = &node.get_long_range_node(1, -1, 0);
							pf::PhaseNode* node_upxupz = &node.get_long_range_node(1, 0, 1);
							pf::PhaseNode* node_downxdownz = &node.get_long_range_node(-1, 0, -1);
							pf::PhaseNode* node_upzdownx = &node.get_long_range_node(-1, 0, 1);
							pf::PhaseNode* node_downzupx = &node.get_long_range_node(1, 0, -1);
							pf::PhaseNode* node_upzupy = &node.get_long_range_node(0, 1, 1);
							pf::PhaseNode* node_downzdowny = &node.get_long_range_node(0, -1, -1);
							pf::PhaseNode* node_upydownz = &node.get_long_range_node(0, 1, -1);
							pf::PhaseNode* node_downyupz = &node.get_long_range_node(0, -1, 1);
							if (node_upxupy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upxupy = &node;
							if (node_downxdowny->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downxdowny = &node;
							if (node_upydownx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upydownx = &node;
							if (node_downyupx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downyupx = &node;
							if (node_upxupz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upxupz = &node;
							if (node_downxdownz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downxdownz = &node;
							if (node_upzdownx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upzdownx = &node;
							if (node_downzupx->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downzupx = &node;
							if (node_upzupy->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upzupy = &node;
							if (node_downzdowny->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downzdowny = &node;
							if (node_upydownz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_upydownz = &node;
							if (node_downyupz->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
								node_downyupz = &node;

							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								p->gradient[0] = (node_upx->potential[p->index].value -
									node_downx->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[1] = (node_upy->potential[p->index].value -
									node_downy->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->gradient[2] = (node_upz->potential[p->index].value -
									node_downz->potential[p->index].value) / 2.0 / phaseMesh->dr;
								p->laplacian = (node_downx->potential[p->index].value
									+ node_upx->potential[p->index].value
									+ node_downy->potential[p->index].value
									+ node_upy->potential[p->index].value
									+ node_downz->potential[p->index].value
									+ node_upz->potential[p->index].value - 6 * p->value) / phaseMesh->dr / phaseMesh->dr;
							}
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
									node.kinetics_coeff(p->index, p2->index).gradient[grad_x] = (node_upx->kinetics_coeff(p->index, p2->index).value -
										node_downx->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_y] = (node_upy->kinetics_coeff(p->index, p2->index).value -
										node_downy->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									node.kinetics_coeff(p->index, p2->index).gradient[grad_z] = (node_upz->kinetics_coeff(p->index, p2->index).value -
										node_downz->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
								}
							}
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					double DIFF_FLUX = 0.0, REAC_FLUX = 0.0, PHASETRANS_FLUX = 0.0;
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh->dr) * (-1.0);
						double s_phi_grad_norm = s_phi_grad.abs();
						for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
							double DIFF_FLUX = 0.0, REAC_FLUX = 0.0, PHASETRANS_FLUX = 0.0;
							// diffusion term
							for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
								Vector3 D_KineticCoef = node.kinetics_coeff.get_gradientVec3(p->index, p2->index);
								double kinetic_coef = node.kinetics_coeff(p->index, p2->index).value;

								DIFF_FLUX += kinetic_coef * (s_phi_grad * p2->gradient)
									+ (D_KineticCoef * p2->gradient + kinetic_coef * p2->laplacian) * (1.0 - node.customValues[ExternalFields::CON_Smooth_Phi]);

							}
#ifdef _DEBUG
							if (_isnan(DIFF_FLUX)) {
								cout << "DEBUG: DIFF_FLUX error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							// reaction term
							if (s_phi_grad_norm > SYS_EPSILON)
								REAC_FLUX = s_phi_grad_norm * int_flux(node, p->index);
							REAC_FLUX += Source(node, p->index) * (1.0 - node.customValues[ExternalFields::CON_Smooth_Phi]);

							// phase transtion term
							for (auto phase = node.begin(); phase < node.end(); phase++) {
									PHASETRANS_FLUX -= (phase->phi - phase->old_phi) / dt * phase->x[p->index].value;
							}
#ifdef _DEBUG
							if (_isnan(PHASETRANS_FLUX) || REAC_FLUX > 100) {
								cout << "DEBUG: PHASETRANS_FLUX error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							// summary
							p->increment = DIFF_FLUX + REAC_FLUX + PHASETRANS_FLUX;
#ifdef _DEBUG
							if (_isnan(p->increment)) {
								cout << "DEBUG: node.potential[x].increment error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							double pre_factor_con_i = 0.0;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								pre_factor_con_i += phase->phi * dphase_x_du(node, *phase, p->index);
#ifdef _DEBUG
							if (Is_Equality(pre_factor_con_i, 0.0)) {
								cout << "DEBUG: pre_factor_con_i error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
							p->increment /= pre_factor_con_i;

#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(DIFF_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX]) // MAX_diffusionFlux
									MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX] = abs(DIFF_FLUX * dt);
								if (abs(REAC_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX]) // MAX_reactionFlux
									MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX] = abs(REAC_FLUX * dt);
								if (abs(PHASETRANS_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX]) // MAX_phaseTransFlux
									MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] = abs(PHASETRANS_FLUX * dt);
							}


						}
					}
				}
		return MAX_VARIATION;
	}

	double CahnHilliardSolver::solve_grand_potential_functional_outside_phis(double dt) {
		double MAX_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
						for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
							double old_p = p->value;
							p->value += p->increment * dt;
							if (p->value < grand_potential_range[0])
								p->value = grand_potential_range[0];
							else if (p->value > grand_potential_range[1])
								p->value = grand_potential_range[1];
							Boundary_Condition_TotalX(node, p->index);
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(old_p - p->value) > MAX_VARIATION) // MAX_x_increment
									MAX_VARIATION = abs(old_p - p->value);
							}
#ifdef _DEBUG
							if (_isnan(p->value)) {
								cout << "DEBUG: node.potential[x].value error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
					node.customValues[ExternalFields::CON_Smooth_Old_Phi] = node.customValues[ExternalFields::CON_Smooth_Phi];
				}
		return MAX_VARIATION;
	}

	void CahnHilliardSolver::summary_phix_to_x() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto x = node.x.begin(); x < node.x.end(); x++)
						x->value = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						for (auto x = phase->x.begin(); x < phase->x.end(); x++)
							node.x[x->index].value += x->value * phase->phi;
				}
	}

	void CahnHilliardSolver::summary_phip_to_p() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto p = node.potential.begin(); p < node.potential.end(); p++)
						p->value = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						for (auto p = phase->potential.begin(); p < phase->potential.end(); p++)
							node.potential[p->index].value += p->value * phase->phi;
				}
	}
};