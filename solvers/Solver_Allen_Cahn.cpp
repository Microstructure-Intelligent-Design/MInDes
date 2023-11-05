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


#include"Solver_Allen_Cahn.h"
namespace pf {
	void AllenCahnSolver::init_phi(bool normalize_phi) {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if(normalize_phi)
						node.normalized_phi();
					for (auto phase = node.begin(); phase < node.end(); phase++)
						phase->old_phi = phase->phi;
				}
	}

	void AllenCahnSolver::init_phi_pair_wise(bool normalize_phi) {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					if (normalize_phi)
						node.normalized_phi();
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						phase->old_phi = phase->phi;
						phase->_flag = phaseMesh->currentFlag(node, phase->index);
					}
				}
	}

	void AllenCahnSolver::pre_calculation_phi_pair_wise_normal(double dt, bool normalize_phi, DifferenceMethod diff_method) {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						//intphase and near intphase
						if (phase->_flag) {
							phase->bulk_increment = 0.0;
							phase->int_increment = 0.0;
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
						else {
							phase->phi_grad[0] = 0.0;
							phase->phi_grad[1] = 0.0;
							phase->phi_grad[2] = 0.0;
							phase->laplacian = 0.0;
						}
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					dfint_dphi(node, normalize_phi);
					// equation //
					for (auto alpha = node.begin(); alpha < node.end() - 1; alpha++)
						for (auto beta = alpha + 1; beta < node.end(); beta++) {
							if (alpha->_flag == pf_INTERFACE && beta->_flag == pf_INTERFACE) {
								double bulk_increment_b_a = Lij(node, *alpha, *beta) * (dfbulk_dphi(node, *beta) - dfbulk_dphi(node, *alpha)) + Source_ij(node, *alpha, *beta);
								if (normalize_phi) {
									if ((bulk_increment_b_a > SYS_EPSILON && (alpha->phi > (1.0 - SYS_EPSILON) || beta->phi < SYS_EPSILON))
										|| (bulk_increment_b_a < -SYS_EPSILON && (alpha->phi < SYS_EPSILON || beta->phi >(1.0 - SYS_EPSILON))))
										bulk_increment_b_a = 0.0;
								}
								alpha->bulk_increment += bulk_increment_b_a;
								beta->bulk_increment -= bulk_increment_b_a;
							}
						}
					// numerical treatment
					if (normalize_phi) {
						double scale = 1.0, increment = 0.0;
						for (auto phase = node.begin(); phase < node.end(); phase++)
							if (phase->_flag) {
								increment = dt * (phase->int_increment + phase->bulk_increment);
								if (Is_Equality(increment, 0.0))
									continue;
								else if ((phase->phi + increment) > 1.0) {
									double p_scale = (1.0 - phase->phi) / increment;
									if (p_scale < scale)
										scale = p_scale;
								}
								else if ((phase->phi + increment) < 0.0) {
									double p_scale = abs(phase->phi / increment);
									if (p_scale < scale)
										scale = p_scale;
								}
							}
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							phase->int_increment *= scale;
							phase->bulk_increment *= scale;
						}
					}
				}
	}

	double AllenCahnSolver::solve_phi_pair_wise_normal(double dt, bool normalize_phi) {
		double MAX_PHI_INCREMENT = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						phase->old_phi = phase->phi;
						if (phase->_flag) {
							phase->phi += dt * (phase->int_increment + phase->bulk_increment);
							Boundary_Condition(node, *phase);
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(phase->old_phi - phase->phi) > MAX_PHI_INCREMENT)
									MAX_PHI_INCREMENT = abs(phase->old_phi - phase->phi);
							}
#ifdef _DEBUG
							if (_isnan(phase->phi)) {
								cout << "DEBUG: phase->phi error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
					// normalize the phi
					if (normalize_phi) {
						node.normalized_phi();
					}
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						// change _flag
						if (phase->phi >= Phi_Num_Cut_Off && phase->phi <= (1.0 - Phi_Num_Cut_Off)) {
							if (phase->_flag != pf_INTERFACE) {
								phase->_flag = pf_INTERFACE;
								for (int direction = 0; direction < 6; direction++) {
									PhaseNode& check_node = node.get_neighbor_node(Direction(direction));
									for (auto check_phase = check_node.begin(); check_phase < check_node.end(); check_phase++)
										if (check_phase->index == phase->index && check_phase->_flag == pf_BULK)
											check_phase->_flag = pf_NEAR_INTERFACE;
								}
							}
						}
						if (phase->phi < Phi_Num_Cut_Off) {
							bool is_near = false;
							for (int direction = 0; direction < 6; direction++) {
								PhaseNode& check_node = node.get_neighbor_node(Direction(direction));
								for (auto check_phase = check_node.begin(); check_phase < check_node.end(); check_phase++)
									if (check_phase->index == phase->index && check_phase->phi >= Phi_Num_Cut_Off)
										is_near = true;
							}
							if(is_near)
								phase->_flag = pf_NEAR_INTERFACE;
							else
								phase->_flag = pf_BULK;
						}
						else if (phase->phi > (1.0 - Phi_Num_Cut_Off)) {
							bool is_near = false;
							for (int direction = 0; direction < 6; direction++) {
								PhaseNode& check_node = node.get_neighbor_node(Direction(direction));
								for (auto check_phase = check_node.begin(); check_phase < check_node.end(); check_phase++)
									if (check_phase->index == phase->index && check_phase->phi <= (1.0 - Phi_Num_Cut_Off))
										is_near = true;
							}
							if (is_near)
								phase->_flag = pf_NEAR_INTERFACE;
							else
								phase->_flag = pf_BULK;
						}
					}
				}
		return MAX_PHI_INCREMENT;
	}

	void AllenCahnSolver::pre_calculation_phi_normal(double dt, bool normalize_phi, DifferenceMethod diff_method) {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						phase->bulk_increment = 0.0;
						phase->int_increment = 0.0;
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
					dfint_dphi(node, normalize_phi);
					// equation //
					for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
						alpha->bulk_increment = Source_i(node, *alpha);
						if (alpha->laplacian > SYS_EPSILON || alpha->laplacian < -SYS_EPSILON)
							for (auto beta = node.begin(); beta < node.end(); beta++)
								if (beta->laplacian > SYS_EPSILON || beta->laplacian < -SYS_EPSILON)
									alpha->bulk_increment += -Lij(node, *alpha, *beta) * dfbulk_dphi(node, *beta);
					}
					/*if (normalize_phi) {
						double scale = 1.0, increment = 0.0;
						for (auto phase = node.begin(); phase < node.end(); phase++)
							if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
								increment = dt * (phase->int_increment + phase->bulk_increment);
								if (Is_Equality(increment, 0.0))
									continue;
								else if ((phase->phi + increment) > 1.0) {
									double p_scale = (1.0 - phase->phi) / increment;
									if (p_scale < scale)
										scale = p_scale;
								}
								else if ((phase->phi + increment) < 0.0) {
									double p_scale = abs(phase->phi / increment);
									if (p_scale < scale)
										scale = p_scale;
								}
							}
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							phase->int_increment *= scale;
							phase->bulk_increment *= scale;
						}
					}*/
				}

	}

	double AllenCahnSolver::solve_phi_normal(double dt, bool normalize_phi) {
		double MAX_PHI_INCREMENT = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						phase->old_phi = phase->phi;
						if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
							phase->phi += dt * (phase->int_increment + phase->bulk_increment);
							Boundary_Condition(node, *phase);
							// normalize the phi
							if (normalize_phi) {
								if (phase->phi > (1.0 - SYS_EPSILON))
									phase->phi = 1.0;
								else if (phase->phi < SYS_EPSILON)
									phase->phi = 0.0;
							}
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(phase->old_phi - phase->phi) > MAX_PHI_INCREMENT)
									MAX_PHI_INCREMENT = abs(phase->old_phi - phase->phi);
							}
#ifdef _DEBUG
							if (_isnan(phase->phi)) {
								cout << "DEBUG: phase->phi error !" << endl;
								SYS_PROGRAM_STOP;
							}
#endif
						}
					}
				}
		return MAX_PHI_INCREMENT;
	}
}