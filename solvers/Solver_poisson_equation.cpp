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


#include"Solver_poisson_equation.h"
namespace pf {
	void PoissonEquationSolver_Explicit::set_RHS_value() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++)
					setRHS_value((*phaseMesh)(x, y, z), RHS_INDEX);
	}
	void PoissonEquationSolver_Explicit::set_LHS_value() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++)
					setLHS_value((*phaseMesh)(x, y, z), LHS_INDEX);
	}
	int PoissonEquationSolver_Explicit::solve_whole_domain(double accuracy, int MAXIterations, bool debug_solver, int output_step) {
		double dx = phaseMesh->dr, MAXVariation = 0.0;
		int iteration_times = 0;
		set_RHS_value();
		set_LHS_value();
		do {
			iteration_times++;
			MAXVariation = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						double old_R_value = node.customValues[R_INDEX];
						Vector3 vec3_R = node.cal_customValues_gradient(R_INDEX, dx), vec3_LHS = node.cal_customValues_gradient(LHS_INDEX, dx);
						node.customValues[R_INDEX] = (
							node.get_neighbor_node(Direction::x_down).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::x_up).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::y_down).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::y_up).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::z_down).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::z_up).customValues[R_INDEX] 
							- dx * dx * (node.customValues[RHS_INDEX] - vec3_R * vec3_LHS)) / 6.0;
						boundary(node, R_INDEX);
						double variation = std::abs(old_R_value - node.customValues[R_INDEX]);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (variation > MAXVariation)
								MAXVariation = variation;
						}
					}
			if ((iteration_times % output_step == 0) && debug_solver) {
				cout << solver_name << " iterate " << iteration_times << " times !" << endl;
				cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
			}
		} while (MAXVariation > accuracy && iteration_times < MAXIterations);
		if (debug_solver) {
			cout << solver_name << " iterate " << iteration_times << " times !" << endl;
			cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
		}
		return iteration_times;
	}
	int PoissonEquationSolver_Explicit::solve_whole_domain_with_average_boundary_condition(double accuracy, int MAXIterations, double average_value, bool debug_solver, int output_step) {
		double dx = phaseMesh->dr, MAXVariation = 0.0, current_average = 0.0;
		int iteration_times = 0;
		set_RHS_value();
		set_LHS_value();
		do {
			iteration_times++;
			MAXVariation = 0.0;
			double age_increment = 0.0;
			if (iteration_times > 1)
				age_increment = average_value - current_average;
			current_average = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						double old_R_value = node.customValues[R_INDEX];
						Vector3 vec3_R = node.cal_customValues_gradient(R_INDEX, dx), vec3_LHS = node.cal_customValues_gradient(LHS_INDEX, dx);
						node.customValues[R_INDEX] = (node.get_neighbor_node(Direction::x_down).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::x_up).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::y_down).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::y_up).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::z_down).customValues[R_INDEX]
							+ node.get_neighbor_node(Direction::z_up).customValues[R_INDEX] - dx * dx * (node.customValues[RHS_INDEX] - vec3_R * vec3_LHS)) / 6.0 + age_increment;
						boundary(node, R_INDEX);
						double variation = std::abs(old_R_value - node.customValues[R_INDEX]);
						current_average += node.customValues[R_INDEX];
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (variation > MAXVariation)
								MAXVariation = variation;
						}
					}
			current_average /= phaseMesh->limit_x * phaseMesh->limit_y * phaseMesh->limit_z;
			if ((iteration_times % output_step == 0) && debug_solver) {
				cout << solver_name << " iterate " << iteration_times << " times !" << endl;
				cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
			}
		} while (MAXVariation > accuracy && iteration_times < MAXIterations);
		if (debug_solver) {
			cout << solver_name << " iterate " << iteration_times << " times !" << endl;
			cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
		}
		return iteration_times;
	}
	int PoissonEquationSolver_Explicit::solve_inside_phases_region(double accuracy, int MAXIterations, vector<int> phaseIndexes, double phi_threshold, bool debug_solver, int output_step) {
		double dx = phaseMesh->dr, MAXVariation = 0.0;
		int iteration_times = 0, Smooth_Phi_index = poison_protect_index + LHS_INDEX;
		set_RHS_value();
		set_LHS_value();
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.customValues.add_double(Smooth_Phi_index, node.cal_phases_fraction_by_index(phaseIndexes));
				}
		do {
			iteration_times++;
			MAXVariation = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						double old_R_value = node.customValues[R_INDEX];
						Vector3 vec3_R = node.cal_customValues_gradient(R_INDEX, dx), vec3_LHS = node.cal_customValues_gradient(LHS_INDEX, dx);
						double r_x_down = old_R_value, r_x_up = old_R_value, r_y_down = old_R_value, r_y_up = old_R_value,
							r_z_down = old_R_value, r_z_up = old_R_value;

						if (node.get_neighbor_node(Direction::x_down).customValues[Smooth_Phi_index] > phi_threshold)
							r_x_down = node.get_neighbor_node(Direction::x_down).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::x_up).customValues[Smooth_Phi_index] > phi_threshold)
							r_x_up = node.get_neighbor_node(Direction::x_up).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::y_down).customValues[Smooth_Phi_index] > phi_threshold)
							r_y_down = node.get_neighbor_node(Direction::y_down).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::y_up).customValues[Smooth_Phi_index] > phi_threshold)
							r_y_up = node.get_neighbor_node(Direction::y_up).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::z_down).customValues[Smooth_Phi_index] > phi_threshold)
							r_z_down = node.get_neighbor_node(Direction::z_down).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::z_up).customValues[Smooth_Phi_index] > phi_threshold)
							r_z_up = node.get_neighbor_node(Direction::z_up).customValues[R_INDEX];

						node.customValues[R_INDEX] = (r_x_down + r_x_up + r_y_down + r_y_up + r_z_down + r_z_up
							- dx * dx * (node.customValues[RHS_INDEX] - vec3_R * vec3_LHS)) / 6.0;
						boundary(node, R_INDEX);
						double variation = std::abs(old_R_value - node.customValues[R_INDEX]);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (variation > MAXVariation)
								MAXVariation = variation;
						}
					}
			if ((iteration_times % output_step == 0) && debug_solver) {
				cout << solver_name << " iterate " << iteration_times << " times !" << endl;
				cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
			}
		} while (MAXVariation > accuracy && iteration_times < MAXIterations);
		if (debug_solver) {
			cout << solver_name << " iterate " << iteration_times << " times !" << endl;
			cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
		}
		return iteration_times;
	}
	int PoissonEquationSolver_Explicit::solve_outside_phases_region(double accuracy, int MAXIterations, vector<int> phaseIndexes, double phi_threshold, bool debug_solver, int output_step) {
		double dx = phaseMesh->dr, MAXVariation = 0.0;
		int iteration_times = 0, Smooth_Phi_index = poison_protect_index + LHS_INDEX;
		set_RHS_value();
		set_LHS_value();
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.customValues.add_double(Smooth_Phi_index, node.cal_phases_fraction_by_index(phaseIndexes));
				}
		do {
			iteration_times++;
			MAXVariation = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						double old_R_value = node.customValues[R_INDEX];
						Vector3 vec3_R = node.cal_customValues_gradient(R_INDEX, dx), vec3_LHS = node.cal_customValues_gradient(LHS_INDEX, dx);
						double r_x_down = old_R_value, r_x_up = old_R_value, r_y_down = old_R_value, r_y_up = old_R_value,
							r_z_down = old_R_value, r_z_up = old_R_value;
						
						if (node.get_neighbor_node(Direction::x_down).customValues[Smooth_Phi_index] < phi_threshold)
							r_x_down = node.get_neighbor_node(Direction::x_down).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::x_up).customValues[Smooth_Phi_index] < phi_threshold)
							r_x_up = node.get_neighbor_node(Direction::x_up).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::y_down).customValues[Smooth_Phi_index] < phi_threshold)
							r_y_down = node.get_neighbor_node(Direction::y_down).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::y_up).customValues[Smooth_Phi_index] < phi_threshold)
							r_y_up = node.get_neighbor_node(Direction::y_up).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::z_down).customValues[Smooth_Phi_index] < phi_threshold)
							r_z_down = node.get_neighbor_node(Direction::z_down).customValues[R_INDEX];
						if (node.get_neighbor_node(Direction::z_up).customValues[Smooth_Phi_index] < phi_threshold)
							r_z_up = node.get_neighbor_node(Direction::z_up).customValues[R_INDEX];

						node.customValues[R_INDEX] = (r_x_down + r_x_up + r_y_down + r_y_up + r_z_down + r_z_up
							- dx * dx * (node.customValues[RHS_INDEX] - vec3_R * vec3_LHS)) / 6.0;
						boundary(node, R_INDEX);
						double variation = std::abs(old_R_value - node.customValues[R_INDEX]);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (variation > MAXVariation)
								MAXVariation = variation;
						}
					}
			if ((iteration_times % output_step == 0) && debug_solver) {
				cout << solver_name << " iterate " << iteration_times << " times !" << endl;
				cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
			}
		} while (MAXVariation > accuracy && iteration_times < MAXIterations);
		if (debug_solver) {
			cout << solver_name << " iterate " << iteration_times << " times !" << endl;
			cout << solver_name << " MAXVariation = " << MAXVariation << " !" << endl;
		}
		return iteration_times;
	}
}