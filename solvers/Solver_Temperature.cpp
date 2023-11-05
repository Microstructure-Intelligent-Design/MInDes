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


#include"Solver_Temperature.h"
namespace pf {

	void TemperatureSolver::pre_calculation_temperature(DifferenceMethod diff_method) {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					node.temperature.D = D_temp(node);
					node.temperature.laplace = node.cal_laplace_temperature(phaseMesh->dr, diff_method);
				}
	}

	double TemperatureSolver::solve_temperature(double dt) {
		double MAX_T_INCREMENT = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					Vector3 T_grand = node.cal_temperature_gradient(phaseMesh->dr);
					Vector3 D_grand;
					D_grand[0] = (node.get_neighbor_node(Direction::x_up).temperature.D - node.get_neighbor_node(Direction::x_down).temperature.D) / 2.0 / phaseMesh->dr;
					D_grand[1] = (node.get_neighbor_node(Direction::y_up).temperature.D - node.get_neighbor_node(Direction::y_down).temperature.D) / 2.0 / phaseMesh->dr;
					D_grand[2] = (node.get_neighbor_node(Direction::z_up).temperature.D - node.get_neighbor_node(Direction::z_down).temperature.D) / 2.0 / phaseMesh->dr;
					node.temperature.increment = T_grand * D_grand + node.temperature.D * node.temperature.laplace + Source(node);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					double old_T = node.temperature.T;
					node.temperature.T += dt * node.temperature.increment;
					BoundaryCondition(node);
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (abs(node.temperature.T - old_T) > MAX_T_INCREMENT)
							MAX_T_INCREMENT = abs(node.temperature.T - old_T);
					}
#ifdef _DEBUG
					if (_isnan(node.temperature.T)) {
						cout << "DEBUG: node.temperature error !" << endl;
						SYS_PROGRAM_STOP;
					}
#endif
				}
		return MAX_T_INCREMENT;
	}
}