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


#include"Solver_custom.h"
namespace pf {
	double CustomAllenCahnSolver::solve_one_step(double dt, bool adjust_phi_0_1) {
		double MAX_PHI_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (int grain = 0; grain < grains_number; grain++)
						node.customValues[-(grain + grains_start_index)] = -L(node, grain, grains_start_index) * dF_dphi(node, grain, grains_start_index) + Source(node, grain, grains_start_index);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (int grain = 0; grain < grains_number; grain++) {
						double old_val = node.customValues[(grain + grains_start_index)];
						node.customValues[(grain + grains_start_index)] += node.customValues[-(grain + grains_start_index)] * dt;
						if (adjust_phi_0_1) {
							if (node.customValues[(grain + grains_start_index)] >= 1.0)
								node.customValues[(grain + grains_start_index)] = 1.0;
							else if (node.customValues[(grain + grains_start_index)] <= 0.0)
								node.customValues[(grain + grains_start_index)] = 0.0;
						}
						BoundaryCondition(node, grain, grains_start_index);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (abs(old_val - node.customValues[(grain + grains_start_index)]) > MAX_PHI_VARIATION)
								MAX_PHI_VARIATION = abs(old_val - node.customValues[(grain + grains_start_index)]);
						}
					}
				}
		return MAX_PHI_VARIATION;
	}

	double CustomCahnHilliardSolver::solve_one_step(double dt, DifferenceMethod diff_method, bool adjust_phi_0_1) {
		double dr = phaseMesh->dr, MAX_COMP_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (int comp = 0; comp < components_number; comp++) {
						node.customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)] = Mobility(node, comp, comps_start_index);
						node.customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)] = dF_dc(node, comp, comps_start_index);
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (int comp = 0; comp < components_number; comp++) {
						Vector3 vec_m((node.get_neighbor_node(Direction::x_down).customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]
							- node.get_neighbor_node(Direction::x_up).customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]) / 2.0 / dr,
							(node.get_neighbor_node(Direction::y_down).customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]
								- node.get_neighbor_node(Direction::y_up).customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]) / 2.0 / dr,
							(node.get_neighbor_node(Direction::z_down).customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]
								- node.get_neighbor_node(Direction::z_up).customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]) / 2.0 / dr)
							, vec_df_dc((node.get_neighbor_node(Direction::x_down).customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]
								- node.get_neighbor_node(Direction::x_up).customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]) / 2.0 / dr,
								(node.get_neighbor_node(Direction::y_down).customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]
									- node.get_neighbor_node(Direction::y_up).customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]) / 2.0 / dr,
								(node.get_neighbor_node(Direction::z_down).customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]
									- node.get_neighbor_node(Direction::z_up).customValues[(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)]) / 2.0 / dr);
						node.customValues[-(comp + comps_start_index)] = vec_m * vec_df_dc
							+ node.customValues[-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD)] * node.cal_customValues_laplace((comp + BUFF_FOR_SOLVER_CAHN_HILLIARD), dr, diff_method)
							+ Source(node, comp, comps_start_index);
					}
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					for (int comp = 0; comp < components_number; comp++) {
						double old_val = node.customValues[(comp + comps_start_index)];
						node.customValues[(comp + comps_start_index)] += node.customValues[-(comp + comps_start_index)] * dt;
						if (adjust_phi_0_1) {
							if (node.customValues[(comp + comps_start_index)] >= 1.0)
								node.customValues[(comp + comps_start_index)] = 1.0;
							else if (node.customValues[(comp + comps_start_index)] <= 0.0)
								node.customValues[(comp + comps_start_index)] = 0.0;
						}
						BoundaryCondition(node, comp, comps_start_index);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (MAX_COMP_VARIATION < abs(old_val - node.customValues[(comp + comps_start_index)]))
								MAX_COMP_VARIATION = abs(old_val - node.customValues[(comp + comps_start_index)]);
						}
					}
				}
		return MAX_COMP_VARIATION;
	}

	void FourierTransformSolver::init_basic_real_space() {
#pragma omp parallel for
		for (int i = 0; i < mirrorMesh.Nx; i++) {
			int ii = 0;
			for (int j = 0; j < mirrorMesh.Ny; j++)
				for (int k = 0; k < mirrorMesh.Nz; k++) {
					PhaseNode& node = get_phaseNode_by_mirrorMesh_XYZ(i, j, k);
					init_real_space(rlDomain[0][k + mirrorMesh.Nz * (j + mirrorMesh.Ny * i)], node, i, j, k);
				}
		}
	}

	void FourierTransformSolver::solve_one_step() {
		// fill real space
#pragma omp parallel for
		for (int i = 0; i < mirrorMesh.Nx; i++)
			for (int j = 0; j < mirrorMesh.Ny; j++)
				for (int k = 0; k < mirrorMesh.Nz; k++) {
					PhaseNode& node = get_phaseNode_by_mirrorMesh_XYZ(i, j, k);
					int XYZ = k + mirrorMesh.Nz * (j + mirrorMesh.Ny * i);
					vector<std::complex<double>> real_space;
					real_space.resize(rlc_space_size, { 0.0, 0.0 });
					fill_node_real_space(real_space, node, i, j, k);
					for (int index = 1; index < rlc_space_size; index++)
						rlDomain[index][XYZ]._Val[FFTW_REAL] = real_space[index]._Val[FFTW_REAL];
				}
		// transform real space to fourier space
#pragma omp parallel sections// OMP BEGIN
		{
			for (int index = 0; index < rlc_space_size; index++) {
#pragma omp section
				{
					fftw_execute(ForwardPlan[index]);
				}
			}
		}
		// calculation in fourier space
#pragma omp parallel for
		for (int i = 0; i < mirrorMesh.Nx; i++)
			for (int j = 0; j < mirrorMesh.Ny; j++)
				for (int k = 0; k < mirrorMesh.Nz; k++) {
					PhaseNode& node = get_phaseNode_by_mirrorMesh_XYZ(i, j, k);
					int XYZ = k + mirrorMesh.Nz * (j + mirrorMesh.Ny * i);
					vector<std::complex<double>> fourier_space;
					for (int index = 0; index < rlc_space_size; index++)
						fourier_space.push_back(rcDomain[index][XYZ]);
					double q2 = Q2(XYZ);
					rcDomain[0][XYZ] = dynamic_equation_fourier_space(fourier_space, node, q2, q2 * q2);
				}
		// transform fourier space to real space
		fftw_execute(BackwardPlan);
		// boundary condition
#pragma omp parallel for
		for (int i = 0; i < mirrorMesh.Nx; i++)
			for (int j = 0; j < mirrorMesh.Ny; j++)
				for (int k = 0; k < mirrorMesh.Nz; k++) {
					int XYZ = k + mirrorMesh.Nz * (j + mirrorMesh.Ny * i);
					std::complex<double>& basic_real_space = rlDomain[0][k + mirrorMesh.Nz * (j + mirrorMesh.Ny * i)];
					PhaseNode& node = get_phaseNode_by_mirrorMesh_XYZ(i, j, k);

					basic_real_space = basic_real_space * Norm;

					boundary_condition_real_space(basic_real_space, node, i, j, k);
				}
	}

}