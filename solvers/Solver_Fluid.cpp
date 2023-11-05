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


#include"Solver_Fluid.h"
using namespace std;

namespace pf {

	void FluidField::set_u_in_velocity_field(int x, int y, int z, double u) {
		U(x, y, z).vals[0] = u;
	}
	void FluidField::set_v_in_velocity_field(int x, int y, int z, double u) {
		V(x, y, z).vals[0] = u;
	}
	void FluidField::set_w_in_velocity_field(int x, int y, int z, double u) {
		W(x, y, z).vals[0] = u;
	}

	void FluidField::cal_parameters_before_calculation() {
		// boundary condition
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					cal_parameters_for_main_domain(node, phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z);
					boundary_condition_for_pressure(node, phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z);
				}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					boundary_condition_for_domain_U(u, node.get_neighbor_node(pf::Direction::x_down), node, U.limit_x, U.limit_y, U.limit_z);
					boundary_condition_for_domain_V(v, node.get_neighbor_node(pf::Direction::y_down), node, V.limit_x, V.limit_y, V.limit_z);
					boundary_condition_for_domain_W(w, node.get_neighbor_node(pf::Direction::z_down), node, W.limit_x, W.limit_y, W.limit_z);
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					boundary_condition_for_domain_U(u, node, node.get_neighbor_node(pf::Direction::x_up), U.limit_x, U.limit_y, U.limit_z);
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					boundary_condition_for_domain_V(v, node, node.get_neighbor_node(pf::Direction::y_up), V.limit_x, V.limit_y, V.limit_z);
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					boundary_condition_for_domain_W(w, node, node.get_neighbor_node(pf::Direction::z_up), W.limit_x, W.limit_y, W.limit_z);
				}
		}
	}

	Vector3 FluidField::evolve_momentum_equation(double dt) {
		Vector3 dVelocity(0.0, 0.0, 0.0);
		double dr = phaseMesh->dr;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					// momentum equation for u
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::x_down).customValues[ExternalFields::FLUID_mass]) / 2.0
							, viscosity = (node.customValues[ExternalFields::FLUID_viscosity] + node.get_neighbor_node(Direction::x_down).customValues[ExternalFields::FLUID_viscosity]) / 2.0;
						double	u_xdown = u.get_neighbor_node(Direction::x_down).vals[0], u_xup = u.get_neighbor_node(Direction::x_up).vals[0],
							u_ydown = u.get_neighbor_node(Direction::y_down).vals[0], u_yup = u.get_neighbor_node(Direction::y_up).vals[0],
							u_zdown = u.get_neighbor_node(Direction::z_down).vals[0], u_zup = u.get_neighbor_node(Direction::z_up).vals[0],
							viscosity_xdown = node.get_neighbor_node(Direction::x_down).customValues[ExternalFields::FLUID_viscosity],
							viscosity_xup = node.customValues[ExternalFields::FLUID_viscosity],
							viscosity_ydown = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, -1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, -1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_yup = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, 1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, 1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zdown = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, 0, -1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, 0, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zup = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, 0, 1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, 0, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							u_laplace = (u_xdown + u_xup + u_ydown + u_yup + u_zdown + u_zup - 6.0 * u.vals[0]) / dr / dr,
							v_ave_down = (v.get_long_range_node(0, 0, 0).vals[0] + v.get_long_range_node(-1, 0, 0).vals[0]) / 2.0,
							v_ave_up = (v.get_long_range_node(0, 1, 0).vals[0] + v.get_long_range_node(-1, 1, 0).vals[0]) / 2.0,
							w_ave_down = (w.get_long_range_node(0, 0, 0).vals[0] + w.get_long_range_node(-1, 0, 0).vals[0]) / 2.0,
							w_ave_up = (w.get_long_range_node(0, 0, 1).vals[0] + w.get_long_range_node(-1, 0, 1).vals[0]) / 2.0;
						Vector3 vec_viscosity, vec_u;
						vec_viscosity[0] = (viscosity_xdown - viscosity_xup) / dr;
						vec_viscosity[1] = (viscosity_ydown - viscosity_yup) / dr;
						vec_viscosity[2] = (viscosity_zdown - viscosity_zup) / dr;
						vec_u[0] = (u_xdown - u_xup) / dr / 2.0;
						vec_u[1] = (u_ydown - u_yup) / dr / 2.0;
						vec_u[2] = (u_zdown - u_zup) / dr / 2.0;
						double A = (viscosity * u_laplace + vec_viscosity * vec_u) / density - ((u_xdown * u_xdown - u_xup * u_xup) + (u_ydown * v_ave_down - u_yup * v_ave_up)
							+ (u_zdown * w_ave_down - u_zup * w_ave_up)) / 2.0 / dr;
						u.vals[1] = A + (node.get_long_range_node(-1, 0, 0).customVec3s[ExternalFields::FLUID_volume_force][0] + node.customVec3s[ExternalFields::FLUID_volume_force][0]) / 2.0;
					}
					// momentum equation for v
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::y_down).customValues[ExternalFields::FLUID_mass]) / 2.0
							, viscosity = (node.customValues[ExternalFields::FLUID_viscosity] + node.get_neighbor_node(Direction::y_down).customValues[ExternalFields::FLUID_viscosity]) / 2.0;
						double v_xdown = v.get_neighbor_node(Direction::x_down).vals[0], v_xup = v.get_neighbor_node(Direction::x_up).vals[0],
							v_ydown = v.get_neighbor_node(Direction::y_down).vals[0], v_yup = v.get_neighbor_node(Direction::y_up).vals[0],
							v_zdown = v.get_neighbor_node(Direction::z_down).vals[0], v_zup = v.get_neighbor_node(Direction::z_up).vals[0],
							viscosity_ydown = node.get_neighbor_node(Direction::y_down).customValues[ExternalFields::FLUID_viscosity],
							viscosity_yup = node.customValues[ExternalFields::FLUID_viscosity],
							viscosity_xdown = (viscosity_ydown + viscosity_yup + node.get_long_range_node(-1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, -1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_xup = (viscosity_ydown + viscosity_yup + node.get_long_range_node(1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, -1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zdown = (viscosity_ydown + viscosity_yup + node.get_long_range_node(0, 0, -1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, -1, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zup = (viscosity_ydown + viscosity_yup + node.get_long_range_node(0, 0, 1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, -1, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							v_laplace = (v_xdown + v_xup + v_ydown + v_yup + v_zdown + v_zup - 6.0 * v.vals[0]) / dr / dr,
							u_ave_down = (u.get_long_range_node(0, 0, 0).vals[0] + u.get_long_range_node(0, -1, 0).vals[0]) / 2.0,
							u_ave_up = (u.get_long_range_node(1, 0, 0).vals[0] + u.get_long_range_node(1, -1, 0).vals[0]) / 2.0,
							w_ave_down = (w.get_long_range_node(0, 0, 0).vals[0] + w.get_long_range_node(0, -1, 0).vals[0]) / 2.0,
							w_ave_up = (w.get_long_range_node(0, 0, 1).vals[0] + w.get_long_range_node(0, -1, 1).vals[0]) / 2.0;
						Vector3 vec_viscosity, vec_v;
						vec_viscosity[0] = (viscosity_xdown - viscosity_xup) / dr;
						vec_viscosity[1] = (viscosity_ydown - viscosity_yup) / dr;
						vec_viscosity[2] = (viscosity_zdown - viscosity_zup) / dr;
						vec_v[0] = (v_xdown - v_xup) / dr / 2.0;
						vec_v[1] = (v_ydown - v_yup) / dr / 2.0;
						vec_v[2] = (v_zdown - v_zup) / dr / 2.0;
						double B = (viscosity * v_laplace + vec_viscosity * vec_v) / density - ((v_xdown * u_ave_down - v_xup * u_ave_up) + (v_ydown * v_ydown - v_yup * v_yup)
							+ (v_zdown * w_ave_down - v_zup * w_ave_up)) / 2.0 / dr;
						v.vals[1] = B + (node.get_long_range_node(0, -1, 0).customVec3s[ExternalFields::FLUID_volume_force][0] + node.customVec3s[ExternalFields::FLUID_volume_force][0]) / 2.0;
					}
					// momentum equation for w
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::z_down).customValues[ExternalFields::FLUID_mass]) / 2.0
							, viscosity = (node.customValues[ExternalFields::FLUID_viscosity] + node.get_neighbor_node(Direction::z_down).customValues[ExternalFields::FLUID_viscosity]) / 2.0;
						double w_xdown = w.get_neighbor_node(Direction::x_down).vals[0], w_xup = w.get_neighbor_node(Direction::x_up).vals[0],
							w_ydown = w.get_neighbor_node(Direction::y_down).vals[0], w_yup = w.get_neighbor_node(Direction::y_up).vals[0],
							w_zdown = w.get_neighbor_node(Direction::z_down).vals[0], w_zup = w.get_neighbor_node(Direction::z_up).vals[0],
							viscosity_zdown = node.get_neighbor_node(Direction::z_down).customValues[ExternalFields::FLUID_viscosity],
							viscosity_zup = node.customValues[ExternalFields::FLUID_viscosity],
							viscosity_xdown = (viscosity_zdown + viscosity_zup + node.get_long_range_node(-1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, 0, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_xup = (viscosity_zdown + viscosity_zup + node.get_long_range_node(1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, 0, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_ydown = (viscosity_zdown + viscosity_zup + node.get_long_range_node(0, -1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, -1, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_yup = (viscosity_zdown + viscosity_zup + node.get_long_range_node(0, 1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, 1, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							w_laplace = (w_xdown + w_xup + w_ydown + w_yup + w_zdown + w_zup - 6.0 * w.vals[0]) / dr / dr,
							u_ave_down = (u.get_long_range_node(0, 0, 0).vals[0] + u.get_long_range_node(0, 0, -1).vals[0]) / 2.0,
							u_ave_up = (u.get_long_range_node(1, 0, 0).vals[0] + u.get_long_range_node(1, 0, -1).vals[0]) / 2.0,
							v_ave_down = (v.get_long_range_node(0, 0, 0).vals[0] + v.get_long_range_node(0, 0, -1).vals[0]) / 2.0,
							v_ave_up = (v.get_long_range_node(0, 1, 0).vals[0] + v.get_long_range_node(0, 1, -1).vals[0]) / 2.0;
						Vector3 vec_viscosity, vec_w;
						vec_viscosity[0] = (viscosity_xdown - viscosity_xup) / dr;
						vec_viscosity[1] = (viscosity_ydown - viscosity_yup) / dr;
						vec_viscosity[2] = (viscosity_zdown - viscosity_zup) / dr;
						vec_w[0] = (w_xdown - w_xup) / dr / 2.0;
						vec_w[1] = (w_ydown - w_yup) / dr / 2.0;
						vec_w[2] = (w_zdown - w_zup) / dr / 2.0;
						double C = (viscosity * w_laplace + vec_viscosity * vec_w) / density - ((w_zdown * u_ave_down - w_zup * u_ave_up) + (w_zdown * v_ave_down - w_zup * v_ave_up)
							+ (w_zdown * w_zdown - w_zup * w_zup)) / 2.0 / dr;
						w.vals[1] = C + (node.get_long_range_node(0, 0, -1).customVec3s[ExternalFields::FLUID_volume_force][0] + node.customVec3s[ExternalFields::FLUID_volume_force][0]) / 2.0;
					}
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::x_up).customValues[ExternalFields::FLUID_mass]) / 2.0
						, viscosity = (node.customValues[ExternalFields::FLUID_viscosity] + node.get_neighbor_node(Direction::x_up).customValues[ExternalFields::FLUID_viscosity]) / 2.0;
					// momentum equation for u
					{
						double u_xdown = u.get_neighbor_node(Direction::x_down).vals[0], u_xup = u.get_neighbor_node(Direction::x_up).vals[0],
							u_ydown = u.get_neighbor_node(Direction::y_down).vals[0], u_yup = u.get_neighbor_node(Direction::y_up).vals[0],
							u_zdown = u.get_neighbor_node(Direction::z_down).vals[0], u_zup = u.get_neighbor_node(Direction::z_up).vals[0],
							viscosity_xdown = node.customValues[ExternalFields::FLUID_viscosity],
							viscosity_xup = node.get_neighbor_node(Direction::x_up).customValues[ExternalFields::FLUID_viscosity],
							viscosity_ydown = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, -1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, -1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_yup = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, 1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, 1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zdown = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, 0, -1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, 0, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zup = (viscosity_xup + viscosity_xdown + node.get_long_range_node(0, 0, 1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, 0, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							u_laplace = (u_xdown + u_xup + u_ydown + u_yup + u_zdown + u_zup - 6.0 * u.vals[0]) / dr / dr,
							v_ave_up = (v.get_long_range_node(0, 0, 0).vals[0] + v.get_long_range_node(1, 0, 0).vals[0]) / 2.0,
							v_ave_down = (v.get_long_range_node(0, -1, 0).vals[0] + v.get_long_range_node(1, -1, 0).vals[0]) / 2.0,
							w_ave_up = (w.get_long_range_node(0, 0, 0).vals[0] + w.get_long_range_node(1, 0, 0).vals[0]) / 2.0,
							w_ave_down = (w.get_long_range_node(0, 0, -1).vals[0] + w.get_long_range_node(1, 0, -1).vals[0]) / 2.0;
						Vector3 vec_viscosity, vec_u;
						vec_viscosity[0] = (viscosity_xdown - viscosity_xup) / dr;
						vec_viscosity[1] = (viscosity_ydown - viscosity_yup) / dr;
						vec_viscosity[2] = (viscosity_zdown - viscosity_zup) / dr;
						vec_u[0] = (u_xdown - u_xup) / dr / 2.0;
						vec_u[1] = (u_ydown - u_yup) / dr / 2.0;
						vec_u[2] = (u_zdown - u_zup) / dr / 2.0;
						double A = (viscosity * u_laplace + vec_viscosity * vec_u) / density - ((u_xdown * u_xdown - u_xup * u_xup) + (u_ydown * v_ave_down - u_yup * v_ave_up)
							+ (u_zdown * w_ave_down - u_zup * w_ave_up)) / 2.0 / dr;
						u.vals[1] = A + (node.get_long_range_node(1, 0, 0).customVec3s[ExternalFields::FLUID_volume_force][0] + node.customVec3s[ExternalFields::FLUID_volume_force][0]) / 2.0;
					}
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::y_up).customValues[ExternalFields::FLUID_mass]) / 2.0
						, viscosity = (node.customValues[ExternalFields::FLUID_viscosity] + node.get_neighbor_node(Direction::y_up).customValues[ExternalFields::FLUID_viscosity]) / 2.0;
					// momentum equation for v
					{
						double v_xdown = v.get_neighbor_node(Direction::x_down).vals[0], v_xup = v.get_neighbor_node(Direction::x_up).vals[0],
							v_ydown = v.get_neighbor_node(Direction::y_down).vals[0], v_yup = v.get_neighbor_node(Direction::y_up).vals[0],
							v_zdown = v.get_neighbor_node(Direction::z_down).vals[0], v_zup = v.get_neighbor_node(Direction::z_up).vals[0],
							viscosity_ydown = node.customValues[ExternalFields::FLUID_viscosity],
							viscosity_yup = node.get_neighbor_node(Direction::y_up).customValues[ExternalFields::FLUID_viscosity],
							viscosity_xdown = (viscosity_ydown + viscosity_yup + node.get_long_range_node(-1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, 1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_xup = (viscosity_ydown + viscosity_yup + node.get_long_range_node(1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, 1, 0).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zdown = (viscosity_ydown + viscosity_yup + node.get_long_range_node(0, 0, -1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, 1, -1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_zup = (viscosity_ydown + viscosity_yup + node.get_long_range_node(0, 0, 1).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, 1, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							v_laplace = (v_xdown + v_xup + v_ydown + v_yup + v_zdown + v_zup - 6.0 * v.vals[0]) / dr / dr,
							u_ave_up = (u.get_long_range_node(0, 0, 0).vals[0] + u.get_long_range_node(0, 1, 0).vals[0]) / 2.0,
							u_ave_down = (u.get_long_range_node(-1, 0, 0).vals[0] + u.get_long_range_node(-1, 1, 0).vals[0]) / 2.0,
							w_ave_up = (w.get_long_range_node(0, 0, 0).vals[0] + w.get_long_range_node(0, 1, 0).vals[0]) / 2.0,
							w_ave_down = (w.get_long_range_node(0, 0, -1).vals[0] + w.get_long_range_node(0, 1, -1).vals[0]) / 2.0;
						Vector3 vec_viscosity, vec_v;
						vec_viscosity[0] = (viscosity_xdown - viscosity_xup) / dr;
						vec_viscosity[1] = (viscosity_ydown - viscosity_yup) / dr;
						vec_viscosity[2] = (viscosity_zdown - viscosity_zup) / dr;
						vec_v[0] = (v_xdown - v_xup) / dr / 2.0;
						vec_v[1] = (v_ydown - v_yup) / dr / 2.0;
						vec_v[2] = (v_zdown - v_zup) / dr / 2.0;
						double B = (viscosity * v_laplace + vec_viscosity * vec_v) / density - ((v_xdown * u_ave_down - v_xup * u_ave_up) + (v_ydown * v_ydown - v_yup * v_yup)
							+ (v_zdown * w_ave_down - v_zup * w_ave_up)) / 2.0 / dr;
						v.vals[1] = B + (node.get_long_range_node(0, 1, 0).customVec3s[ExternalFields::FLUID_volume_force][0] + node.customVec3s[ExternalFields::FLUID_volume_force][0]) / 2.0;
					}
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::z_up).customValues[ExternalFields::FLUID_mass]) / 2.0
						, viscosity = (node.customValues[ExternalFields::FLUID_viscosity] + node.get_neighbor_node(Direction::z_up).customValues[ExternalFields::FLUID_viscosity]) / 2.0;
					// momentum equation for w
					{
						double w_xdown = w.get_neighbor_node(Direction::x_down).vals[0], w_xup = w.get_neighbor_node(Direction::x_up).vals[0],
							w_ydown = w.get_neighbor_node(Direction::y_down).vals[0], w_yup = w.get_neighbor_node(Direction::y_up).vals[0],
							w_zdown = w.get_neighbor_node(Direction::z_down).vals[0], w_zup = w.get_neighbor_node(Direction::z_up).vals[0],
							viscosity_zdown = node.customValues[ExternalFields::FLUID_viscosity],
							viscosity_zup = node.get_neighbor_node(Direction::z_up).customValues[ExternalFields::FLUID_viscosity],
							viscosity_xdown = (viscosity_zdown + viscosity_zup + node.get_long_range_node(-1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(-1, 0, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_xup = (viscosity_zdown + viscosity_zup + node.get_long_range_node(1, 0, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(1, 0, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_ydown = (viscosity_zdown + viscosity_zup + node.get_long_range_node(0, -1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, -1, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							viscosity_yup = (viscosity_zdown + viscosity_zup + node.get_long_range_node(0, 1, 0).customValues[ExternalFields::FLUID_viscosity]
								+ node.get_long_range_node(0, 1, 1).customValues[ExternalFields::FLUID_viscosity]) / 4.0,
							w_laplace = (w_xdown + w_xup + w_ydown + w_yup + w_zdown + w_zup - 6.0 * w.vals[0]) / dr / dr,
							u_ave_up = (u.get_long_range_node(0, 0, 0).vals[0] + u.get_long_range_node(0, 0, 1).vals[0]) / 2.0,
							u_ave_down = (u.get_long_range_node(-1, 0, 0).vals[0] + u.get_long_range_node(-1, 0, 1).vals[0]) / 2.0,
							v_ave_up = (v.get_long_range_node(0, 0, 0).vals[0] + v.get_long_range_node(0, 0, 1).vals[0]) / 2.0,
							v_ave_down = (v.get_long_range_node(0, -1, 0).vals[0] + v.get_long_range_node(0, -1, 1).vals[0]) / 2.0;
						Vector3 vec_viscosity, vec_w;
						vec_viscosity[0] = (viscosity_xdown - viscosity_xup) / dr;
						vec_viscosity[1] = (viscosity_ydown - viscosity_yup) / dr;
						vec_viscosity[2] = (viscosity_zdown - viscosity_zup) / dr;
						vec_w[0] = (w_xdown - w_xup) / dr / 2.0;
						vec_w[1] = (w_ydown - w_yup) / dr / 2.0;
						vec_w[2] = (w_zdown - w_zup) / dr / 2.0;
						double C = (viscosity * w_laplace + vec_viscosity * vec_w) / density - ((w_zdown * u_ave_down - w_zup * u_ave_up) + (w_zdown * v_ave_down - w_zup * v_ave_up)
							+ (w_zdown * w_zdown - w_zup * w_zup)) / 2.0 / dr;
						w.vals[1] = C + (node.get_long_range_node(0, 0, 1).customVec3s[ExternalFields::FLUID_volume_force][0] + node.customVec3s[ExternalFields::FLUID_volume_force][0]) / 2.0;
					}
				}
		}
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					double old_u = u.vals[0], old_v = v.vals[0], old_w = w.vals[0];
					u.vals[0] += dt * u.vals[1];
					v.vals[0] += dt * v.vals[1];
					w.vals[0] += dt * w.vals[1];
					boundary_condition_for_domain_U(u, node.get_neighbor_node(pf::Direction::x_down), node, U.limit_x, U.limit_y, U.limit_z);
					boundary_condition_for_domain_V(v, node.get_neighbor_node(pf::Direction::y_down), node, V.limit_x, V.limit_y, V.limit_z);
					boundary_condition_for_domain_W(w, node.get_neighbor_node(pf::Direction::z_down), node, W.limit_x, W.limit_y, W.limit_z);
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (abs(u.vals[0] - old_u) > dVelocity[0])
							dVelocity[0] = abs(u.vals[0] - old_u);
						if (abs(v.vals[0] - old_v) > dVelocity[1])
							dVelocity[1] = abs(v.vals[0] - old_v);
						if (abs(w.vals[0] - old_w) > dVelocity[2])
							dVelocity[2] = abs(w.vals[0] - old_w);
					}
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					double old_u = u.vals[0];
					u.vals[0] += dt * u.vals[1];
					boundary_condition_for_domain_U(u, node, node.get_neighbor_node(pf::Direction::x_up), U.limit_x, U.limit_y, U.limit_z);
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (abs(u.vals[0] - old_u) > dVelocity[0])
							dVelocity[0] = abs(u.vals[0] - old_u);
					}
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					double old_v = v.vals[0];
					v.vals[0] += dt * v.vals[1];
					boundary_condition_for_domain_V(v, node, node.get_neighbor_node(pf::Direction::y_up), V.limit_x, V.limit_y, V.limit_z);
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (abs(v.vals[0] - old_v) > dVelocity[1])
							dVelocity[1] = abs(v.vals[0] - old_v);
					}
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					double old_w = w.vals[0];
					w.vals[0] += dt * w.vals[1];
					boundary_condition_for_domain_W(w, node, node.get_neighbor_node(pf::Direction::z_up), W.limit_x, W.limit_y, W.limit_z);
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (abs(w.vals[0] - old_w) > dVelocity[2])
							dVelocity[2] = abs(w.vals[0] - old_w);
					}
				}
		}
		return dVelocity;
	}

	vector<double> FluidField::do_pressure_correction(double fluid_threshold, double accuracy, int max_iterate_steps, bool debug_solver, int output_step) {
		vector<double> re_value;
		re_value.push_back(0.0); // MAX_abs_dPressure
		re_value.push_back(SYS_EPSILON); // pressure_iterate_times
		double dr = phaseMesh->dr;
		do {
			re_value[0] = 0.0;
			re_value[1] += 1.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						if (node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_threshold) {
							node.customValues[ExternalFields::FLUID_pressure] = 0.0;
							continue;
						}
						else {
							VectorNode& u = U(x + Minus_0_5, y, z);
							VectorNode& v = V(x, y + Minus_0_5, z);
							VectorNode& w = W(x, y, z + Minus_0_5);
							double old_lhs_value = node.customValues[ExternalFields::FLUID_pressure], density = node.customValues[ExternalFields::FLUID_mass];
							double lhs_x_down = old_lhs_value, lhs_x_up = old_lhs_value, lhs_y_down = old_lhs_value, lhs_y_up = old_lhs_value,
								lhs_z_down = old_lhs_value, lhs_z_up = old_lhs_value;
							double dispersion = (u.vals[0] + v.vals[0] + w.vals[0] - u.get_neighbor_node(Direction::x_up).vals[0]
								- v.get_neighbor_node(Direction::y_up).vals[0] - w.get_neighbor_node(Direction::z_up).vals[0]) / dr;
							// - solve for fluid phases
							if (node.get_neighbor_node(Direction::x_down).customValues[ExternalFields::FLUID_Fluid_Domain] > fluid_threshold)
								lhs_x_down = node.get_neighbor_node(Direction::x_down).customValues[ExternalFields::FLUID_pressure];
							if (node.get_neighbor_node(Direction::x_up).customValues[ExternalFields::FLUID_Fluid_Domain] > fluid_threshold)
								lhs_x_up = node.get_neighbor_node(Direction::x_up).customValues[ExternalFields::FLUID_pressure];
							if (node.get_neighbor_node(Direction::y_down).customValues[ExternalFields::FLUID_Fluid_Domain] > fluid_threshold)
								lhs_y_down = node.get_neighbor_node(Direction::y_down).customValues[ExternalFields::FLUID_pressure];
							if (node.get_neighbor_node(Direction::y_up).customValues[ExternalFields::FLUID_Fluid_Domain] > fluid_threshold)
								lhs_y_up = node.get_neighbor_node(Direction::y_up).customValues[ExternalFields::FLUID_pressure];
							if (node.get_neighbor_node(Direction::z_down).customValues[ExternalFields::FLUID_Fluid_Domain] > fluid_threshold)
								lhs_z_down = node.get_neighbor_node(Direction::z_down).customValues[ExternalFields::FLUID_pressure];
							if (node.get_neighbor_node(Direction::z_up).customValues[ExternalFields::FLUID_Fluid_Domain] > fluid_threshold)
								lhs_z_up = node.get_neighbor_node(Direction::z_up).customValues[ExternalFields::FLUID_pressure];
							// -
							node.customValues[ExternalFields::FLUID_pressure] = (lhs_x_down + lhs_x_up + lhs_y_down + lhs_y_up + lhs_z_down + lhs_z_up
								- dr * dr * dispersion * density) / 6.0;
							boundary_condition_for_pressure(node, phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z);
							double variation = std::abs(old_lhs_value - node.customValues[ExternalFields::FLUID_pressure]);
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (variation > re_value[0])
									re_value[0] = variation;
							}
						}
					}
			if ((int(re_value[1]) % output_step == 0) && debug_solver) {
				cout << "> pressure solver iterate " << int(re_value[1]) << " times !" << endl;
				cout << "> pressure MAXVariation = " << re_value[0] << " !" << endl;
			}
		} while (re_value[0] > accuracy && int(re_value[1]) < max_iterate_steps);
		if (debug_solver) {
			cout << "> pressure solver iterate " << int(re_value[1]) << " times !" << endl;
			cout << "> pressure MAXVariation = " << re_value[0] << " !" << endl;
		}
		return re_value;
	}

	void FluidField::correcting_velocity_field(double dt) {
		double dr = phaseMesh->dr;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					// equation for u
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::x_down).customValues[ExternalFields::FLUID_mass]) / 2.0;
						u.vals[0] += -dt * (node.get_long_range_node(-1, 0, 0).customValues[ExternalFields::FLUID_pressure] - node.customValues[ExternalFields::FLUID_pressure]) / dr / density;
					}
					// equation for v
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::y_down).customValues[ExternalFields::FLUID_mass]) / 2.0;
						v.vals[0] += -dt * (node.get_long_range_node(0, -1, 0).customValues[ExternalFields::FLUID_pressure] - node.customValues[ExternalFields::FLUID_pressure]) / dr / density;
					}
					// equation for w
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::z_down).customValues[ExternalFields::FLUID_mass]) / 2.0;
						w.vals[0] += -dt * (node.get_long_range_node(0, 0, -1).customValues[ExternalFields::FLUID_pressure] - node.customValues[ExternalFields::FLUID_pressure]) / dr / density;
					}
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					// momentum equation for u
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::x_up).customValues[ExternalFields::FLUID_mass]) / 2.0;
						u.vals[0] += -dt * (node.customValues[ExternalFields::FLUID_pressure] - node.get_long_range_node(1, 0, 0).customValues[ExternalFields::FLUID_pressure]) / dr / density;
					}
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					// momentum equation for v
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::y_up).customValues[ExternalFields::FLUID_mass]) / 2.0;
						v.vals[0] += -dt * (node.customValues[ExternalFields::FLUID_pressure] - node.get_long_range_node(0, 1, 0).customValues[ExternalFields::FLUID_pressure]) / dr / density;
					}
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					// momentum equation for w
					{
						double  density = (node.customValues[ExternalFields::FLUID_mass] + node.get_neighbor_node(Direction::z_up).customValues[ExternalFields::FLUID_mass]) / 2.0;
						w.vals[0] += -dt * (node.customValues[ExternalFields::FLUID_pressure] - node.get_long_range_node(0, 0, 1).customValues[ExternalFields::FLUID_pressure]) / dr / density;
					}
				}
		}
	}

	void FluidField::assign_velocity_to_main_domain() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					node.customVec3s[ExternalFields::FLUID_velocity][0] = (u.vals[0] + u.get_neighbor_node(Direction::x_up).vals[0]) / 2.0;
					node.customVec3s[ExternalFields::FLUID_velocity][1] = (v.vals[0] + v.get_neighbor_node(Direction::y_up).vals[0]) / 2.0;
					node.customVec3s[ExternalFields::FLUID_velocity][2] = (w.vals[0] + w.get_neighbor_node(Direction::z_up).vals[0]) / 2.0;
				}
	}

	void FluidField::assign_velocity_to_fluid_domain() {
		bool is_fluid_field_exist = false;
		for (auto cusVec3 = (*phaseMesh)(0, 0, 0).customVec3s.begin(); cusVec3 < (*phaseMesh)(0, 0, 0).customVec3s.end(); cusVec3++)
			if (cusVec3->index == ExternalFields::FLUID_velocity)
				is_fluid_field_exist = true;
		if (!is_fluid_field_exist)
			return;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					if (x == 0)
						u.vals[0] = node.customVec3s[ExternalFields::FLUID_velocity][0];
					else
						u.vals[0] = (node.customVec3s[ExternalFields::FLUID_velocity][0] + node.get_neighbor_node(Direction::x_down).customVec3s[ExternalFields::FLUID_velocity][0]) / 2.0;

					if (y == 0)
						v.vals[0] = node.customVec3s[ExternalFields::FLUID_velocity][1];
					else
						v.vals[0] = (node.customVec3s[ExternalFields::FLUID_velocity][1] + node.get_neighbor_node(Direction::y_down).customVec3s[ExternalFields::FLUID_velocity][1]) / 2.0;

					if (z == 0)
						w.vals[0] = node.customVec3s[ExternalFields::FLUID_velocity][2];
					else
						w.vals[0] = (node.customVec3s[ExternalFields::FLUID_velocity][2] + node.get_neighbor_node(Direction::z_down).customVec3s[ExternalFields::FLUID_velocity][2]) / 2.0;
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					u.vals[0] = node.customVec3s[ExternalFields::FLUID_velocity][0];
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					v.vals[0] = node.customVec3s[ExternalFields::FLUID_velocity][1];
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					w.vals[0] = node.customVec3s[ExternalFields::FLUID_velocity][2];
				}
		}
	}

	void LBM::init_distribution_functions() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					_init_distribution_functions(node, LBM_F_INDEX);
				}
	}

	void LBM::collision() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					_collision(node, LBM_F_INDEX);
				}
		if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++)
							node.customValues[LBM_F_INDEX + index] += node.customValues[LBM_F_INDEX + index + 19];
					}
		}
		else if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++)
							node.customValues[LBM_F_INDEX + index] += node.customValues[LBM_F_INDEX + index + 19];
					}
		}
	}

	void LBM::streaming() {
		if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
			int z = 0;
#pragma omp parallel sections// OMP BEGIN
			{
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = 0; y < phaseMesh->limit_y; y++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_1] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_1];
						}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_2] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_2];
						}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_3] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_3];
						}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_4] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_4];
						}
				}
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_5] = node.get_long_range_node(-1, -1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_5];
						}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_6] = node.get_long_range_node(1, -1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_6];
						}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_7] = node.get_long_range_node(1, 1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_7];
						}
				}
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_8] = node.get_long_range_node(-1, 1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_8];
						}
				}
			}//OMP END
		}
		else {
#pragma omp parallel sections// OMP BEGIN
			{
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_1] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_1];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_2] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_2];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_3] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_3];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_4] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_4];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = phaseMesh->limit_z - 1; z > 0; z--) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_5] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_5];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z - 1; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_6] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_6];
							}
				}
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_7] = node.get_long_range_node(-1, -1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_7];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_8] = node.get_long_range_node(1, -1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_8];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_9] = node.get_long_range_node(1, 1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_9];
							}
				}
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_10] = node.get_long_range_node(-1, 1, 0).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_10];
							}
				}
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = phaseMesh->limit_z - 1; z > 0; z--) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_11] = node.get_long_range_node(-1, 0, -1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_11];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = phaseMesh->limit_z - 1; z > 0; z--) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_12] = node.get_long_range_node(1, 0, -1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_12];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x - 1; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z - 1; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_13] = node.get_long_range_node(1, 0, 1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_13];
							}
				}
#pragma omp section
				{
					for (int x = phaseMesh->limit_x - 1; x > 0; x--)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z - 1; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_14] = node.get_long_range_node(-1, 0, 1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_14];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--)
							for (int z = phaseMesh->limit_z - 1; z > 0; z--) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_15] = node.get_long_range_node(0, -1, -1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_15];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++)
							for (int z = phaseMesh->limit_z - 1; z > 0; z--) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_16] = node.get_long_range_node(0, 1, -1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_16];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y - 1; y++)
							for (int z = 0; z < phaseMesh->limit_z - 1; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_17] = node.get_long_range_node(0, 1, 1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_17];
							}
				}
#pragma omp section
				{
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = phaseMesh->limit_y - 1; y > 0; y--)
							for (int z = 0; z < phaseMesh->limit_z - 1; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_18] = node.get_long_range_node(0, -1, 1).customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_18];
							}
				}
			}//OMP END
		}
	}

	void LBM::boundary_condition() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					_boundary_condition(node, LBM_F_INDEX);
				}
	}

	Vector3 LBM::cal_macro_variables() {
		Vector3 MAX_MOMENTUM_VARIATION;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					Vector3 old_momuntum = node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro];
					_cal_macro_variables(node, LBM_F_INDEX);
					Vector3 MOMENTUM_VARIATION = old_momuntum - node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro];
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (MAX_MOMENTUM_VARIATION[0] < abs(MOMENTUM_VARIATION[0]))
							MAX_MOMENTUM_VARIATION[0] = abs(MOMENTUM_VARIATION[0]);
						if (MAX_MOMENTUM_VARIATION[1] < abs(MOMENTUM_VARIATION[1]))
							MAX_MOMENTUM_VARIATION[1] = abs(MOMENTUM_VARIATION[1]);
						if (MAX_MOMENTUM_VARIATION[2] < abs(MOMENTUM_VARIATION[2]))
							MAX_MOMENTUM_VARIATION[2] = abs(MOMENTUM_VARIATION[2]);
					}
				}
		return MAX_MOMENTUM_VARIATION;
	}

	void LBM::cal_macro_variables(double& F_MACRO_MAX_VARIATION) {
		F_MACRO_MAX_VARIATION = 0.0;
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					double old_f_macro = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					_cal_macro_variables(node, LBM_F_INDEX);
					double VARIATION = old_f_macro - node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (F_MACRO_MAX_VARIATION < abs(VARIATION))
							F_MACRO_MAX_VARIATION = abs(VARIATION);
					}
				}
	}
}
