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
#include "../../../Base.h"

namespace pf {
	namespace lbm_boundary_condition {
		enum Fluid_Domain_Boundary_Condition { FDBC_Wall, FDBC_Period, FDBC_Free, FDBC_Pressure, FDBC_Normal_Flow };
		enum Fluid_Boundary_Property { FBP_WallRoughness, FBP_WallSpeed, FBP_DensityValue, FBP_NormalFlowSpeed };
		// boundary condition
		static vector<double_box> fluid_boundary_condition; // fluid_boundary_condition[Direction]
		// smooth boundary
		const double solid_liquid_interface_threshold = 0.5;
		const double q_c = 0.75;
		static vector<int> solid_phases;
		static int Nx = 0;
		static int Ny = 0;
		static int Nz = 0;
		static double PCT_dt = 0.0;
		static double Cs2 = 1.0 / 3.0;
		static double Cs4 = 1.0 / 9.0;
		static vector<Vector3> d2q9_w;
		static vector<Vector3> d3q19_w;
		static vector<double> w;
		// viscosity
		static double viscosity_liquid = 0.0;
		static double viscosity_gas = 0.0;
		static double viscosity_two_phase(PhaseNode& node) {
			return node.customValues[ExternalFieldsPlus::DF_Macro_TwoPhase + LBM_Symbols::LBM_f_macro] * (viscosity_liquid - viscosity_gas) + viscosity_gas;
		}
		static double viscosity_one_phase(PhaseNode& node) {
			return viscosity_liquid;
		}
		static double (*viscosity)(PhaseNode& node);
		// density
		double density_liquid = 1.0;
		double density_gas = 1.0;
		static double density_two_phase(PhaseNode& node) {
			return node.customValues[ExternalFieldsPlus::DF_Macro_TwoPhase + LBM_Symbols::LBM_f_macro] * (density_liquid - density_gas) + density_gas;
		}
		static double density_one_phase(PhaseNode& node) {
			return density_liquid;
		}
		static double (*density)(PhaseNode& node);
		// tau
		static double mobility_two_phase = 0.0;
		static double _tau_const = DBL_MAX;
		static double _tau_two_phase = DBL_MAX;
		static double tau_const(PhaseNode& node) {
			return _tau_const;
		}
		static double tau_standard(PhaseNode& node) {
			return viscosity(node) / PCT_dt / Cs2 + 0.5;
		}
		static double tau_two_phase_const(PhaseNode& node) {
			return _tau_two_phase;
		}
		static double (*tau)(PhaseNode& node);
		static double (*tau_two_phase)(PhaseNode& node);

		namespace bc_funcs {
			static void default_domain_boundary_condition(pf::PhaseNode& node, int LBM_F_INDEX) {};
			static double f_eq_a_virtual_d2q9(int INDEX_A, double f_macro, Vector3& U) {
				double CU = d2q9_w[INDEX_A] * U;
				return w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double f_eq_a_virtual_d3q19(int INDEX_A, double f_macro, Vector3& U) {
				double CU = d3q19_w[INDEX_A] * U;
				return w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static void d2q9_fluid_solid_boundary_Guo2002(pf::PhaseNode& node, int LBM_F_INDEX) {
				double f_eq = 0.0, f_neq = 0.0, f_macro = 0.0, _tau = tau(node);
				Vector3 U;
				// LBM_f_0
				node.customValues[LBM_F_INDEX + LBM_f_0] = 0.0;
				PhaseNode* near_node;
				// LBM_f_1
				near_node = &node.get_neighbor_node(Direction::x_up);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_1, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_1] - f_eq_a_virtual_d2q9(LBM_f_1, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::x_up);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_1, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_1] - f_eq_a_virtual_d2q9(LBM_f_1, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_1] - f_eq_a_virtual_d2q9(LBM_f_1, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_1] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_1] = 0.0;
				}
				// LBM_f_2
				near_node = &node.get_neighbor_node(Direction::y_up);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_2, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_2] - f_eq_a_virtual_d2q9(LBM_f_2, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::y_up);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_2, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_2] - f_eq_a_virtual_d2q9(LBM_f_2, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_2] - f_eq_a_virtual_d2q9(LBM_f_2, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_2] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_2] = 0.0;
				}
				// LBM_f_3
				near_node = &node.get_neighbor_node(Direction::x_down);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_3, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_3] - f_eq_a_virtual_d2q9(LBM_f_3, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::x_down);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_3, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_3] - f_eq_a_virtual_d2q9(LBM_f_3, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_3] - f_eq_a_virtual_d2q9(LBM_f_3, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_3] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_3] = 0.0;
				}
				// LBM_f_4
				near_node = &node.get_neighbor_node(Direction::y_down);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_4, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_4] - f_eq_a_virtual_d2q9(LBM_f_4, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::y_down);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_4, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_4] - f_eq_a_virtual_d2q9(LBM_f_4, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_4] - f_eq_a_virtual_d2q9(LBM_f_4, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_4] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_4] = 0.0;
				}
				// LBM_f_5
				near_node = &node.get_long_range_node(1, 1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_5, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_5] - f_eq_a_virtual_d2q9(LBM_f_5, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(1, 1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_5, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_5] - f_eq_a_virtual_d2q9(LBM_f_5, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_5] - f_eq_a_virtual_d2q9(LBM_f_5, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_5] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_5] = 0.0;
				}
				// LBM_f_6
				near_node = &node.get_long_range_node(-1, 1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_6, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_6] - f_eq_a_virtual_d2q9(LBM_f_6, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(-1, 1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_6, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_6] - f_eq_a_virtual_d2q9(LBM_f_6, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_6] - f_eq_a_virtual_d2q9(LBM_f_6, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_6] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_6] = 0.0;
				}
				// LBM_f_7
				near_node = &node.get_long_range_node(-1, -1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_7, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_7] - f_eq_a_virtual_d2q9(LBM_f_7, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(-1, -1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_7, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_7] - f_eq_a_virtual_d2q9(LBM_f_7, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_7] - f_eq_a_virtual_d2q9(LBM_f_7, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_7] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_7] = 0.0;
				}
				// LBM_f_8
				near_node = &node.get_long_range_node(1, -1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_f_8, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_8] - f_eq_a_virtual_d2q9(LBM_f_8, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(1, -1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_f_8, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_8] - f_eq_a_virtual_d2q9(LBM_f_8, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_8] - f_eq_a_virtual_d2q9(LBM_f_8, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_8] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_8] = 0.0;
				}
			}
			namespace lbm_bc_d2q9 {
				// FDBC_Wall_No_Slip
				static void wall_no_slip_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double roughness = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] * roughness + node.customValues[LBM_F_INDEX + LBM_f_6] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] * roughness + node.customValues[LBM_F_INDEX + LBM_f_7] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double roughness = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] * roughness + node.customValues[LBM_F_INDEX + LBM_f_5] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] * roughness + node.customValues[LBM_F_INDEX + LBM_f_8] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double roughness = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] * roughness + node.customValues[LBM_F_INDEX + LBM_f_8] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] * roughness + node.customValues[LBM_F_INDEX + LBM_f_7] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double roughness = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] * roughness + node.customValues[LBM_F_INDEX + LBM_f_6] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] * roughness + node.customValues[LBM_F_INDEX + LBM_f_5] * (1.0 - roughness);
					}
				}
				// FDBC_Wall_Slip
				static void wall_slip_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double loc_density = node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_4]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_6] + node.customValues[LBM_F_INDEX + LBM_f_7]),
							tail = loc_density * fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] * roughness + node.customValues[LBM_F_INDEX + LBM_f_6] * (1.0 - roughness) + tail;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] * roughness + node.customValues[LBM_F_INDEX + LBM_f_7] * (1.0 - roughness) - tail;
					}
				}
				static void wall_slip_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double loc_density = node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_4]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_8]),
							tail = loc_density * fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] * roughness + node.customValues[LBM_F_INDEX + LBM_f_5] * (1.0 - roughness) + tail;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] * roughness + node.customValues[LBM_F_INDEX + LBM_f_8] * (1.0 - roughness) - tail;
					}
				}
				static void wall_slip_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double loc_density = node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_3]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8]),
							tail = loc_density * fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] * roughness + node.customValues[LBM_F_INDEX + LBM_f_8] * (1.0 - roughness) + tail;
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] * roughness + node.customValues[LBM_F_INDEX + LBM_f_7] * (1.0 - roughness) - tail;
					}
				}
				static void wall_slip_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double loc_density = node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_3]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_6] + node.customValues[LBM_F_INDEX + LBM_f_5]),
							tail = loc_density * fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] * roughness + node.customValues[LBM_F_INDEX + LBM_f_6] * (1.0 - roughness) - tail;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] * roughness + node.customValues[LBM_F_INDEX + LBM_f_5] * (1.0 - roughness) + tail;
					}
				}
				// FDBC_Period
				static void period_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_8];
					}
				}
				static void period_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_6];
					}
				}
				static void period_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_6];
					}
				}
				static void period_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_8];
					}
				}
				// FDBC_Free
				static void free_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_8];
						// node.customValues[LBM_F_INDEX + LBM_f_1] = 2.0 * node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_1] - node.get_neighbor_node(Direction::x_up).get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_1];
						// node.customValues[LBM_F_INDEX + LBM_f_5] = 2.0 * node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_5] - node.get_neighbor_node(Direction::x_up).get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_5];
						// node.customValues[LBM_F_INDEX + LBM_f_8] = 2.0 * node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_8] - node.get_neighbor_node(Direction::x_up).get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_8];
					}
				}
				static void free_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_6];
						// node.customValues[LBM_F_INDEX + LBM_f_3] = 2.0 * node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_3] - node.get_neighbor_node(Direction::x_down).get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_3];
						// node.customValues[LBM_F_INDEX + LBM_f_7] = 2.0 * node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_7] - node.get_neighbor_node(Direction::x_down).get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_7];
						// node.customValues[LBM_F_INDEX + LBM_f_6] = 2.0 * node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_6] - node.get_neighbor_node(Direction::x_down).get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_6];
					}
				}
				static void free_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_6];
						//node.customValues[LBM_F_INDEX + LBM_f_2] = 2.0 * node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_2] - node.get_neighbor_node(Direction::y_up).get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_2];
						//node.customValues[LBM_F_INDEX + LBM_f_5] = 2.0 * node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_5] - node.get_neighbor_node(Direction::y_up).get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_5];
						//node.customValues[LBM_F_INDEX + LBM_f_6] = 2.0 * node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_6] - node.get_neighbor_node(Direction::y_up).get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_6];
					}
				}
				static void free_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_8];
						//node.customValues[LBM_F_INDEX + LBM_f_4] = 2.0 * node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_4] - node.get_neighbor_node(Direction::y_down).get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_4];
						//node.customValues[LBM_F_INDEX + LBM_f_7] = 2.0 * node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_7] - node.get_neighbor_node(Direction::y_down).get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_7];
						//node.customValues[LBM_F_INDEX + LBM_f_8] = 2.0 * node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_8] - node.get_neighbor_node(Direction::y_down).get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_8];
					}
				}
				// FDBC_Pressure
				static void pressure_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double pu = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] - (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_2] +
							node.customValues[LBM_F_INDEX + LBM_f_4] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_6] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_7]),
							diff = (node.customValues[LBM_F_INDEX + LBM_f_4] - node.customValues[LBM_F_INDEX + LBM_f_2]) / 2.0;
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_3] + 2.0 / 3.0 * pu;
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] + diff + pu / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] - diff + pu / 6.0;
					}
				}
				static void pressure_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double pu = node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_2] +
							node.customValues[LBM_F_INDEX + LBM_f_4] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_5] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_8] - fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue],
							diff = (node.customValues[LBM_F_INDEX + LBM_f_4] - node.customValues[LBM_F_INDEX + LBM_f_2]) / 2.0;
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_1] - 2.0 / 3.0 * pu;
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] + diff - pu / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] - diff - pu / 6.0;
					}
				}
				static void pressure_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double pu = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] - (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_1] +
							node.customValues[LBM_F_INDEX + LBM_f_3] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_7] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_8]),
							diff = (node.customValues[LBM_F_INDEX + LBM_f_3] - node.customValues[LBM_F_INDEX + LBM_f_1]) / 2.0;
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_4] + 2.0 / 3.0 * pu;
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] + diff + pu / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] - diff + pu / 6.0;
					}
				}
				static void pressure_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double pu = (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_1] +
							node.customValues[LBM_F_INDEX + LBM_f_3] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_5] + 2.0 * node.customValues[LBM_F_INDEX + LBM_f_6]) - fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue],
							diff = (node.customValues[LBM_F_INDEX + LBM_f_3] - node.customValues[LBM_F_INDEX + LBM_f_1]) / 2.0;
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_2] - 2.0 / 3.0 * pu;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] - diff - pu / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] + diff - pu / 6.0;
					}
				}
				// FDBC_Normal_Micro_Flow
				static void normal_micro_flow_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double local_density = (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_4]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_6] + node.customValues[LBM_F_INDEX + LBM_f_7])) / (1.0 - fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_3] + 2.0 / 3.0 * tail;
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] + tail / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] + tail / 6.0;
					}
				}
				static void normal_micro_flow_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double local_density = (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_4]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_8])) / (1.0 - fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_1] - 2.0 / 3.0 * tail;
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] - tail / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] - tail / 6.0;
					}
				}
				static void normal_micro_flow_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double local_density = (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_3]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8])) / (1.0 - fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_4] + 2.0 / 3.0 * tail;
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_7] + tail / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_8] + tail / 6.0;
					}
				}
				static void normal_micro_flow_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double local_density = (node.customValues[LBM_F_INDEX + LBM_f_0] + node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_3]
							+ 2 * (node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6])) / (1.0 - fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_2] - 2.0 / 3.0 * tail;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_5] - tail / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_6] - tail / 6.0;
					}
				}
			}

			static void d3q19_fluid_solid_boundary_Guo2002(pf::PhaseNode& node, int LBM_F_INDEX) {
				double f_eq = 0.0, f_neq = 0.0, f_macro = 0.0, _tau = tau(node);
				Vector3 U;
				// LBM_f_0
				node.customValues[LBM_F_INDEX + LBM_f_0] = 0.0;
				PhaseNode* near_node;
				// LBM_f_1
				near_node = &node.get_neighbor_node(Direction::x_up);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_1, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_1] - f_eq_a_virtual_d3q19(LBM_f_1, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::x_up);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_1, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_1] - f_eq_a_virtual_d3q19(LBM_f_1, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_1] - f_eq_a_virtual_d3q19(LBM_f_1, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_1] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_1] = 0.0;
				}
				// LBM_f_2
				near_node = &node.get_neighbor_node(Direction::x_down);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_2, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_2] - f_eq_a_virtual_d3q19(LBM_f_2, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::x_down);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_2, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_2] - f_eq_a_virtual_d3q19(LBM_f_2, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_2] - f_eq_a_virtual_d3q19(LBM_f_2, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_2] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_2] = 0.0;
				}
				// LBM_f_3
				near_node = &node.get_neighbor_node(Direction::y_up);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_3, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_3] - f_eq_a_virtual_d3q19(LBM_f_3, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::y_up);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_3, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_3] - f_eq_a_virtual_d3q19(LBM_f_3, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_3] - f_eq_a_virtual_d3q19(LBM_f_3, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_3] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_3] = 0.0;
				}
				// LBM_f_4
				near_node = &node.get_neighbor_node(Direction::y_down);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_4, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_4] - f_eq_a_virtual_d3q19(LBM_f_4, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::y_down);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_4, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_4] - f_eq_a_virtual_d3q19(LBM_f_4, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_4] - f_eq_a_virtual_d3q19(LBM_f_4, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_4] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_4] = 0.0;
				}
				// LBM_f_5
				near_node = &node.get_neighbor_node(Direction::z_up);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_5, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_5] - f_eq_a_virtual_d3q19(LBM_f_5, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::z_up);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_5, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_5] - f_eq_a_virtual_d3q19(LBM_f_5, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_5] - f_eq_a_virtual_d3q19(LBM_f_5, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_5] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_5] = 0.0;
				}
				// LBM_f_6
				near_node = &node.get_neighbor_node(Direction::z_down);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_6, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_6] - f_eq_a_virtual_d3q19(LBM_f_6, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_neighbor_node(Direction::z_down);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_6, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_6] - f_eq_a_virtual_d3q19(LBM_f_6, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_6] - f_eq_a_virtual_d3q19(LBM_f_6, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_6] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_6] = 0.0;
				}
				// LBM_f_7
				near_node = &node.get_long_range_node(1, 1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_7, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_7] - f_eq_a_virtual_d3q19(LBM_f_7, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(1, 1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_7, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_7] - f_eq_a_virtual_d3q19(LBM_f_7, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_7] - f_eq_a_virtual_d3q19(LBM_f_7, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_7] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_7] = 0.0;
				}
				// LBM_f_8
				near_node = &node.get_long_range_node(-1, 1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_8, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_8] - f_eq_a_virtual_d3q19(LBM_f_8, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(-1, 1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_8, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_8] - f_eq_a_virtual_d3q19(LBM_f_8, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_8] - f_eq_a_virtual_d3q19(LBM_f_8, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_8] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_8] = 0.0;
				}
				// LBM_f_9
				near_node = &node.get_long_range_node(-1, -1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_9, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_9] - f_eq_a_virtual_d3q19(LBM_f_9, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(-1, -1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_9, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_9] - f_eq_a_virtual_d3q19(LBM_f_9, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_9] - f_eq_a_virtual_d3q19(LBM_f_9, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_9] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_9] = 0.0;
				}
				// LBM_f_10
				near_node = &node.get_long_range_node(1, -1, 0);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_10, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_10] - f_eq_a_virtual_d3q19(LBM_f_10, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(1, -1, 0);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_10, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_10] - f_eq_a_virtual_d3q19(LBM_f_10, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_10] - f_eq_a_virtual_d3q19(LBM_f_10, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_10] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_10] = 0.0;
				}
				// LBM_f_11
				near_node = &node.get_long_range_node(1, 0, 1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_11, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_11] - f_eq_a_virtual_d3q19(LBM_f_11, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(1, 0, 1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_11, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_11] - f_eq_a_virtual_d3q19(LBM_f_11, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_11] - f_eq_a_virtual_d3q19(LBM_f_11, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_11] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_11] = 0.0;
				}
				// LBM_f_12
				near_node = &node.get_long_range_node(-1, 0, 1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_12, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_12] - f_eq_a_virtual_d3q19(LBM_f_12, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(-1, 0, 1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_12, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_12] - f_eq_a_virtual_d3q19(LBM_f_12, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_12] - f_eq_a_virtual_d3q19(LBM_f_12, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_12] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_12] = 0.0;
				}
				// LBM_f_13
				near_node = &node.get_long_range_node(-1, 0, -1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_13, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_13] - f_eq_a_virtual_d3q19(LBM_f_13, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(-1, 0, -1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_13, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_13] - f_eq_a_virtual_d3q19(LBM_f_13, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_13] - f_eq_a_virtual_d3q19(LBM_f_13, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_13] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_13] = 0.0;
				}
				// LBM_f_14
				near_node = &node.get_long_range_node(1, 0, -1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_14, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_14] - f_eq_a_virtual_d3q19(LBM_f_14, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(1, 0, -1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_14, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_14] - f_eq_a_virtual_d3q19(LBM_f_14, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_14] - f_eq_a_virtual_d3q19(LBM_f_14, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_14] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_14] = 0.0;
				}
				// LBM_f_15
				near_node = &node.get_long_range_node(0, 1, 1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_15, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_15] - f_eq_a_virtual_d3q19(LBM_f_15, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(0, 1, 1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_15, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_15] - f_eq_a_virtual_d3q19(LBM_f_15, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_15] - f_eq_a_virtual_d3q19(LBM_f_15, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_15] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_15] = 0.0;
				}
				// LBM_f_16
				near_node = &node.get_long_range_node(0, -1, 1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_16, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_16] - f_eq_a_virtual_d3q19(LBM_f_16, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(0, -1, 1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_16, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_16] - f_eq_a_virtual_d3q19(LBM_f_16, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_16] - f_eq_a_virtual_d3q19(LBM_f_16, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_16] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_16] = 0.0;
				}
				// LBM_f_17
				near_node = &node.get_long_range_node(0, -1, -1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_17, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_17] - f_eq_a_virtual_d3q19(LBM_f_17, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(0, -1, -1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_17, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_17] - f_eq_a_virtual_d3q19(LBM_f_17, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_17] - f_eq_a_virtual_d3q19(LBM_f_17, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_17] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_17] = 0.0;
				}
				// LBM_f_18
				near_node = &node.get_long_range_node(0, 1, -1);
				if (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_node->customValues[ExternalFields::FLUID_solid_fraction]) / (node.customValues[ExternalFields::FLUID_solid_fraction] - near_node->customValues[ExternalFields::FLUID_solid_fraction]);
					double qqq = 1.0 - (solid_liquid_interface_threshold - node.customValues[ExternalFields::FLUID_Fluid_Domain]) / (near_node->customValues[ExternalFields::FLUID_Fluid_Domain] - node.customValues[ExternalFields::FLUID_Fluid_Domain]);
					f_macro = near_node->customValues[LBM_F_INDEX + LBM_f_macro];
					if (qqq >= q_c) {
						U = (node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d3q19(LBM_f_18, f_macro, U);
						f_neq = near_node->customValues[LBM_F_INDEX + LBM_f_18] - f_eq_a_virtual_d3q19(LBM_f_18, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity]);
					}
					else {
						PhaseNode& near_near_node = near_node->get_long_range_node(0, 1, -1);
						U = node.customVec3s[ExternalFields::FLUID_velocity] + near_node->customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)
							+ (node.customVec3s[ExternalFields::FLUID_velocity] * 2.0 + near_near_node.customVec3s[ExternalFields::FLUID_velocity] * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d3q19(LBM_f_18, f_macro, U);
						f_neq = (near_node->customValues[LBM_F_INDEX + LBM_f_18] - f_eq_a_virtual_d3q19(LBM_f_18, f_macro, near_node->customVec3s[ExternalFields::FLUID_velocity])) * qqq
							+ (near_near_node.customValues[LBM_F_INDEX + LBM_f_18] - f_eq_a_virtual_d3q19(LBM_f_18, near_near_node.customValues[LBM_F_INDEX + LBM_f_macro], near_near_node.customVec3s[ExternalFields::FLUID_velocity])) * (1.0 - qqq);
					}
					node.customValues[LBM_F_INDEX + LBM_f_18] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_f_18] = 0.0;
				}
			}
			namespace lbm_bc_d3q19 {
				// FDBC_Wall_No_Slip
				static void wall_no_slip_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double roughness = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_9] * roughness + node.customValues[LBM_F_INDEX + LBM_f_8] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.customValues[LBM_F_INDEX + LBM_f_8] * roughness + node.customValues[LBM_F_INDEX + LBM_f_9] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.customValues[LBM_F_INDEX + LBM_f_13] * roughness + node.customValues[LBM_F_INDEX + LBM_f_12] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.customValues[LBM_F_INDEX + LBM_f_12] * roughness + node.customValues[LBM_F_INDEX + LBM_f_13] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double roughness = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_10] * roughness + node.customValues[LBM_F_INDEX + LBM_f_7] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.customValues[LBM_F_INDEX + LBM_f_7] * roughness + node.customValues[LBM_F_INDEX + LBM_f_10] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.customValues[LBM_F_INDEX + LBM_f_14] * roughness + node.customValues[LBM_F_INDEX + LBM_f_11] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.customValues[LBM_F_INDEX + LBM_f_11] * roughness + node.customValues[LBM_F_INDEX + LBM_f_14] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double roughness = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_9] * roughness + node.customValues[LBM_F_INDEX + LBM_f_10] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_10] * roughness + node.customValues[LBM_F_INDEX + LBM_f_9] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.customValues[LBM_F_INDEX + LBM_f_17] * roughness + node.customValues[LBM_F_INDEX + LBM_f_16] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.customValues[LBM_F_INDEX + LBM_f_16] * roughness + node.customValues[LBM_F_INDEX + LBM_f_17] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double roughness = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.customValues[LBM_F_INDEX + LBM_f_7] * roughness + node.customValues[LBM_F_INDEX + LBM_f_8] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.customValues[LBM_F_INDEX + LBM_f_8] * roughness + node.customValues[LBM_F_INDEX + LBM_f_7] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.customValues[LBM_F_INDEX + LBM_f_18] * roughness + node.customValues[LBM_F_INDEX + LBM_f_15] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.customValues[LBM_F_INDEX + LBM_f_15] * roughness + node.customValues[LBM_F_INDEX + LBM_f_18] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_z_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == 0) {
						double roughness = fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_6];
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.customValues[LBM_F_INDEX + LBM_f_13] * roughness + node.customValues[LBM_F_INDEX + LBM_f_14] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.customValues[LBM_F_INDEX + LBM_f_14] * roughness + node.customValues[LBM_F_INDEX + LBM_f_13] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.customValues[LBM_F_INDEX + LBM_f_17] * roughness + node.customValues[LBM_F_INDEX + LBM_f_18] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.customValues[LBM_F_INDEX + LBM_f_18] * roughness + node.customValues[LBM_F_INDEX + LBM_f_17] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_z_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == Nz - 1) {
						double roughness = fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_WallRoughness];
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.customValues[LBM_F_INDEX + LBM_f_11] * roughness + node.customValues[LBM_F_INDEX + LBM_f_12] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.customValues[LBM_F_INDEX + LBM_f_12] * roughness + node.customValues[LBM_F_INDEX + LBM_f_11] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.customValues[LBM_F_INDEX + LBM_f_15] * roughness + node.customValues[LBM_F_INDEX + LBM_f_16] * (1.0 - roughness);
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.customValues[LBM_F_INDEX + LBM_f_16] * roughness + node.customValues[LBM_F_INDEX + LBM_f_15] * (1.0 - roughness);
					}
				}
				// FDBC_Period
				static void period_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_10];
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_11];
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_14];
					}
				}
				static void period_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_8];
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_9];
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_12];
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_13];
					}
				}
				static void period_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_8];
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_15];
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_18];
					}
				}
				static void period_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_9];
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_10];
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_16];
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_17];
					}
				}
				static void period_z_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_11];
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_12];
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_15];
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_16];
					}
				}
				static void period_z_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == Nz - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_6];
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_13];
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_14];
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_17];
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_18];
					}
				}
				// FDBC_Free
				static void free_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_1];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_10];
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_11];
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.get_neighbor_node(Direction::x_up).customValues[LBM_F_INDEX + LBM_f_14];
					}
				}
				static void free_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_2];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_8];
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_9];
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_12];
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.get_neighbor_node(Direction::x_down).customValues[LBM_F_INDEX + LBM_f_13];
					}
				}
				static void free_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_3];
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_7];
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_8];
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_15];
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.get_neighbor_node(Direction::y_up).customValues[LBM_F_INDEX + LBM_f_18];
					}
				}
				static void free_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_4];
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_9];
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_10];
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_16];
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.get_neighbor_node(Direction::y_down).customValues[LBM_F_INDEX + LBM_f_17];
					}
				}
				static void free_z_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == 0) {
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_5];
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_11];
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_12];
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_15];
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.get_neighbor_node(Direction::z_up).customValues[LBM_F_INDEX + LBM_f_16];
					}
				}
				static void free_z_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == Nz - 1) {
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_6];
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_13];
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_14];
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_17];
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.get_neighbor_node(Direction::z_down).customValues[LBM_F_INDEX + LBM_f_18];
					}
				}
				// FDBC_Pressure
				static void pressure_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double locDensity = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_15] + node.customValues[LBM_F_INDEX + LBM_f_16] + node.customValues[LBM_F_INDEX + LBM_f_17] + node.customValues[LBM_F_INDEX + LBM_f_18]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_12]
									+ node.customValues[LBM_F_INDEX + LBM_f_13])) / locDensity;
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_2] + locDensity * u / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_9] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.customValues[LBM_F_INDEX + LBM_f_8] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.customValues[LBM_F_INDEX + LBM_f_13] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.customValues[LBM_F_INDEX + LBM_f_12] + locDensity * u / 6.0;
					}
				}
				static void pressure_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double locDensity = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_15] + node.customValues[LBM_F_INDEX + LBM_f_16] + node.customValues[LBM_F_INDEX + LBM_f_17] + node.customValues[LBM_F_INDEX + LBM_f_18]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_10] + node.customValues[LBM_F_INDEX + LBM_f_11]
									+ node.customValues[LBM_F_INDEX + LBM_f_14])) / locDensity;
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_1] + locDensity * u / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_10] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.customValues[LBM_F_INDEX + LBM_f_7] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.customValues[LBM_F_INDEX + LBM_f_14] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.customValues[LBM_F_INDEX + LBM_f_11] + locDensity * u / 6.0;
					}
				}
				static void pressure_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double locDensity = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_11] + node.customValues[LBM_F_INDEX + LBM_f_12] + node.customValues[LBM_F_INDEX + LBM_f_13] + node.customValues[LBM_F_INDEX + LBM_f_14]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_10] + node.customValues[LBM_F_INDEX + LBM_f_16]
									+ node.customValues[LBM_F_INDEX + LBM_f_17])) / locDensity;
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_4] + locDensity * u / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_9] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_10] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.customValues[LBM_F_INDEX + LBM_f_17] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.customValues[LBM_F_INDEX + LBM_f_16] + locDensity * u / 6.0;
					}
				}
				static void pressure_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double locDensity = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_11] + node.customValues[LBM_F_INDEX + LBM_f_12] + node.customValues[LBM_F_INDEX + LBM_f_13] + node.customValues[LBM_F_INDEX + LBM_f_14]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_15]
									+ node.customValues[LBM_F_INDEX + LBM_f_18])) / locDensity;
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_3] + locDensity * u / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.customValues[LBM_F_INDEX + LBM_f_7] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.customValues[LBM_F_INDEX + LBM_f_8] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.customValues[LBM_F_INDEX + LBM_f_18] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.customValues[LBM_F_INDEX + LBM_f_15] + locDensity * u / 6.0;
					}
				}
				static void pressure_z_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == 0) {
						double locDensity = fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4]
								+ node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_10]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_6] + node.customValues[LBM_F_INDEX + LBM_f_13] + node.customValues[LBM_F_INDEX + LBM_f_14] + node.customValues[LBM_F_INDEX + LBM_f_17]
									+ node.customValues[LBM_F_INDEX + LBM_f_18])) / locDensity;
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_6] + locDensity * u / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.customValues[LBM_F_INDEX + LBM_f_13] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.customValues[LBM_F_INDEX + LBM_f_14] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.customValues[LBM_F_INDEX + LBM_f_17] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.customValues[LBM_F_INDEX + LBM_f_18] + locDensity * u / 6.0;
					}
				}
				static void pressure_z_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == Nz - 1) {
						double locDensity = fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4]
								+ node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_10]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_11] + node.customValues[LBM_F_INDEX + LBM_f_12] + node.customValues[LBM_F_INDEX + LBM_f_15]
									+ node.customValues[LBM_F_INDEX + LBM_f_16])) / locDensity;
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_5] + locDensity * u / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.customValues[LBM_F_INDEX + LBM_f_11] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.customValues[LBM_F_INDEX + LBM_f_12] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.customValues[LBM_F_INDEX + LBM_f_15] + locDensity * u / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.customValues[LBM_F_INDEX + LBM_f_16] + locDensity * u / 6.0;
					}
				}
				// FDBC_Normal_Micro_Flow
				static void normal_micro_flow_x_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == 0) {
						double locVelocity = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_15] + node.customValues[LBM_F_INDEX + LBM_f_16] + node.customValues[LBM_F_INDEX + LBM_f_17] + node.customValues[LBM_F_INDEX + LBM_f_18]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_12]
									+ node.customValues[LBM_F_INDEX + LBM_f_13])) / (1.0 - locVelocity);
						node.customValues[LBM_F_INDEX + LBM_f_1] = node.customValues[LBM_F_INDEX + LBM_f_2] + locVelocity * locDensity / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_9] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.customValues[LBM_F_INDEX + LBM_f_8] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.customValues[LBM_F_INDEX + LBM_f_13] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.customValues[LBM_F_INDEX + LBM_f_12] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_x_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._x == Nx - 1) {
						double locVelocity = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_15] + node.customValues[LBM_F_INDEX + LBM_f_16] + node.customValues[LBM_F_INDEX + LBM_f_17] + node.customValues[LBM_F_INDEX + LBM_f_18]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_10] + node.customValues[LBM_F_INDEX + LBM_f_11]
									+ node.customValues[LBM_F_INDEX + LBM_f_14])) / (1.0 - locVelocity);
						node.customValues[LBM_F_INDEX + LBM_f_2] = node.customValues[LBM_F_INDEX + LBM_f_1] + locVelocity * locDensity / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_10] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.customValues[LBM_F_INDEX + LBM_f_7] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.customValues[LBM_F_INDEX + LBM_f_14] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.customValues[LBM_F_INDEX + LBM_f_11] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_y_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == 0) {
						double locVelocity = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_11] + node.customValues[LBM_F_INDEX + LBM_f_12] + node.customValues[LBM_F_INDEX + LBM_f_13] + node.customValues[LBM_F_INDEX + LBM_f_14]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_4] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_10] + node.customValues[LBM_F_INDEX + LBM_f_16]
									+ node.customValues[LBM_F_INDEX + LBM_f_17])) / (1.0 - locVelocity);
						node.customValues[LBM_F_INDEX + LBM_f_3] = node.customValues[LBM_F_INDEX + LBM_f_4] + locVelocity * locDensity / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_7] = node.customValues[LBM_F_INDEX + LBM_f_9] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_8] = node.customValues[LBM_F_INDEX + LBM_f_10] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.customValues[LBM_F_INDEX + LBM_f_17] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.customValues[LBM_F_INDEX + LBM_f_16] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_y_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._y == Ny - 1) {
						double locVelocity = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_6]
								+ node.customValues[LBM_F_INDEX + LBM_f_11] + node.customValues[LBM_F_INDEX + LBM_f_12] + node.customValues[LBM_F_INDEX + LBM_f_13] + node.customValues[LBM_F_INDEX + LBM_f_14]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_15]
									+ node.customValues[LBM_F_INDEX + LBM_f_18])) / (1.0 - locVelocity);
						node.customValues[LBM_F_INDEX + LBM_f_4] = node.customValues[LBM_F_INDEX + LBM_f_3] + locVelocity * locDensity / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_9] = node.customValues[LBM_F_INDEX + LBM_f_7] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_10] = node.customValues[LBM_F_INDEX + LBM_f_8] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.customValues[LBM_F_INDEX + LBM_f_18] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.customValues[LBM_F_INDEX + LBM_f_15] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_z_down(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == 0) {
						double locVelocity = fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4]
								+ node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_10]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_6] + node.customValues[LBM_F_INDEX + LBM_f_13] + node.customValues[LBM_F_INDEX + LBM_f_14] + node.customValues[LBM_F_INDEX + LBM_f_17]
									+ node.customValues[LBM_F_INDEX + LBM_f_18])) / (1.0 - locVelocity);
						node.customValues[LBM_F_INDEX + LBM_f_5] = node.customValues[LBM_F_INDEX + LBM_f_6] + locVelocity * locDensity / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_11] = node.customValues[LBM_F_INDEX + LBM_f_13] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_12] = node.customValues[LBM_F_INDEX + LBM_f_14] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_15] = node.customValues[LBM_F_INDEX + LBM_f_17] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_16] = node.customValues[LBM_F_INDEX + LBM_f_18] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_z_up(pf::PhaseNode& node, int LBM_F_INDEX) {
					if (node._z == Nz - 1) {
						double locVelocity = fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (node.customValues[LBM_F_INDEX + LBM_f_0]
								+ node.customValues[LBM_F_INDEX + LBM_f_1] + node.customValues[LBM_F_INDEX + LBM_f_2] + node.customValues[LBM_F_INDEX + LBM_f_3] + node.customValues[LBM_F_INDEX + LBM_f_4]
								+ node.customValues[LBM_F_INDEX + LBM_f_7] + node.customValues[LBM_F_INDEX + LBM_f_8] + node.customValues[LBM_F_INDEX + LBM_f_9] + node.customValues[LBM_F_INDEX + LBM_f_10]
								+ 2.0 * (node.customValues[LBM_F_INDEX + LBM_f_5] + node.customValues[LBM_F_INDEX + LBM_f_11] + node.customValues[LBM_F_INDEX + LBM_f_12] + node.customValues[LBM_F_INDEX + LBM_f_15]
									+ node.customValues[LBM_F_INDEX + LBM_f_16])) / (1.0 - locVelocity);
						node.customValues[LBM_F_INDEX + LBM_f_6] = node.customValues[LBM_F_INDEX + LBM_f_5] + locVelocity * locDensity / 3.0;
						node.customValues[LBM_F_INDEX + LBM_f_13] = node.customValues[LBM_F_INDEX + LBM_f_11] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_14] = node.customValues[LBM_F_INDEX + LBM_f_12] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_17] = node.customValues[LBM_F_INDEX + LBM_f_15] + locVelocity * locDensity / 6.0;
						node.customValues[LBM_F_INDEX + LBM_f_18] = node.customValues[LBM_F_INDEX + LBM_f_16] + locVelocity * locDensity / 6.0;
					}
				}
			}
		}

		static void (*d2q9_domain_boundary_x_down)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d2q9_domain_boundary_x_up)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d2q9_domain_boundary_y_down)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d2q9_domain_boundary_y_up)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void(*d2q9_fluid_solid_boundary)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void boundary_condition_d2q9(pf::PhaseNode& node, int LBM_F_INDEX) {
			d2q9_domain_boundary_x_down(node, LBM_F_INDEX);
			d2q9_domain_boundary_x_up(node, LBM_F_INDEX);
			d2q9_domain_boundary_y_down(node, LBM_F_INDEX);
			d2q9_domain_boundary_y_up(node, LBM_F_INDEX);
			if (node.customValues[ExternalFields::FLUID_Fluid_Domain] < solid_liquid_interface_threshold)
				d2q9_fluid_solid_boundary(node, LBM_F_INDEX);
		}

		static void (*d3q19_domain_boundary_x_down)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d3q19_domain_boundary_x_up)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d3q19_domain_boundary_y_down)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d3q19_domain_boundary_y_up)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d3q19_domain_boundary_z_down)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void (*d3q19_domain_boundary_z_up)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void(*d3q19_fluid_solid_boundary)(pf::PhaseNode& node, int LBM_F_INDEX);
		static void boundary_condition_d3q19(pf::PhaseNode& node, int LBM_F_INDEX) {
			d3q19_domain_boundary_x_down(node, LBM_F_INDEX);
			d3q19_domain_boundary_x_up(node, LBM_F_INDEX);
			d3q19_domain_boundary_y_down(node, LBM_F_INDEX);
			d3q19_domain_boundary_y_up(node, LBM_F_INDEX);
			d3q19_domain_boundary_z_down(node, LBM_F_INDEX);
			d3q19_domain_boundary_z_up(node, LBM_F_INDEX);
			if (node.customValues[ExternalFields::FLUID_Fluid_Domain] < solid_liquid_interface_threshold)
				d3q19_fluid_solid_boundary(node, LBM_F_INDEX);
		}

		static void cal_fluid_domain(FieldStorage_forPhaseNode& phaseMesh) {
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						double solid_phases_volume_fraction = 0.0;
						for (auto phi = node.begin(); phi < node.end(); phi++)
							for (auto phase = solid_phases.begin(); phase < solid_phases.end(); phase++)
								if (phi->property == *phase) {
									solid_phases_volume_fraction += phi->phi;
								}
						if (solid_phases_volume_fraction > 1.0)
							solid_phases_volume_fraction = 1.0;
						else if (solid_phases_volume_fraction < 0.0)
							solid_phases_volume_fraction = 0.0;
						node.customValues.add_double(ExternalFields::FLUID_Fluid_Domain, 1.0 - solid_phases_volume_fraction);
					}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh, LBM& fluid_lbm_solver) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			Nx = phaseMesh.limit_x;
			Ny = phaseMesh.limit_y;
			Nz = phaseMesh.limit_z;
			PCT_dt = Solvers::get_instance()->parameters.dt;
			tau = tau_const;
			viscosity = viscosity_one_phase;
			density = density_one_phase;
			double cc = phaseMesh.dr / PCT_dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			fluid_boundary_condition.resize(6);
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.FluidDynamics.LatticeBoltzmann.solid_phases = (phase_name, ... ) \n", InputFileReader::get_instance()->debug_file);
			string fluid_phase_key = "Postprocess.FluidDynamics.LatticeBoltzmann.solid_phases", fluid_phase_input = "()";
			InputFileReader::get_instance()->read_string_value(fluid_phase_key, fluid_phase_input, infile_debug);
			vector<input_value> fluid_phase_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, fluid_phase_key, fluid_phase_input, infile_debug);
			for (auto fluid_name = fluid_phase_value.begin(); fluid_name < fluid_phase_value.end(); fluid_name++) {
				solid_phases.push_back(Solvers::get_instance()->parameters.Phases[fluid_name->string_value].phi_property);
			}
			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# tau = viscosity / fluid_dt / Cs2 + 0.5 \n", InputFileReader::get_instance()->debug_file);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.liquid_viscosity", viscosity_liquid, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.liquid_density", density_liquid, infile_debug);
			_tau_const = viscosity_liquid / PCT_dt / Cs2 + 0.5;

			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				w.resize(9);
				w[0] = 4.0 / 9.0;
				w[1] = 1.0 / 9.0; w[2] = 1.0 / 9.0; w[3] = 1.0 / 9.0; w[4] = 1.0 / 9.0;
				w[5] = 1.0 / 36.0; w[6] = 1.0 / 36.0; w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0;
				fluid_lbm_solver._boundary_condition = boundary_condition_d2q9;
				// load boundary condition
				if (solid_phases.size() != 0)
					d2q9_fluid_solid_boundary = bc_funcs::d2q9_fluid_solid_boundary_Guo2002;
				else
					d2q9_fluid_solid_boundary = bc_funcs::default_domain_boundary_condition;
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .LatticeBoltzmann.boundary_condition = (down_x,up_x,down_y,up_y) \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                        0 - Wall, 1 - Period, 2 - Free, 3 - Pressure, 4 - Normal_Flow \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                            .pressure = p0 , density0 = p0 / Cs^2 , Cs = 1 / sqrt(3) \n", InputFileReader::get_instance()->debug_file);
				string bc_key = "Postprocess.FluidDynamics.LatticeBoltzmann.boundary_condition", bc_input = "(0,0,0,0)";
				InputFileReader::get_instance()->read_string_value(bc_key, bc_input, infile_debug);
				vector<input_value> bc_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, bc_key, bc_input, infile_debug);
				switch (Fluid_Domain_Boundary_Condition(bc_value[0].int_value)) // down_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_WallSpeed, 0.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.wall_roughness", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.wall_speed", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed], infile_debug);
					if (Is_Equality(fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::wall_no_slip_x_down;
					else
						d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::wall_slip_x_down;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::period_x_down;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::free_x_down;
					break;
				case FDBC_Pressure:
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::pressure_x_down;
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.pressure", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::normal_micro_flow_x_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.normal_velocity", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[1].int_value)) // up_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_WallSpeed, 0.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.wall_roughness", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.wall_speed", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallSpeed], infile_debug);
					if (Is_Equality(fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::wall_no_slip_x_up;
					else
						d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::wall_slip_x_up;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::period_x_up;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::free_x_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::pressure_x_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.pressure", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::normal_micro_flow_x_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.normal_velocity", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[2].int_value)) // down_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_WallSpeed, 0.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.wall_roughness", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.wall_speed", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed], infile_debug);
					if (Is_Equality(fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::wall_no_slip_y_down;
					else
						d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::wall_slip_y_down;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::period_y_down;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::free_y_down;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::pressure_y_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.pressure", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::normal_micro_flow_y_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.normal_velocity", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[3].int_value)) // up_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_WallSpeed, 0.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.wall_roughness", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.wall_speed", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallSpeed], infile_debug);
					if (Is_Equality(fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::wall_no_slip_y_up;
					else
						d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::wall_slip_y_up;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::period_y_up;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::free_y_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::pressure_y_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.pressure", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::normal_micro_flow_y_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.normal_velocity", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				d2q9_w.push_back(Vector3(0.0, 0.0, 0.0));	 // f0  c( 0,  0)
				d2q9_w.push_back(Vector3(1.0, 0.0, 0.0));	 // f1  c( 1,  0)
				d2q9_w.push_back(Vector3(0.0, 1.0, 0.0));	 // f2  c( 0,  1)
				d2q9_w.push_back(Vector3(-1.0, 0.0, 0.0));	 // f3  c(-1,  0)
				d2q9_w.push_back(Vector3(0.0, -1.0, 0.0));	 // f4  c( 0, -1)
				d2q9_w.push_back(Vector3(1.0, 1.0, 0.0));	 // f5  c( 1,  1)
				d2q9_w.push_back(Vector3(-1.0, 1.0, 0.0));	 // f6  c(-1,  1)
				d2q9_w.push_back(Vector3(-1.0, -1.0, 0.0));	 // f7  c(-1, -1)
				d2q9_w.push_back(Vector3(1.0, -1.0, 0.0));	 // f8  c( 1, -1)
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				w.resize(19);
				w[0] = 12.0 / 36.0;
				w[1] = 2.0 / 36.0; w[2] = 2.0 / 36.0; w[3] = 2.0 / 36.0; w[4] = 2.0 / 36.0; w[5] = 2.0 / 36.0; w[6] = 2.0 / 36.0;
				w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0; w[9] = 1.0 / 36.0; w[10] = 1.0 / 36.0; w[11] = 1.0 / 36.0; w[12] = 1.0 / 36.0;
				w[13] = 1.0 / 36.0; w[14] = 1.0 / 36.0; w[15] = 1.0 / 36.0; w[16] = 1.0 / 36.0; w[17] = 1.0 / 36.0; w[18] = 1.0 / 36.0;
				fluid_lbm_solver._boundary_condition = boundary_condition_d3q19;
				if (solid_phases.size() != 0)
					d3q19_fluid_solid_boundary = bc_funcs::d3q19_fluid_solid_boundary_Guo2002;
				else
					d3q19_fluid_solid_boundary = bc_funcs::default_domain_boundary_condition;
				// load boundary condition
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .LatticeBoltzmann.boundary_condition = (down_x,up_x,down_y,up_y,down_z,up_z) \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                        0 - Wall, 1 - Period, 2 - Free, 3 - Pressure, 4 - Normal_Flow \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                            .pressure = p0 , density0 = p0 / Cs^2 , Cs = 1 / sqrt(3) \n", InputFileReader::get_instance()->debug_file);
				string bc_key = "Postprocess.FluidDynamics.LatticeBoltzmann.boundary_condition", bc_input = "(0,0,0,0,0,0)";
				InputFileReader::get_instance()->read_string_value(bc_key, bc_input, infile_debug);
				vector<input_value> bc_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, bc_key, bc_input, infile_debug);
				switch (Fluid_Domain_Boundary_Condition(bc_value[0].int_value)) // down_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.wall_roughness", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::wall_no_slip_x_down;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::period_x_down;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::free_x_down;
					break;
				case FDBC_Pressure:
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::pressure_x_down;
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.pressure", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::DOWN_X].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::normal_micro_flow_x_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.normal_velocity", fluid_boundary_condition[Boundary::DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[1].int_value)) // up_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.wall_roughness", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::wall_no_slip_x_up;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::period_x_up;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::free_x_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::pressure_x_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.pressure", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::UP_X].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::normal_micro_flow_x_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.normal_velocity", fluid_boundary_condition[Boundary::UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[2].int_value)) // down_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.wall_roughness", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::wall_no_slip_y_down;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::period_y_down;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::free_y_down;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::pressure_y_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.pressure", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::DOWN_Y].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::normal_micro_flow_y_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.normal_velocity", fluid_boundary_condition[Boundary::DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[3].int_value)) // up_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.wall_roughness", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::wall_no_slip_y_up;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::period_y_up;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::free_y_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::pressure_y_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.pressure", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::UP_Y].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::normal_micro_flow_y_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.normal_velocity", fluid_boundary_condition[Boundary::UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[4].int_value)) // down_z
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::DOWN_Z].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Z.wall_roughness", fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::wall_no_slip_z_down;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::period_z_down;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::free_z_down;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::DOWN_Z].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::pressure_z_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Z.pressure", fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::DOWN_Z].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::normal_micro_flow_z_down;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Z.normal_velocity", fluid_boundary_condition[Boundary::DOWN_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[5].int_value)) // up_z
				{
				case FDBC_Wall:
					fluid_boundary_condition[Boundary::UP_Z].add_double(Fluid_Boundary_Property::FBP_WallRoughness, 1.0);
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Z.wall_roughness", fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_WallRoughness], infile_debug);
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::wall_no_slip_z_up;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::period_z_up;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::free_z_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Boundary::UP_Z].add_double(Fluid_Boundary_Property::FBP_DensityValue, 1.0);
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::pressure_z_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Z.pressure", fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_DensityValue], infile_debug);
					fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Boundary::UP_Z].add_double(Fluid_Boundary_Property::FBP_NormalFlowSpeed, 0.0);
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::normal_micro_flow_z_up;
					InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Z.normal_velocity", fluid_boundary_condition[Boundary::UP_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed], infile_debug);
					break;
				default:
					break;
				}
				d3q19_w.push_back(Vector3(0.0, 0.0, 0.0));	 // f0   c( 0,  0,  0)
				d3q19_w.push_back(Vector3(1.0, 0.0, 0.0));	 // f1   c( 1,  0,  0)
				d3q19_w.push_back(Vector3(-1.0, 0.0, 0.0));	 // f2   c(-1,  0,  0)
				d3q19_w.push_back(Vector3(0.0, 1.0, 0.0));	 // f3   c( 0,  1,  0)
				d3q19_w.push_back(Vector3(0.0, -1.0, 0.0));	 // f4   c( 0, -1,  0)
				d3q19_w.push_back(Vector3(0.0, 0.0, 1.0));	 // f5   c( 0,  0,  1)
				d3q19_w.push_back(Vector3(0.0, 0.0, -1.0));	 // f6   c( 0,  0, -1)
				d3q19_w.push_back(Vector3(1.0, 1.0, 0.0));	 // f7   c( 1,  1,  0)
				d3q19_w.push_back(Vector3(-1.0, 1.0, 0.0));	 // f8   c(-1,  1,  0)
				d3q19_w.push_back(Vector3(-1.0, -1.0, 0.0)); // f9   c(-1, -1,  0)
				d3q19_w.push_back(Vector3(1.0, -1.0, 0.0));	 // f10  c( 1, -1,  0)
				d3q19_w.push_back(Vector3(1.0, 0.0, 1.0));	 // f11  c( 1,  0,  1)
				d3q19_w.push_back(Vector3(-1.0, 0.0, 1.0));	 // f12  c(-1,  0,  1)
				d3q19_w.push_back(Vector3(-1.0, 0.0, -1.0)); // f13  c(-1,  0, -1)
				d3q19_w.push_back(Vector3(1.0, 0.0, -1.0));	 // f14  c( 1,  0, -1)
				d3q19_w.push_back(Vector3(0.0, 1.0, 1.0));	 // f15  c( 0,  1,  1)
				d3q19_w.push_back(Vector3(0.0, -1.0, 1.0));	 // f16  c( 0, -1,  1)
				d3q19_w.push_back(Vector3(0.0, -1.0, -1.0)); // f17  c( 0, -1, -1)
				d3q19_w.push_back(Vector3(0.0, 1.0, -1.0));	 // f18  c( 0,  1, -1)
			}
		}

		static void init_two_phase_solver(FieldStorage_forPhaseNode& phaseMesh, LBM& field_lbm_two_phase_solver) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			tau = tau_standard;
			tau_two_phase = tau_two_phase_const;
			viscosity = viscosity_two_phase;
			density = density_two_phase;
			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# tau_two_phase = Mobility / fluid_dt / Cs2 + 0.5 \n", InputFileReader::get_instance()->debug_file);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.Mobility", mobility_two_phase, infile_debug);
			_tau_two_phase = mobility_two_phase / PCT_dt / Cs2 + 0.5;
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.gas_viscosity", viscosity_gas, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.gas_density", density_gas, infile_debug);

			if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				field_lbm_two_phase_solver._boundary_condition = boundary_condition_d2q9;
			}
			else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				field_lbm_two_phase_solver._boundary_condition = boundary_condition_d3q19;
			}
		}

		static void lbm_properties_automatically_change(FieldStorage_forPhaseNode& phaseMesh) {
			PCT_dt = Solvers::get_instance()->parameters.dt;
			double cc = phaseMesh.dr / PCT_dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			_tau_const = viscosity_liquid / PCT_dt / Cs2 + 0.5;
			_tau_two_phase = mobility_two_phase / PCT_dt / Cs2 + 0.5;
		}

		void deinit() {
			fluid_boundary_condition.clear();
			solid_phases.clear();
			w.clear();
			d2q9_w.clear();
			d3q19_w.clear();
			tau = nullptr;
			viscosity = nullptr;
			density = nullptr;
			d2q9_domain_boundary_x_down = nullptr;
			d2q9_domain_boundary_x_up = nullptr;
			d2q9_domain_boundary_y_down = nullptr;
			d2q9_domain_boundary_y_up = nullptr;
			d3q19_domain_boundary_x_down = nullptr;
			d3q19_domain_boundary_x_up = nullptr;
			d3q19_domain_boundary_y_down = nullptr;
			d3q19_domain_boundary_y_up = nullptr;
			d3q19_domain_boundary_z_down = nullptr;
			d3q19_domain_boundary_z_up = nullptr;
		}
	}
}