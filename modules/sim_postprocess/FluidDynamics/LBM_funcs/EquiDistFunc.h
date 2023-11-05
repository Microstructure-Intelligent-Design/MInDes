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
	namespace lbm_equilibrium_distribution_function {
		static vector<Vector3> d2q9_w;
		static vector<Vector3> d3q19_w;
		static vector<double> w;
		static double Cs2 = 1.0 / 3.0;
		static double Cs4 = 1.0 / 9.0;
		static int index_h_liang;
		namespace equilibrium_distribution_functions {
			static double f_eq_a_d2q9_standard(int INDEX_A, double p_macro, double f_macro, Vector3& U) {
				double CU = d2q9_w[INDEX_A] * U;
				return w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double f_eq_a_d3q19_standard(int INDEX_A, double p_macro, double f_macro, Vector3& U) {
				double CU = d3q19_w[INDEX_A] * U;
				return w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double s_a_d2q9_two_phase_flow(int INDEX_A, Vector3& U) {
				double CU = d2q9_w[INDEX_A] * U;
				return w[INDEX_A] * (CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double s_a_d3q19_two_phase_flow(int INDEX_A, Vector3& U) {
				double CU = d3q19_w[INDEX_A] * U;
				return w[INDEX_A] * (CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double f_eq_a_d2q9_two_phase_flow(int INDEX_A, double p_macro, double f_macro, Vector3& U) {
				if(INDEX_A == 0)
					return p_macro / Cs2 * (w[INDEX_A] - 1.0) + f_macro * s_a_d2q9_two_phase_flow(INDEX_A, U);
				else
					return p_macro / Cs2 * w[INDEX_A] + f_macro * s_a_d2q9_two_phase_flow(INDEX_A, U);
			}
			static double f_eq_a_d3q19_two_phase_flow(int INDEX_A, double p_macro, double f_macro, Vector3& U) {
				if (INDEX_A == 0)
					return p_macro / Cs2 * (w[INDEX_A] - 1.0) + f_macro * s_a_d3q19_two_phase_flow(INDEX_A, U);
				else
					return p_macro / Cs2 * w[INDEX_A] + f_macro * s_a_d3q19_two_phase_flow(INDEX_A, U);
			}
			static double f_eq_d2q9_two_phase(int INDEX_A, double p_macro, double f_macro, Vector3& U) {
				double CU = d2q9_w[INDEX_A] * U;
				return w[INDEX_A] * f_macro * (1 + CU / Cs2);
			}
			static double f_eq_d3q19_two_phase(int INDEX_A, double p_macro, double f_macro, Vector3& U) {
				double CU = d3q19_w[INDEX_A] * U;
				return w[INDEX_A] * f_macro * (1 + CU / Cs2);
			}

		}

		static double (*f_eq_i)(int INDEX_i, double p_macro, double f_macro, Vector3& U);

		static double (*f_eq_two_phase_i)(int INDEX_i, double p_macro, double f_macro, Vector3& U);

		static void lbm_properties_automatically_change(FieldStorage_forPhaseNode& phaseMesh) {
			double cc = phaseMesh.dr / Solvers::get_instance()->parameters.dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh, LBM& fluid_lbm_solver) {
			double cc = phaseMesh.dr / Solvers::get_instance()->parameters.dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				f_eq_i = equilibrium_distribution_functions::f_eq_a_d2q9_standard;
				w.resize(9);
				w[0] = 4.0 / 9.0;
				w[1] = 1.0 / 9.0; w[2] = 1.0 / 9.0; w[3] = 1.0 / 9.0; w[4] = 1.0 / 9.0;
				w[5] = 1.0 / 36.0; w[6] = 1.0 / 36.0; w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0;
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
				f_eq_i = equilibrium_distribution_functions::f_eq_a_d3q19_standard;
				w.resize(19);
				w[0] = 12.0 / 36.0;
				w[1] = 2.0 / 36.0; w[2] = 2.0 / 36.0; w[3] = 2.0 / 36.0; w[4] = 2.0 / 36.0; w[5] = 2.0 / 36.0; w[6] = 2.0 / 36.0;
				w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0; w[9] = 1.0 / 36.0; w[10] = 1.0 / 36.0; w[11] = 1.0 / 36.0; w[12] = 1.0 / 36.0;
				w[13] = 1.0 / 36.0; w[14] = 1.0 / 36.0; w[15] = 1.0 / 36.0; w[16] = 1.0 / 36.0; w[17] = 1.0 / 36.0; w[18] = 1.0 / 36.0;
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
		static void init_two_phase(FieldStorage_forPhaseNode& phaseMesh, LBM& field_lbm_two_phase_solver) {
			if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				f_eq_two_phase_i = equilibrium_distribution_functions::f_eq_d2q9_two_phase;
				f_eq_i = equilibrium_distribution_functions::f_eq_a_d2q9_two_phase_flow;
			}
			else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				f_eq_two_phase_i = equilibrium_distribution_functions::f_eq_d3q19_two_phase;
				f_eq_i = equilibrium_distribution_functions::f_eq_a_d3q19_two_phase_flow;
			}
		}
	}
}