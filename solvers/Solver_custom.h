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
#include "base.h"
using namespace std;
namespace pf {
	namespace custom_solver {
		enum CUSTOM_SOLVER { SOLVER_ALLEN_CAHN = 10000, SOLVER_CAHN_HILLIARD = 11000, BUFF_FOR_SOLVER_CAHN_HILLIARD = 12000 };
		namespace allenCahnEquationSolver {
			// grain_index from 0 to grain_number - 1
			static double dF_dphi_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
				return 0.0;
			}
			static double L_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
				return 0.0;
			}
			static double Source_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
				return 0.0;
			}
			static void BoundaryCondition(pf::PhaseNode& node, int comp_index, int comp_start_index) {
				return;
			}
		}

		// partial(phi_i) / partial(t) = - Mobility_i * variation(F) / variation(phi_i) + Source
		class CustomAllenCahnSolver {
		public:
			CustomAllenCahnSolver() {};
			CustomAllenCahnSolver(FieldStorage_forPhaseNode& _phaseMesh, int _grains_number = 0, int _grains_start_index = SOLVER_ALLEN_CAHN, string _solver_name = "AllenCahnSolver") {
				init(_phaseMesh, _grains_number, _grains_start_index, _solver_name);
			}
			~CustomAllenCahnSolver() {
				clear();
			};
			void init(FieldStorage_forPhaseNode& _phaseMesh, int _grains_number = 0, int _grains_start_index = SOLVER_ALLEN_CAHN, string _solver_name = "AllenCahnSolver") {
				solver_name = _solver_name;
				phaseMesh = &_phaseMesh;
				grains_number = _grains_number;
				grains_start_index = _grains_start_index;
				dF_dphi = allenCahnEquationSolver::dF_dphi_cal;
				L = allenCahnEquationSolver::L_cal;
				Source = allenCahnEquationSolver::Source_cal;
				BoundaryCondition = allenCahnEquationSolver::BoundaryCondition;
			}
			void define_funcs_for_AC_solver(double(*dF_dphi_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::dF_dphi_cal,
				double(*L_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::L_cal,
				double(*Source_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::Source_cal,
				void(*BoundaryCondition_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::BoundaryCondition) {
				dF_dphi = dF_dphi_cal;
				L = L_cal;
				Source = Source_cal;
				BoundaryCondition = BoundaryCondition_cal;
			}
			void init_mesh() {
				for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
					for (int grain = 0; grain < grains_number; grain++) {
						node->customValues.add_double((grain + grains_start_index), 0.0);
						node->customValues.add_double(-(grain + grains_start_index), 0.0);
					}
				}
			}
			string print_model() {
				stringstream rep;
				rep << "External Allen Cahn Model : dphi_dt = - L * df_dphi + S " << endl;
				return rep.str();
			}
			// retuan MAX_PHI_VARIATION;
			double solve_one_step(double dt, bool adjust_phi_0_1 = false);
			string solver_name;
			int grains_number;
			int grains_start_index;
			double(*dF_dphi)(pf::PhaseNode&, int, int);
			double(*L)(pf::PhaseNode&, int, int);
			double(*Source)(pf::PhaseNode&, int, int);
			void(*BoundaryCondition)(pf::PhaseNode&, int, int);
		private:
			FieldStorage_forPhaseNode* phaseMesh;
			void clear() {
				phaseMesh = nullptr;
				dF_dphi = nullptr;
				L = nullptr;
				Source = nullptr;
				BoundaryCondition = nullptr;
			}
		};

		namespace cahnHilliardEquationSolver {
			// grain_index from 0 to grain_number - 1
			static double dF_dc_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
				return 0.0;
			}
			static double Mobility_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
				return 0.0;
			}
			static double Source_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
				return 0.0;
			}
			static void BoundaryCondition(pf::PhaseNode& node, int comp_index, int comp_start_index) {
				return;
			}
		}

		// partial(c_i) / partial(t) = delt(M_i * delt(variation(F) / variation(c_i))) + Source
		class CustomCahnHilliardSolver {
		public:
			CustomCahnHilliardSolver() { };
			CustomCahnHilliardSolver(FieldStorage_forPhaseNode& _phaseMesh, int _components_number = 0, int _comps_start_index = SOLVER_CAHN_HILLIARD, string _solver_name = "CahnHilliardSolver") {
				init(_phaseMesh, _components_number, _comps_start_index, _solver_name);
			}
			~CustomCahnHilliardSolver() {
				clear();
			};
			void init(FieldStorage_forPhaseNode& _phaseMesh, int _components_number = 0, int _comps_start_index = SOLVER_CAHN_HILLIARD, string _solver_name = "CahnHilliardSolver") {
				solver_name = _solver_name;
				phaseMesh = &_phaseMesh;
				components_number = _components_number;
				comps_start_index = _comps_start_index;
				dF_dc = cahnHilliardEquationSolver::dF_dc_cal;
				Mobility = cahnHilliardEquationSolver::Mobility_cal;
				Source = cahnHilliardEquationSolver::Source_cal;
				BoundaryCondition = cahnHilliardEquationSolver::BoundaryCondition;
			}
			void init_mesh() {
				for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
					for (int comp = 0; comp < components_number; comp++) {
						node->customValues.add_double(comp + comps_start_index, 0.0);
						node->customValues.add_double(-(comp + comps_start_index), 0.0);
						node->customValues.add_double(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD, 0.0);
						node->customValues.add_double(-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD), 0.0);
					}
				}
			}
			void define_funcs_for_AC_solver(double(*dF_dc_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::dF_dc_cal,
				double(*Mobility_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::Mobility_cal,
				double(*Source_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::Source_cal,
				void(*BoundaryCondition_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::BoundaryCondition) {
				dF_dc = dF_dc_cal;
				Mobility = Mobility_cal;
				Source = Source_cal;
				BoundaryCondition = BoundaryCondition_cal;
			}
			string print_model() {
				stringstream rep;
				rep << "External Cahn Hilliard Model : dc_dt = \\labla{ M * \\labla( df_dc ) } + S " << endl;
				return rep.str();
			}
			// return MAX_COMP_VARIATION
			double solve_one_step(double dt, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT, bool adjust_phi_0_1 = false);
			string solver_name;
			int components_number;
			int comps_start_index;
		private:
			FieldStorage_forPhaseNode* phaseMesh;
			double(*dF_dc)(pf::PhaseNode&, int, int);
			double(*Mobility)(pf::PhaseNode&, int, int);
			double(*Source)(pf::PhaseNode&, int, int);
			void(*BoundaryCondition)(pf::PhaseNode&, int, int);
			void clear() {
				phaseMesh = nullptr;
				dF_dc = nullptr;
				Mobility = nullptr;
				Source = nullptr;
				BoundaryCondition = nullptr;
			}
		};
	}
}