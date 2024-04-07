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
	static int poison_protect_index = -3333;
	namespace poissonEquationSolver {
		// - explicit
		static void rhs_cal(pf::PhaseNode& node, int rhs_index) {
			node.customValues[rhs_index] = 0.0;
		}
		static void lhs_cal(pf::PhaseNode& node, int lhs_index) {
			node.customValues[lhs_index] = 0.0;
		}
		static void boundary(pf::PhaseNode& node, int r_index) {
			return;
		}
	}

	class PoissonEquationSolver_Explicit
	{
	public:
		PoissonEquationSolver_Explicit() {};
		PoissonEquationSolver_Explicit(FieldStorage_forPhaseNode& _phaseMesh, int LHS_index, int R_index, int RHS_index, string _solver_name = "PoissonEquationSolver_Explicit") {
			init_field(_phaseMesh, LHS_index, R_index, RHS_index, _solver_name);
		}
		~PoissonEquationSolver_Explicit() {
			clear();
		};
		void init_field(FieldStorage_forPhaseNode& _phaseMesh, int LHS_index, int R_index, int RHS_index, string _solver_name = "PoissonEquationSolver_Explicit") {
			phaseMesh = &_phaseMesh;
			solver_name = _solver_name;
			LHS_INDEX = LHS_index;
			RHS_INDEX = RHS_index;
			R_INDEX = R_index;
			setRHS_value = poissonEquationSolver::rhs_cal;
			setLHS_value = poissonEquationSolver::lhs_cal;
			boundary = poissonEquationSolver::boundary;
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				node->customValues.add_double(LHS_INDEX, 0.0);
				node->customValues.add_double(RHS_INDEX, 0.0);
				node->customValues.add_double(R_INDEX, 0.0);
			}
		}
		void set_field_variable(double LHS_value, double R_value, double RHS_value) {
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				node->customValues.add_double(LHS_INDEX, LHS_value);
				node->customValues.add_double(R_INDEX, R_value);
				node->customValues.add_double(RHS_INDEX, RHS_value);
			}
		}
		int solve_whole_domain(double accuracy, int MAXIterations, bool debug_solver = false, int output_step = 1000);
		int solve_whole_domain_with_average_boundary_condition(double accuracy, int MAXIterations, double average_value = 0.0, bool debug_solver = false, int output_step = 1000);
		int solve_inside_phases_region(double accuracy, int MAXIterations, vector<int> phaseIndexes, double phi_threshold = Phi_Num_Cut_Off, bool debug_solver = false, int output_step = 1000);
		int solve_outside_phases_region(double accuracy, int MAXIterations, vector<int> phaseIndexes, double phi_threshold = Phi_Num_Cut_Off, bool debug_solver = false, int output_step = 1000);
		void set_RHS_calfunc(void(*RHS_cal)(pf::PhaseNode&, int)) {
			setRHS_value = RHS_cal;
		}
		void set_LHS_calfunc(void(*LHS_cal)(pf::PhaseNode&, int)) {
			setLHS_value = LHS_cal;
		}
		void set_BoundaryCondition_calfunc(void(*Boundary)(pf::PhaseNode&, int)) {
			boundary = Boundary;
		}

		string solver_name;
	private:
		int LHS_INDEX, R_INDEX, RHS_INDEX;
		FieldStorage_forPhaseNode* phaseMesh;
		void(*setRHS_value)(pf::PhaseNode&, int);
		void(*setLHS_value)(pf::PhaseNode&, int);
		void(*boundary)(pf::PhaseNode&, int);
		void set_RHS_value();
		void set_LHS_value();
		void clear() {
			phaseMesh = nullptr;
			setRHS_value = nullptr;
			setLHS_value = nullptr;
			boundary = nullptr;
		}
	};

}