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

	namespace Funcs_AllenCahnSolver {
		//
		static void dfint_dphi(pf::PhaseNode& node, bool normalize_phi) {
			return;
		}
		static double dfbulk_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		//
		static double Source_i(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		static double Lij(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return 0.0;
		}
		static double Source_ij(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return 0.0;
		}
		// boundary condition
		static void boundary_condition(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
	}

	// partial(A_i) / partial(t) = sum_j { - L_ij * variation(F_int) / variation(A_j) + Source_ij} + Source_i
	class AllenCahnSolver {
	public:
		AllenCahnSolver() {};
		AllenCahnSolver(FieldStorage_forPhaseNode& _phaseMesh) {
			init(_phaseMesh);
		}
		~AllenCahnSolver() {
			clear();
		};
		void init(FieldStorage_forPhaseNode& _phaseMesh) {
			phaseMesh = &_phaseMesh;
			dfint_dphi = Funcs_AllenCahnSolver::dfint_dphi;
			dfbulk_dphi = Funcs_AllenCahnSolver::dfbulk_dphi;
			Lij = Funcs_AllenCahnSolver::Lij;
			Source_ij = Funcs_AllenCahnSolver::Source_ij;
			Source_i = Funcs_AllenCahnSolver::Source_i;
			Boundary_Condition = Funcs_AllenCahnSolver::boundary_condition;
		}
		void define_funcs_for_pair_wise_equation(double(*_Lij)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::Lij,
			void(*_dfint_dphi)(pf::PhaseNode&, bool) = Funcs_AllenCahnSolver::dfint_dphi,
			double(*_dfbulk_dphi)(pf::PhaseNode&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::dfbulk_dphi,
			double(*_Source_ij)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::Source_ij,
			void(*_Boundary_Condition)(pf::PhaseNode&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::boundary_condition) {
			Lij = _Lij;
			dfint_dphi = _dfint_dphi;
			dfbulk_dphi = _dfbulk_dphi;
			Source_ij = _Source_ij;
			Boundary_Condition = _Boundary_Condition;
		}
		void define_funcs_for_normal_equation(double(*_Lij)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::Lij,
			void(*_dfint_dphi)(pf::PhaseNode&, bool) = Funcs_AllenCahnSolver::dfint_dphi,
			double(*_dfbulk_dphi)(pf::PhaseNode&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::dfbulk_dphi,
			double(*_Source_i)(pf::PhaseNode&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::Source_i,
			void(*_Boundary_Condition)(pf::PhaseNode&, pf::PhaseEntry&) = Funcs_AllenCahnSolver::boundary_condition) {
			Lij = _Lij;
			dfint_dphi = _dfint_dphi;
			dfbulk_dphi = _dfbulk_dphi;
			Source_i = _Source_i;
			Boundary_Condition = _Boundary_Condition;
		}

		void init_phi(bool normalize_phi = false);

		void init_phi_pair_wise(bool normalize_phi = false);

		void pre_calculation_phi_normal(double dt, bool normalize_phi = false, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT);

		void pre_calculation_phi_pair_wise_normal(double dt, bool normalize_phi = false, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT);

		// MAX_PHI_INCREMENT
		double solve_phi_normal(double dt, bool normalize_phi = false);

		double solve_phi_pair_wise_normal(double dt, bool normalize_phi = false);

		string print_phi_model_normal() {
			stringstream rep;
			rep << "Phi Standard Model : dPa_dt = - La * (dfint_dPa + dfbulk_dPa) + Sa" << endl;
			return rep.str();
		}

		string print_phi_model_pairwise() {
			stringstream rep;
			rep << "Phi Pair-WIse Model : dPa_dt = \\sum_b[ Lab * ( dfint_dPb - dfint_dPa + dfbulk_dPb - dfbulk_dPa ) + Sab ]" << endl;
			return rep.str();
		}

		void(*dfint_dphi)(pf::PhaseNode&, bool);
		double(*dfbulk_dphi)(pf::PhaseNode&, pf::PhaseEntry&);
		double(*Source_i)(pf::PhaseNode&, pf::PhaseEntry&);
		double(*Lij)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&);
		double(*Source_ij)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&);
		void(*Boundary_Condition)(pf::PhaseNode&, pf::PhaseEntry&);
	private:
		FieldStorage_forPhaseNode* phaseMesh;
		void clear() {
			dfint_dphi = nullptr;
			dfbulk_dphi = nullptr;
			Source_i = nullptr;
			Lij = nullptr;
			Source_ij = nullptr;
			Boundary_Condition = nullptr;
		}
	};

}