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

	namespace temperatureSolver {
		static double D_temp(pf::PhaseNode& node) {
			return 0.0;
		}
		static double Source(pf::PhaseNode& node) {
			return 0.0;
		}
		static void BoundaryCondition(pf::PhaseNode& node) {
			return;
		}
	}

	class TemperatureSolver
	{
	public:
		TemperatureSolver() {};
		TemperatureSolver(FieldStorage_forPhaseNode& _phaseMesh) {
			init(_phaseMesh);
		};
		~TemperatureSolver() {
			free();
		};

		void init(FieldStorage_forPhaseNode& _phaseMesh) {
			phaseMesh = &_phaseMesh;
			D_temp = temperatureSolver::D_temp;
			Source = temperatureSolver::Source;
			BoundaryCondition = temperatureSolver::BoundaryCondition;
		}

		void define_funcs_for_pair_wise_equation(double(*_D_temp)(pf::PhaseNode&) = temperatureSolver::D_temp,
			double(*_Source)(pf::PhaseNode&) = temperatureSolver::Source,
			void(*_BoundaryCondition)(pf::PhaseNode&) = temperatureSolver::BoundaryCondition) {
			D_temp = _D_temp;
			Source = _Source;
			BoundaryCondition = _BoundaryCondition;
		}

		void pre_calculation_temperature(DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT);

		// MAX_T_INCREMENT
		double solve_temperature(double dt);

		string print_temperature_model() {
			stringstream rep;
			rep << "Temperature Model : dT_dt = \\labla{ D * \\labla( T ) } + S" << endl;
			return rep.str();
		}

		void free() {
			phaseMesh = nullptr;
			D_temp = nullptr;
			Source = nullptr;
			BoundaryCondition = nullptr;
		}
		double(*D_temp)(pf::PhaseNode&);
		double(*Source)(pf::PhaseNode&);
		void(*BoundaryCondition)(pf::PhaseNode&);
	private:
		FieldStorage_forPhaseNode* phaseMesh;
	};
}