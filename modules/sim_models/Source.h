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
#include "../Base.h"
#include "Source/BulkReaction.h"
#include "Source/Convection.h"
#include "Source/InterfaceReaction.h"

namespace pf {
	namespace extern_source {
		// Phi
		static double (*convection_a)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*convection_ab)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta);
		static double (*reaction_a)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*interface_reaction_ab)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta);
		static double Source_a(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double source = 0.0;
			source += convection_a(node, phase);
			source += reaction_a(node, phase);
			return source;
		}
		static double Source_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			double source = 0.0;
			source += convection_ab(node, alpha, beta);
			source += interface_reaction_ab(node, alpha, beta);
			return source;
		}
		// Concentration
		static void (*convection_A)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static void (*reaction_A)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static void (*interface_reaction_AB)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double abs_dPaPb);
		static double (*convection_i)(pf::PhaseNode& node, int con_i);
		static double (*reaction_i)(pf::PhaseNode& node, int con_i);
		static double (*interface_reaction_i)(pf::PhaseNode& node, int con_i);
		static void Source_A(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			convection_A(node, phase);
			reaction_A(node, phase);
			return;
		}
		static void Source_AB(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double abs_dPaPb) {
			// AbsGradPhi_Phi * ( standard_reaction_flux_on_interface )
			interface_reaction_AB(node, alpha, beta, abs_dPaPb);
			return;
		}
		static double Source_i(pf::PhaseNode& node, int con_i) {
			double source = 0.0;
			source += convection_i(node, con_i);
			source += reaction_i(node, con_i);
			return source;
		}

		// Temperature
		static double (*convection_T)(pf::PhaseNode& node);
		static double (*reaction_T)(pf::PhaseNode& node);
		static double Source_T(pf::PhaseNode& node) {
			double source = 0.0;
			source += convection_T(node);
			source += reaction_T(node);
			return source;
		}
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bulk_reaction::init(phaseMesh);
			reaction_a = bulk_reaction::reaction_a;
			reaction_A = bulk_reaction::reaction_A;
			reaction_i = bulk_reaction::reaction_i;
			reaction_T = bulk_reaction::reaction_T;

			convection::init(phaseMesh);
			convection_a = convection::convection_a;
			convection_ab = convection::convection_ab;
			convection_A = convection::convection_A;
			convection_i = convection::convection_i;
			convection_T = convection::convection_T;

			interface_reaction::init(phaseMesh);
			interface_reaction_ab = interface_reaction::interface_reaction_ab;
			interface_reaction_AB = interface_reaction::interface_reaction_AB;
			interface_reaction_i = interface_reaction::interface_reaction_i;

			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				Solvers::get_instance()->C_Solver.Source_A = Source_A;
				Solvers::get_instance()->C_Solver.Source_AB = Source_AB;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				Solvers::get_instance()->C_Solver.int_flux = interface_reaction_i;
				Solvers::get_instance()->C_Solver.Source = Source_i;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				Solvers::get_instance()->C_Solver.int_flux = interface_reaction_i;
				Solvers::get_instance()->C_Solver.Source = Source_i;
			}

			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				// calculated when alpha._flag and beta._flag are pf_INTERFACE
				Solvers::get_instance()->Phi_Solver_AC.Source_ij = Source_ab;
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
				Solvers::get_instance()->Phi_Solver_AC.Source_i = Source_a;
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				Solvers::get_instance()->Phi_Solver_CH.Source_a = Source_a;
			}

			if (Solvers::get_instance()->parameters.TempEType == pf::TemperatureEquationType::TType_Standard) {
				Solvers::get_instance()->T_Solver.Source = Source_T;
			}
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			bulk_reaction::exec_pre(phaseMesh);
			convection::exec_pre(phaseMesh);
			interface_reaction::exec_pre(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			report += bulk_reaction::exec_loop(phaseMesh);
			report += convection::exec_loop(phaseMesh);
			report += interface_reaction::exec_loop(phaseMesh);
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			bulk_reaction::deinit(phaseMesh);
			convection::deinit(phaseMesh);
			interface_reaction::deinit(phaseMesh);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			bulk_reaction::write_scalar(fout, phaseMesh);
			convection::write_scalar(fout, phaseMesh);
			interface_reaction::write_scalar(fout, phaseMesh);
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			bulk_reaction::write_vec3(fout, phaseMesh);
			convection::write_vec3(fout, phaseMesh);
			interface_reaction::write_vec3(fout, phaseMesh);
		}
	}
}