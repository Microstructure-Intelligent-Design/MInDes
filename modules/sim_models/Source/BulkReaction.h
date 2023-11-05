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
#include "../../Base.h"

namespace pf {
	namespace bulk_reaction {
		// Phi
		static double reaction_a_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		static double (*reaction_a)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function

		// Concentration
		static void reaction_A_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static void (*reaction_A)(pf::PhaseNode& node, pf::PhaseEntry& phase);   // main function

		static double reaction_i_none(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		static double (*reaction_i)(pf::PhaseNode& node, int con_i);   // main function

		// Temperature
		static double reaction_T_none(pf::PhaseNode& node) {
			return 0.0;
		}
		static double (*reaction_T)(pf::PhaseNode& node);   // main function

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			reaction_a = reaction_a_none;
			reaction_A = reaction_A_none;
			reaction_i = reaction_i_none;
			reaction_T = reaction_T_none;
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Bulk Reaction !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}