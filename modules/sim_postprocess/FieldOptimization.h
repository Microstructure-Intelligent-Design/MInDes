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
using namespace std;
namespace pf {

	namespace field_optimization {
		// - optimization grand-potential field
		static double is_optimize_grand = false;
		static int optimize_grand_steps = 1000;

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			//Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Pretreatment !\n", LOG_FILE_NAME);
		}

		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_optimize_grand) {
				//optimization_grand_potential_field(phaseMesh);
			}
		}

		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string str = "";
			return str;
		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			
		}

		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}

		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}