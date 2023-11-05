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

namespace pf {
	namespace auto_input {
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh, string input_file_path, int input_file_lines_buff = 1000) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			InputFileReader::get_instance()->auto_read_infile(input_file_path, input_file_lines_buff);
			InputFileReader::get_instance()->auto_read_bool_value("Solver.Loop.stop", Solvers::get_instance()->is_stop_loop, infile_debug);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh, string input_file_path, int input_file_lines_buff = 1000) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			InputFileReader::get_instance()->auto_read_infile(input_file_path, input_file_lines_buff);
			stringstream report;
			InputFileReader::get_instance()->auto_read_bool_value("Solver.Loop.stop", Solvers::get_instance()->is_stop_loop, false);
			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			if (Solvers::get_instance()->is_stop_loop) {
				stringstream out;
				out << "> The main loop is terminated by key Solver.Loop.stop at simulation step: " << Solvers::get_instance()->current_istep << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(out.str(), LOG_FILE_NAME);
			}
		}
	}
}