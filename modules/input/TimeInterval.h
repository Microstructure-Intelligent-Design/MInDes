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
	namespace time_interval {
		static double time_interval = 1.0;
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			InputFileReader::get_instance()->read_double_value("Solver.PCT.RealTime.init", Solvers::get_instance()->real_time, infile_debug);

			Solvers::get_instance()->parameters.dt = 1.0; // DEFAULT
			InputFileReader::get_instance()->read_double_value("Solver.PCT.TimeInterval.dt", time_interval, infile_debug);
			Solvers::get_instance()->parameters.dt = time_interval;

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Time interval automatically adjust !\n", LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream log, report;
			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			
		}
	}
}