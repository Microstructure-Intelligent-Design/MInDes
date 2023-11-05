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
		static bool is_time_auto_change = false;
		static double dt_scale = 1.0;
		static int delt_step = 100;
		static double max_scale = 1e3;
		static bool is_reduce_output = true;
		static double phi_increment_limit = 1e-3;
		static double con_increment_limit = 1e-3;
		static double temp_increment_limit = 1e-3;
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			InputFileReader::get_instance()->read_double_value("Solver.PCT.RealTime.init", Solvers::get_instance()->real_time, infile_debug);

			Solvers::get_instance()->parameters.dt = 1.0; // DEFAULT
			InputFileReader::get_instance()->read_double_value("Solver.PCT.TimeInterval.dt", time_interval, infile_debug);
			Solvers::get_instance()->parameters.dt = time_interval;

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.PCT.TimeInterval.auto_adjust = (delt_step,max_scale,is_reduce_output)\n", InputFileReader::get_instance()->debug_file);
			string time_interval_key = "Solver.PCT.TimeInterval.auto_adjust", time_interval_input = "(100,1e3,true)";
			if (InputFileReader::get_instance()->read_string_value(time_interval_key, time_interval_input, infile_debug)) {
				is_time_auto_change = true;
				vector<InputValueType> time_interval_structure;
				time_interval_structure.push_back(InputValueType::IVType_INT);
				time_interval_structure.push_back(InputValueType::IVType_DOUBLE);
				time_interval_structure.push_back(InputValueType::IVType_BOOL);
				vector<input_value> time_interval_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(time_interval_structure, time_interval_key, time_interval_input, infile_debug);
				delt_step = time_interval_value[0].int_value;
				max_scale = time_interval_value[1].double_value;
				is_reduce_output = time_interval_value[2].bool_value;
			}

			InputFileReader::get_instance()->read_double_value("Solver.PCT.phi_increment_limit", phi_increment_limit, infile_debug);
			InputFileReader::get_instance()->read_double_value("Solver.PCT.con_increment_limit", con_increment_limit, infile_debug);
			InputFileReader::get_instance()->read_double_value("Solver.PCT.temp_increment_limit", temp_increment_limit, infile_debug);

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Time interval automatically adjust !\n", LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream log, report;
			if (!is_time_auto_change)
				return report.str();

			double MAX_dphi_versus_limitPhi = Solvers::get_instance()->MAX_PHI_INCREMENT / phi_increment_limit;
			double MAX_dcon_versus_limitCon = Solvers::get_instance()->MAX_CON_INCREMENT / con_increment_limit;
			double MAX_dT_versus_limitT = Solvers::get_instance()->MAX_TEMP_INCREMENT / temp_increment_limit;
			double delt_scale = 0.05;

			if (MAX_dcon_versus_limitCon > 1.0) {
				dt_scale /= MAX_dcon_versus_limitCon + delt_scale;
				MAX_dcon_versus_limitCon = 1.0;
				if (!is_reduce_output)
					log << "> Adjust time interval for concentration evolving stability, dt_scale = " << dt_scale << endl;
			}
			else if (MAX_dphi_versus_limitPhi > 1.0) {
				dt_scale /= MAX_dphi_versus_limitPhi + delt_scale;
				MAX_dphi_versus_limitPhi = 1.0;
				if (!is_reduce_output)
					log << "> Adjust time interval for phi evolving stability, dt_scale = " << dt_scale << endl;
			}
			else if (MAX_dT_versus_limitT > 1.0) {
				dt_scale /= MAX_dT_versus_limitT + delt_scale;
				MAX_dT_versus_limitT = 1.0;
				if (!is_reduce_output)
					log << "> Adjust time interval for temperature evolving stability, dt_scale = " << dt_scale << endl;
			}
			else if (Solvers::get_instance()->current_istep % delt_step == 0 && (dt_scale * 1.1) <= max_scale) {
				dt_scale *= 1.1;
				if (!is_reduce_output)
					log << "> Adjust time interval for fields evolving quickly, dt_scale = " << dt_scale << endl;
			}
			else if (Solvers::get_instance()->current_istep % delt_step == 0 && (dt_scale * 1.1) > max_scale) {
				dt_scale = max_scale;
				if (!is_reduce_output)
					log << "> Adjust time interval for fields evolving quickly, dt_scale = " << dt_scale << endl;
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(log.str(), LOG_FILE_NAME);
			Solvers::get_instance()->parameters.dt = time_interval * dt_scale;

			if (Solvers::get_instance()->current_istep % Solvers::get_instance()->parameters.screen_output_step == 0)
				report << "> Time interval is automatically adjusting, time interval = " << Solvers::get_instance()->parameters.dt << endl;

			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			
		}
	}
}