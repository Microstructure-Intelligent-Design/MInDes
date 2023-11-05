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
	namespace solver_init {
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			InputFileReader::get_instance()->read_int_value("Solver.Loop.begin_step", Solvers::get_instance()->parameters.begin_step, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Loop.end_step", Solvers::get_instance()->parameters.end_step, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Loop.screen_loop_step", Solvers::get_instance()->parameters.screen_loop_step, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Loop.screen_output_step", Solvers::get_instance()->parameters.screen_output_step, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Loop.vts_output_step", Solvers::get_instance()->parameters.vts_output_step, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Loop.data_output_step", Solvers::get_instance()->parameters.data_output_step, infile_debug);

			int difference_method = pf::DifferenceMethod::FIVE_POINT;
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.difference_method : 0 - FIVE_POINT , 1 - NINE_POINT \n", InputFileReader::get_instance()->debug_file);
			InputFileReader::get_instance()->read_int_value("Solver.difference_method", difference_method, infile_debug);
			Solvers::get_instance()->parameters.Difference_Method = pf::DifferenceMethod(difference_method);

			InputFileReader::get_instance()->read_bool_value("Solver.Phi.is_normalize", Solvers::get_instance()->parameters.is_Normalize_Phi, infile_debug);

			InputFileReader::get_instance()->read_bool_value("Solver.Con.is_normalize", Solvers::get_instance()->parameters.is_Normalize_Con, infile_debug);

			InputFileReader::get_instance()->read_int_value("Solver.Parallel.openmp_thread_counts", Solvers::get_instance()->parameters.OpenMP_Thread_Counts, infile_debug);
			stringstream log;
			// protect computer system
			if (Solvers::get_instance()->parameters.OpenMP_Thread_Counts >= omp_get_num_procs())
				Solvers::get_instance()->parameters.OpenMP_Thread_Counts = omp_get_num_procs() - 1;
			else if (Solvers::get_instance()->parameters.OpenMP_Thread_Counts < 1)
				Solvers::get_instance()->parameters.OpenMP_Thread_Counts = 1;
			log << "> Simulation OMP threads number is " << Solvers::get_instance()->parameters.OpenMP_Thread_Counts << endl;
			omp_set_num_threads(Solvers::get_instance()->parameters.OpenMP_Thread_Counts);
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(log.str(), LOG_FILE_NAME);

			Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_CON] = 1;
			Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP] = 1;
			InputFileReader::get_instance()->read_int_value("Solver.Loop.Iterate.Con", Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_CON], infile_debug);
			if (Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_CON] < 1)
				Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_CON] = 1;
			InputFileReader::get_instance()->read_int_value("Solver.Loop.Iterate.Temp", Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP], infile_debug);
			if (Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP] < 1)
				Solvers::get_instance()->parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP] = 1;

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.Comps = (c0, c1, c2, ... ) \n", InputFileReader::get_instance()->debug_file);
			string component_key = "Solver.Comps", component_input = "()";
			InputFileReader::get_instance()->read_string_value(component_key, component_input, infile_debug);
			vector<input_value> component_value;
			component_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, component_key, component_input, infile_debug);
			for (int i = 0; i < component_value.size(); i++)
				Solvers::get_instance()->parameters.Components.add_nodeEntry(i, component_value[i].string_value);

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.Phases = {[(phase0),(c0, c1, ... )], [(phase1),(c0, c1, ... )], ... } \n", InputFileReader::get_instance()->debug_file);
			string phase_key = "Solver.Phases", phase_input = "{}";
			InputFileReader::get_instance()->read_string_value(phase_key, phase_input, infile_debug);
			vector<InputValueType> phase_structure; phase_structure.push_back(InputValueType::IVType_STRING); phase_structure.push_back(InputValueType::IVType_STRING);
			vector<vector<vector<input_value>>> phase_value;
			phase_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(phase_structure, phase_key, phase_input, infile_debug);
			for (int i = 0; i < phase_value.size(); i++) {
				Solvers::get_instance()->parameters.Phases.add_Phase(i, phase_value[i][0][0].string_value);
				for (int j = 0; j < phase_value[i][1].size(); j++)
					Solvers::get_instance()->parameters.Phases[i].x.add_nodeEntry(Solvers::get_instance()->parameters.Components[phase_value[i][1][j].string_value].index, phase_value[i][1][j].string_value);
			}


			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.GrainsOrientations = {[(phi_index_0, phi_index_2, ... ),(rotation_angle_1, rotation_angle_2, rotation_angle_3)],  ... } \n", InputFileReader::get_instance()->debug_file);
			string grains_key = "Solver.GrainsOrientations", grains_input = "{}";
			InputFileReader::get_instance()->read_string_value(grains_key, grains_input, infile_debug);
			vector<InputValueType> grains_structure; grains_structure.push_back(InputValueType::IVType_INT); grains_structure.push_back(InputValueType::IVType_DOUBLE);
			vector<vector<vector<input_value>>> grains_value;
			grains_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(grains_structure, grains_key, grains_input, infile_debug);
			for (int i = 0; i < grains_value.size(); i++) {
				vector<int> phis;
				for (int j = 0; j < grains_value[i][0].size(); j++)
					phis.push_back(grains_value[i][0][j].int_value);
				Solvers::get_instance()->parameters.Grains.set_phi_orientation(phis, Vector3(AngleToRadians(grains_value[i][1][0].double_value),
					AngleToRadians(grains_value[i][1][1].double_value), AngleToRadians(grains_value[i][1][2].double_value)));
			}
			if (infile_debug) {
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.GrainsOrientations.rotation_gauge = 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                            6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n", InputFileReader::get_instance()->debug_file);
			}
			int rotation_gauge = RotationGauge::RG_ZXZ;
			InputFileReader::get_instance()->read_int_value("Solver.GrainsOrientations.rotation_gauge", rotation_gauge, infile_debug);

			Solvers::get_instance()->Phi_Solver_CH.diff_method = Solvers::get_instance()->parameters.Difference_Method;
			Solvers::get_instance()->C_Solver.diff_method = Solvers::get_instance()->parameters.Difference_Method;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			Solvers::get_instance()->parameters.Components.clear();
			Solvers::get_instance()->parameters.Phases.clear();
		}
	}
}