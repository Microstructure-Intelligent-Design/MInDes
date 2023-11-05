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
#include "Base.h"
#include "input/SolverInit.h"
#include "input/TimeInterval.h"
#include "input/MeshStructure.h"
#include "input/AutoInput.h"

namespace pf {
	namespace input_manager {
		static string _input_file_path = "";
		static int _input_file_lines_buff = 1000;
		static char _split = ' ';
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			// debug input file
			bool input_file_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", input_file_debug, true);
			if (input_file_debug) {
				InputFileReader::get_instance()->debug_infile_and_valid_words();
				InputFileReader::get_instance()->debug_custom_variavle_and_funcs();
			}
			if (input_file_debug) {
				stringstream out;
				out << "######################### define custom variables and functions #########################" << endl;
				out << "# Define.Var = A,0.1" << endl;
				out << "# Define.Func = ABC@{[(A*PHI<1>)]}@" << endl;
//				out << "# default field variables: \"PHI\", \"dPHI_dt\", \"lap_PHI\", \"PHI_X\", \"dPHI_X_dt\", \"X\", \"dX_dt\"" << endl;
//				out << "#                          \"T\", \"dT_dt\", \"lap_T\", \"P\", \"dP_dt\", \"lap_P\", \"PHI_P\", \"dPHI_P_dt\", \"lap_PHI_P\"" << endl;
				out << "# default functions      : \"pow\", \"sqrt\", \"abs\", \"exp\", \"ln\", \"log\", \"sin\", \"cos\", \"tan\", \"asin\", \"acos\", \"atan\"" << endl;
				out << "#########################################################################################" << endl;
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(out.str(), InputFileReader::get_instance()->debug_file);
			}
			pf::solver_init::init(phaseMesh);
			pf::time_interval::init(phaseMesh);
			pf::mesh_structure::init(phaseMesh);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			pf::auto_input::exec_pre(phaseMesh, _input_file_path, _input_file_lines_buff);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			report += pf::time_interval::exec_loop(phaseMesh);
			report += pf::auto_input::exec_loop(phaseMesh, _input_file_path, _input_file_lines_buff);
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			pf::time_interval::deinit(phaseMesh);
			pf::solver_init::deinit(phaseMesh);
			pf::mesh_structure::deinit(phaseMesh);
			pf::auto_input::deinit(phaseMesh);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void load_module(string input_file_path, int input_file_lines_buff = 1000, char split = ' ') {
			_input_file_path = input_file_path;
			_input_file_lines_buff = input_file_lines_buff;
			_split = split;
			// init input file reader
			InputFileReader::get_instance()->init(input_file_path, Solvers::get_instance()->writer, false, input_file_lines_buff, split);
			Solvers::get_instance()->create_a_new_module(init, exec_pre, exec_loop, deinit, write_scalar, write_vec3);
		};
	}
}