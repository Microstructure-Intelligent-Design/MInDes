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
#include "Mechanics/ElasticSolver.h"
#include "Mechanics/PlasticSolver.h"

namespace pf {
	namespace mechanical_field {

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			elastic_solver::init(phaseMesh);
			plastic_solver::init(phaseMesh);
			crack_propagation::init(phaseMesh);
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : MechanicalField !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			stiffness_eigenstrain::exec_pre(phaseMesh);
			if (crack_propagation::crack_propagation_model() != crack_propagation::CPM_None) {
				bool is_iterate = true;
				int iterate_times = 0;
				elastic_solver::exec_pre(phaseMesh);
				crack_propagation::crack_pre(phaseMesh);
				Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Do crack propagation:\n", LOG_FILE_NAME);
				while (is_iterate && iterate_times < crack_propagation::crack_solver_max_iterate_steps())
				{
					iterate_times++;
					stringstream output;
					output << elastic_solver::exec_loop(phaseMesh);
					double variation = crack_propagation::crack_loop(phaseMesh);
					output << "> iterate step: " << iterate_times << ", crack variation = " << variation << endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
					if (variation < crack_propagation::crack_solver_epsilon())
						is_iterate = false;
				}
			}
			else {
				elastic_solver::exec_pre(phaseMesh);
			}
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string mech_report = "";
			stiffness_eigenstrain::exec_pre(phaseMesh);
			if (pf::crack_propagation::crack_propagation_model() != pf::crack_propagation::CPM_None) {
				bool is_iterate = true;
				int iterate_times = 0;
				double max_variation = 0.0;
				while (is_iterate && iterate_times < pf::crack_propagation::crack_solver_max_iterate_steps()) {
					iterate_times++;
					mech_report = pf::elastic_solver::exec_loop(phaseMesh);
					double variation = pf::crack_propagation::crack_loop(phaseMesh);
					if (max_variation < variation)
						max_variation = variation;
					if (variation < pf::crack_propagation::crack_solver_epsilon())
						is_iterate = false;
				}
				mech_report += "> crack solver iterate " + to_string(iterate_times) + " times, max crack variation = " + to_string(max_variation) + "\n";
			}
			else {
				mech_report = pf::elastic_solver::exec_loop(phaseMesh);
			}
			return mech_report;
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			elastic_solver::write_vec3(fout, phaseMesh);
			plastic_solver::write_vec3(fout, phaseMesh);
			crack_propagation::write_vec3(fout, phaseMesh);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			elastic_solver::write_scalar(fout, phaseMesh);
			plastic_solver::write_scalar(fout, phaseMesh);
			crack_propagation::write_scalar(fout, phaseMesh);
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			elastic_solver::deinit(phaseMesh);
			plastic_solver::deinit(phaseMesh);
			crack_propagation::deinit(phaseMesh);
		}
	}
}