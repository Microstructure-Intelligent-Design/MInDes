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
#include "Reaction/ElectrodeReaction.h"

namespace pf {
	namespace interface_reaction {
		// Phi
		
		// main function
		static vector<double(*)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&)> interface_reaction_ab_list;
		static double interface_reaction_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			double buff = 0.0;
			for (auto func = interface_reaction_ab_list.begin(); func < interface_reaction_ab_list.end(); func++)
				buff += (*func)(node, alpha, beta);
			return buff;
		}

		// Concentration
		static double_box k_const;
		static tensor3_double k_matrix;
		static double K_AB_i_const(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, int con_i) {
			for(auto k = k_const.begin(); k < k_const.end(); k++)
				if(k->index == con_i)
					return k->value;
			return 0.0;
		}
		static double K_AB_i_matrix(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, int con_i) {
			for (auto A_index = k_matrix.begin(); A_index < k_matrix.end(); A_index++)
				if (A_index->index == alpha.index)
					for (auto B_index = A_index->begin(); B_index < A_index->end(); B_index++)
						if (B_index->index == beta.index)
							for (auto con = B_index->begin(); con < B_index->end(); con++)
								if (con->index == con_i)
									return con->val;
			return 0.0;
		}
		static double (*K_AB_i)(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, int con_i);
		static void interface_reaction_AB_dissipation(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double abs_dPaPb) {
			// AbsGradPhi_Phi * ( standard_reaction_flux_on_interface )
			for (auto comp1 = alpha.x.begin(); comp1 < alpha.x.end(); comp1++)
				for (auto comp2 = beta.x.begin(); comp2 < beta.x.end(); comp2++)
					if (comp1->index == comp2->index) {
						double reaction_rate = K_AB_i(alpha, beta, comp1->index);
						double pre = reaction_rate * abs_dPaPb * (beta.potential[comp2->index].value - alpha.potential[comp1->index].value);
						if (pre > 0 && comp2->value > SYS_EPSILON && comp1->value < (1.0 - SYS_EPSILON)) {
							comp1->ChemicalReactionFlux += pre;
						}
						else if (pre < 0 && comp1->value > SYS_EPSILON && comp2->value < (1.0 - SYS_EPSILON)) {
							comp1->ChemicalReactionFlux += pre;
						}
					}
		}
		static vector<void(*)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, double)> interface_reaction_AB_list;

		// main function
		static void interface_reaction_AB(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double abs_dPaPb) {
			for (auto func = interface_reaction_AB_list.begin(); func < interface_reaction_AB_list.end(); func++)
				(*func)(node, alpha, beta, abs_dPaPb);
		}

		
		static vector<double(*)(pf::PhaseNode&, int)> interface_reaction_i_list;

		// main function
		static double interface_reaction_i(pf::PhaseNode& node, int con_i) {
			double buff = 0.0;
			for (auto func = interface_reaction_i_list.begin(); func < interface_reaction_i_list.end(); func++)
				buff += (*func)(node, con_i);
			return buff;
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Con.DissipationRate.const  = [( con_name , k_i ) ... ] \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                  .matrix = [(phi_name_a, phi_name_b, con_name, k_AB_i), ...] \n", InputFileReader::get_instance()->debug_file);
				string matrix_key = "ModelsManager.Con.DissipationRate.matrix", matrix_input = "[()]",
				       const_key = "ModelsManager.Con.DissipationRate.const", const_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(const_key, const_input, infile_debug)) {
					K_AB_i = K_AB_i_const;
					interface_reaction_AB_list.push_back(interface_reaction_AB_dissipation);
					vector<InputValueType> const_structure; const_structure.push_back(InputValueType::IVType_STRING); const_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> const_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(const_structure, const_key, const_input, infile_debug);
					for (int index = 0; index < const_value.size(); index++)
						k_const.add_double(Solvers::get_instance()->parameters.Components[const_value[index][0].string_value].index, const_value[index][1].double_value);
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_input, infile_debug)) {
					K_AB_i = K_AB_i_matrix;
					interface_reaction_AB_list.push_back(interface_reaction_AB_dissipation);
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); 
					matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key, matrix_input, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						k_matrix.add_double(matrix_value[index][0].int_value, matrix_value[index][1].int_value, Solvers::get_instance()->parameters.Components[matrix_value[index][2].string_value].index, matrix_value[index][3].double_value);
				}
				int electrolyte_phase_index = electrode_reaction::electrolyte_phase_index;
				if (InputFileReader::get_instance()->read_int_value("ModelsManager.Con.ElectrodeReaction.electrolyte_index", electrolyte_phase_index, infile_debug)) {
					electrode_reaction::init(phaseMesh);
					interface_reaction_AB_list.push_back(electrode_reaction::calc_battery_charge_smoothed_boundary);
				}
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP ||
				Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				int electrolyte_phase_index = electrode_reaction::electrolyte_phase_index;
				if (InputFileReader::get_instance()->read_int_value("ModelsManager.Con.ElectrodeReaction.electrolyte_index", electrolyte_phase_index, infile_debug)) {
					electrode_reaction::init(phaseMesh);
					interface_reaction_i_list.push_back(electrode_reaction::int_flux);
				}
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Interface Reaction !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			electrode_reaction::exec_pre(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			report += electrode_reaction::exec_loop(phaseMesh);
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