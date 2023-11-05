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
	namespace interface_mobility {
		// unit = m/s
		static double const_Lij = 0.0;
		static double const_Mij = 0.0;
		static vector<tensor1_int> const_block_Lij;
		static tensor2_double matirx_Lij;
		static tensor2_double matirx_Mij;
		static PairValue block_Lij;
		static PairValue block_Mij;
		//-
		static tensor2_double get_matirx_Lij() {
			return matirx_Lij;
		}
		static PairValue get_block_Lij() {
			return block_Lij;
		}
		//-
		static double Lij_pair_wise_const_block(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			int check = 0;
			for (auto block = const_block_Lij.begin(); block < const_block_Lij.end(); block++) {
				check = 0;
				for (auto elem = block->begin(); elem < block->end(); elem++)
					if (elem->index == alpha.index || elem->index == beta.index)
						check++;
				if(check == 2)
					return const_Lij;
			}
			return 0.0;
		};
		static double Lij_pair_wise_const(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return const_Lij;
		};
		static double Lij_standard_const(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			if (alpha.index == beta.index)
				return const_Lij;
			else
				return 0.0;
		};
		static double Mij_const(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			if(alpha.index == beta.index)
				return const_Mij;
			else
				return 0.0;
		};
		static double Mij_matrix(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			// to be defined
			double _M = 0.0;
			for (auto Mi = matirx_Mij.begin(); Mi < matirx_Mij.end(); Mi++)
				for (auto Mj = Mi->begin(); Mj < Mi->end(); Mj++)
					if ((Mi->index == alpha.index && Mj->index == beta.index) || (Mi->index == beta.index && Mj->index == alpha.index))
						_M = Mj->val;
			return _M;
		};
		static double Lij_matrix(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			// to be defined
			double _L = 0.0;
			for (auto Li = matirx_Lij.begin(); Li < matirx_Lij.end(); Li++)
				for (auto Lj = Li->begin(); Lj < Li->end(); Lj++)
					if ((Li->index == alpha.index && Lj->index == beta.index) || (Li->index == beta.index && Lj->index == alpha.index))
						_L = Lj->val;
			return _L;
		};
		static double Mij_block(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			// to be defined
			double _M = 0.0;
			for (auto block = block_Mij.begin(); block < block_Mij.end(); block++)
				if ((block->pairIndex_1 <= alpha.index && alpha.index <= block->pairIndex_2) && (block->pairIndex_1 <= beta.index && beta.index <= block->pairIndex_2))
					_M = block->value;
			return _M;
		};
		static double Lij_block(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			// to be defined
			double _L = 0.0;
			for (auto block = block_Lij.begin(); block < block_Lij.end(); block++)
				if ((block->pairIndex_1 <= alpha.index && alpha.index <= block->pairIndex_2) && (block->pairIndex_1 <= beta.index && beta.index <= block->pairIndex_2))
					_L = block->value;
			return _L;
		};

		static void load_mobility(PhiEquationType _phi_type, bool infile_debug) {
			if (_phi_type == PhiEquationType::PEType_AC_Pairwise) {
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.Lij.const  = Lij_value \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                      .matrix = [(phi_i, phi_j, Lij_value), ... ] \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                      .block = [(phi_begin, phi_end, Lij_value), ... ] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key = "ModelsManager.Phi.Lij.matrix", matrix_input = "[()]", block_key = "ModelsManager.Phi.Lij.block", block_input = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Lij.const", const_Lij, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_pair_wise_const;
					if (infile_debug) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.Lij.Const.block = [(phi_index1, phi_index2, ... ), ... ] \n", InputFileReader::get_instance()->debug_file);
					}
					string const_block_key = "ModelsManager.Phi.Lij.Const.block", const_block_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(const_block_key, const_block_input, infile_debug)) {
						Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_pair_wise_const_block;
						vector<vector<input_value>> const_block_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, infile_debug);
						for (int index = 0; index < const_block_value.size(); index++) {
							tensor1_int const_block;
							for (int phi_index = 0; phi_index < const_block_value[index].size(); phi_index++)
								const_block.add_int(const_block_value[index][phi_index].int_value, 0);
							const_block_Lij.push_back(const_block);
						}
					}
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_input, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key, matrix_input, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_Lij.add_double(matrix_value[index][0].int_value, matrix_value[index][1].int_value, matrix_value[index][2].double_value);
				}
				else if (InputFileReader::get_instance()->read_string_value(block_key, block_input, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_block;
					vector<InputValueType> block_structure; block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> block_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(block_structure, block_key, block_input, infile_debug);
					for (int index = 0; index < block_value.size(); index++)
						block_Lij.add(block_value[index][0].int_value, block_value[index][1].int_value, block_value[index][2].double_value);
				}
			}
			else if (_phi_type == PhiEquationType::PEType_AC_Standard) {
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.Lij.const  = Lii_value \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                      .matrix = [(phi_i, phi_j, Lij_value), ... ] \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                      .block = [(phi_begin, phi_end, Lij_value), ... ] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key = "ModelsManager.Phi.Lij.matrix", matrix_input = "[()]", block_key = "ModelsManager.Phi.Lij.block", block_input = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Lij.const", const_Lij, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_standard_const;
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_input, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key, matrix_input, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_Lij.add_double(matrix_value[index][0].int_value, matrix_value[index][1].int_value, matrix_value[index][2].double_value);
				}
				else if (InputFileReader::get_instance()->read_string_value(block_key, block_input, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_block;
					vector<InputValueType> block_structure; block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> block_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(block_structure, block_key, block_input, infile_debug);
					for (int index = 0; index < block_value.size(); index++)
						block_Lij.add(block_value[index][0].int_value, block_value[index][1].int_value, block_value[index][2].double_value);
				}
			}
			else if (_phi_type == PhiEquationType::PEType_CH_Standard) {
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.Mij.const  = Mij_value \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                      .matrix = [(phi_i, phi_j, Mij_value), ... ] \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                      .block = [(phi_begin, phi_end, Mij_value), ... ] \n", InputFileReader::get_instance()->debug_file);
				}
				string matrix_key = "ModelsManager.Phi.Mij.matrix", matrix_input = "[()]", block_key = "ModelsManager.Phi.Mij.block", block_input = "[()]";
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Mij.const", const_Mij, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_CH.M_ab = Mij_const;
				}
				else if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_input, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_CH.M_ab = Mij_matrix;
					vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key, matrix_input, infile_debug);
					for (int index = 0; index < matrix_value.size(); index++)
						matirx_Mij.add_double(matrix_value[index][0].int_value, matrix_value[index][1].int_value, matrix_value[index][2].double_value);
				}
				else if (InputFileReader::get_instance()->read_string_value(block_key, block_input, infile_debug)) {
					Solvers::get_instance()->Phi_Solver_CH.M_ab = Mij_block;
					vector<InputValueType> block_structure; block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> block_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(block_structure, matrix_key, matrix_input, infile_debug);
					for (int index = 0; index < block_value.size(); index++)
						block_Mij.add(block_value[index][0].int_value, block_value[index][1].int_value, block_value[index][2].double_value);
				}
			}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			block_Lij.set(pf::PairValueProperty::pf_STORAGE);
			block_Mij.set(pf::PairValueProperty::pf_STORAGE);
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_pair_wise_const;
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
				Solvers::get_instance()->Phi_Solver_AC.Lij = Lij_standard_const;
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				Solvers::get_instance()->Phi_Solver_CH.M_ab = Mij_const;
			}
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			load_mobility(Solvers::get_instance()->parameters.PhiEType, infile_debug);
			
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Mobility !\n", LOG_FILE_NAME);
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