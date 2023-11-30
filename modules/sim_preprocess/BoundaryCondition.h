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
	namespace boundary_condition {
		static tensor3_double temperature_source_point;
		static double_box temperature_source_boundary;
		static int_box phi_source_boundary;
		static tensor2_double con_source_boundary;
		static tensor3_double phase_con_source_boundary;
		static tensor2_double grand_potential_source_boundary;

		static int Nx = 1;
		static int Ny = 1;
		static int Nz = 1;

		static double (*df_dx)(pf::PhaseNode&, int);

		namespace boundary_condition_funcs {
			// Allen Cahn
			static void boundary_condition_phi(pf::PhaseNode& node, pf::PhaseEntry& phase) {
				for (auto bc = phi_source_boundary.begin(); bc < phi_source_boundary.end(); bc++) {
					switch (bc->index)
					{
					case 0:  // down_x
						if (node._x == 0 && phase.index == bc->value)
							phase.phi = 1.0;
						break;
					case 1:  // up_x
						if (node._x == Nx - 1 && phase.index == bc->value)
							phase.phi = 1.0;
						break;
					case 2:  // down_y
						if (node._y == 0 && phase.index == bc->value)
							phase.phi = 1.0;
						break;
					case 3:  // up_y
						if (node._y == Ny - 1 && phase.index == bc->value)
							phase.phi = 1.0;
						break;
					case 4:  // down_z
						if (node._z == 0 && phase.index == bc->value)
							phase.phi = 1.0;
						break;
					case 5:  // up_z
						if (node._z == Nz - 1 && phase.index == bc->value)
							phase.phi = 1.0;
						break;
					default:
						break;
					}
				}
			}
			static void boundary_condition_total_x(pf::PhaseNode& node, int comp_index) {
				for (auto bc = con_source_boundary.begin(); bc < con_source_boundary.end(); bc++) {
					switch (bc->index)
					{
					case 0:  // down_x
						if (node._x == 0) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 1:  // up_x
						if (node._x == Nx - 1) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 2:  // down_y
						if (node._y == 0) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 3:  // up_y
						if (node._y == Ny - 1) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 4:  // down_z
						if (node._z == 0) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 5:  // up_z
						if (node._z == Nz - 1) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					default:
						break;
					}
				}
			}
			static void boundary_condition_grand_potential(pf::PhaseNode& node, int comp_index) {
				for (auto bc = grand_potential_source_boundary.begin(); bc < grand_potential_source_boundary.end(); bc++) {
					switch (bc->index)
					{
					case 0:  // down_x
						if (node._x == 0) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.potential.begin(); comp < node.potential.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 1:  // up_x
						if (node._x == Nx - 1) {
							for (auto bc_x = bc->begin(); bc_x < bc->end(); bc_x++)
								for (auto comp = node.potential.begin(); comp < node.potential.end(); comp++)
									if (bc_x->index == comp->index)
										comp->value = bc_x->val;
						}
						break;
					case 2:  // down_y
						if (node._y == 0) {
							for (auto bc_y = bc->begin(); bc_y < bc->end(); bc_y++)
								for (auto comp = node.potential.begin(); comp < node.potential.end(); comp++)
									if (bc_y->index == comp->index)
										comp->value = bc_y->val;
						}
						break;
					case 3:  // up_y
						if (node._y == Ny - 1) {
							for (auto bc_y = bc->begin(); bc_y < bc->end(); bc_y++)
								for (auto comp = node.potential.begin(); comp < node.potential.end(); comp++)
									if (bc_y->index == comp->index)
										comp->value = bc_y->val;
						}
						break;
					case 4:  // down_z
						if (node._z == 0) {
							for (auto bc_z = bc->begin(); bc_z < bc->end(); bc_z++)
								for (auto comp = node.potential.begin(); comp < node.potential.end(); comp++)
									if (bc_z->index == comp->index)
										comp->value = bc_z->val;
						}
						break;
					case 5:  // up_z
						if (node._z == Nz - 1) {
							for (auto bc_z = bc->begin(); bc_z < bc->end(); bc_z++)
								for (auto comp = node.potential.begin(); comp < node.potential.end(); comp++)
									if (bc_z->index == comp->index)
										comp->value = bc_z->val;
						}
						break;
					default:
						break;
					}
				}
			}
			static void boundary_condition_phase_x(pf::PhaseNode& node, pf::PhaseEntry& phase, int comp_index) {
				for (auto bc = phase_con_source_boundary.begin(); bc < phase_con_source_boundary.end(); bc++) {
					switch (bc->index)
					{
					case 0:  // down_x
						for (auto bc_phi = bc->begin(); bc_phi < bc->end(); bc_phi++)
							if (node._x == 0 && bc_phi->index == phase.index && phase.phi > Phi_Num_Cut_Off) {
								for (auto bc_x = bc_phi->begin(); bc_x < bc_phi->end(); bc_x++)
									for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
										if (bc_x->index == comp->index)
											comp->value = bc_x->val;
							}
						break;
					case 1:  // up_x
						for (auto bc_phi = bc->begin(); bc_phi < bc->end(); bc_phi++)
							if (node._x == Nx - 1 && bc_phi->index == phase.index && phase.phi > Phi_Num_Cut_Off) {
								for (auto bc_x = bc_phi->begin(); bc_x < bc_phi->end(); bc_x++)
									for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
										if (bc_x->index == comp->index)
											comp->value = bc_x->val;
							}
						break;
					case 2:  // down_y
						for (auto bc_phi = bc->begin(); bc_phi < bc->end(); bc_phi++)
							if (node._y == 0 && bc_phi->index == phase.index && phase.phi > Phi_Num_Cut_Off) {
								for (auto bc_x = bc_phi->begin(); bc_x < bc_phi->end(); bc_x++)
									for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
										if (bc_x->index == comp->index)
											comp->value = bc_x->val;
							}
						break;
					case 3:  // up_y
						for (auto bc_phi = bc->begin(); bc_phi < bc->end(); bc_phi++)
							if (node._y == Ny - 1 && bc_phi->index == phase.index && phase.phi > Phi_Num_Cut_Off) {
								for (auto bc_x = bc_phi->begin(); bc_x < bc_phi->end(); bc_x++)
									for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
										if (bc_x->index == comp->index)
											comp->value = bc_x->val;
							}
						break;
					case 4:  // down_z
						for (auto bc_phi = bc->begin(); bc_phi < bc->end(); bc_phi++)
							if (node._z == 0 && bc_phi->index == phase.index && phase.phi > Phi_Num_Cut_Off) {
								for (auto bc_x = bc_phi->begin(); bc_x < bc_phi->end(); bc_x++)
									for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
										if (bc_x->index == comp->index)
											comp->value = bc_x->val;
							}
						break;
					case 5:  // up_z
						for (auto bc_phi = bc->begin(); bc_phi < bc->end(); bc_phi++)
							if (node._z == Nz - 1 && bc_phi->index == phase.index && phase.phi > Phi_Num_Cut_Off) {
								for (auto bc_x = bc_phi->begin(); bc_x < bc_phi->end(); bc_x++)
									for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
										if (bc_x->index == comp->index)
											comp->value = bc_x->val;
							}
						break;
					default:
						break;
					}
				}
			}
			// Temperature
			static void BoundaryCondition_temp(pf::PhaseNode& node) {
				for (auto bc = temperature_source_boundary.begin(); bc < temperature_source_boundary.end(); bc++) {
					switch (bc->index)
					{
					case 0:  // down_x
						if (node._x == 0)
							node.temperature.T = bc->value;
						break;
					case 1:  // up_x
						if (node._x == Nx - 1)
							node.temperature.T = bc->value;
						break;
					case 2:  // down_y
						if (node._y == 0)
							node.temperature.T = bc->value;
						break;
					case 3:  // up_y
						if (node._y == Ny - 1)
							node.temperature.T = bc->value;
						break;
					case 4:  // down_z
						if (node._z == 0)
							node.temperature.T = bc->value;
						break;
					case 5:  // up_z
						if (node._z == Nz - 1)
							node.temperature.T = bc->value;
						break;
					default:
						break;
					}
				}
				for (auto x_pos = temperature_source_point.begin(); x_pos < temperature_source_point.end(); x_pos++)
					for (auto y_pos = x_pos->begin(); y_pos < x_pos->end(); y_pos++)
						for (auto z_pos = y_pos->begin(); z_pos < y_pos->end(); z_pos++)
							if (node._x == x_pos->index && node._y == y_pos->index && node._z == z_pos->index)
								node.temperature.T = z_pos->val;
			}
		}
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			Nx = phaseMesh.limit_x;
			Ny = phaseMesh.limit_y;
			Nz = phaseMesh.limit_z;
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard)
				Solvers::get_instance()->Phi_Solver_AC.Boundary_Condition = boundary_condition_funcs::boundary_condition_phi;
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard)
				Solvers::get_instance()->Phi_Solver_CH.Boundary_Condition_Phi = boundary_condition_funcs::boundary_condition_phi;

			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX)
				Solvers::get_instance()->C_Solver.Boundary_Condition_TotalX = boundary_condition_funcs::boundary_condition_total_x;
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP)
				Solvers::get_instance()->C_Solver.Boundary_Condition_TotalX = boundary_condition_funcs::boundary_condition_grand_potential;
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX)
				Solvers::get_instance()->C_Solver.Boundary_Condition_PhaseX = boundary_condition_funcs::boundary_condition_phase_x;

			if (Solvers::get_instance()->parameters.TempEType == TemperatureEquationType::TType_Standard)
				Solvers::get_instance()->T_Solver.BoundaryCondition = boundary_condition_funcs::BoundaryCondition_temp;

			if (Solvers::get_instance()->parameters.TempEType == TemperatureEquationType::TType_Standard) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .point = [(point_0_x, point_0_y, point_0_z, temperature), ... ] \n", InputFileReader::get_instance()->debug_file);
				string temp_point_key = "Solver.Mesh.BoundaryCondition.Temperature.Fix.point", temp_point_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(temp_point_key, temp_point_input, infile_debug)) {
					vector<InputValueType> temp_point_structure; temp_point_structure.push_back(InputValueType::IVType_INT); temp_point_structure.push_back(InputValueType::IVType_INT);
					temp_point_structure.push_back(InputValueType::IVType_INT); temp_point_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> temp_point_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(temp_point_structure, temp_point_key, temp_point_input, infile_debug);
					for (auto temp_point = temp_point_value.begin(); temp_point < temp_point_value.end(); temp_point++) {
						temperature_source_point.add_double((*temp_point)[0].int_value, (*temp_point)[1].int_value, (*temp_point)[2].int_value, (*temp_point)[3].double_value);
					}
				}
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .boundary = [(boundary, temperature), ... ] \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              boundary : 0 - x_down , 1 - x_up , 2 - y_down , 3 - y_up , 4 - z_down , 5 - z_up \n", InputFileReader::get_instance()->debug_file);
				string temp_boundary_key = "Solver.Mesh.BoundaryCondition.Temperature.Fix.boundary", temp_boundary_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(temp_boundary_key, temp_boundary_input, infile_debug)) {
					vector<InputValueType> temp_boundary_structure; temp_boundary_structure.push_back(InputValueType::IVType_INT); temp_boundary_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> temp_boundary_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(temp_boundary_structure, temp_boundary_key, temp_boundary_input, infile_debug);
					for (auto temp_boundary = temp_boundary_value.begin(); temp_boundary < temp_boundary_value.end(); temp_boundary++) {
						temperature_source_boundary.add_double((*temp_boundary)[0].int_value, (*temp_boundary)[1].double_value);
					}
				}
			}
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .boundary = [(boundary, phi_index), ... ] \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              boundary : 0 - x_down , 1 - x_up , 2 - y_down , 3 - y_up , 4 - z_down , 5 - z_up \n", InputFileReader::get_instance()->debug_file);
				string phi_boundary_key = "Solver.Mesh.BoundaryCondition.Phi.Fix.boundary", phi_boundary_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(phi_boundary_key, phi_boundary_input, infile_debug)) {
					vector<InputValueType> phi_boundary_structure; phi_boundary_structure.push_back(InputValueType::IVType_INT); phi_boundary_structure.push_back(InputValueType::IVType_INT);
					vector<vector<input_value>> phi_boundary_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(phi_boundary_structure, phi_boundary_key, phi_boundary_input, infile_debug);
					for (auto phi_boundary = phi_boundary_value.begin(); phi_boundary < phi_boundary_value.end(); phi_boundary++) {
						phi_source_boundary.add_int((*phi_boundary)[0].int_value, (*phi_boundary)[1].int_value);
					}
				}
			}
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .boundary = [(boundary, con_name, value), ... ] \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              boundary : 0 - x_down , 1 - x_up , 2 - y_down , 3 - y_up , 4 - z_down , 5 - z_up \n", InputFileReader::get_instance()->debug_file);
				string con_boundary_key = "Solver.Mesh.BoundaryCondition.Con.Fix.boundary", con_boundary_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(con_boundary_key, con_boundary_input, infile_debug)) {
					vector<InputValueType> con_boundary_structure; con_boundary_structure.push_back(InputValueType::IVType_INT); con_boundary_structure.push_back(InputValueType::IVType_STRING); con_boundary_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> con_boundary_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(con_boundary_structure, con_boundary_key, con_boundary_input, infile_debug);
					for (auto con_boundary = con_boundary_value.begin(); con_boundary < con_boundary_value.end(); con_boundary++) {
						con_source_boundary.add_double((*con_boundary)[0].int_value, Solvers::get_instance()->parameters.Components[(*con_boundary)[1].string_value].index, (*con_boundary)[2].double_value);
					}
				}
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				df_dx = Solvers::get_instance()->C_Solver.df_dx;
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .boundary = [(boundary, con_name, value), ... ] \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              boundary : 0 - x_down , 1 - x_up , 2 - y_down , 3 - y_up , 4 - z_down , 5 - z_up \n", InputFileReader::get_instance()->debug_file);
				string potential_boundary_key = "Solver.Mesh.BoundaryCondition.Potential.Fix.boundary", potential_boundary_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(potential_boundary_key, potential_boundary_input, infile_debug)) {
					vector<InputValueType> potential_boundary_structure; potential_boundary_structure.push_back(InputValueType::IVType_INT); 
					potential_boundary_structure.push_back(InputValueType::IVType_STRING); potential_boundary_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> potential_boundary_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(potential_boundary_structure, potential_boundary_key, potential_boundary_input, infile_debug);
					for (auto potential_boundary = potential_boundary_value.begin(); potential_boundary < potential_boundary_value.end(); potential_boundary++) {
						grand_potential_source_boundary.add_double((*potential_boundary)[0].int_value, Solvers::get_instance()->parameters.Components[(*potential_boundary)[1].string_value].index, (*potential_boundary)[2].double_value);
					}
				}
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .boundary = [(boundary, phi_index, con_name, value), ... ] \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              boundary : 0 - x_down , 1 - x_up , 2 - y_down , 3 - y_up , 4 - z_down , 5 - z_up \n", InputFileReader::get_instance()->debug_file);
				string phi_con_boundary_key = "Solver.Mesh.BoundaryCondition.Con.Fix.boundary", phi_con_boundary_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(phi_con_boundary_key, phi_con_boundary_input, infile_debug)) {
					vector<InputValueType> phi_con_boundary_structure; phi_con_boundary_structure.push_back(InputValueType::IVType_INT); phi_con_boundary_structure.push_back(InputValueType::IVType_INT);
					phi_con_boundary_structure.push_back(InputValueType::IVType_STRING); phi_con_boundary_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> phi_con_boundary_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(phi_con_boundary_structure, phi_con_boundary_key, phi_con_boundary_input, infile_debug);
					for (auto phi_con_boundary = phi_con_boundary_value.begin(); phi_con_boundary < phi_con_boundary_value.end(); phi_con_boundary++) {
						phase_con_source_boundary.add_double((*phi_con_boundary)[0].int_value, (*phi_con_boundary)[1].int_value, Solvers::get_instance()->parameters.Components[(*phi_con_boundary)[2].string_value].index, (*phi_con_boundary)[3].double_value);
					}
				}
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : BoundaryCondition for Phi C T !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP)
				df_dx = Solvers::get_instance()->C_Solver.df_dx;
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

		static void init_module() {
			Solvers::get_instance()->create_a_new_module(init, exec_pre, exec_loop, deinit, write_scalar, write_vec3);
		};
	}
}