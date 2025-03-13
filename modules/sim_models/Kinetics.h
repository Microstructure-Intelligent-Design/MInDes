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
	namespace kinetics_funcs {
		static double D_default(pf::PhaseNode& node) {
			return 0.0;
		}
		static void M_phase_default(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}

		static double M_bulk_default(pf::PhaseNode& node, int con_i, int con_j) {
			return 0.0;
		}

	}
	namespace kinetics {
		// phase x
		static tensor2_double phase_Mii;
		static tensor3_double phase_Mij;
		// for electrodeposition
		static int active_comp_index = 0;
		static double Mii_elec = 0.0;
		static double Dconst{};
		// temp
		static tensor1_double phase_Dtemp;

		static double D_temp(pf::PhaseNode& node) {
			double D_temp = 0.0;
			for (auto phi = node.begin(); phi < node.end(); phi++)
				for (auto D = phase_Dtemp.begin(); D < phase_Dtemp.end(); D++)
					if (phi->property == D->index)
						D_temp += phi->phi * D->val;
			return D_temp;
		}

		static double D_const(pf::PhaseNode& node) {
			return Dconst;
		}

		static void Mphase_ii(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			tensor1_double M_phase = phase_Mii(phase.property);
			for (auto xi = phase.x.begin(); xi < phase.x.end(); xi++)
				for (auto xj = phase.x.begin(); xj < phase.x.end(); xj++)
					if (xi->index == xj->index)
						phase.kinetics_coeff.set(xi->index, xj->index, M_phase(xi->index));
		}
		static void Mphase_ij(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			tensor2_double M_phase = phase_Mij(phase.property);
			for (auto xi = phase.x.begin(); xi < phase.x.end(); xi++)
				for (auto xj = phase.x.begin(); xj < phase.x.end(); xj++)
					phase.kinetics_coeff.set(xi->index, xj->index, M_phase(xi->index, xj->index));
		}
		static double Mtotal_ii(pf::PhaseNode& node, int con_i, int con_j) {
			double Mii = 0.0;
			if (con_i == con_j) {
				for (auto phase = node.begin(); phase < node.end(); phase++)
					if (phase->phi > SYS_EPSILON)
						for (auto phi_M = phase_Mii.begin(); phi_M < phase_Mii.end(); phi_M++)
							if (phase->property == phi_M->index) {
								for (auto M = phi_M->begin(); M < phi_M->end(); M++)
									if (con_i == M->index)
										Mii += phase->phi * M->val;
							}
			}
			return Mii;
		}
		static double Mtotal_electrodeposition_ii(pf::PhaseNode& node, int con_i, int con_j) {
			double Mii = 0.0;
			if (con_i == con_j && con_i == active_comp_index) {
				double solid_phi = 0.0;
				for (auto phase = node.begin(); phase < node.end(); phase++)
					if (phase->phi > SYS_EPSILON) {
						solid_phi += phase->phi;
						for (auto phi_M = phase_Mii.begin(); phi_M < phase_Mii.end(); phi_M++)
							if (phase->property == phi_M->index) {
								for (auto M = phi_M->begin(); M < phi_M->end(); M++)
									if (con_i == M->index)
										Mii += phase->phi * M->val;
							}
					}
				if (solid_phi > SYS_EPSILON)
					Mii = Mii / solid_phi;
				double h_func = interpolation_func(solid_phi);
				Mii = Mii * h_func + Mii_elec * (1.0 - h_func);
			}
			return Mii;
		}
		static double Mtotal_ij(pf::PhaseNode& node, int con_i, int con_j) {
			double Mij = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++)
				if (phase->phi > SYS_EPSILON)
					for (auto phi_M = phase_Mij.begin(); phi_M < phase_Mij.end(); phi_M++)
						if (phase->property == phi_M->index) {
							for (auto Mi = phi_M->begin(); Mi < phi_M->end(); Mi++)
								if (con_i == Mi->index)
									for (auto Mj = Mi->begin(); Mj < Mi->end(); Mj++)
										Mij += phase->phi * Mj->val;
						}
			return 0.0;
		}
		static double Mgrand_ii(pf::PhaseNode& node, int con_i, int con_j) {
			double Mii = 0.0;
			if (con_i == con_j) {
				for (auto phi = node.begin(); phi < node.end(); phi++)
					for (auto phi_M = phase_Mii.begin(); phi_M < phase_Mii.end(); phi_M++)
						if (phi->property == phi_M->index) {
							for (auto M = phi_M->begin(); M < phi_M->end(); M++)
								if (con_i == M->index)
									Mii += phi->phi * M->val;
						}
			}
			return Mii;
		}

		static void (*M_phase_ij)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*M_bulk_ij)(pf::PhaseNode& node, int con_i, int con_j);

		static double (*D_phase_temp)(pf::PhaseNode& node);

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			M_phase_ij = kinetics_funcs::M_phase_default;
			M_bulk_ij = kinetics_funcs::M_bulk_default;
			D_phase_temp = kinetics_funcs::D_default;
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Con.Mii = [(phase_0_M_00 , phase_0_M_11, ...) , ... ]\n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                  .Mij = {[(phase_0_M_00, phase_0_M_01, ...), (phase_0_M_10, phase_0_M_11, ...), ... ], ... }\n", InputFileReader::get_instance()->debug_file);
				string Mii_key = "ModelsManager.Con.Mii", Mii_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(Mii_key, Mii_input, infile_debug)) {
					M_phase_ij = Mphase_ii;
					vector<vector<input_value>> Mii_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, Mii_key, Mii_input, infile_debug);
					if (Mii_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mii[] isn't equal to Phases \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					int index_1 = 0;
					for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
						phase_Mii.add_tensor(phi->phi_property);
						if (Mii_value[index_1].size() != phi->x.size()) {
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mii[()] isn't equal to Phase Cons \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
						int index_2 = 0;
						for (auto phi_con = phi->x.begin(); phi_con < phi->x.end(); phi_con++) {
							phase_Mii(phi->phi_property).add_double(phi_con->index, Mii_value[index_1][index_2].double_value);
							index_2++;
						}
						index_1++;
					}
				}
				string Mij_key = "ModelsManager.Con.Mij", Mij_input = "{[()]}";
				if (InputFileReader::get_instance()->read_string_value(Mij_key, Mij_input, infile_debug)) {
					M_phase_ij = Mphase_ij;
					vector<vector<vector<input_value>>> Mij_value = InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_DOUBLE, Mij_key, Mij_input, infile_debug);
					if (Mij_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mij{} isn't equal to Phases \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					int a_index = 0, i_index = 0, j_index = 0;
					for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
						i_index = 0;
						phase_Mij.add_tensor(phi->phi_property);
						if (Mij_value[a_index].size() != phi->x.size()) {
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mij{[]} isn't equal to the number of components in aim phase \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
						for (auto xi = phi->x.begin(); xi < phi->x.end(); xi++) {
							j_index = 0;
							phase_Mij(phi->phi_property).add_tensor(xi->index);
							if (Mij_value[a_index][i_index].size() != phi->x.size()) {
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mij{[()]} isn't equal to the number of components in aim phase \n", InputFileReader::get_instance()->debug_file);
								exit(0);
							}
							for (auto xj = phi->x.begin(); xj < phi->x.end(); xj++) {
								phase_Mij(phi->phi_property)(xi->index).add_double(xj->index, Mij_value[a_index][i_index][j_index].double_value);
								j_index++;
							}
							i_index++;
						}
						a_index++;
					}
				}
				Solvers::get_instance()->C_Solver.Mbulk_ij = M_phase_ij;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Con.Mii = [(phase_0_M_00 , phase_0_M_11, ...) , ... ]\n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                  .Mij = {[(phase_0_M_00, phase_0_M_01, ...), (phase_0_M_10, phase_0_M_11, ...), ... ], ... }\n", InputFileReader::get_instance()->debug_file);
				string Mii_key = "ModelsManager.Con.Mii", Mii_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(Mii_key, Mii_input, infile_debug)) {
					M_bulk_ij = Mtotal_ii;
					vector<vector<input_value>> Mii_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, Mii_key, Mii_input, infile_debug);
					if (Mii_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mii[] isn't equal to Phases \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					int index_1 = 0;
					for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
						phase_Mii.add_tensor(phi->phi_property);
						if (Mii_value[index_1].size() != phi->x.size()) {
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mii[()] isn't equal to Phase Cons \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
						int index_2 = 0;
						for (auto phi_con = phi->x.begin(); phi_con < phi->x.end(); phi_con++) {
							phase_Mii(phi->phi_property).add_double(phi_con->index, Mii_value[index_1][index_2].double_value);
							index_2++;
						}
						index_1++;
					}
				}
				string Mij_key = "ModelsManager.Con.Mij", Mij_input = "{[()]}";
				if (InputFileReader::get_instance()->read_string_value(Mij_key, Mij_input, infile_debug)) {
					M_bulk_ij = Mtotal_ij;
					vector<vector<vector<input_value>>> Mij_value = InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_DOUBLE, Mij_key, Mij_input, infile_debug);
					if (Mij_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mij{} isn't equal to Phases \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					int a_index = 0, i_index = 0, j_index = 0;
					for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
						i_index = 0;
						phase_Mij.add_tensor(phi->phi_property);
						if (Mij_value[a_index].size() != phi->x.size()) {
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mij{[]} isn't equal to the number of components in aim phase \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
						for (auto xi = phi->x.begin(); xi < phi->x.end(); xi++) {
							j_index = 0;
							phase_Mij(phi->phi_property).add_tensor(xi->index);
							if (Mij_value[a_index][i_index].size() != phi->x.size()) {
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mij{[()]} isn't equal to the number of components in aim phase \n", InputFileReader::get_instance()->debug_file);
								exit(0);
							}
							for (auto xj = phi->x.begin(); xj < phi->x.end(); xj++) {
								phase_Mij(phi->phi_property)(xi->index).add_double(xj->index, Mij_value[a_index][i_index][j_index].double_value);
								j_index++;
							}
							i_index++;
						}
						a_index++;
					}
				}
				// - for electrodeposition
				bool is_electric_field_on = false;
				InputFileReader::get_instance()->read_bool_value("Postprocess.PhysicalFields.electric", is_electric_field_on, false);
				string active_comp_name = "";
				if (InputFileReader::get_instance()->read_string_value("ModelsManager.PhiCon.ElectroDeposition.active_component", active_comp_name, false) && is_electric_field_on) {
					active_comp_index = Solvers::get_instance()->parameters.Components[active_comp_name].index;
					InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectroDeposition.Electrolyte_Mii", Mii_elec, infile_debug);
					M_bulk_ij = Mtotal_electrodeposition_ii;
				}
				Solvers::get_instance()->C_Solver.M_ij = M_bulk_ij;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Con.Mii = [(phase_0_M_00 , phase_0_M_11, ...) , ... ]\n", InputFileReader::get_instance()->debug_file);
				string Mii_key = "ModelsManager.Con.Mii", Mii_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(Mii_key, Mii_input, infile_debug)) {
					M_bulk_ij = Mgrand_ii;
					vector<vector<input_value>> Mii_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, Mii_key, Mii_input, infile_debug);
					if (Mii_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mii[] isn't equal to Phases \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					int index_1 = 0;
					for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
						phase_Mii.add_tensor(phi->phi_property);
						if (Mii_value[index_1].size() != phi->x.size()) {
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Con.Mii[()] isn't equal to Phase Cons \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
						int index_2 = 0;
						for (auto phi_con = phi->x.begin(); phi_con < phi->x.end(); phi_con++) {
							phase_Mii(phi->phi_property).add_double(phi_con->index, Mii_value[index_1][index_2].double_value);
							index_2++;
						}
						index_1++;
					}
				}
				Solvers::get_instance()->C_Solver.M_ij = M_bulk_ij;
			}
			if (Solvers::get_instance()->parameters.TempEType == TemperatureEquationType::TType_Standard) {
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Temp.Dtemp = ( phase_0_D , phase_1_D , ... ) \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                   .Dconst= 0.0 \n", InputFileReader::get_instance()->debug_file);
				}
				string D_a_key = "ModelsManager.Temp.Dtemp", D_a_input = "()";
				if (InputFileReader::get_instance()->read_string_value(D_a_key, D_a_input, infile_debug)) {
					D_phase_temp = D_temp;
					vector<input_value> D_a_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, D_a_key, D_a_input, infile_debug);
					if (D_a_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of ModelsManager.Temp.Dtemp isn't equal to Phases \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					int a_index = 0;
					for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
						phase_Dtemp.add_double(phi->phi_property, D_a_value[a_index].double_value);
						a_index++;
					}
				}
				if (InputFileReader::get_instance()->read_double_value("ModelsManager.Temp.Dconst", Dconst, infile_debug)) {
					D_phase_temp = D_const;
				}
				Solvers::get_instance()->T_Solver.D_temp = D_phase_temp;
			}
			Solvers::get_instance()->Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Kinetics !\n", LOG_FILE_NAME);
		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			M_phase_ij = nullptr;
			M_bulk_ij = nullptr;
			D_phase_temp = nullptr;
		}

	}
}