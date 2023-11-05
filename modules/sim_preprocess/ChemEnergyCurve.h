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
#include "../sim_models/BulkEnergy/ChemicalEnergy.h"

namespace pf {
	namespace chemical_energy_curve {
		enum OutputFileFormat { OFF_Matlab };
		//-
		static PhaseNode sample_node;
		//-
		static tensor2_double phases_fix_comp;
		static double con_epsilon = 0.1;
		//-
		static OutputFileFormat outFile_format = OFF_Matlab;
		static string file_name = "chemical_energy_curve";
		static bool is_plot_energy_curve = false;
		static void write_chemical_energy_hlfuncs_to_m_one_comp_change(PhaseEntry& phase) {
			PhaseNode node;
			tensor1_double& phase_fix_comp = phases_fix_comp(phase.property);
			Info_Phases& Phases = Solvers::get_instance()->parameters.Phases;
			int con_index = 0;
			for (auto x = phase.x.begin(); x < phase.x.end(); x++) {
				bool is_fixed = false;
				for (auto fix_con = phase_fix_comp.begin(); fix_con < phase_fix_comp.end(); fix_con++)
					if (x->index == fix_con->index)
						is_fixed = true;
				if (!is_fixed) {
					con_index = x->index;
					break;
				}
			}
			stringstream comp_i, fchem, dfchem_dcomp_i;
			comp_i << "con_" << Phases[phase.property].x[con_index].name << " = [ ";
			fchem << "fchem_" << Phases[phase.property].phi_name << " = [ ";
			dfchem_dcomp_i << "df_" << Phases[phase.property].phi_name << "_dcon_" << Phases[phase.property].x[con_index].name << " = [ ";
			double other_con = 0.0;
			for (auto x = phase.x.begin(); x < phase.x.end(); x++)
				for (auto fix_con = phase_fix_comp.begin(); fix_con < phase_fix_comp.end(); fix_con++)
					if (x->index == fix_con->index) {
						x->value = fix_con->val;
						other_con += x->value;
					}
			for (double con = 0.0; con < 1.0; con += con_epsilon) {
				phase.x[con_index].value = con;
				comp_i << con << ", ";
				if (con + other_con <= 1.0) {
					double energy_density = chemical_energy::fchem_hlfuncs(node, phase);
					fchem << energy_density << ", ";
					double diffusion_potential = chemical_energy::dfchem_dcon_hlfuncs(phase, preset_function::comp2hlf(phase.x), con_index);
					dfchem_dcomp_i << diffusion_potential << ", ";
				}
				else {
					fchem << "NaN, ";
					dfchem_dcomp_i<< "NaN, ";
				}
			}
			comp_i << "];" << endl;
			fchem << "];" << endl;
			dfchem_dcomp_i << "];" << endl;
			Solvers::get_instance()->writer.add_string_to_m(comp_i.str(), file_name);
			Solvers::get_instance()->writer.add_string_to_m(fchem.str(), file_name);
			Solvers::get_instance()->writer.add_string_to_m(dfchem_dcomp_i.str(), file_name);
		}
		static void write_chemical_energy_hlfuncs_to_m_two_comp_change(PhaseEntry& phase) {
			PhaseNode node;
			tensor1_double& phase_fix_comp = phases_fix_comp(phase.property);
			Info_Phases& Phases = Solvers::get_instance()->parameters.Phases;
			int con_index_i = 0, con_index_j = 0;
			bool is_i_set = false;
			for (auto x = phase.x.begin(); x < phase.x.end(); x++) {
				bool is_fixed = false;
				for (auto fix_con = phase_fix_comp.begin(); fix_con < phase_fix_comp.end(); fix_con++)
					if (x->index == fix_con->index)
						is_fixed = true;
				if (!is_fixed && !is_i_set) {
					con_index_i = x->index;
				}
				else if (!is_fixed && is_i_set) {
					con_index_j = x->index;
					break;
				}
			}
			stringstream comp_i, comp_j, fchem, dfchem_dcomp_i, dfchem_dcomp_j;
			comp_i << "con_" << Phases[phase.property].x[con_index_i].name << " = [ ";
			comp_j << "con_" << Phases[phase.property].x[con_index_j].name << " = [ ";
			fchem << "fchem_" << Phases[phase.property].phi_name << " = [ ";
			dfchem_dcomp_i << "df_" << Phases[phase.property].phi_name << "_dcon_" << Phases[phase.property].x[con_index_i].name << " = [ ";
			dfchem_dcomp_j << "df_" << Phases[phase.property].phi_name << "_dcon_" << Phases[phase.property].x[con_index_j].name << " = [ ";
			double other_con = 0.0;
			for (auto x = phase.x.begin(); x < phase.x.end(); x++)
				for (auto fix_con = phase_fix_comp.begin(); fix_con < phase_fix_comp.end(); fix_con++)
					if (x->index == fix_con->index) {
						x->value = fix_con->val;
						other_con += x->value;
					}
			for (double con_i = 0.0; con_i < 1.0; con_i += con_epsilon) {
				for (double con_j = 0.0; con_j < 1.0; con_j += con_epsilon) {
					phase.x[con_index_i].value = con_i;
					phase.x[con_index_j].value = con_j;
					comp_i << con_i << ", ";
					comp_j << con_j << ", ";
					if (con_i + con_j + other_con <= 1.0) {
						double energy_density = chemical_energy::fchem_hlfuncs(node, phase);
						fchem << energy_density << ", ";
						vector<double> concentration = preset_function::comp2hlf(phase.x);
						double diffusion_potential = chemical_energy::dfchem_dcon_hlfuncs(phase, concentration, con_index_i);
						dfchem_dcomp_i << diffusion_potential << ", ";
						diffusion_potential = chemical_energy::dfchem_dcon_hlfuncs(phase, concentration, con_index_j);
						dfchem_dcomp_j << diffusion_potential << ", ";
					}
					else {
						fchem << "NaN, ";
						dfchem_dcomp_i << "NaN, ";
						dfchem_dcomp_j << "NaN, ";
					}
				}
				comp_i << "; ";
				comp_j << "; ";
				fchem << "; ";
				dfchem_dcomp_i << "; ";
				dfchem_dcomp_j << "; ";
			}
			comp_i << "];" << endl;
			comp_j << "];" << endl;
			fchem << "];" << endl;
			dfchem_dcomp_i << "];" << endl;
			dfchem_dcomp_j << "];" << endl;
			Solvers::get_instance()->writer.add_string_to_m(comp_i.str(), file_name);
			Solvers::get_instance()->writer.add_string_to_m(comp_j.str(), file_name);
			Solvers::get_instance()->writer.add_string_to_m(fchem.str(), file_name);
			Solvers::get_instance()->writer.add_string_to_m(dfchem_dcomp_i.str(), file_name);
			Solvers::get_instance()->writer.add_string_to_m(dfchem_dcomp_j.str(), file_name);
		}
		static void write_chemical_energy_hlfuncs_to_m() {
			Solvers::get_instance()->writer.init_m_file(file_name);
			stringstream out;
			Info_Phases& Phases = Solvers::get_instance()->parameters.Phases;
			for (auto phase = sample_node.begin(); phase < sample_node.end(); phase++) {
				out.str("");
				out << endl;
				out << "% following are energy curve for phase: " << Phases[phase->property].phi_name << endl;
				out << "% there are " << Phases[phase->property].x.size() << " components in this phase: ";
				for (auto comp = Phases[phase->property].x.begin(); comp < Phases[phase->property].x.end(); comp++)
					out << comp->name << ", ";
				out << endl;
				out << "% where " << phases_fix_comp(phase->property).size() << " components are fixed: ";
				for (auto comp = phases_fix_comp(phase->property).begin(); comp < phases_fix_comp(phase->property).end(); comp++)
					out << Phases[phase->property].x[comp->index].name << " - " << comp->val << ", ";
				out << endl;
				Solvers::get_instance()->writer.add_string_to_m(out.str(), file_name);
				if (phase->x.size() - phases_fix_comp(phase->property).size() == 1) {
					write_chemical_energy_hlfuncs_to_m_one_comp_change(*phase);
				}
				else if (phase->x.size() - phases_fix_comp(phase->property).size() == 2) {
					write_chemical_energy_hlfuncs_to_m_two_comp_change(*phase);
				}
				else {
					out.str("");
					out << "% the energy curve plotting for this phase fails and you should set one or two independent components!" << endl;
					Solvers::get_instance()->writer.add_string_to_m(out.str(), file_name);
				}
			}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			//- init sample node
			Info_Phases& Phases = Solvers::get_instance()->parameters.Phases;
			sample_node.clear();
			int model_type = 0;
			if (InputFileReader::get_instance()->read_int_value("ModelsManager.PhiCon.BulkEnergy.type", model_type, false)) {
				if (model_type == dfdphiType::dfdphi_HighOrder) {
					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Preprocess.plotChemEnergy = {[(plot_phase_name),(fix_con_name, ... ),(fix_value, ...)], ... }\n", InputFileReader::get_instance()->debug_file);
					string chem_energy_key = "Preprocess.plotChemEnergy", chem_energy_input = "{[()]}";
					if (InputFileReader::get_instance()->read_string_value(chem_energy_key, chem_energy_input, infile_debug)) {
						is_plot_energy_curve = true;
						vector<InputValueType> chem_energy_structure; chem_energy_structure.push_back(InputValueType::IVType_STRING);
						chem_energy_structure.push_back(InputValueType::IVType_STRING); chem_energy_structure.push_back(InputValueType::IVType_DOUBLE);
						vector<vector<vector<input_value>>> chem_energy_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(chem_energy_structure, chem_energy_key, chem_energy_input, infile_debug);
						for (int phi_index = 0; phi_index < chem_energy_value.size(); phi_index++) {
							Info_Phase& defined_phase = Phases[chem_energy_value[phi_index][0][0].string_value];
							phases_fix_comp.add_tensor(defined_phase.phi_property);
							sample_node.add_phase(defined_phase.phi_property, defined_phase.phi_property, 0, 0.0);
							for (auto comp = defined_phase.x.begin(); comp < defined_phase.x.end(); comp++)
								sample_node[defined_phase.phi_property].x.add_con(comp->index, 0.0);
							if (chem_energy_value[phi_index][1].size() != chem_energy_value[phi_index][2].size()) {
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR, Preprocess.ChemEnergyCurve, number of fix cons and their values mismatch! \n", InputFileReader::get_instance()->debug_file);
								exit(0);
							}
							for (int con_index = 0; con_index < chem_energy_value[phi_index][1].size(); con_index++)
								for (auto comp = defined_phase.x.begin(); comp < defined_phase.x.end(); comp++)
									if (chem_energy_value[phi_index][1][con_index].string_value.compare(comp->name) == 0)
										phases_fix_comp(defined_phase.phi_property).add_double(comp->index, chem_energy_value[phi_index][2][con_index].double_value);
						}
						InputFileReader::get_instance()->read_double_value("Preprocess.PlotChemEnergy.con_epsilon", con_epsilon, infile_debug);
						if (con_epsilon < SYS_EPSILON)
							con_epsilon = SYS_EPSILON;
						InputFileReader::get_instance()->read_string_value("Preprocess.PlotChemEnergy.file_name", file_name, infile_debug);
					}
				}
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Chemical Energy Curve !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_plot_energy_curve) {
				Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Do chemcia energy curve plot:\n", LOG_FILE_NAME);
				write_chemical_energy_hlfuncs_to_m();
			}
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