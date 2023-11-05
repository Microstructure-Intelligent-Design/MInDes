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
#include "../sim_models/BulkEnergy.h"

namespace pf {
	namespace WritePhiCT {
		static bool phi = false;
		static bool con = false;
		static bool potential = false;
		static bool energy_density = false;
		static bool temperature = false;
		static bool phase_con = false;
		static bool phase_potential = false;
		static bool grains_rev = false;
		static bool phi_index = false;
		static bool phi_grad = false;
		static bool phi_name = false;
		static bool phi_summary = false;
		static bool interface_flag = false;
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false); 
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phi", phi, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.con", con, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.potential", potential, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.energy_density", energy_density, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.temperature", temperature, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phase_con", phase_con, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phase_potential", phase_potential, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.grains_rev", grains_rev, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phi_index", phi_index, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phi_gradient", phi_grad, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phi_name", phi_name, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.phi_summary", phi_summary, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.VTS.interface_flag", interface_flag, infile_debug);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (phi) {
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					string name;
					name = "\"Phi_" + to_string(phase->index) + "_" + Solvers::get_instance()->parameters.Phases[phase->property].phi_name + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << name <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; ++k)
						for (int j = 0; j < phaseMesh.limit_y; ++j)
							for (int i = 0; i < phaseMesh.limit_x; ++i) {
								PhaseNode& node = phaseMesh(i, j, k);
								double pp = 0.0;
								for (auto p = node.begin(); p < node.end(); p++)
									if (p->index == phase->index)
										pp = p->phi;
								if (IS_NAN(pp))
									fout << NaN << endl;
								else
									fout << pp << endl;
							}
					fout << "</DataArray>" << endl;
				}
			}
			if (con) {
				for (auto comp = phaseMesh.info_node.x.begin(); comp < phaseMesh.info_node.x.end(); comp++) {
					string name = "\"Con_" + to_string(comp->index) + "_" + Solvers::get_instance()->parameters.Components[comp->index].name + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << name <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; ++k)
						for (int j = 0; j < phaseMesh.limit_y; ++j)
							for (int i = 0; i < phaseMesh.limit_x; ++i) {
								PhaseNode& node = phaseMesh(i, j, k);
								double con = 0.0;
								for (auto x = node.x.begin(); x < node.x.end(); x++)
									if (x->index == comp->index)
										con = x->value;
								if (IS_NAN(con))
									fout << NaN << endl;
								else
									fout << con << endl;
							}
					fout << "</DataArray>" << endl;
				}
			}
			if (potential) {
				for (auto comp = phaseMesh.info_node.x.begin(); comp < phaseMesh.info_node.x.end(); comp++) {
					string name = "\"Potential_" + to_string(comp->index) + "_" + Solvers::get_instance()->parameters.Components[comp->index].name + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << name <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++) {
								PhaseNode& node = phaseMesh(i, j, k);
								double pot = 0.0;
								for (auto p = node.potential.begin(); p < node.potential.end(); p++)
									if (p->index == comp->index)
										pot = p->value;
								if (IS_NAN(pot))
									fout << NaN << endl;
								else
									fout << pot << endl;
							}
					fout << "</DataArray>" << endl;
				}
			}
			if (energy_density) {
				pf::ConEquationType _type = Solvers::get_instance()->parameters.ConEType;
				pf::ConEquationDomain _domain = Solvers::get_instance()->parameters.ConEDomain;
				double threshold = Solvers::get_instance()->C_Solver.threshold;
				string name = "\"energy_density\" ";
				fout << "<DataArray type = \"Float64\" Name = " << name <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double energy_density = 0.0;
							if (_type == ConEquationType::CEType_PhaseX) {
								for (auto p = node.begin(); p < node.end(); p++)
									if (p->phi > Phi_Num_Cut_Off)
										energy_density += p->phi * bulk_energy::bulk_energy_density(node, *p);
							}
							else if (_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
								for (auto p = node.begin(); p < node.end(); p++)
									if (p->phi > Phi_Num_Cut_Off)
										energy_density += p->phi * bulk_energy::bulk_energy_density(node, *p);
							}
							else if (_domain == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
								for (auto p = node.begin(); p < node.end(); p++)
									if (p->phi > Phi_Num_Cut_Off)
										energy_density += p->phi * bulk_energy::bulk_energy_density(node, *p);
							}
							if (IS_NAN(energy_density))
								fout << NaN << endl;
							else
								fout << energy_density << endl;
						}
				fout << "</DataArray>" << endl;
			}
			if (temperature) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "temperature" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							double temp = phaseMesh(i, j, k).temperature.T;
							if (IS_NAN(temp))
								fout << NaN << endl;
							else
								fout << temp << endl;
						}
				fout << "</DataArray>" << endl;
			}
			if (phase_con) {
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
						string name = "\"Phi_" + to_string(phase->index) + "_Con_" + to_string(comp->index) + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << name <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < phaseMesh.limit_z; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i) {
									PhaseNode& node = phaseMesh(i, j, k);
									if (node[phase->index].phi > Phi_Num_Cut_Off) {
										double con = node[phase->index].x[comp->index].value;
										if (IS_NAN(con))
											fout << NaN << endl;
										else
											fout << con << endl;
									}
									else {
										fout << 0.0 << endl;
									}
								}
						fout << "</DataArray>" << endl;
					}
				}
			}
			if (phase_potential) {
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
						string name = "\"Phi_" + to_string(phase->index) + "_Potential_" + to_string(comp->index) + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << name <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < phaseMesh.limit_z; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i) {
									PhaseNode& node = phaseMesh(i, j, k);
									if (node[phase->index].phi > Phi_Num_Cut_Off) {
										double con = node[phase->index].potential[comp->index].value;
										if (IS_NAN(con))
											fout << NaN << endl;
										else
											fout << con << endl;
									}
									else {
										fout << 0.0 << endl;
									}
								}
						fout << "</DataArray>" << endl;
					}
				}
			}
			if (grains_rev) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "Grains_Rev" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double fix = 0.0;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								fix += phase->phi * phase->phi;
							if (IS_NAN(fix))
								fout << NaN << endl;
							else
								fout << 1.0 - fix << endl;
						}
				fout << "</DataArray>" << endl;
			}
			if (phi_index) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "Phi_Indexs" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double fix = 0.0;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								fix += phase->index * phase->phi;
							fout << fix << endl;
						}
				fout << "</DataArray>" << endl;
			}
			if (phi_name) {
				Info_Phases& phases = Solvers::get_instance()->parameters.Phases;
				fout << "<DataArray type = \"Float64\" Name = \"" << "Phi_Properties" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double fix = 0.0;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								fix += phase->property * phase->phi;
							fout << fix << endl;
						}
				fout << "</DataArray>" << endl;
				for (auto phi_property = phases.begin(); phi_property < phases.end(); phi_property++) {
					fout << "<DataArray type = \"Float64\" Name = \"" << "Phi_" + phi_property->phi_name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++) {
								PhaseNode& node = phaseMesh(i, j, k);
								double fix = 0.0;
								for (auto phase = node.begin(); phase < node.end(); phase++)
									if (phase->property == phi_property->phi_property)
										fix += phase->phi;
								fout << fix << endl;
							}
					fout << "</DataArray>" << endl;
				}
			}
			if (phi_summary) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "Phi_Summary" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double fix = 0.0;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								fix += phase->phi;
							fout << fix << endl;
						}
				fout << "</DataArray>" << endl;
			}
			if (interface_flag) {
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					string name = "\"Phi" + to_string(phase->index) + "_Flag\" ";
					fout << "<DataArray type = \"Float64\" Name = " << name <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; ++k)
						for (int j = 0; j < phaseMesh.limit_y; ++j)
							for (int i = 0; i < phaseMesh.limit_x; ++i) {
								PhaseNode& node = phaseMesh(i, j, k);
								int flag = pf_BULK;
								for (auto p = node.begin(); p < node.end(); p++)
									if (p->index == phase->index)
										flag = p->_flag;
								fout << flag << endl;
							}
					fout << "</DataArray>" << endl;
				}
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (phi_grad) {
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					string name;
					name = "\"Phi_" + to_string(phase->index) + "_Gradient\" ";
					fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++) {
								PhaseNode& node = phaseMesh(i, j, k);
								for (auto p = node.begin(); p < node.end(); p++) {
									if (p->index == phase->index)
										fout << p->phi_grad[0] << " " << p->phi_grad[1] << " " << p->phi_grad[2] << endl;
									else
										fout << 0.0 << " " << 0.0 << " " << 0.0 << endl;
								}
							}
					fout << "</DataArray>" << endl;
				}
			}
		}
	}
}