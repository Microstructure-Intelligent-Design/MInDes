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
#include "../sim_preprocess/MicroStructureInit.h"
#include "../sim_models/BulkEnergy.h"

namespace pf {
	namespace noise {
		enum NoiseType { NT_UNIFORM, NT_INTERFACE };
		static NoiseType noise_type = NoiseType::NT_UNIFORM;
		// - Phi
		static int phi_noise_steps = 10000;
		static int phi_noise_min_steps = 0;
		static int phi_noise_max_steps = 100000000;
		static double phi_noise_amplitude = 0.0;
		static tensor2_int noise_phi_properties;   // val[applied_phi_index][generated_phi_index] = generated_phi_property
		static bool is_noise_step_0 = false;
		static vector<double> noise_range;
		static bool is_noise_on_phi = true;
		// -
		static PhaseNode buff_node;
		// -
		static double driving_force_threshold = SYS_EPSILON;

		static void default_noise_func(FieldStorage_forPhaseNode& phaseMesh) {
			return;
		}

		static void init_noise_phi_property(FieldStorage_forPhaseNode& phaseMesh) {
			pf::Info_Phases& phases = Solvers::get_instance()->parameters.Phases;
			for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++)
				for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++) {
					//- generated_phi_property  generated_phi_index  should be checked in mesh
					for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++)
						if (phase->index == gene_phi->index && phase->property != gene_phi->val) {
							stringstream _report;
							_report << "> Error : Postprocess.Phi.noise, phase index mismatches with its name, index = " << gene_phi->index;
							_report << ", phase name = " << phases[gene_phi->val].phi_name << " ! " << endl;
							Solvers::get_instance()->writer.add_string_to_txt_and_screen(_report.str(), LOG_FILE_NAME);
							exit(0);
						}
					//- init this phi
					micro_structure_init::add_new_phi_by_index_with_index(gene_phi->index, gene_phi->val, app_phi->index, 0.0);
				}
		}

		static double cal_abs_distribution_grand_potential(pf::PhaseNode& node, pf::PhaseEntry& app_phase, int gene_phase_index) {
			if(!is_noise_on_phi && node[gene_phase_index].phi > SYS_EPSILON)
				return 0.0;
			PhaseEntry& buff_phase = buff_node[gene_phase_index];
			buff_phase.potential = app_phase.potential;
			chemical_energy::phase_con(node, buff_phase);
			double drivingForce = bulk_energy::dfbulk_dphi(node, app_phase) - bulk_energy::dfbulk_dphi(node, buff_phase);
			if (drivingForce > driving_force_threshold) {
				//- generate phase will grows
				return RAND_0_1 * (noise_range[1] - noise_range[0]) + noise_range[0];
			}
			else {
				return 0.0;
			}
		}

		static void cal_noise_grand_potential_multiphase_uniform(FieldStorage_forPhaseNode& phaseMesh) {
			pf::Info_Phases& phases = Solvers::get_instance()->parameters.Phases;
			pf::ConEquationDomain _domain = Solvers::get_instance()->parameters.ConEDomain;
			double threshold = Solvers::get_instance()->C_Solver.threshold;
			bool is_normalized = Solvers::get_instance()->parameters.is_Normalize_Phi;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						if ((_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) ||
							(_domain == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] > threshold)) {
							continue;
						}
						else {
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								if (phase->phi < SYS_EPSILON)
									continue;
								for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++)
									if (phase->index == app_phi->index) {
										for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++) {
											// calculate distribution
											double distribution = cal_abs_distribution_grand_potential(node, *phase, gene_phi->index);
											// calculate interpolation
											double interpolation = phase->phi;
											// calculate noise
											double noise = interpolation * distribution * phi_noise_amplitude;
											// cout << "> interpolation = " << interpolation << ", distribution = " << distribution << ", phi_noise_amplitude = " << phi_noise_amplitude << endl;
											// adjust phi for generate phase and applied phase
											node[gene_phi->index].phi += noise;
											phase->phi -= noise;
										}
										if (is_normalized) {
											double others_phi = 0.0, aim_phi = 0.0;
											if (phase->phi < 0.0)
												phase->phi = 0.0;
											else if (phase->phi > 1.0)
												phase->phi = 1.0;
											aim_phi = phase->phi;
											for (auto phase2 = node.begin(); phase2 < node.end(); phase2++) {
												bool is_gene_phase = false;
												for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
													if (phase2->index == gene_phi->index)
														is_gene_phase = true;
												if (is_gene_phase) {
													if (phase2->phi < 0.0)
														phase2->phi = 0.0;
													else if (phase2->phi > 1.0)
														phase2->phi = 1.0;
													aim_phi += phase2->phi;
												}
												else if (phase2->index != phase->index) {
													others_phi += phase2->phi;
												}
											}
											if (aim_phi < SYS_EPSILON) {
												int phi_num = app_phi->size() + 1;
												phase->phi = (1.0 - others_phi) / phi_num;
												for (auto phase2 = node.begin(); phase2 < node.end(); phase2++) {
													for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
														if (phase2->index == gene_phi->index)
															phase2->phi = (1.0 - others_phi) / phi_num;
												}
											}
											else {
												phase->phi *= (1.0 - others_phi) / aim_phi;
												for (auto phase2 = node.begin(); phase2 < node.end(); phase2++) {
													for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
														if (phase2->index == gene_phi->index)
															phase2->phi *= (1.0 - others_phi) / aim_phi;
												}
											}
										}
									}
							}
						}
					}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						for (auto phase = node.begin(); phase < node.end(); phase++)
							phase->_flag = phaseMesh.currentFlag(node, phase->index);
					}
		}

		static void cal_noise_grand_potential_multiphase_interface(FieldStorage_forPhaseNode& phaseMesh) {
			pf::Info_Phases& phases = Solvers::get_instance()->parameters.Phases;
			pf::ConEquationDomain _domain = Solvers::get_instance()->parameters.ConEDomain;
			double threshold = Solvers::get_instance()->C_Solver.threshold;
			bool is_normalized = Solvers::get_instance()->parameters.is_Normalize_Phi;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						if ((_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) ||
							(_domain == ConEquationDomain::CEDomain_Reverse  && node.customValues[ExternalFields::CON_Smooth_Phi] > threshold)) {
							continue;
						}
						else {
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								if (phase->phi < SYS_EPSILON)
									continue;
								for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++)
									if (phase->index == app_phi->index) {
										for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++) {
											// calculate distribution
											double distribution = cal_abs_distribution_grand_potential(node, *phase, gene_phi->index);
											// calculate interpolation
											double interpolation = 0.0;
											for (auto phase2 = node.begin(); phase2 < node.end(); phase2++)
												if (phase2->index != app_phi->index && phase2->index != gene_phi->index && phase2->phi > SYS_EPSILON) {
													interpolation += phase->phi * phase2->phi;
												}
											// calculate noise
											double noise = interpolation * distribution * phi_noise_amplitude;
											// cout << "> interpolation = " << interpolation << ", distribution = " << distribution << ", phi_noise_amplitude = " << phi_noise_amplitude << endl;
											// adjust phi for generate phase and applied phase
											node[gene_phi->index].phi += noise;
											phase->phi -= noise;
										}
										if (is_normalized) {
											double others_phi = 0.0, aim_phi = 0.0;
											if (phase->phi < 0.0)
												phase->phi = 0.0;
											else if (phase->phi > 1.0)
												phase->phi = 1.0;
											aim_phi = phase->phi;
											for (auto phase2 = node.begin(); phase2 < node.end(); phase2++) {
												bool is_gene_phase = false;
												for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
													if (phase2->index == gene_phi->index)
														is_gene_phase = true;
												if (is_gene_phase) {
													if (phase2->phi < 0.0)
														phase2->phi = 0.0;
													else if (phase2->phi > 1.0)
														phase2->phi = 1.0;
													aim_phi += phase2->phi;
												}
												else if(phase2->index != phase->index){
													others_phi += phase2->phi;
												}
											}
											if (aim_phi < SYS_EPSILON) {
												int phi_num = app_phi->size() + 1;
												phase->phi = (1.0 - others_phi) / phi_num;
												for (auto phase2 = node.begin(); phase2 < node.end(); phase2++) {
													for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
														if (phase2->index == gene_phi->index)
															phase2->phi = (1.0 - others_phi) / phi_num;
												}
											}
											else {
												phase->phi *= (1.0 - others_phi) / aim_phi;
												for (auto phase2 = node.begin(); phase2 < node.end(); phase2++) {
													for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
														if (phase2->index == gene_phi->index)
															phase2->phi *= (1.0 - others_phi) / aim_phi;
												}
											}
										}
									}
							}
						}
					}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						for (auto phase = node.begin(); phase < node.end(); phase++)
							phase->_flag = phaseMesh.currentFlag(node, phase->index);
					}
		}

		static void (*init_noise)(FieldStorage_forPhaseNode&) ;

		static void(*generate_noise)(FieldStorage_forPhaseNode&);

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			init_noise = default_noise_func;
			generate_noise = default_noise_func;
			noise_range.resize(2);
			noise_range[0] = 0.0;
			noise_range[1] = 1.0;
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise &&
				Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				bool infile_debug = false;
				InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
				if(infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Phi.interface_noise = {[(applied_phi_index),(generate_phi_index_0, ... ),(generate_phi_name_0, ... )], ... } \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                .uniform_noise   = {[(applied_phi_index),(generate_phi_index_0, ... ),(generate_phi_name_0, ... )], ... } \n", InputFileReader::get_instance()->debug_file);
				string int_noise_key = "Postprocess.Phi.interface_noise", uni_noise_key = "Postprocess.Phi.uniform_noise", int_noise_input = "{[()]}", uni_noise_input = "{[()]}";
				if (InputFileReader::get_instance()->read_string_value(int_noise_key, int_noise_input, infile_debug)) {
					vector<InputValueType> int_noise_structure; int_noise_structure.push_back(InputValueType::IVType_INT);
					int_noise_structure.push_back(InputValueType::IVType_INT); int_noise_structure.push_back(InputValueType::IVType_STRING);
					vector<vector<vector<input_value>>> int_noise_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(int_noise_structure, int_noise_key, int_noise_input, infile_debug);
					for (int app_index = 0; app_index < int_noise_value.size(); app_index++) {
						if (int_noise_value[app_index][1].size() != int_noise_value[app_index][2].size()) {
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Phi.interface_noise error, generate phi and its name mismatch ! \n", InputFileReader::get_instance()->debug_file);
							Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Postprocess.Phi.interface_noise error, generate phi and its name mismatch ! \n", LOG_FILE_NAME);
							exit(0);
						}
						for (int gene_index = 0; gene_index < int_noise_value[app_index][1].size(); gene_index++) {
							noise_phi_properties.add_int(int_noise_value[app_index][0][0].int_value, int_noise_value[app_index][1][gene_index].int_value,
								Solvers::get_instance()->parameters.Phases[int_noise_value[app_index][2][gene_index].string_value].phi_property);
						}
					}
					noise_type = NoiseType::NT_INTERFACE;
					init_noise = init_noise_phi_property;
					generate_noise = cal_noise_grand_potential_multiphase_interface;
				}
				if (InputFileReader::get_instance()->read_string_value(uni_noise_key, uni_noise_input, infile_debug)) {
					vector<InputValueType> uni_noise_structure; uni_noise_structure.push_back(InputValueType::IVType_INT);
					uni_noise_structure.push_back(InputValueType::IVType_INT); uni_noise_structure.push_back(InputValueType::IVType_STRING);
					vector<vector<vector<input_value>>> uni_noise_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(uni_noise_structure, uni_noise_key, uni_noise_input, infile_debug);
					for (int app_index = 0; app_index < uni_noise_value.size(); app_index++) {
						if (uni_noise_value[app_index][1].size() != uni_noise_value[app_index][2].size()) {
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Phi.uniform_noise error, generate phi and its name mismatch ! \n", InputFileReader::get_instance()->debug_file);
							Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Postprocess.Phi.uniform_noise error, generate phi and its name mismatch ! \n", LOG_FILE_NAME);
							exit(0);
						}
						for (int gene_index = 0; gene_index < uni_noise_value[app_index][1].size(); gene_index++) {
							noise_phi_properties.add_int(uni_noise_value[app_index][0][0].int_value, uni_noise_value[app_index][1][gene_index].int_value,
								Solvers::get_instance()->parameters.Phases[uni_noise_value[app_index][2][gene_index].string_value].phi_property);
						}
					}
					noise_type = NoiseType::NT_UNIFORM;
					init_noise = init_noise_phi_property;
					generate_noise = cal_noise_grand_potential_multiphase_uniform;
				}

				InputFileReader::get_instance()->read_bool_value("Postprocess.Phi.noise.is_start", is_noise_step_0, infile_debug);

				InputFileReader::get_instance()->read_double_value("Postprocess.Phi.noise.driving_force_threshold", driving_force_threshold, infile_debug);

				InputFileReader::get_instance()->read_bool_value("Postprocess.Phi.noise.is_noise_on_phi", is_noise_on_phi, infile_debug);

				InputFileReader::get_instance()->read_int_value("Postprocess.Phi.noise.frequency", phi_noise_steps, infile_debug);
				if (phi_noise_steps < 1)
					phi_noise_steps = 1;

				InputFileReader::get_instance()->read_int_value("Postprocess.Phi.noise.min_step", phi_noise_min_steps, infile_debug);

				InputFileReader::get_instance()->read_int_value("Postprocess.Phi.noise.max_step", phi_noise_max_steps, infile_debug);

				InputFileReader::get_instance()->read_double_value("Postprocess.Phi.noise.amplitude", phi_noise_amplitude, infile_debug);
				if (phi_noise_amplitude < 0.0)
					phi_noise_amplitude = 0.0;

				string noise_range_input = "(0.0,1.0)", noise_range_key = "Postprocess.Phi.noise.random";
				if (InputFileReader::get_instance()->read_string_value(noise_range_key, noise_range_input, infile_debug)) {
					vector<input_value> noise_range_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, noise_range_key, noise_range_input, infile_debug);
					noise_range[0] = noise_range_value[0].double_value;
					noise_range[1] = noise_range_value[1].double_value;
				}
				if (noise_range[1] < noise_range[0])
					noise_range[1] = noise_range[0];

			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Noise !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			init_noise(phaseMesh);
			buff_node = phaseMesh.info_node;
			if (is_noise_step_0)
				generate_noise(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			if (Solvers::get_instance()->current_istep % phi_noise_steps == 0 && (Solvers::get_instance()->current_istep >= phi_noise_min_steps && Solvers::get_instance()->current_istep <= phi_noise_max_steps)) {
				generate_noise(phaseMesh);
				report = "> Add noise to phase-field, current step = " + to_string(Solvers::get_instance()->current_istep) + "\n";
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(report, LOG_FILE_NAME);
			}
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			pf::ConEquationDomain _domain = Solvers::get_instance()->parameters.ConEDomain;
			double threshold = Solvers::get_instance()->C_Solver.threshold;
			pf::PhaseNode& sample_node = phaseMesh.info_node;
			pf::Info_Phases& phases = Solvers::get_instance()->parameters.Phases;
			int_box output_phis;
			for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++) {
				output_phis.add_int(sample_node[app_phi->index].property, 1);
			}
			for (auto out_phi = output_phis.begin(); out_phi < output_phis.end(); out_phi++) {
				string distribution_name = "\"Phi_" + phases[out_phi->index].phi_name + "_noise_distribution\"";
				fout << "<DataArray type = \"Float64\" Name = " << distribution_name <<
					" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i) {
							PhaseNode& node = phaseMesh(i, j, k);
							if ((_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) ||
								(_domain == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] > threshold)) {
								fout << 0.0 << endl;
							}
							else {
								double merge_distribution = 0.0, merge_num = 0.0;
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									if (phase->phi > SYS_EPSILON && phase->property == out_phi->index) {
										for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++)
											if (phase->index == app_phi->index) {
												for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++) {
													// calculate distribution
													merge_distribution += cal_abs_distribution_grand_potential(node, *phase, gene_phi->index);
													merge_num += 1.0;
												}
											}
									}
								}
								if (merge_num > 0.5)
									fout << merge_distribution / merge_num << endl;
								else
									fout << 0.0 << endl;
							}
						}
				fout << "</DataArray>" << endl;
			}
			output_phis.clear();
			for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++) {
				for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++) {
					output_phis.add_int(gene_phi->val, 1); // save property
				}
			}
			for (auto out_phi = output_phis.begin(); out_phi < output_phis.end(); out_phi++) {
				string noise_name = "\"Phi_" + phases[out_phi->index].phi_name + "_noise\"";
				fout << "<DataArray type = \"Float64\" Name = " << noise_name <<
					" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i) {
							PhaseNode& node = phaseMesh(i, j, k);
							if ((_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) ||
								(_domain == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] > threshold)) {
								fout << 0.0 << endl;
							}
							else if (Solvers::get_instance()->current_istep % phi_noise_steps != 0) {
								fout << 0.0 << endl;
							}
							else {
								double merge_noise = 0.0, merge_num = 0.0;
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									if (phase->phi > SYS_EPSILON) {
										for (auto app_phi = noise_phi_properties.begin(); app_phi < noise_phi_properties.end(); app_phi++)
											if (phase->index == app_phi->index)
												for (auto gene_phi = app_phi->begin(); gene_phi < app_phi->end(); gene_phi++)
													if (gene_phi->val == out_phi->index) {
														// calculate distribution
														double distribution = cal_abs_distribution_grand_potential(node, *phase, gene_phi->index);
														// calculate interpolation
														double interpolation = 0.0;
														if (noise_type == NoiseType::NT_UNIFORM)
															interpolation = phase->phi;
														else if (noise_type == NoiseType::NT_INTERFACE) {
															for (auto phase2 = node.begin(); phase2 < node.end(); phase2++)
																if (phase2->index != app_phi->index && phase2->index != gene_phi->index && phase2->phi > SYS_EPSILON) {
																	interpolation += phase->phi * phase2->phi;
																}
														}
														// calculate noise
														merge_noise += interpolation * distribution * phi_noise_amplitude;
														merge_num += 1.0;
													}
									}
								}
								if (merge_num > 0.5)
									fout << merge_noise / merge_num << endl;
								else
									fout << 0.0 << endl;
							}
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}