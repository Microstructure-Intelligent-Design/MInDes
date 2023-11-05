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

namespace pf {
	namespace stiffness_eigenstrain {
		enum StiffnessEigenStrainType { ESType_None, ESType_PhaseDependent, ESType_GrainDependent };
		static tensor1_matrix6 phase_stiffness;
		static tensor1_matrix6 phi_stiffness;
		static tensor1_strain phase_eigen_strain;
		static Matrix6x6 get_phi_stiffness(PhaseEntry& phase) {
			return phi_stiffness(phase.index);
		}
		static tensor1_matrix6 get_stiffness() {
			return phi_stiffness;
		}
		static void do_grain_orientation(FieldStorage_forPhaseNode& phaseMesh) {
			for (auto phi = phaseMesh.info_node.begin(); phi < phaseMesh.info_node.end(); phi++)
				phi_stiffness.add_matrix6(phi->index, phase_stiffness(phi->property).get_rotated_matrix(Solvers::get_instance()->parameters.Grains.get_phi_rotationMatrix(phi->index)));
		}
		static void cal_eigenstrain_phi_dependent(PhaseNode& node, Vector6& eigenstrain) {
			for (auto phase = node.begin(); phase < node.end(); phase++) {
				if (phase->phi < SYS_EPSILON)
					continue;
				for (auto pES = phase_eigen_strain.begin(); pES < phase_eigen_strain.end(); pES++)
					if (pES->index == phase->property)
						eigenstrain += pES->val * phase->phi;
			}
		}
		static void cal_stiffness_phi_dependent(PhaseNode& node, Matrix6x6& stiffness) {
			for (auto phase = node.begin(); phase < node.end(); phase++) {
				if (phase->phi < SYS_EPSILON)
					continue;
				for (auto pC = phi_stiffness.begin(); pC < phi_stiffness.end(); pC++)
					if (pC->index == phase->index)
						stiffness += pC->val * phase->phi;
			}
		}

		namespace molar_volume {
			static double ref_molar_volume = 1e-5;
			static double_box phase_molar_volume;
			//-
			static StiffnessEigenStrainType SESType = StiffnessEigenStrainType::ESType_None;
			static bool debug_grain_eigenstrain = false;
			static int region_num = 0;
			static vector<vector<int>> region_phi_index_init;  // .[region_index] = phi_index_list
			static vector<double> region_ave_vm_init;  // .[region_index] = ave_vm_init
			static vector<vector<int>> region_phi_index;  // .[region_index] = phi_index_list
			static vector<double> region_ave_vm_current;  // .[region_index] = ave_vm_current
			//static vector<Matrix6x6> region_ave_stiffness;  // .[region_index] = ave_vm_current

			static void eigenstrain_node_dependent_molarvolume(PhaseNode& node, Vector6& eigenstrain) {
				for (auto phase = node.begin(); phase < node.end(); phase++) {
					if (phase->phi < SYS_EPSILON)
						continue;
					for (auto pVm = phase_molar_volume.begin(); pVm < phase_molar_volume.end(); pVm++)
						if (pVm->index == phase->property) {
							double eigenstrain_0 = pow(pVm->value / ref_molar_volume, 1.0 / 3.0) - 1.0;
							eigenstrain[0] += eigenstrain_0 * phase->phi;
							eigenstrain[1] += eigenstrain_0 * phase->phi;
							eigenstrain[2] += eigenstrain_0 * phase->phi;
						}
				}
			}
			static void eigenstrain_region_dependent_molarvolume(PhaseNode& node, Vector6& eigenstrain) {
				for (int region_index = 0; region_index < region_num; region_index++) {
					double eigenstrain_0 = pow(region_ave_vm_current[region_index] / region_ave_vm_init[region_index], 1.0 / 3.0) - 1.0;
					eigenstrain[0] += eigenstrain_0 * node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index];
					eigenstrain[1] += eigenstrain_0 * node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index];
					eigenstrain[2] += eigenstrain_0 * node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index];
				}
			}

			/*static void stiffness_region_dependent(PhaseNode& node, Matrix6x6& stiffness) {
				for (int region_index = 0; region_index < region_num; region_index++)
					stiffness += region_ave_stiffness[region_index] * node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index];
			}*/

			static void check_eigenstrain_regions(PhaseNode& example_node) {
				//- 1
				vector<int> phis;
				for (auto region = region_phi_index_init.begin(); region < region_phi_index_init.end(); region++) {
					for (auto phi = region->begin(); phi < region->end(); phi++) {
						for (auto save_phi = phis.begin(); save_phi < phis.end(); save_phi++)
							if (*phi == *save_phi) {
								Solvers::get_instance()->writer.add_string_to_txt_and_screen("ERROR, eigenstrain_regions, one phi used two times !\n", LOG_FILE_NAME);
								exit(0);
							}
						phis.push_back(*phi);
					}
				}
				if (region_phi_index_init.size() != region_phi_index.size() || region_phi_index_init.size() != region_num) {
					string error = "ERROR, eigenstrain_regions, region of init and loop size mismatch ! size of region init = " + to_string(region_phi_index_init.size()) + ", size of region loop = " + to_string(region_phi_index.size()) + "\n";
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(error, LOG_FILE_NAME);
					exit(0);
				}
				for (int region_index = 0; region_index < region_num; region_index++) {
					if (region_phi_index_init[region_index].size() == 0) {
						string error = "ERROR, eigenstrain_regions, region of init hasn't been defined ! region index = " + to_string(region_index) + "\n";
						Solvers::get_instance()->writer.add_string_to_txt_and_screen(error, LOG_FILE_NAME);
						exit(0);
					}
					if (region_phi_index[region_index].size() == 0) {
						string error = "ERROR, eigenstrain_regions, region of loop hasn't been defined ! region index = " + to_string(region_index) + "\n";
						Solvers::get_instance()->writer.add_string_to_txt_and_screen(error, LOG_FILE_NAME);
						exit(0);
					}
				}
			}
			static void init_regions(FieldStorage_forPhaseNode& phaseMesh, FieldStorage_forPhaseNode& buff_mesh) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							PhaseNode& buff_node = buff_mesh(x, y, z);
							for (int region_index = 0; region_index < region_num; region_index++) {
								node.customValues.add_double(ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index, 0.0);
								for (auto r_phi = region_phi_index_init[region_index].begin(); r_phi < region_phi_index_init[region_index].end(); r_phi++)
									for (auto phase = buff_node.begin(); phase < buff_node.end(); phase++)
										if (*r_phi == phase->index)
											node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index] += phase->phi;
							}
						}
			}
			static void init_reference_region_molar_volume(FieldStorage_forPhaseNode& phaseMesh, FieldStorage_forPhaseNode& buff_mesh) {
				region_ave_vm_init.clear();
				vector<double_box> region_phi_volume;
				for (int region_index = 0; region_index < region_num; region_index++) {
					region_ave_vm_init.push_back(0.0);
					double_box rphi;
					for (auto r_phi = region_phi_index[region_index].begin(); r_phi < region_phi_index[region_index].end(); r_phi++)
						rphi.add_double(*r_phi, 0.0);
					region_phi_volume.push_back(rphi);
				}
#pragma omp parallel for
				for (int x = 0; x < buff_mesh.limit_x; x++)
					for (int y = 0; y < buff_mesh.limit_y; y++)
						for (int z = 0; z < buff_mesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							PhaseNode& buff_node = buff_mesh(x, y, z);
							for (int region_index = 0; region_index < region_num; region_index++) {
								if (node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index] > SYS_EPSILON)
									for (auto phase = buff_node.begin(); phase < buff_node.end(); phase++)
										for (auto s_phi = region_phi_volume[region_index].begin(); s_phi < region_phi_volume[region_index].end(); s_phi++)
											if(s_phi->index == phase->index)
												s_phi->value += phase->phi;
							}
						}
				for (int region_index = 0; region_index < region_num; region_index++) {
					double sum_volume = 0.0, sum_molar_volume = 0.0;
					for (auto rphi = region_phi_volume[region_index].begin(); rphi < region_phi_volume[region_index].end(); rphi++) {
						sum_molar_volume += phase_molar_volume[phaseMesh.info_node[rphi->index].property] * rphi->value;
						sum_volume += rphi->value;
					}
					if (sum_volume < Phi_Num_Cut_Off) {
						cout << "> ERROR ! Region : " << region_index << " , phis don't exist !" << endl;
						SYS_PROGRAM_STOP;
					}
					else {
						region_ave_vm_init[region_index] = sum_molar_volume / sum_volume;
					}
				}
			}
			static void loop_grain_eigenstrain_model(FieldStorage_forPhaseNode& phaseMesh) {
				region_ave_vm_current.clear();
				//region_ave_stiffness.clear();
				vector<double_box> region_phi_volume;
				for (int region_index = 0; region_index < region_num; region_index++) {
					region_ave_vm_current.push_back(0.0);
					//Matrix6x6 C_0;
					//region_ave_stiffness.push_back(C_0);
					double_box rphi;
					for (auto r_phi = region_phi_index[region_index].begin(); r_phi < region_phi_index[region_index].end(); r_phi++)
						rphi.add_double(*r_phi, 0.0);
					region_phi_volume.push_back(rphi);
				}
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							for (int region_index = 0; region_index < region_num; region_index++) {
								if (node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index] > SYS_EPSILON)
									for (auto phase = node.begin(); phase < node.end(); phase++)
										for (auto s_phi = region_phi_volume[region_index].begin(); s_phi < region_phi_volume[region_index].end(); s_phi++)
											if (s_phi->index == phase->index)
												s_phi->value += phase->phi;
							}
						}
				for (int region_index = 0; region_index < region_num; region_index++) {
					double sum_volume = 0.0, sum_molar_volume = 0.0;
					for (auto rphi = region_phi_volume[region_index].begin(); rphi < region_phi_volume[region_index].end(); rphi++) {
						sum_molar_volume += phase_molar_volume[phaseMesh.info_node[rphi->index].property] * rphi->value;
						//region_ave_stiffness[region_index] += phi_stiffness(rphi->index) * rphi->value;
						sum_volume += rphi->value;
					}
					if (sum_volume < Phi_Num_Cut_Off) {
						cout << "> ERROR ! Region : " << region_index << " , phis don't exist !" << endl;
						SYS_PROGRAM_STOP;
					}
					else {
						region_ave_vm_current[region_index] = sum_molar_volume / sum_volume;
						//region_ave_stiffness[region_index] /= sum_volume;
					}
				}
			}
		}

		// function list, which is called in elastic solver
		static vector<void(*)(PhaseNode&, Vector6&)> eigenstrain_list;
		static void(*stiffness)(PhaseNode&, Matrix6x6&);

		static void cal_eigenstrain_physics(PhaseNode& node, Vector6& eigenstrain) {
			for (auto func = eigenstrain_list.begin(); func < eigenstrain_list.end(); func++)
				(*func)(node, eigenstrain);
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh, vector<int> solid_phases) {
			bool infile_debug = false, is_es_phi_dependent = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			stiffness = cal_stiffness_phi_dependent;
			for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
				bool is_solid = false;
				for (auto sphi = solid_phases.begin(); sphi < solid_phases.end(); sphi++)
					if (phi->phi_property == *sphi)
						is_solid = true;
				if (is_solid) {
					string stiffness_key = "Postprocess.SolidMechanics.Stiffness." + phi->phi_name, stiffness_input = "[(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0),(0,0,0,0,0,0)]";
					InputFileReader::get_instance()->read_string_value(stiffness_key, stiffness_input, infile_debug);
					vector<vector<input_value>> stiffness_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, stiffness_key, stiffness_input, infile_debug);
					Matrix6x6 C;
					for (int xx = 0; xx < 6; xx++)
						for (int yy = 0; yy < 6; yy++)
							C(xx, yy) = stiffness_value[xx][yy].double_value;
					phase_stiffness.add_matrix6(phi->phi_property, C);
					// (Exx, Eyy, Ezz, Eyz, Exz, Exy)
					string eigenstrain_key = "Postprocess.SolidMechanics.EigenStrain." + phi->phi_name, eigenstrain_input = "(0,0,0,0,0,0)";
					if (InputFileReader::get_instance()->read_string_value(eigenstrain_key, eigenstrain_input, infile_debug))
						is_es_phi_dependent = true;
					vector<input_value> eigenstrain_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, eigenstrain_key, eigenstrain_input, infile_debug);
					Vector6 E;
					for (int xx = 0; xx < 6; xx++)
						E[xx] = eigenstrain_value[xx].double_value;
					phase_eigen_strain.add_strain(phi->phi_property, E);
				}
				else {
					Matrix6x6 C;
					C.set_to_zero();
					phase_stiffness.add_matrix6(phi->phi_property, C);
					Vector6 E;
					E.set_to_zero();
					phase_eigen_strain.add_strain(phi->phi_property, E);
				}
			}
			if (is_es_phi_dependent)
				eigenstrain_list.push_back(cal_eigenstrain_phi_dependent);

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.StiffnessEigenStrain.model = (model_1, model_2, ...)   0 - Normal, 1 - PhaseDependent_MolarVolume , 2 - RegionDependent \n", InputFileReader::get_instance()->debug_file);
			string stiffness_eigenstrain_key = "Postprocess.SolidMechanics.StiffnessEigenStrain.model", stiffness_eigenstrain_input = "(0)";
			InputFileReader::get_instance()->read_string_value(stiffness_eigenstrain_key, stiffness_eigenstrain_input, infile_debug);
			vector<input_value> stiffness_eigenstrain_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, stiffness_eigenstrain_key, stiffness_eigenstrain_input, infile_debug);
			for (int index = 0; index < stiffness_eigenstrain_value.size(); index++) {
				if (stiffness_eigenstrain_value[index].int_value == StiffnessEigenStrainType::ESType_PhaseDependent && molar_volume::SESType == StiffnessEigenStrainType::ESType_None) {
					molar_volume::SESType = StiffnessEigenStrainType::ESType_PhaseDependent;
					eigenstrain_list.push_back(molar_volume::eigenstrain_node_dependent_molarvolume);
					for (int index = 0; index < solid_phases.size(); index++) {
						double Vm = 1e-5;
						InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.EigenStrain.MolarVolume." + Solvers::get_instance()->parameters.Phases[solid_phases[index]].phi_name, Vm, infile_debug);
						molar_volume::phase_molar_volume.add_double(solid_phases[index], Vm);
					}
					InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.EigenStrain.MolarVolume.reference", molar_volume::ref_molar_volume, infile_debug);
				}
				else if (stiffness_eigenstrain_value[index].int_value == StiffnessEigenStrainType::ESType_GrainDependent && molar_volume::SESType == StiffnessEigenStrainType::ESType_None) {
					molar_volume::SESType = StiffnessEigenStrainType::ESType_GrainDependent;
					eigenstrain_list.push_back(molar_volume::eigenstrain_region_dependent_molarvolume);
					//stiffness = molar_volume::stiffness_region_dependent;
					for (int index = 0; index < solid_phases.size(); index++) {
						double Vm = 1e-5;
						InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.EigenStrain.MolarVolume." + Solvers::get_instance()->parameters.Phases[solid_phases[index]].phi_name, Vm, infile_debug);
						molar_volume::phase_molar_volume.add_double(solid_phases[index], Vm);
					}
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.EigenStrain.MolarVolume.regions = [(region_1),(region_2), ... ], (region) = (phi_index, ... ) \n", InputFileReader::get_instance()->debug_file);
					string eigenstrain_region_key = "Postprocess.SolidMechanics.EigenStrain.MolarVolume.regions", eigenstrain_region_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(eigenstrain_region_key, eigenstrain_region_input, infile_debug)) {
						vector<vector<input_value>> eigenstrain_region_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, eigenstrain_region_key, eigenstrain_region_input, infile_debug);
						molar_volume::region_num = int(eigenstrain_region_value.size());
						for (int region_index = 0; region_index < molar_volume::region_num; region_index++) {
							vector<int> def_region;
							for (int phi_index = 0; phi_index < eigenstrain_region_value[region_index].size(); phi_index++)
								def_region.push_back(eigenstrain_region_value[region_index][phi_index].int_value);
							molar_volume::region_phi_index.push_back(def_region);
						}
					}
					bool is_read_datafile_by_path = false;
					string region_datafile_path = "DATA.dat", reference_datafile_path = "DATA.dat";
					if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.EigenStrain.MolarVolume.RegionsInit.datafile_path : relative path from infile folder.\n", InputFileReader::get_instance()->debug_file);
					if (InputFileReader::get_instance()->read_string_value("Postprocess.SolidMechanics.EigenStrain.MolarVolume.RegionsInit.datafile_path", region_datafile_path, infile_debug))
						is_read_datafile_by_path = true;
					region_datafile_path = Solvers::get_instance()->Infile_Folder_Path + dirSeparator + region_datafile_path;
					if (is_read_datafile_by_path) {
						if (infile_debug)
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.EigenStrain.MolarVolume.Reference.datafile_path : relative path from infile folder.\n", InputFileReader::get_instance()->debug_file);
						InputFileReader::get_instance()->read_string_value("Postprocess.SolidMechanics.EigenStrain.MolarVolume.Reference.datafile_path", reference_datafile_path, infile_debug);
						reference_datafile_path = Solvers::get_instance()->Infile_Folder_Path + dirSeparator + reference_datafile_path;
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.EigenStrain.MolarVolume.regions_init = [(region_1),(region_2), ... ], (region) = (phi_index, ... ) \n", InputFileReader::get_instance()->debug_file);
						string eigenstrain_region_init_key = "Postprocess.SolidMechanics.EigenStrain.MolarVolume.regions_init", eigenstrain_region_init_input = "[()]";
						if (InputFileReader::get_instance()->read_string_value(eigenstrain_region_init_key, eigenstrain_region_init_input, infile_debug)) {
							vector<vector<input_value>> eigenstrain_region_init_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, eigenstrain_region_init_key, eigenstrain_region_init_input, infile_debug);
							for (int region_index = 0; region_index < eigenstrain_region_init_value.size(); region_index++) {
								vector<int> def_region;
								for (int phi_index = 0; phi_index < eigenstrain_region_init_value[region_index].size(); phi_index++)
									def_region.push_back(eigenstrain_region_init_value[region_index][phi_index].int_value);
								molar_volume::region_phi_index_init.push_back(def_region);
							}
						}
						molar_volume::check_eigenstrain_regions(phaseMesh.info_node);
						FieldStorage_forPhaseNode buff_mesh;
						buff_mesh.init(phaseMesh.limit_x, phaseMesh.limit_y, phaseMesh.limit_z, phaseMesh.dr, phaseMesh._bc_x_up, phaseMesh._bc_y_up, phaseMesh._bc_z_up,
							phaseMesh._bc_x_down, phaseMesh._bc_y_down, phaseMesh._bc_z_down);
						Data_report buff_report;
						Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Stiffness & EigenStrain RegionDependent init Regions: \n", LOG_FILE_NAME);
						micro_structure_init::init_mesh_with_datafile(buff_mesh, buff_report, region_datafile_path, infile_debug);
						molar_volume::init_regions(phaseMesh, buff_mesh);
						Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Stiffness & EigenStrain RegionDependent init reference state: \n", LOG_FILE_NAME);
						micro_structure_init::init_mesh_with_datafile(buff_mesh, buff_report, reference_datafile_path, infile_debug);
						molar_volume::init_reference_region_molar_volume(phaseMesh, buff_mesh);
						buff_mesh.free();
					}
					else {
						molar_volume::region_phi_index_init = molar_volume::region_phi_index;
						molar_volume::check_eigenstrain_regions(phaseMesh.info_node);
						molar_volume::init_regions(phaseMesh, phaseMesh);
						molar_volume::init_reference_region_molar_volume(phaseMesh, phaseMesh);
					}
				}
			}

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Stiffness & EigenStrain \n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			do_grain_orientation(phaseMesh);
			if (molar_volume::SESType == StiffnessEigenStrainType::ESType_GrainDependent)
				molar_volume::loop_grain_eigenstrain_model(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			do_grain_orientation(phaseMesh);
			if (molar_volume::SESType == StiffnessEigenStrainType::ESType_GrainDependent)
				molar_volume::loop_grain_eigenstrain_model(phaseMesh);
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (molar_volume::debug_grain_eigenstrain) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "grain_eigenstrain_region" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i) {
							PhaseNode& node = phaseMesh(i, j, k);
							double region_index_sum = 0.0;
							for (int region_index = 0; region_index < molar_volume::region_num; region_index++)
								region_index_sum += region_index * node.customValues[ExternalFieldsPlus::EFP_GrainEigenStrain_Region + region_index];
							fout << region_index_sum << endl;
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}