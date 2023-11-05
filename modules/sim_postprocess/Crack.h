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
#include "../sim_postprocess/Mechanics/PlasticSolver.h"

namespace pf {
	namespace crack_propagation {
		enum CrackPropagationModel { CPM_None, CPM_Single_Order_Parameter, CPM_Multiple_Order_Parameter };
		static CrackPropagationModel crack_model = CrackPropagationModel::CPM_None;
		static DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT;
		//-
		static int		solve_crack_steps = 1;
		static int		max_iterate_steps = 0;
		static double		damage_threshold = 0.9;
		static double		crack_epsilon = 1e-3;
		static double		crack_dt = 1.0;
		static double    crack_int_width = 2.0;
		static double    crack_dr = 1.0;
		static double_box	pGc;
		//- noise
		enum CrackNoiseType { CNT_UNIFORM, CNT_INTERFACE };
		static CrackNoiseType crack_noise_type = CrackNoiseType::CNT_UNIFORM;
		static int crack_noise_steps = 10000;
		static double crack_noise_amplitude = 0.0;
		static bool is_noise_start = false;
		static vector<double> crack_noise_range;
		static bool is_noise_on_crack = true;
		static bool is_dfcrack_dphi_output = false;
		//- initial
		static bool is_crack_init_from_datafile = false;
		static string crack_init_datafile = "";

		static CrackPropagationModel crack_propagation_model() {
			return crack_model;
		}

		static int crack_solver_max_iterate_steps() {
			return max_iterate_steps;
		}

		static double crack_solver_epsilon() {
			return crack_epsilon;
		}

		static double crack_fraction_single(pf::PhaseNode& node) {
			return node.customValues[ExternalFieldsPlus::EFP_Crack];
		}

		static double crack_fraction_multiple(pf::PhaseNode& node, pf::PhaseEntry phase) {
			return node.customValues[ExternalFieldsPlus::EFP_Crack + phase.index];
		}

		static void normalize_crack(double& crack) {
			if (crack > 1.0)
				crack = 1.0;
			else if (crack < 0.0)
				crack = 0.0;
		}

		static void init_crack_field_default(FieldStorage_forPhaseNode& phaseMesh) {
			bool is_crack_init_by_datafile = false;
			if (micro_structure_init::is_init_by_datefile()) {
				pf::Data_report data_report = micro_structure_init::get_datafile_info();
				for (auto custom_var = data_report.custom_vars.begin(); custom_var < data_report.custom_vars.end(); custom_var++)
					if (custom_var->index == ExternalFieldsPlus::EFP_Crack)
						is_crack_init_by_datafile = true;
			}
			if (is_crack_init_by_datafile) {
				if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								node.customValues.add_double(ExternalFieldsPlus::EFP_Crack_Incre, 0.0);
							}
				}
				else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									node.customValues.add_double(ExternalFieldsPlus::EFP_Crack_Incre + phase->index, 0.0);
								}
							}
				}
			}
			else {
				if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								node.customValues.add_double(ExternalFieldsPlus::EFP_Crack, 0.0);
								node.customValues.add_double(ExternalFieldsPlus::EFP_Crack_Incre, 0.0);
							}
				}
				else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									node.customValues.add_double(ExternalFieldsPlus::EFP_Crack + phase->index, 0.0);
									node.customValues.add_double(ExternalFieldsPlus::EFP_Crack_Incre + phase->index, 0.0);
								}
							}
				}
			}
		}

		static void init_crack_field_from_datafile(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "> Init crack with datafile: " + crack_init_datafile + "\n";
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(report, LOG_FILE_NAME);
			string datafile_path = Solvers::get_instance()->Infile_Folder_Path + dirSeparator + crack_init_datafile;
			Data_report buff_report;
			FieldStorage_forPhaseNode buff_mesh;
			buff_mesh.init(phaseMesh.limit_x, phaseMesh.limit_y, phaseMesh.limit_z, phaseMesh.dr, phaseMesh._bc_x_up, phaseMesh._bc_y_up, phaseMesh._bc_z_up, phaseMesh._bc_x_down, phaseMesh._bc_y_down, phaseMesh._bc_z_down);
			micro_structure_init::init_mesh_with_datafile(buff_mesh, buff_report, datafile_path, false);
			if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							PhaseNode& buff_node = buff_mesh(x, y, z);
							node.customValues.add_double(ExternalFieldsPlus::EFP_Crack, buff_node.customValues[ExternalFieldsPlus::EFP_Crack]);
							node.customValues.add_double(ExternalFieldsPlus::EFP_Crack_Incre, 0.0);
						}
			}
			else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							PhaseNode& buff_node = buff_mesh(x, y, z);
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								node.customValues.add_double(ExternalFieldsPlus::EFP_Crack + phase->index, buff_node.customValues[ExternalFieldsPlus::EFP_Crack + phase->index]);
								node.customValues.add_double(ExternalFieldsPlus::EFP_Crack_Incre + phase->index, 0.0);
							}
						}
			}
		}

		static double Gc(pf::PhaseNode& node) {
			double gc = 0.0, sum_h = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++) {
				gc += phase->phi * pGc[phase->property];
			}
			return gc;
		}

		namespace crack_propagation_models {
			const double K = 9.0 / 64.0;
			static double df_dcrack_default(pf::PhaseNode& node) {
				return 0.0;
			}
			static double dfint_dcrack_single(pf::PhaseNode& node) {
				// obstacle
				//return Gc(node) * (K / crack_int_width - crack_int_width * 2.0 * node.cal_customValues_laplace(ExternalFieldsPlus::EFP_Crack, crack_dr, diff_method));
				// well
				return Gc(node) * (node.customValues[ExternalFieldsPlus::EFP_Crack] / crack_int_width * 2.0
					- crack_int_width * 2.0 * node.cal_customValues_laplace(ExternalFieldsPlus::EFP_Crack, crack_dr, diff_method));
			}

			static double dfelas_dcrack_single(pf::PhaseNode& node) {
				double phi_crack = node.customValues[ExternalFieldsPlus::EFP_Crack];
				return - /*2.0 * (1.0 - phi_crack) * */(node.customVec6s[ExternalFields::MECH_stress]
					* (node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain]));
			}

			static double dfmech_dcrack_single(pf::PhaseNode& node) {
				double phi_crack = node.customValues[ExternalFieldsPlus::EFP_Crack], ave_plastic_strain = node.customValues[ExternalFields::MECH_ave_plastic_strain];
				return - /*2.0 * (1.0 - phi_crack) * */(node.customVec6s[ExternalFields::MECH_stress] * (node.customVec6s[ExternalFields::MECH_strain]
					- node.customVec6s[ExternalFields::MECH_eigen_strain] - node.customVec6s[ExternalFields::MECH_plastic_strain]) +
					plastic_solver::get_hardening_modulus(node) * ave_plastic_strain * ave_plastic_strain);
			}
		}

		static double (*dfint_dcrack)(pf::PhaseNode& node);

		static double (*dfmech_dcrack)(pf::PhaseNode& node);

		static void noise_crack(FieldStorage_forPhaseNode& phaseMesh) {
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						// calculate distribution
						double distribution = RAND_0_1 * (crack_noise_range[1] - crack_noise_range[0]) + crack_noise_range[0];
						// calculate interpolation
						double interpolation = 0.0;
						if (crack_noise_type == CrackNoiseType::CNT_UNIFORM) {
							interpolation = 1.0;
						}
						else if (crack_noise_type == CrackNoiseType::CNT_INTERFACE) {
							for (auto alpha = node.begin(); alpha < node.end() - 1; alpha++)
								for (auto beta = alpha + 1; beta < node.end(); beta++)
									interpolation += alpha->phi * beta->phi;
						}
						// calculate noise
						double noise = interpolation * distribution * crack_noise_amplitude;

						if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
							if (!is_noise_on_crack && node.customValues[ExternalFieldsPlus::EFP_Crack] > SYS_EPSILON)
								noise = 0.0;
							node.customValues[ExternalFieldsPlus::EFP_Crack] += noise;
							normalize_crack(node.customValues[ExternalFieldsPlus::EFP_Crack]);
						}
						else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {
							for (auto phase = node.begin(); phase < node.end(); phase++)
								if (phase->phi > SYS_EPSILON) {
									if (!is_noise_on_crack && node.customValues[ExternalFieldsPlus::EFP_Crack + phase->index] > SYS_EPSILON)
										noise = 0.0;
									node.customValues[ExternalFieldsPlus::EFP_Crack + phase->index] += noise;
									normalize_crack(node.customValues[ExternalFieldsPlus::EFP_Crack + phase->index]);
								}
						}
					}
		}

		static double crack_loop_single(FieldStorage_forPhaseNode& phaseMesh) {
			double MAX_VARIATION = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.customValues[ExternalFieldsPlus::EFP_Crack_Incre] = -(dfint_dcrack(node) + dfmech_dcrack(node));
						if (node.customValues[ExternalFieldsPlus::EFP_Crack] > damage_threshold)
							if (node.customValues[ExternalFieldsPlus::EFP_Crack_Incre] < 0.0)
								node.customValues[ExternalFieldsPlus::EFP_Crack_Incre] = 0.0;
					}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						double old_crack = node.customValues[ExternalFieldsPlus::EFP_Crack];
						node.customValues[ExternalFieldsPlus::EFP_Crack] += node.customValues[ExternalFieldsPlus::EFP_Crack_Incre] * crack_dt;
						normalize_crack(node.customValues[ExternalFieldsPlus::EFP_Crack]);
						double abs_vari = abs(node.customValues[ExternalFieldsPlus::EFP_Crack] - old_crack);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (abs_vari > MAX_VARIATION)
								MAX_VARIATION = abs_vari;
						}
					}
			return MAX_VARIATION;
		}

		static double crack_loop_multiple(FieldStorage_forPhaseNode& phaseMesh) {
			return 0.0;
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Crack.model = 1 - Single \n", InputFileReader::get_instance()->debug_file);
			int input_crack_model = 0;
			dfint_dcrack = crack_propagation_models::df_dcrack_default;
			dfmech_dcrack = crack_propagation_models::df_dcrack_default;
			if (InputFileReader::get_instance()->read_int_value("Postprocess.Crack.model", input_crack_model, infile_debug)) {
				crack_model = CrackPropagationModel(input_crack_model);
				diff_method = Solvers::get_instance()->parameters.Difference_Method;
				crack_dr = phaseMesh.dr;
				crack_noise_range.resize(2);
				crack_noise_range[0] = 0.0;
				crack_noise_range[1] = 1.0;
				if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
					dfint_dcrack = crack_propagation_models::dfint_dcrack_single;
					if (plastic_solver::is_plasticity_on())
						dfmech_dcrack = crack_propagation_models::dfmech_dcrack_single;
					else
						dfmech_dcrack = crack_propagation_models::dfelas_dcrack_single;

					if (InputFileReader::get_instance()->read_string_value("Postprocess.Crack.Init.datafile_path", crack_init_datafile, infile_debug)) {
						is_crack_init_from_datafile = true;
						init_crack_field_from_datafile(phaseMesh);
					}
					else {
						is_crack_init_from_datafile = false;
						init_crack_field_default(phaseMesh);
					}

					Solvers::get_instance()->data_writer.register_custom_value(ExternalFieldsPlus::EFP_Crack, Data_Custom_Flag::DCW_Value);

					InputFileReader::get_instance()->read_int_value("Postprocess.Crack.Solver.frequency", solve_crack_steps, infile_debug);

					InputFileReader::get_instance()->read_int_value("Postprocess.Crack.Solver.max_iterate_steps", max_iterate_steps, infile_debug);

					InputFileReader::get_instance()->read_double_value("Postprocess.Crack.Solver.damage_threshold", damage_threshold, infile_debug);

					InputFileReader::get_instance()->read_double_value("Postprocess.Crack.Solver.crack_epsilon", crack_epsilon, infile_debug);

					InputFileReader::get_instance()->read_double_value("Postprocess.Crack.Solver.dt", crack_dt, infile_debug);

					InputFileReader::get_instance()->read_double_value("Postprocess.Crack.Solver.dr", crack_dr, infile_debug);

					InputFileReader::get_instance()->read_double_value("Postprocess.Crack.Solver.int_width", crack_int_width, infile_debug);

					InputFileReader::get_instance()->read_bool_value("Postprocess.Crack.VTS.debug", is_dfcrack_dphi_output, infile_debug);

					if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Crack.resistance = (Phase_0_G, Phase_1_G, ... ) \n", InputFileReader::get_instance()->debug_file);
					string resistance_key = "Postprocess.Crack.resistance", resistance_input = "()";
					if (InputFileReader::get_instance()->read_string_value(resistance_key, resistance_input, infile_debug)) {
						vector<input_value> resistance_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(pf::InputValueType::IVType_DOUBLE, resistance_key, resistance_input, infile_debug);
						if (Solvers::get_instance()->parameters.Phases.size() != resistance_value.size()) {
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Error, size of .Crack.resistance mismatch with Phases \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
						pGc.clear();
						for (int index = 0; index < resistance_value.size(); index++) {
							pGc.add_double(index, resistance_value[index].double_value);
						}
					}

					if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Crack.noise = 0 - uniform, 1 - interface \n", InputFileReader::get_instance()->debug_file);
					int input_crack_noise = 0;
					if (InputFileReader::get_instance()->read_int_value("Postprocess.Crack.noise", input_crack_noise, infile_debug)) {
						crack_noise_type = CrackNoiseType(input_crack_noise);

						InputFileReader::get_instance()->read_int_value("Postprocess.Crack.Noise.frequency", crack_noise_steps, infile_debug);

						InputFileReader::get_instance()->read_double_value("Postprocess.Crack.Noise.amplitude", crack_noise_amplitude, infile_debug);

						InputFileReader::get_instance()->read_bool_value("Postprocess.Crack.Noise.is_start", is_noise_start, infile_debug);

						InputFileReader::get_instance()->read_bool_value("Postprocess.Crack.Noise.on_crack", is_noise_on_crack, infile_debug);

						string noise_range_input = "(0.0,1.0)", noise_range_key = "Postprocess.Crack.Noise.random";
						if (InputFileReader::get_instance()->read_string_value(noise_range_key, noise_range_input, infile_debug)) {
							vector<input_value> noise_range_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, noise_range_key, noise_range_input, infile_debug);
							crack_noise_range[0] = noise_range_value[0].double_value;
							crack_noise_range[1] = noise_range_value[1].double_value;

						}
						if (crack_noise_range[1] < crack_noise_range[0])
							crack_noise_range[1] = crack_noise_range[0];
					}
				}
				else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {
					if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Error, crack model for Multiple_Order_Parameters hasn't been defined ! \n", InputFileReader::get_instance()->debug_file);
					exit(0);
				}
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Crack propagation !\n", LOG_FILE_NAME);
		}
		// is crack change smaller than threshold
		static int noise_gene_step = -10000000;
		static void crack_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_noise_start) {
				noise_crack(phaseMesh);
				noise_gene_step = Solvers::get_instance()->parameters.begin_step;
			}
		}
		// 
		static double crack_loop(FieldStorage_forPhaseNode& phaseMesh) {
			if (Solvers::get_instance()->current_istep != noise_gene_step
				&& Solvers::get_instance()->current_istep % crack_noise_steps == 0
				&& Solvers::get_instance()->current_istep != Solvers::get_instance()->parameters.end_step) {
				noise_crack(phaseMesh);
				noise_gene_step = Solvers::get_instance()->current_istep;
			}
			if (Solvers::get_instance()->current_istep % solve_crack_steps == 0) {
				if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
					diff_method = Solvers::get_instance()->parameters.Difference_Method;
					double VARIATION = crack_loop_single(phaseMesh);
					return VARIATION;
				}
				else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {
					return 0.0;
				}
			}
			return 0.0;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (crack_model == CrackPropagationModel::CPM_Single_Order_Parameter) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "crack" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							fout << node.customValues[ExternalFieldsPlus::EFP_Crack] << endl;
						}
				fout << "</DataArray>" << endl;
				if (is_dfcrack_dphi_output) {
					fout << "<DataArray type = \"Float64\" Name = \"" << "dfint_dcrack" <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++)
								fout << -dfint_dcrack(phaseMesh(i, j, k)) << endl;
					fout << "</DataArray>" << endl;
					fout << "<DataArray type = \"Float64\" Name = \"" << "dfmech_dcrack" <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++)
								fout << -dfmech_dcrack(phaseMesh(i, j, k)) << endl;
					fout << "</DataArray>" << endl;
					fout << "<DataArray type = \"Float64\" Name = \"" << "crack_resistance" <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++)
								fout << Gc(phaseMesh(i, j, k)) << endl;
					fout << "</DataArray>" << endl;
				}
			}
			else if (crack_model == CrackPropagationModel::CPM_Multiple_Order_Parameter) {

			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}