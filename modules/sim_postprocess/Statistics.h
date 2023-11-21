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
#include "MechanicalField.h"
#include "../sim_models/Source/Reaction/ElectrodeReaction.h"
#include "../sim_models/BulkEnergy.h"

namespace pf {
	namespace statistics {
		//-
		static bool is_mechanical_field_on = false;
		//-
		static bool is_phi_c_t_statistics = false;
		static bool is_mechanical_statistics = false;
		static bool is_electirc_statistics = false;
		static ConEquationType _con_type = ConEquationType::CEType_Const;
		static ConEquationDomain _domain = ConEquationDomain::CEDomain_Standard;
		static string file_name = "data_statistics";
		static string statistics_separator = "     ";
		static void init(FieldStorage_forPhaseNode& phaseMesh, bool is_mechanicas) {
			bool infile_debug = false;
			is_mechanical_field_on = is_mechanicas;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			InputFileReader::get_instance()->read_string_value("Postprocess.Statistics.file_name", file_name, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.Statistics.is_phi_c_t", is_phi_c_t_statistics, infile_debug);
			if (is_mechanicas)
				InputFileReader::get_instance()->read_bool_value("Postprocess.Statistics.is_mechanics", is_mechanical_statistics, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.Statistics.is_electricity", is_electirc_statistics, infile_debug);

			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.Statistics.datafiles = (datafile1, ... ) \n", InputFileReader::get_instance()->debug_file);

			_con_type = Solvers::get_instance()->parameters.ConEType;
			_domain = Solvers::get_instance()->parameters.ConEDomain;

			Solvers::get_instance()->writer.init_txt_file(file_name);
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Statistics !\n", LOG_FILE_NAME);
		}

		static void statistics_info(FieldStorage_forPhaseNode& phaseMesh, Vector6& ave_stress, Vector6& ave_strain, double& ave_plas_strain, double& int_potential, double& int_sum_phi, double& active_phi) {
			int active_component = electrode_reaction::get_active_component(), electrolyte_phase_index = electrode_reaction::get_electrolyte_phase_index();
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (is_phi_c_t_statistics) {
								;
							}
							if (is_mechanical_statistics) {
								ave_stress += node.customVec6s[ExternalFields::MECH_stress];
								ave_strain += node.customVec6s[ExternalFields::MECH_strain];
								if (plastic_solver::is_plasticity_on())
									ave_plas_strain += node.customValues[ExternalFields::MECH_ave_plastic_strain];
							}
							if (is_electirc_statistics) {
								if (_con_type == ConEquationType::CEType_GrandP || _con_type == ConEquationType::CEType_TotalX) {
									if (_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] > Solvers::get_instance()->C_Solver.threshold/* && node.customValues[ExternalFields::CON_Smooth_Phi] < 1.0 - SYS_EPSILON*/) {
										double local_potential = 0.0, local_phis = 0.0;
										for (auto phase = node.begin(); phase < node.end(); phase++)
											if (phase->phi > SYS_EPSILON && phase->index != electrolyte_phase_index) {
												double phase_potential = bulk_energy::bulk_energy_density(node, *phase);
												for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
													if (comp->index == active_component) {
														phase_potential += node.potential[comp->index].value * (1.0 - comp->value);
													}
													else {
														phase_potential -= node.potential[comp->index].value * comp->value;
													}
												}
												local_potential += phase_potential * phase->phi;
												local_phis += phase->phi;
											}
										if (local_phis > SYS_EPSILON) {
											int_potential += local_potential;
											int_sum_phi += local_phis;
										}
									}
									else if (_domain == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] < Solvers::get_instance()->C_Solver.threshold/* && node.customValues[ExternalFields::CON_Smooth_Phi] > SYS_EPSILON*/) {
										double local_potential = 0.0, local_phis = 0.0;
										for (auto phase = node.begin(); phase < node.end(); phase++)
											if (phase->phi > SYS_EPSILON && phase->index != electrolyte_phase_index) {
												double phase_potential = bulk_energy::bulk_energy_density(node, *phase);
												for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
													if (comp->index == active_component) {
														phase_potential += node.potential[comp->index].value * (1.0 - comp->value);
													}
													else {
														phase_potential -= node.potential[comp->index].value * comp->value;
													}
												}
												local_potential += phase_potential * phase->phi;
												local_phis += phase->phi;
											}
										if (local_phis > SYS_EPSILON) {
											int_potential += local_potential;
											int_sum_phi += local_phis;
										}
									}
								}
							}
							if (_con_type == ConEquationType::CEType_GrandP || _con_type == ConEquationType::CEType_TotalX) {
								if (_domain == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] > Solvers::get_instance()->C_Solver.threshold) {
									active_phi += node.customValues[ExternalFields::CON_Smooth_Phi];
								}
								else if (_domain == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] < Solvers::get_instance()->C_Solver.threshold) {
									active_phi += 1.0 - node.customValues[ExternalFields::CON_Smooth_Phi];
								}
							}
						}
					}
		}

		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			//- other should be exe

			// - 
			PhaseNode& info_node = Solvers::get_instance()->statistics_information_in_phaseMesh();
			Info_Phases& phases = Solvers::get_instance()->parameters.Phases;
			Info_Node& components = Solvers::get_instance()->parameters.Components;
			stringstream str;
			str << "istep" << statistics_separator << "real_time" << statistics_separator << "dtime";
			if (is_phi_c_t_statistics) {
				for (auto phi = info_node.begin(); phi < info_node.end(); phi++)
					str << statistics_separator << "phi_" << to_string(phi->index) << "_" << phases[phi->property].phi_name; // for each phi
				for (auto con = info_node.x.begin(); con < info_node.x.end(); con++)
					str << statistics_separator << "con_" << to_string(con->index) << "_" << components[con->index].name; // for each con
				for (auto pot = info_node.potential.begin(); pot < info_node.potential.end(); pot++)
					str << statistics_separator << "pot_" << to_string(pot->index) << "_" << components[pot->index].name; // for each potential
				if (_con_type == ConEquationType::CEType_GrandP || _con_type == ConEquationType::CEType_TotalX)
					str << statistics_separator << "active_phi";
			}
			if (is_mechanical_statistics) {
				str << statistics_separator << "app_strain_x" << statistics_separator << "app_strain_y" << statistics_separator << "app_strain_z";
				str << statistics_separator << "app_stress_x" << statistics_separator << "app_stress_y" << statistics_separator << "app_stress_z";
				str << statistics_separator << "ave_strain_x" << statistics_separator << "ave_strain_y" << statistics_separator << "ave_strain_z";
				str << statistics_separator << "ave_stress_x" << statistics_separator << "ave_stress_y" << statistics_separator << "ave_stress_z";
				str << statistics_separator << "ave_plas_strain";
			}
			if (is_electirc_statistics) {
				str << statistics_separator << "electrode_int_potential";
			}
			str << endl;
			double active_phi = 0.0, int_sum_phi = 0.0, int_potential = 0.0;
			Solvers::get_instance()->writer.add_string_to_txt(str.str(), file_name);
			Vector6 ave_stress, ave_strain; double ave_plas_strain = 0.0;
			ave_stress.set_to_zero(); ave_strain.set_to_zero();

			statistics_info(phaseMesh, ave_stress, ave_strain, ave_plas_strain, int_potential, int_sum_phi, active_phi);

			stringstream info_values;
			info_values << Solvers::get_instance()->current_istep << statistics_separator << Solvers::get_instance()->real_time << statistics_separator << Solvers::get_instance()->parameters.dt;
			if (is_phi_c_t_statistics) {
				for (auto phi = info_node.begin(); phi < info_node.end(); phi++)
					info_values << statistics_separator << phi->phi; // for each phi
				for (auto con = info_node.x.begin(); con < info_node.x.end(); con++)
					info_values << statistics_separator << con->value; // for each con
				for (auto pot = info_node.potential.begin(); pot < info_node.potential.end(); pot++)
					info_values << statistics_separator << pot->value; // for each potential
				if (_con_type == ConEquationType::CEType_GrandP || _con_type == ConEquationType::CEType_TotalX)
					info_values << statistics_separator << active_phi / phaseMesh.limit_x / phaseMesh.limit_y / phaseMesh.limit_z;
			}
			if (is_mechanical_statistics) {
				ave_stress /= phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z;
				ave_strain /= phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z;
				ave_plas_strain /= phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z;
				vStrain app_strain = elastic_solver::get_applied_strain();
				vStress app_stress = elastic_solver::get_applied_stress();
				info_values << statistics_separator << app_strain[0];
				info_values << statistics_separator << app_strain[1];
				info_values << statistics_separator << app_strain[2];
				info_values << statistics_separator << app_stress[0];
				info_values << statistics_separator << app_stress[1];
				info_values << statistics_separator << app_stress[2];
				info_values << statistics_separator << ave_strain[0];
				info_values << statistics_separator << ave_strain[1];
				info_values << statistics_separator << ave_strain[2];
				info_values << statistics_separator << ave_stress[0];
				info_values << statistics_separator << ave_stress[1];
				info_values << statistics_separator << ave_stress[2];
				info_values << statistics_separator << ave_plas_strain;
			}
			if (is_electirc_statistics) {
				if (int_sum_phi > SYS_EPSILON) {
					int_potential /= int_sum_phi;
				}
				else {
					int_potential = 0.0;
				}
				info_values << statistics_separator << int_potential;
			}
			info_values << endl;
			Solvers::get_instance()->writer.add_string_to_txt(info_values.str(), file_name);

		}

		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			//- other should be exe

			// - 
			int current_step = Solvers::get_instance()->current_istep;
			if (current_step % Solvers::get_instance()->parameters.screen_output_step != 0)
				return "";
			double active_phi = 0.0, int_sum_phi = 0.0, int_potential = 0.0;
			PhaseNode& info_node = Solvers::get_instance()->statistics_information_in_phaseMesh();
			Vector6 ave_stress, ave_strain; double ave_plas_strain = 0.0;
			ave_stress.set_to_zero(); ave_strain.set_to_zero();

			statistics_info(phaseMesh, ave_stress, ave_strain, ave_plas_strain, int_potential, int_sum_phi, active_phi);

			stringstream info_values;
			info_values << current_step << statistics_separator << Solvers::get_instance()->real_time << statistics_separator << Solvers::get_instance()->parameters.dt;
			if (is_phi_c_t_statistics) {
				for (auto phi = info_node.begin(); phi < info_node.end(); phi++)
					info_values << statistics_separator << phi->phi; // for each phi
				for (auto con = info_node.x.begin(); con < info_node.x.end(); con++)
					info_values << statistics_separator << con->value; // for each con
				for (auto pot = info_node.potential.begin(); pot < info_node.potential.end(); pot++)
					info_values << statistics_separator << pot->value; // for each potential
				if (_con_type == ConEquationType::CEType_GrandP || _con_type == ConEquationType::CEType_TotalX)
					info_values << statistics_separator << active_phi / phaseMesh.limit_x / phaseMesh.limit_y / phaseMesh.limit_z;
			}
			if (is_mechanical_statistics) {
				ave_stress /= phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z;
				ave_strain /= phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z;
				ave_plas_strain /= phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z;
				vStrain app_strain = elastic_solver::get_applied_strain();
				vStress app_stress = elastic_solver::get_applied_stress();
				info_values << statistics_separator << app_strain[0];
				info_values << statistics_separator << app_strain[1];
				info_values << statistics_separator << app_strain[2];
				info_values << statistics_separator << app_stress[0];
				info_values << statistics_separator << app_stress[1];
				info_values << statistics_separator << app_stress[2];
				info_values << statistics_separator << ave_strain[0];
				info_values << statistics_separator << ave_strain[1];
				info_values << statistics_separator << ave_strain[2];
				info_values << statistics_separator << ave_stress[0];
				info_values << statistics_separator << ave_stress[1];
				info_values << statistics_separator << ave_stress[2];
				info_values << statistics_separator << ave_plas_strain;
			}
			if (is_electirc_statistics) {
				if (int_sum_phi > SYS_EPSILON) {
					int_potential /= int_sum_phi;
				}
				else {
					int_potential = 0.0;
				}
				info_values << statistics_separator << int_potential;
			}
			info_values << endl;
			Solvers::get_instance()->writer.add_string_to_txt(info_values.str(), file_name);
			return "";
		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}