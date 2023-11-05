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
	namespace plastic_solver {
		// Plastic Parameters
		static int plastic_max_iterate_times = 0;
		static tensor1_double phi_yield_stress;
		static tensor1_double phi_hardening_modulus;
		static tensor1_double phi_shear_modulus;
		static bool is_plastic_on = false;
		static bool is_plasticity_on() {
			return is_plastic_on;
		}
		static double get_phase_hardening_modulus(PhaseEntry& phase) {
			return phi_hardening_modulus(phase.property);
		}
		static double get_hardening_modulus(PhaseNode& node) {
			double H = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++)
				H += phase->phi * phi_hardening_modulus(phase->property);
			return H;
		}

		static void init_plastic_field(FieldStorage_forPhaseNode& phaseMesh) {
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						bool plastic_elem_exists = false;
						for (auto elem = node.customValues.begin(); elem < node.customValues.end(); elem++)
							if (elem->index == ExternalFields::MECH_ave_plastic_strain)
								plastic_elem_exists = true;
						if (!plastic_elem_exists)
							node.customValues.add_double(ExternalFields::MECH_ave_plastic_strain, 0.0);
						plastic_elem_exists = false;
						Vector6 vec; vec.set_to_zero();
						for (auto elem = node.customVec6s.begin(); elem < node.customVec6s.end(); elem++)
							if (elem->index == ExternalFields::MECH_plastic_strain)
								plastic_elem_exists = true;
						if (!plastic_elem_exists)
							node.customVec6s.add_vec(ExternalFields::MECH_plastic_strain, vec);
					}
		}

		static double Prandtl_Reuss_formula(double deviatoric_part_of_mises_stress_norm, double yield_stress
			, double hardening_modulus, double ave_plastic_strain) {
			// sqrt(2.0 / 3.0) = 0.8164965809
			return deviatoric_part_of_mises_stress_norm - 0.8164965809 * (yield_stress + hardening_modulus * ave_plastic_strain);
		}

		static void solve_plastic_flow(FieldStorage_forPhaseNode& phaseMesh, int& max_iteration_steps, double& max_dplastic_strain, double& max_dave_plastic_strain) {
			// sqrt(2.0 / 3.0) = 0.8164965809
			max_iteration_steps = 0;
			max_dplastic_strain = 0.0;
			max_dave_plastic_strain = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.customVec6s[ExternalFields::MECH_stress] = node.customMatrix6x6s[ExternalFields::MECH_stiffness] *
							(node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain] - node.customVec6s[ExternalFields::MECH_plastic_strain]);
						Vector6 deviatoric_stress = node.customVec6s[ExternalFields::MECH_stress];
						double mean_stress = (deviatoric_stress[0] + deviatoric_stress[1] + deviatoric_stress[2]) / 3.0, yield_stress = 0.0, hardening_modulus = 0.0, shear_modulus = 0.0;
						deviatoric_stress[0] -= mean_stress;
						deviatoric_stress[1] -= mean_stress;
						deviatoric_stress[2] -= mean_stress;
						double deviatoric_stress_norm = deviatoric_stress.norm();
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							yield_stress += phase->phi * phi_yield_stress(phase->property);
							hardening_modulus += phase->phi * phi_hardening_modulus(phase->property);
							shear_modulus += phase->phi * phi_shear_modulus(phase->property);
						}
						double plastic_judgement = Prandtl_Reuss_formula(deviatoric_stress_norm, yield_stress, hardening_modulus, node.customValues[ExternalFields::MECH_ave_plastic_strain]);
						int iterate_steps = 0;
						if (plastic_max_iterate_times < 1)
							plastic_judgement = -1.0;
						Vector6 old_plastic_strain(node.customVec6s[ExternalFields::MECH_plastic_strain]);
						double old_ave_plastic_strain = node.customValues[ExternalFields::MECH_ave_plastic_strain];
						while (plastic_judgement > 0.0) {
							iterate_steps++;
							double magnitude_plastic_flow = plastic_judgement / (2.0 * shear_modulus + 2.0 / 3.0 * hardening_modulus);
							node.customVec6s[ExternalFields::MECH_plastic_strain] += deviatoric_stress / deviatoric_stress_norm * magnitude_plastic_flow;
							node.customValues[ExternalFields::MECH_ave_plastic_strain] += 0.8164965809 * magnitude_plastic_flow;
							node.customVec6s[ExternalFields::MECH_stress] = node.customMatrix6x6s[ExternalFields::MECH_stiffness] *
								(node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain] - node.customVec6s[ExternalFields::MECH_plastic_strain]);
							if (iterate_steps >= plastic_max_iterate_times) {
								plastic_judgement = -1.0;
							}
							else {
								deviatoric_stress = node.customVec6s[ExternalFields::MECH_stress];
								mean_stress = (deviatoric_stress[0] + deviatoric_stress[1] + deviatoric_stress[2]) / 3.0;
								deviatoric_stress[0] -= mean_stress;
								deviatoric_stress[1] -= mean_stress;
								deviatoric_stress[2] -= mean_stress;
								deviatoric_stress_norm = deviatoric_stress.norm();
								plastic_judgement = Prandtl_Reuss_formula(deviatoric_stress_norm, yield_stress, hardening_modulus, node.customValues[ExternalFields::MECH_ave_plastic_strain]);
							}
						}
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							if (iterate_steps > max_iteration_steps)
								max_iteration_steps = iterate_steps;
							for (int i = 0; i < 6; i++) {
								double delt_strain = abs(old_plastic_strain[i] - node.customVec6s[ExternalFields::MECH_plastic_strain][i]);
								if (delt_strain > max_dplastic_strain)
									max_dplastic_strain = delt_strain;
							}
							double delt_ave_strain = abs(old_ave_plastic_strain - node.customValues[ExternalFields::MECH_ave_plastic_strain]);
							if (delt_ave_strain > max_dave_plastic_strain)
								max_dave_plastic_strain = delt_ave_strain;
						}
					}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.plasticity", is_plastic_on, false);
			if (is_plastic_on) {
				InputFileReader::get_instance()->read_int_value("Postprocess.SolidMechanics.Plasticity.max_iteration_steps", plastic_max_iterate_times, infile_debug);
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.Plasticity.yield_stress = ( phase_0_stress, ... ) \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                      .hardening_modulus = (phase_0_modulus, ...) \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                      .shear_modulus = (phase_0_modulus, ...) \n", InputFileReader::get_instance()->debug_file);
				}
				string yield_stress_key = "Postprocess.SolidMechanics.Plasticity.yield_stress", yield_stress_input = "()";
				InputFileReader::get_instance()->read_string_value(yield_stress_key, yield_stress_input, infile_debug);
				vector<input_value> yield_stress_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, yield_stress_key, yield_stress_input, infile_debug);
				string hardening_modulus_key = "Postprocess.SolidMechanics.Plasticity.hardening_modulus", hardening_modulus_input = "()";
				InputFileReader::get_instance()->read_string_value(hardening_modulus_key, hardening_modulus_input, infile_debug);
				vector<input_value> hardening_modulus_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, hardening_modulus_key, hardening_modulus_input, infile_debug);
				string shear_modulus_key = "Postprocess.SolidMechanics.Plasticity.shear_modulus", shear_modulus_input = "()";
				InputFileReader::get_instance()->read_string_value(shear_modulus_key, shear_modulus_input, infile_debug);
				vector<input_value> shear_modulus_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, shear_modulus_key, shear_modulus_input, infile_debug);
				int index = 0;
				for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
					phi_yield_stress.add_double(phi->phi_property, yield_stress_value[index].double_value);
					phi_hardening_modulus.add_double(phi->phi_property, hardening_modulus_value[index].double_value);
					phi_shear_modulus.add_double(phi->phi_property, shear_modulus_value[index].double_value);
					index++;
				}
				init_plastic_field(phaseMesh);
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Plastic Explicit \n", LOG_FILE_NAME);
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (is_plastic_on) {
				vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
				for (int ele_index = 0; ele_index < 6; ele_index++)
				{
					string compname = "\"plastic_strain_" + compNameV[ele_index] + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << compname <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < phaseMesh.limit_z; ++k)
						for (int j = 0; j < phaseMesh.limit_y; ++j)
							for (int i = 0; i < phaseMesh.limit_x; ++i)
							{
								fout << phaseMesh(i, j, k).customVec6s[ExternalFields::MECH_plastic_strain][ele_index] << endl;
							}
					fout << "</DataArray>" << endl;
				}
				string compname = "\"ave_plastic_strain\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i)
						{
							fout << phaseMesh(i, j, k).customValues[ExternalFields::MECH_ave_plastic_strain] << endl;
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}