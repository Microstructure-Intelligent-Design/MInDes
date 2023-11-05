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
#include "../../sim_preprocess/MicroStructureInit.h"
#include "../Crack.h"
#include "StiffnessEigenStrain.h"
#include "PlasticSolver.h"

namespace pf {
	enum MechanicalFieldType { MFType_None, MFType_Explicit, MFType_Implicit_Steinbach, MFType_Implicit_Khachaturyan };
	enum FixBoundaryCondition { FixBC_Average, FixBC_Strain, FixBC_Stress };
	namespace elastic_solver {
		MechanicalField_Explicit mechanical_field_solver_ex;
		MechanicalField_Implicit mechanical_field_solver_im;
		// General Parameters
		static int MFType = MechanicalFieldType::MFType_None;
		static double bc_incre_rate = 1.0;
		//-
		static vector<int> fix_domain_boundary;
		static int solver_max_iterate_times = 0;
		static double solver_strain_accuracy = 1e-3;
		static bool solver_debug = false;
		static vector<Vector3> fix_boundary_x_change_rate;
		static vector<Vector3> fix_boundary_y_change_rate;
		static vector<Vector3> fix_boundary_z_change_rate;
		// MFType_Explicit Parameters
		static double mass_density = 1.0;
		static double mechanical_dt = 0.0;
		// MFType_Implicit Parameters
		static bool output_displacement_field = false;
		static double virtual_strain_iterate_rate = 1.0;
		static bool restart_iterator_in_loop = false;
		// plastic solver
		static bool is_plastic_on = false;
		static int mechanic_map_steps = 1;
		// boundary condition
		static void change_fix_boundaty_condition_implicity() {
			if (mechanical_field_solver_im.AppStrainMask[0]) {
				for (auto vec = fix_boundary_x_change_rate.begin(); vec < fix_boundary_x_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_im.applied_strain[0] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			else if (mechanical_field_solver_im.LoadStressMask[0]) {
				for (auto vec = fix_boundary_x_change_rate.begin(); vec < fix_boundary_x_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_im.applied_stress[0] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			if (mechanical_field_solver_im.AppStrainMask[1]) {
				for (auto vec = fix_boundary_y_change_rate.begin(); vec < fix_boundary_y_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_im.applied_strain[1] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			else if (mechanical_field_solver_im.LoadStressMask[1]) {
				for (auto vec = fix_boundary_y_change_rate.begin(); vec < fix_boundary_y_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_im.applied_stress[1] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			if (mechanical_field_solver_im.AppStrainMask[2]) {
				for (auto vec = fix_boundary_z_change_rate.begin(); vec < fix_boundary_z_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_im.applied_strain[2] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			else if (mechanical_field_solver_im.LoadStressMask[2]) {
				for (auto vec = fix_boundary_z_change_rate.begin(); vec < fix_boundary_z_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_im.applied_stress[2] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
		}
		static void change_fix_boundaty_condition_explicity() {
			if (mechanical_field_solver_ex.AppStrainMask[0]) {
				for (auto vec = fix_boundary_x_change_rate.begin(); vec < fix_boundary_x_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_ex.applied_strain[0] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			else if (mechanical_field_solver_ex.LoadStressMask[0]) {
				for (auto vec = fix_boundary_x_change_rate.begin(); vec < fix_boundary_x_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_ex.applied_stress[0] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			if (mechanical_field_solver_ex.AppStrainMask[1]) {
				for (auto vec = fix_boundary_y_change_rate.begin(); vec < fix_boundary_y_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_ex.applied_strain[1] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			else if (mechanical_field_solver_ex.LoadStressMask[1]) {
				for (auto vec = fix_boundary_y_change_rate.begin(); vec < fix_boundary_y_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_ex.applied_stress[1] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			if (mechanical_field_solver_ex.AppStrainMask[2]) {
				for (auto vec = fix_boundary_z_change_rate.begin(); vec < fix_boundary_z_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_ex.applied_strain[2] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
			else if (mechanical_field_solver_ex.LoadStressMask[2]) {
				for (auto vec = fix_boundary_z_change_rate.begin(); vec < fix_boundary_z_change_rate.end(); vec++)
					if (Solvers::get_instance()->real_time >= (*vec)[0] && Solvers::get_instance()->real_time < (*vec)[1]) {
						mechanical_field_solver_ex.applied_stress[2] += (*vec)[2] * Solvers::get_instance()->parameters.dt;
						break;
					}
			}
		}
		// get infomation
		static vStrain get_applied_strain() {
			if (MFType == MechanicalFieldType::MFType_Explicit)
				return mechanical_field_solver_ex.applied_strain;
			else if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach || MFType == MechanicalFieldType::MFType_Implicit_Khachaturyan)
				return mechanical_field_solver_im.applied_strain;
			else
				return vStrain();
		}
		static vStress get_applied_stress() {
			if (MFType == MechanicalFieldType::MFType_Explicit)
				return mechanical_field_solver_ex.applied_stress;
			else if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach || MFType == MechanicalFieldType::MFType_Implicit_Khachaturyan)
				return mechanical_field_solver_im.applied_stress;
			else
				return vStress();
		}
		static MechanicalFieldType _MFType() {
			return MechanicalFieldType(MFType);
		}
		// eigenstrain and stiffness
		static void cal_parameters_for_main_domain_phi_dependent(PhaseNode& node, int Nx, int Ny, int Nz) {
			node.customMatrix6x6s[ExternalFields::MECH_stiffness].set_to_zero();
			node.customVec6s[ExternalFields::MECH_eigen_strain].set_to_zero();
			stiffness_eigenstrain::stiffness(node, node.customMatrix6x6s[ExternalFields::MECH_stiffness]);
			for (auto func = stiffness_eigenstrain::eigenstrain_list.begin(); func < stiffness_eigenstrain::eigenstrain_list.end(); func++)
				(*func)(node, node.customVec6s[ExternalFields::MECH_eigen_strain]);
			return;
		}
		static void cal_parameters_for_main_domain_phi_dependent_single_crack(PhaseNode& node, int Nx, int Ny, int Nz) {
			node.customMatrix6x6s[ExternalFields::MECH_stiffness].set_to_zero();
			node.customVec6s[ExternalFields::MECH_eigen_strain].set_to_zero();
			stiffness_eigenstrain::cal_stiffness_phi_dependent(node, node.customMatrix6x6s[ExternalFields::MECH_stiffness]);
			stiffness_eigenstrain::cal_eigenstrain_phi_dependent(node, node.customVec6s[ExternalFields::MECH_eigen_strain]);
			stiffness_eigenstrain::cal_eigenstrain_physics(node, node.customVec6s[ExternalFields::MECH_eigen_strain]);
			node.customMatrix6x6s[ExternalFields::MECH_stiffness] *= (1.0 - crack_propagation::crack_fraction_single(node));
			node.customVec6s[ExternalFields::MECH_eigen_strain] *= (1.0 - crack_propagation::crack_fraction_single(node));
			return;
		}
		static void cal_parameters_for_main_domain_phi_dependent_multiple_crack(PhaseNode& node, int Nx, int Ny, int Nz) {
			node.customMatrix6x6s[ExternalFields::MECH_stiffness].set_to_zero();
			node.customVec6s[ExternalFields::MECH_eigen_strain].set_to_zero();
			double crack_fraction = 0.0;
			stiffness_eigenstrain::cal_stiffness_phi_dependent(node, node.customMatrix6x6s[ExternalFields::MECH_stiffness]);
			stiffness_eigenstrain::cal_eigenstrain_phi_dependent(node, node.customVec6s[ExternalFields::MECH_eigen_strain]);
			stiffness_eigenstrain::cal_eigenstrain_physics(node, node.customVec6s[ExternalFields::MECH_eigen_strain]);
			node.customMatrix6x6s[ExternalFields::MECH_stiffness] *= (1.0 - crack_fraction);
			node.customVec6s[ExternalFields::MECH_eigen_strain] *= (1.0 - crack_fraction);
			return;
		}
		// solver loop for MFType_Explicit
		static void exec_pre_ex(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			if (solver_debug) {
				output << "> Mechanical field solver debug:" << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
			double MAX_ABS_STRAIN = 0.0, MAX_ABS_TARGET_STRAIN = 0.0;
			int plastic_iteration = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			vStress average_stress; vStrain average_strain; average_stress.set_to_zero(); average_strain.set_to_zero();
			mechanical_field_solver_ex.cal_parameters_before_calculation();
			mechanical_field_solver_ex.cal_stress(average_stress, is_plastic_on);
			mechanical_field_solver_ex.boundary_condition(average_stress, average_strain, bc_incre_rate);
			int ITERATE_STEPS = 0;
			for (int map_step = 1; map_step <= mechanic_map_steps; map_step++) {
				for (int istep = 1; istep <= solver_max_iterate_times; istep++) {
					ITERATE_STEPS++;
					MAX_ABS_STRAIN = mechanical_field_solver_ex.evolve_momentum_equation(mass_density, mechanical_dt, average_strain);
					mechanical_field_solver_ex.cal_stress(average_stress, is_plastic_on);
					MAX_ABS_TARGET_STRAIN = mechanical_field_solver_ex.boundary_condition(average_stress, average_strain, bc_incre_rate);
					if (solver_debug) {
						output.str("");
						output << "(Elastic Solver) iterate step:   " << istep << endl
							<< "                 dstrain:        " << MAX_ABS_STRAIN << endl
							<< "                 dtargetstrain:  " << MAX_ABS_TARGET_STRAIN << endl
							<< "                 Average Strain: " << "( " << average_strain[0] << ", "
							<< average_strain[1] << ", "
							<< average_strain[2] << ", "
							<< average_strain[3] << ", "
							<< average_strain[4] << ", "
							<< average_strain[5] << " )" << endl
							<< "                 Average Stress: " << "( " << average_stress[0] << ", "
							<< average_stress[1] << ", "
							<< average_stress[2] << ", "
							<< average_stress[3] << ", "
							<< average_stress[4] << ", "
							<< average_stress[5] << " )" << endl;
						Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
					}
					if (MAX_ABS_STRAIN < solver_strain_accuracy && MAX_ABS_TARGET_STRAIN < solver_strain_accuracy)
						break;
				}
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(phaseMesh, plastic_iteration, dplastic_strain, dave_plastic_strain);
					if (solver_debug) {
						output.str("");
						output << "(plastic Solver) iterate step:         " << plastic_iteration << endl
							<< "                 dplastic strain:     " << dplastic_strain << endl
							<< "                 dave plastic strain: " << dave_plastic_strain << endl;
						Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
					}
					if (dplastic_strain < solver_strain_accuracy && dave_plastic_strain < solver_strain_accuracy)
						break;
				}
				else {
					break;
				}
			}
			output.str("");
			output << "> Elastic Solver" << endl
				<< "  iterate step:		" << ITERATE_STEPS << endl
				<< "  MAX_ABS_dStrain:	" << MAX_ABS_STRAIN << endl
				<< "  MAX_ABS_dTarStrain: " << MAX_ABS_TARGET_STRAIN << endl
				<< "  Average Strain:     " << "( " << average_strain[0] << ", "
				<< average_strain[1] << ", "
				<< average_strain[2] << ", "
				<< average_strain[3] << ", "
				<< average_strain[4] << ", "
				<< average_strain[5] << " )" << endl
				<< "  Average Stress:	    " << "( " << average_stress[0] << ", "
				<< average_stress[1] << ", "
				<< average_stress[2] << ", "
				<< average_stress[3] << ", "
				<< average_stress[4] << ", "
				<< average_stress[5] << " )" << endl;
			if (is_plastic_on) {
				output << "  plastic steps:       " << plastic_iteration << endl
					<< "  dplastic strain:     " << dplastic_strain << endl
					<< "  dave plastic strain: " << dave_plastic_strain << endl;
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
		}
		static string exec_loop_ex(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			double MAX_ABS_STRAIN = 0.0, MAX_ABS_TARGET_STRAIN = 0.0;
			int plastic_iteration = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			vStress average_stress; vStrain average_strain; average_stress.set_to_zero(); average_strain.set_to_zero();
			change_fix_boundaty_condition_explicity();
			mechanical_field_solver_ex.cal_parameters_before_calculation();
			mechanical_field_solver_ex.cal_stress(average_stress, is_plastic_on);
			mechanical_field_solver_ex.boundary_condition(average_stress, average_strain, bc_incre_rate);
			int ITERATE_STEPS = 0;
			for (int map_step = 1; map_step <= mechanic_map_steps; map_step++) {
				for (int istep = 1; istep <= solver_max_iterate_times; istep++) {
					ITERATE_STEPS++;
					MAX_ABS_STRAIN = mechanical_field_solver_ex.evolve_momentum_equation(mass_density, mechanical_dt, average_strain);
					mechanical_field_solver_ex.cal_stress(average_stress, is_plastic_on);
					MAX_ABS_TARGET_STRAIN = mechanical_field_solver_ex.boundary_condition(average_stress, average_strain, bc_incre_rate);
					if (MAX_ABS_STRAIN < solver_strain_accuracy && MAX_ABS_TARGET_STRAIN < solver_strain_accuracy)
						break;
				}
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(phaseMesh, plastic_iteration, dplastic_strain, dave_plastic_strain);
					if (dplastic_strain < solver_strain_accuracy && dave_plastic_strain < solver_strain_accuracy)
						break;
				}
				else {
					break;
				}
			}
			mechanical_field_solver_ex.cal_stress(average_stress, is_plastic_on);
			output.str("");
			output << "> Elastic Solver" << endl
				<< "  iterate step:		" << ITERATE_STEPS << endl
				<< "  MAX_ABS_dStrain:	" << MAX_ABS_STRAIN << endl
				<< "  MAX_ABS_dTarStrain: " << MAX_ABS_TARGET_STRAIN << endl
				<< "  Average Strain:     " << "( " << average_strain[0] << ", "
				<< average_strain[1] << ", "
				<< average_strain[2] << ", "
				<< average_strain[3] << ", "
				<< average_strain[4] << ", "
				<< average_strain[5] << " )" << endl
				<< "  Average Stress:	    " << "( " << average_stress[0] << ", "
				<< average_stress[1] << ", "
				<< average_stress[2] << ", "
				<< average_stress[3] << ", "
				<< average_stress[4] << ", "
				<< average_stress[5] << " )" << endl;
			if (is_plastic_on) {
				output << "  plastic steps:       " << plastic_iteration << endl
					<< "  dplastic strain:     " << dplastic_strain << endl
					<< "  dave plastic strain: " << dave_plastic_strain << endl;
			}
			return output.str();
		}
		static void write_vec3_ex(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			string name;
			name = "\"mech_velocity\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << mechanical_field_solver_ex.get_u_main_node(i, j, k) << " "
							<< mechanical_field_solver_ex.get_v_main_node(i, j, k) << " "
							<< mechanical_field_solver_ex.get_w_main_node(i, j, k) << endl;
					}
			fout << "</DataArray>" << endl;
		}
		// solver loop for MFType_Implicit
		static void exec_pre_im_steinbach(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			if (solver_debug) {
				output << "> Mechanical field solver debug:" << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
			mechanical_field_solver_im.SetMAXElasticConstants(stiffness_eigenstrain::get_stiffness());
			mechanical_field_solver_im.cal_parameters_before_calculation(is_plastic_on);
			mechanical_field_solver_im.initStrainIncrements();
			int plastic_iteration = 0, map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				output.str("");
				output << "> Map: " << map_steps << endl;
				elastic_solver_output = mechanical_field_solver_im.Solve(solver_strain_accuracy, solver_max_iterate_times, Solvers::get_instance()->writer, bc_incre_rate, solver_debug, output_displacement_field);
				output << elastic_solver_output;
				if (is_plastic_on) {					plastic_solver::solve_plastic_flow(phaseMesh, plastic_iteration, dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain_with_plasticity();
					output  << "  plastic steps:       " << plastic_iteration << endl
						    << "  dplastic strain:     " << dplastic_strain << endl
						    << "  dave plastic strain: " << dave_plastic_strain << endl;
				}
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy || plastic_iteration > 1);
		}
		static string exec_loop_im_steinbach(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			change_fix_boundaty_condition_implicity();
			mechanical_field_solver_im.cal_parameters_before_calculation(is_plastic_on);
			if (restart_iterator_in_loop)
				mechanical_field_solver_im.initStrainIncrements();
			int plastic_iteration = 0, map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				elastic_solver_output = mechanical_field_solver_im.Solve(solver_strain_accuracy, solver_max_iterate_times, Solvers::get_instance()->writer, bc_incre_rate, false, output_displacement_field);
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(phaseMesh, plastic_iteration, dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain_with_plasticity();
				}
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy || plastic_iteration > 1);

			output << "  Mapping steps:       " << map_steps << endl
				<< elastic_solver_output
				<< "  plastic steps:       " << plastic_iteration << endl
				<< "  dplastic strain:     " << dplastic_strain << endl
				<< "  dave plastic strain: " << dave_plastic_strain << endl;
			return output.str();
		}
		static void exec_pre_im_khachaturyan(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			if (solver_debug) {
				output << "> Mechanical field solver debug:" << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
			mechanical_field_solver_im.SetMAXElasticConstants(stiffness_eigenstrain::get_stiffness());
			mechanical_field_solver_im.cal_parameters_before_calculation(is_plastic_on);
			mechanical_field_solver_im.initVirtualEigenstrain();
			int plastic_iteration = 0, map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				output.str("");
				output << "> Map: " << map_steps << endl;
				elastic_solver_output = mechanical_field_solver_im.Solve2(solver_strain_accuracy, solver_max_iterate_times, virtual_strain_iterate_rate, Solvers::get_instance()->writer, solver_debug, output_displacement_field);
				output << elastic_solver_output;
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(phaseMesh, plastic_iteration, dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain_with_plasticity();
					output << "  plastic steps:       " << plastic_iteration << endl
						<< "  dplastic strain:     " << dplastic_strain << endl
						<< "  dave plastic strain: " << dave_plastic_strain << endl;
				}
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy || plastic_iteration > 1);
		}
		static string exec_loop_im_khachaturyan(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			change_fix_boundaty_condition_implicity();
			mechanical_field_solver_im.cal_parameters_before_calculation(is_plastic_on);
			if (restart_iterator_in_loop)
				mechanical_field_solver_im.initVirtualEigenstrain();
			int plastic_iteration = 0, map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				elastic_solver_output = mechanical_field_solver_im.Solve2(solver_strain_accuracy, solver_max_iterate_times, virtual_strain_iterate_rate, Solvers::get_instance()->writer, false, output_displacement_field);
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(phaseMesh, plastic_iteration, dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain_with_plasticity();
				}
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy || plastic_iteration > 1);

			output << "  Mapping steps:       " << map_steps << endl
				<< elastic_solver_output
				<< "  plastic steps:       " << plastic_iteration << endl
				<< "  dplastic strain:     " << dplastic_strain << endl
				<< "  dave plastic strain: " << dave_plastic_strain << endl;
			return output.str();
		}
		static void write_vec3_im(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			string name;
			name = "\"mech_velocity\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << mechanical_field_solver_im.get_u_main_node(i, j, k) << " "
							 << mechanical_field_solver_im.get_v_main_node(i, j, k) << " "
							 << mechanical_field_solver_im.get_w_main_node(i, j, k) << endl;
					}
			fout << "</DataArray>" << endl;
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.momentum_balance = 0 - None , 1 - Explicit , 2 - Implicit (Ingo Steinbach) , 3 - Implicit (Armen G. Khachaturyan) \n", InputFileReader::get_instance()->debug_file);
			InputFileReader::get_instance()->read_int_value("Postprocess.SolidMechanics.momentum_balance", MFType, infile_debug);
			if (MechanicalFieldType(MFType) == MechanicalFieldType::MFType_Explicit) {
				mechanical_field_solver_ex.init(phaseMesh, phaseMesh._bc_x_up, phaseMesh._bc_x_down, phaseMesh._bc_y_up, phaseMesh._bc_y_down, phaseMesh._bc_z_up, phaseMesh._bc_z_down);
				mechanical_field_solver_ex.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent;
				fix_domain_boundary.resize(3);
				InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.Explicit.mass_density", mass_density, infile_debug);
				InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.Explicit.mechanical_dt", mechanical_dt, infile_debug);
				InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.Explicit.bc_ite_rate", bc_incre_rate, infile_debug);
				// (BC_X, BC_Y, BC_Z)
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.type = (BC_X, BC_Y, BC_Z) , 0 - Average , 1 - Strain , 2 - Stress \n", InputFileReader::get_instance()->debug_file);
				string fix_boundary_key = "Postprocess.SolidMechanics.fix_boundary.type", fix_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(fix_boundary_key, fix_boundary_input, infile_debug)) {
					vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, fix_boundary_key, fix_boundary_input, infile_debug);
					fix_domain_boundary[Axis::AXIS_X] = FixBoundaryCondition(fix_boundary_value[0].int_value);
					fix_domain_boundary[Axis::AXIS_Y] = FixBoundaryCondition(fix_boundary_value[1].int_value);
					fix_domain_boundary[Axis::AXIS_Z] = FixBoundaryCondition(fix_boundary_value[2].int_value);
					for (int direction = 0; direction < 3; direction++) {
						if (fix_domain_boundary[direction] == 1) {
							mechanical_field_solver_ex.AppStrainMask[direction] = true;
							mechanical_field_solver_ex.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (Axis(direction))
							{
							case pf::AXIS_X:
								_d = "x";
								break;
							case pf::AXIS_Y:
								_d = "y";
								break;
							case pf::AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d;
							double strain = 0.0;
							InputFileReader::get_instance()->read_double_value(fix_boundary_val_key, strain, infile_debug);
							mechanical_field_solver_ex.applied_strain[direction] = strain;
							// [(time_begin,time_end,change rate), ... ]
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.strain_? = [(real_time_begin, real_time_end, dstrain_dt), ... ] \n", InputFileReader::get_instance()->debug_file);
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (InputFileReader::get_instance()->read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, infile_debug)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_rate_key, fix_boundary_rate_input, infile_debug);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].double_value;
									vec[1] = fix_boundary_rate_value[index][1].double_value;
									vec[2] = fix_boundary_rate_value[index][2].double_value;
									switch (Axis(direction))
									{
									case pf::AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case pf::AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case pf::AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
						else if (fix_domain_boundary[direction] == 2) {
							mechanical_field_solver_im.LoadStressMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (Axis(direction))
							{
							case pf::AXIS_X:
								_d = "x";
								break;
							case pf::AXIS_Y:
								_d = "y";
								break;
							case pf::AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d;
							double stress = 0.0;
							InputFileReader::get_instance()->read_double_value(fix_boundary_val_key, stress, infile_debug);
							mechanical_field_solver_im.applied_stress[direction] = stress;
							// [(time_begin,time_end,change rate), ... ]
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.stress_? = [(real_time_begin, real_time_end, dstress_dt), ... ] \n", InputFileReader::get_instance()->debug_file);
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (InputFileReader::get_instance()->read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, infile_debug)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_rate_key, fix_boundary_rate_input, infile_debug);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].double_value;
									vec[1] = fix_boundary_rate_value[index][1].double_value;
									vec[2] = fix_boundary_rate_value[index][2].double_value;
									switch (Axis(direction))
									{
									case pf::AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case pf::AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case pf::AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
					}
				}
			}
			else if (MechanicalFieldType(MFType) == MechanicalFieldType::MFType_Implicit_Steinbach) {
				pf::BoundaryCondition bc_x = pf::BoundaryCondition::PERIODIC, bc_y = pf::BoundaryCondition::PERIODIC, bc_z = pf::BoundaryCondition::PERIODIC;
				if (Solvers::get_instance()->phaseMesh._bc_x_down != pf::BoundaryCondition::PERIODIC
					|| Solvers::get_instance()->phaseMesh._bc_x_up != pf::BoundaryCondition::PERIODIC)
					bc_x = pf::BoundaryCondition::ADIABATIC;
				if (Solvers::get_instance()->phaseMesh._bc_y_down != pf::BoundaryCondition::PERIODIC
					|| Solvers::get_instance()->phaseMesh._bc_y_up != pf::BoundaryCondition::PERIODIC)
					bc_y = pf::BoundaryCondition::ADIABATIC;
				if (Solvers::get_instance()->phaseMesh._bc_z_down != pf::BoundaryCondition::PERIODIC
					|| Solvers::get_instance()->phaseMesh._bc_z_up != pf::BoundaryCondition::PERIODIC)
					bc_z = pf::BoundaryCondition::ADIABATIC;
				mechanical_field_solver_im.init(Solvers::get_instance()->phaseMesh, bc_x, bc_y, bc_z);
				mechanical_field_solver_im.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent;
				InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.restart_iterator_in_loop", restart_iterator_in_loop, infile_debug);
				InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.Implicit.bc_ite_rate", bc_incre_rate, infile_debug);
				fix_domain_boundary.resize(3);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.type = (BC_X, BC_Y, BC_Z) , 0 - Average , 1 - Strain , 2 - Stress \n", InputFileReader::get_instance()->debug_file);
				string fix_boundary_key = "Postprocess.SolidMechanics.fix_boundary.type", fix_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(fix_boundary_key, fix_boundary_input, infile_debug)) {
					vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, fix_boundary_key, fix_boundary_input, infile_debug);
					fix_domain_boundary[Axis::AXIS_X] = FixBoundaryCondition(fix_boundary_value[0].int_value);
					fix_domain_boundary[Axis::AXIS_Y] = FixBoundaryCondition(fix_boundary_value[1].int_value);
					fix_domain_boundary[Axis::AXIS_Z] = FixBoundaryCondition(fix_boundary_value[2].int_value);
					for (int direction = 0; direction < 3; direction++) {
						if (fix_domain_boundary[direction] == 1) {
							mechanical_field_solver_im.AppStrainMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (Axis(direction))
							{
							case pf::AXIS_X:
								_d = "x";
								break;
							case pf::AXIS_Y:
								_d = "y";
								break;
							case pf::AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d;
							double strain = 0.0;
							InputFileReader::get_instance()->read_double_value(fix_boundary_val_key, strain, infile_debug);
							mechanical_field_solver_im.applied_strain[direction] = strain;
							// [(time_begin,time_end,change rate), ... ]
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.strain_? = [(real_time_begin, real_time_end, dstrain_dt), ... ] \n", InputFileReader::get_instance()->debug_file);
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (InputFileReader::get_instance()->read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, infile_debug)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_rate_key, fix_boundary_rate_input, infile_debug);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].double_value;
									vec[1] = fix_boundary_rate_value[index][1].double_value;
									vec[2] = fix_boundary_rate_value[index][2].double_value;
									switch (Axis(direction))
									{
									case pf::AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case pf::AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case pf::AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
						else if (fix_domain_boundary[direction] == 2) {
							mechanical_field_solver_im.LoadStressMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (Axis(direction))
							{
							case pf::AXIS_X:
								_d = "x";
								break;
							case pf::AXIS_Y:
								_d = "y";
								break;
							case pf::AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d;
							double stress = 0.0;
							InputFileReader::get_instance()->read_double_value(fix_boundary_val_key, stress, infile_debug);
							mechanical_field_solver_im.applied_stress[direction] = stress;
							// [(time_begin,time_end,change rate), ... ]
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.stress_? = [(real_time_begin, real_time_end, dstress_dt), ... ] \n", InputFileReader::get_instance()->debug_file);
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (InputFileReader::get_instance()->read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, infile_debug)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_rate_key, fix_boundary_rate_input, infile_debug);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].double_value;
									vec[1] = fix_boundary_rate_value[index][1].double_value;
									vec[2] = fix_boundary_rate_value[index][2].double_value;
									switch (Axis(direction))
									{
									case pf::AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case pf::AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case pf::AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
					}
				}
			}
			else if (MechanicalFieldType(MFType) == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				pf::BoundaryCondition bc_x = pf::BoundaryCondition::PERIODIC, bc_y = pf::BoundaryCondition::PERIODIC, bc_z = pf::BoundaryCondition::PERIODIC;
				if (Solvers::get_instance()->phaseMesh._bc_x_down != pf::BoundaryCondition::PERIODIC
					|| Solvers::get_instance()->phaseMesh._bc_x_up != pf::BoundaryCondition::PERIODIC)
					bc_x = pf::BoundaryCondition::ADIABATIC;
				if (Solvers::get_instance()->phaseMesh._bc_y_down != pf::BoundaryCondition::PERIODIC
					|| Solvers::get_instance()->phaseMesh._bc_y_up != pf::BoundaryCondition::PERIODIC)
					bc_y = pf::BoundaryCondition::ADIABATIC;
				if (Solvers::get_instance()->phaseMesh._bc_z_down != pf::BoundaryCondition::PERIODIC
					|| Solvers::get_instance()->phaseMesh._bc_z_up != pf::BoundaryCondition::PERIODIC)
					bc_z = pf::BoundaryCondition::ADIABATIC;
				mechanical_field_solver_im.init(Solvers::get_instance()->phaseMesh, bc_x, bc_y, bc_z);
				mechanical_field_solver_im.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent;
				InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.VisualStrain.L_ijkl", virtual_strain_iterate_rate, infile_debug);
				InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.restart_iterator_in_loop", restart_iterator_in_loop, infile_debug);
				fix_domain_boundary.resize(3);
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.type = (BC_X, BC_Y, BC_Z) , 0 - Average , 1 - Strain , 2 - Stress \n", InputFileReader::get_instance()->debug_file);
				string fix_boundary_key = "Postprocess.SolidMechanics.fix_boundary.type", fix_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(fix_boundary_key, fix_boundary_input, infile_debug)) {
					vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, fix_boundary_key, fix_boundary_input, infile_debug);
					fix_domain_boundary[Axis::AXIS_X] = FixBoundaryCondition(fix_boundary_value[0].int_value);
					fix_domain_boundary[Axis::AXIS_Y] = FixBoundaryCondition(fix_boundary_value[1].int_value);
					fix_domain_boundary[Axis::AXIS_Z] = FixBoundaryCondition(fix_boundary_value[2].int_value);
					for (int direction = 0; direction < 3; direction++) {
						if (fix_domain_boundary[direction] == 1) {
							mechanical_field_solver_im.AppStrainMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (Axis(direction))
							{
							case pf::AXIS_X:
								_d = "x";
								break;
							case pf::AXIS_Y:
								_d = "y";
								break;
							case pf::AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d;
							double strain = 0.0;
							InputFileReader::get_instance()->read_double_value(fix_boundary_val_key, strain, infile_debug);
							mechanical_field_solver_im.applied_strain[direction] = strain;
							// [(time_begin,time_end,change rate), ... ]
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.strain_? = [(real_time_begin, real_time_end, dstrain_dt), ... ] \n", InputFileReader::get_instance()->debug_file);
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (InputFileReader::get_instance()->read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, infile_debug)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_rate_key, fix_boundary_rate_input, infile_debug);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].double_value;
									vec[1] = fix_boundary_rate_value[index][1].double_value;
									vec[2] = fix_boundary_rate_value[index][2].double_value;
									switch (Axis(direction))
									{
									case pf::AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case pf::AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case pf::AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
						else if (fix_domain_boundary[direction] == 2) {
							mechanical_field_solver_im.LoadStressMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (Axis(direction))
							{
							case pf::AXIS_X:
								_d = "x";
								break;
							case pf::AXIS_Y:
								_d = "y";
								break;
							case pf::AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d;
							double stress = 0.0;
							InputFileReader::get_instance()->read_double_value(fix_boundary_val_key, stress, infile_debug);
							mechanical_field_solver_im.applied_stress[direction] = stress;
							// [(time_begin,time_end,change rate), ... ]
							if (infile_debug)
								InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.SolidMechanics.fix_boundary.stress_? = [(real_time_begin, real_time_end, dstress_dt), ... ] \n", InputFileReader::get_instance()->debug_file);
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (InputFileReader::get_instance()->read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, infile_debug)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_rate_key, fix_boundary_rate_input, infile_debug);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].double_value;
									vec[1] = fix_boundary_rate_value[index][1].double_value;
									vec[2] = fix_boundary_rate_value[index][2].double_value;
									switch (Axis(direction))
									{
									case pf::AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case pf::AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case pf::AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
					}
				}
			}
			InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.write_displacement_field", output_displacement_field, infile_debug);
			InputFileReader::get_instance()->read_int_value("Postprocess.SolidMechanics.max_iteration_steps", solver_max_iterate_times, infile_debug);
			InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.debug", solver_debug, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.SolidMechanics.strain_accuracy", solver_strain_accuracy, infile_debug);
			string solid_phi_key = "Postprocess.SolidMechanics.solid_phases", solid_phi_input = "()";
			InputFileReader::get_instance()->read_string_value(solid_phi_key, solid_phi_input, infile_debug);
			vector<input_value> solid_phi_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, solid_phi_key, solid_phi_input, infile_debug);
			vector<int> solid_phases;
			for (auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++) {
				bool is_solid = false;
				for (auto in_phi = solid_phi_value.begin(); in_phi < solid_phi_value.end(); in_phi++)
					if (phi->phi_name.compare((*in_phi).string_value) == 0)
						is_solid = true;
				if (is_solid) {
					solid_phases.push_back(phi->phi_property);
				}
			}

			InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.plasticity", is_plastic_on, infile_debug);
			if (is_plastic_on) {
				InputFileReader::get_instance()->read_int_value("Postprocess.SolidMechanics.Elastoplasticity.mapping_steps", mechanic_map_steps, infile_debug);
				if (mechanic_map_steps < 1)
					mechanic_map_steps = 1;
			}

			stiffness_eigenstrain::init(phaseMesh, solid_phases);

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Elastic Solver \n", LOG_FILE_NAME);
		}

		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (_MFType() == MechanicalFieldType::MFType_Explicit) {
				if (crack_propagation::crack_propagation_model() == crack_propagation::CPM_Single_Order_Parameter)
					mechanical_field_solver_ex.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent_single_crack;
				else if (crack_propagation::crack_propagation_model() == crack_propagation::CPM_Multiple_Order_Parameter)
					mechanical_field_solver_ex.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent_multiple_crack;
				exec_pre_ex(phaseMesh);
			}
			else if (_MFType() == MechanicalFieldType::MFType_Implicit_Steinbach) {
				if (crack_propagation::crack_propagation_model() == crack_propagation::CPM_Single_Order_Parameter)
					mechanical_field_solver_im.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent_single_crack;
				else if (crack_propagation::crack_propagation_model() == crack_propagation::CPM_Multiple_Order_Parameter)
					mechanical_field_solver_im.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent_multiple_crack;
				exec_pre_im_steinbach(phaseMesh);
			}
			else if (_MFType() == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				if (crack_propagation::crack_propagation_model() == crack_propagation::CPM_Single_Order_Parameter)
					mechanical_field_solver_im.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent_single_crack;
				else if (crack_propagation::crack_propagation_model() == crack_propagation::CPM_Multiple_Order_Parameter)
					mechanical_field_solver_im.cal_parameters_for_main_domain = cal_parameters_for_main_domain_phi_dependent_multiple_crack;
				exec_pre_im_khachaturyan(phaseMesh);
			}
		}

		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			if (_MFType() == MechanicalFieldType::MFType_Explicit) {
				return exec_loop_ex(phaseMesh);
			}
			else if (_MFType() == MechanicalFieldType::MFType_Implicit_Steinbach) {
				return exec_loop_im_steinbach(phaseMesh);
			}
			else if (_MFType() == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				return exec_loop_im_khachaturyan(phaseMesh);
			}
			else {
				return "";
			}
		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			stiffness_eigenstrain::deinit(phaseMesh);
			mechanical_field_solver_ex.free();
			mechanical_field_solver_im.free();
		}

		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			stiffness_eigenstrain::write_scalar(fout, phaseMesh);
			vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
			for (int ele_index = 0; ele_index < 6; ele_index++)
			{
				string compname = "\"stress_" + compNameV[ele_index] + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i)
						{
							fout << phaseMesh(i, j, k).customVec6s[ExternalFields::MECH_stress][ele_index] << endl;
						}
				fout << "</DataArray>" << endl;
			}
			for (int ele_index = 0; ele_index < 6; ele_index++)
			{
				string compname = "\"strain_" + compNameV[ele_index] + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i)
						{
							fout << phaseMesh(i, j, k).customVec6s[ExternalFields::MECH_strain][ele_index] << endl;
						}
				fout << "</DataArray>" << endl;
			}
			for (int ele_index = 0; ele_index < 6; ele_index++)
			{
				string compname = "\"eigenStrain_" + compNameV[ele_index] + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i)
						{
							fout << phaseMesh(i, j, k).customVec6s[ExternalFields::MECH_eigen_strain][ele_index] << endl;
						}
				fout << "</DataArray>" << endl;
			}

			fout << "<DataArray type = \"Float64\" Name = \"" << "stress_J1" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; ++k)
				for (int j = 0; j < phaseMesh.limit_y; ++j)
					for (int i = 0; i < phaseMesh.limit_x; ++i)
					{
						vStress stress(phaseMesh(i, j, k).customVec6s[ExternalFields::MECH_stress]);
						fout << stress.J1() << endl;
					}
			fout << "</DataArray>" << endl;

			fout << "<DataArray type = \"Float64\" Name = \"" << "stress_vMises" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; ++k)
				for (int j = 0; j < phaseMesh.limit_y; ++j)
					for (int i = 0; i < phaseMesh.limit_x; ++i)
					{
						vStress stress(phaseMesh(i, j, k).customVec6s[ExternalFields::MECH_stress]);
						fout << stress.Mises() << endl;
					}
			fout << "</DataArray>" << endl;
		}

		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			stiffness_eigenstrain::write_vec3(fout, phaseMesh);
			if (!output_displacement_field)
				return;
			if (elastic_solver::_MFType() == MechanicalFieldType::MFType_Explicit) {
				write_vec3_ex(fout, phaseMesh);
			}
			else if (elastic_solver::_MFType() == MechanicalFieldType::MFType_Implicit_Steinbach
				|| elastic_solver::_MFType() == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				write_vec3_im(fout, phaseMesh);
			}
		}

	}
}