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
#include "LBM_funcs/BoundaryCondition.h"
#include "LBM_funcs/Source.h"
#include "LBM_funcs/EquiDistFunc.h"
#include "LBM_funcs/MacroVar.h"

namespace pf {
	enum Fluid_Solver_Type { FST_NONE, FST_Pressure_Correction, FST_Lattice_Boltzmann };
	namespace lattice_boltzmann {
		namespace fluid_lbm_solver_funcs {
			// FST_LBM_Difference
			static LBM fluid_lbm_solver;
			double momentum_accuracy = 1e-4;
			static double_box field_variable_init;
			static double PCT_dt = 0.0;
			static bool debug_solver = false;
			static int debug_output_step = 1000;
			static int max_iterate_steps = 0;
			// normal parameters
			static vector<double> w;
			static double Cs2 = 1.0 / 3.0;
			static double Cs4 = 1.0 / 9.0;
			static double dr = 1.0;
			static int Nx = 0;
			static int Ny = 0;
			static int Nz = 0;
			static vector<Vector3> d2q9_w;
			static vector<Vector3> d3q19_w;
			// two phase flow
			static bool is_two_phase_flow = false;
			static LBM field_lbm_two_phase_solver;
			
			// Source
			static double (*fluid_source_i)(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i);
			static double (*fluid_two_phase_source_i)(pf::PhaseNode& node, double tau, Vector3 prefactor, int LBM_F_INDEX, int LBM_F_i);

			static double (*f_eq_i)(int INDEX_i, double p_macro, double f_macro, Vector3& U);
			static double (*f_eq_two_phase_i)(int INDEX_i, double p_macro, double f_macro, Vector3& U);

			// fluid flow
			// fi = fi + mi
			static void init_distribution_functions_d2q9(pf::PhaseNode& node, int LBM_F_INDEX) {
				double density = lbm_boundary_condition::density(node);
				node.customValues[ExternalFields::FLUID_mass] = density;
				node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = density;
				for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++)
					node.customValues[LBM_F_INDEX + index] = w[index] * density;
			}
			static void init_two_phase_d2q9(pf::PhaseNode& node, int LBM_F_INDEX) {
				double phi = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
				node.customValues.add_double(DF_Macro_TwoPhase_f_macro_old, phi);
				node.customVec3s.add_vec(DF_Macro_TwoPhase_velocity_old, node.customVec3s[ExternalFields::FLUID_velocity]);
				for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++)
					node.customValues[LBM_F_INDEX + index] = w[index] * phi;
			}
			static void init_distribution_functions_d3q19(pf::PhaseNode& node, int LBM_F_INDEX) {
				double density = lbm_boundary_condition::density(node);
				node.customValues[ExternalFields::FLUID_mass] = density;
				node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = density;
				for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++)
					node.customValues[LBM_F_INDEX + index] = w[index] * density;
			}
			static void init_two_phase_d3q19(pf::PhaseNode& node, int LBM_F_INDEX) {
				double phi = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
				node.customValues.add_double(DF_Macro_TwoPhase_f_macro_old, phi);
				node.customVec3s.add_vec(DF_Macro_TwoPhase_velocity_old, node.customVec3s[ExternalFields::FLUID_velocity]);
				for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++)
					node.customValues[LBM_F_INDEX + index] = w[index] * phi;
			}
			static void collision_SRT_d2q9(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					Vector3& U = node.customVec3s[ExternalFields::FLUID_velocity];
					double f_eq = 1.0, f_macro = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro], 
						p_macro = node.customValues[ExternalFields::FLUID_pressure], tau = lbm_boundary_condition::tau(node);
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++) {
						f_eq = f_eq_i(index, p_macro, f_macro, U);
						node.customValues[LBM_F_INDEX + index + delt_index] = (f_eq - node.customValues[LBM_F_INDEX + index]) / tau
							+ PCT_dt * fluid_source_i(node, tau, LBM_F_INDEX, index);
					}
				}
				else {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++) {
						node.customValues[LBM_F_INDEX + index + delt_index] = 0.0;
					}
				}
			}
			static void collision_two_phase_d2q9(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					Vector3& U = node.customVec3s[ExternalFields::FLUID_velocity];
					Vector3 phi_grad = node.cal_customValues_gradient(LBM_F_INDEX + LBM_Symbols::LBM_f_macro, dr);
					double f_eq = 1.0, f_macro = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro], p_macro = 0.0, tau = lbm_boundary_condition::tau_two_phase(node);
					Vector3 prefactor = U * (f_macro - node.customValues[DF_Macro_TwoPhase_f_macro_old]) / PCT_dt + (U - node.customVec3s[DF_Macro_TwoPhase_velocity_old]) / PCT_dt * f_macro
						+ phi_grad.normalize() * Cs2 * 4.0 * f_macro * (1.0 - f_macro) / lbm_source::get_interface_thickness();
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++) {
						f_eq = f_eq_two_phase_i(index, p_macro, f_macro, U);
						node.customValues[LBM_F_INDEX + index + delt_index] = (f_eq - node.customValues[LBM_F_INDEX + index]) / tau
							+ PCT_dt * fluid_two_phase_source_i(node, tau, prefactor, LBM_F_INDEX, index);
					}
					node.customValues[DF_Macro_TwoPhase_f_macro_old] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					node.customVec3s[DF_Macro_TwoPhase_velocity_old] = node.customVec3s[ExternalFields::FLUID_velocity];
				}
				else {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++) {
						node.customValues[LBM_F_INDEX + index + delt_index] = 0.0;
					}
					node.customValues[DF_Macro_TwoPhase_f_macro_old] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					node.customVec3s[DF_Macro_TwoPhase_velocity_old] = node.customVec3s[ExternalFields::FLUID_velocity];
				}
			}
			static void collision_SRT_d3q19(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					Vector3& U = node.customVec3s[ExternalFields::FLUID_velocity];
					double f_eq = 1.0, density = node.customValues[ExternalFields::FLUID_mass],
						p_macro = node.customValues[ExternalFields::FLUID_pressure], tau = lbm_boundary_condition::tau(node);
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++) {
						f_eq = f_eq_i(index, p_macro, density, U);
						node.customValues[LBM_F_INDEX + index + delt_index] = (f_eq - node.customValues[LBM_F_INDEX + index]) / tau
							+ PCT_dt * fluid_source_i(node, tau, LBM_F_INDEX, index);
					}
				}
				else {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++) {
						node.customValues[LBM_F_INDEX + index + delt_index] = 0.0;
					}
				}
			}
			static void collision_two_phase_d3q19(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					Vector3& U = node.customVec3s[ExternalFields::FLUID_velocity];
					Vector3 phi_grad = node.cal_customValues_gradient(LBM_F_INDEX + LBM_Symbols::LBM_f_macro, dr);
					double f_eq = 1.0, f_macro = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro],
						p_macro = 0.0, tau = lbm_boundary_condition::tau_two_phase(node);
					Vector3 prefactor = U * (f_macro - node.customValues[DF_Macro_TwoPhase_f_macro_old]) / PCT_dt + (U - node.customVec3s[DF_Macro_TwoPhase_velocity_old]) / PCT_dt * f_macro
						+ phi_grad.normalize() * Cs2 * 4.0 * f_macro * (1.0 - f_macro) / lbm_source::get_interface_thickness();
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++) {
						f_eq = f_eq_two_phase_i(index, p_macro, f_macro, U);
						node.customValues[LBM_F_INDEX + index + delt_index] = (f_eq - node.customValues[LBM_F_INDEX + index]) / tau
							+ PCT_dt * fluid_two_phase_source_i(node, tau, prefactor, LBM_F_INDEX, index);
					}
					node.customValues[DF_Macro_TwoPhase_f_macro_old] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					node.customVec3s[DF_Macro_TwoPhase_velocity_old] = node.customVec3s[ExternalFields::FLUID_velocity];
				}
				else {
					int delt_index = LBM_Symbols::LBM_m_0 - LBM_Symbols::LBM_f_0;
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++) {
						node.customValues[LBM_F_INDEX + index + delt_index] = 0.0;
					}
					node.customValues[DF_Macro_TwoPhase_f_macro_old] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					node.customVec3s[DF_Macro_TwoPhase_velocity_old] = node.customVec3s[ExternalFields::FLUID_velocity];
				}
			}
		}
		static void lbm_properties_automatically_change(FieldStorage_forPhaseNode& phaseMesh) {
			using namespace fluid_lbm_solver_funcs;
			PCT_dt = Solvers::get_instance()->parameters.dt;
			double cc = phaseMesh.dr / PCT_dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			lbm_boundary_condition::lbm_properties_automatically_change(phaseMesh);
			lbm_equilibrium_distribution_function::lbm_properties_automatically_change(phaseMesh);
			lbm_source::lbm_properties_automatically_change(phaseMesh);
			lbm_macro_variable::lbm_properties_automatically_change(phaseMesh);
		}
		static void init_two_phasae(FieldStorage_forPhaseNode& phaseMesh) {
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.customValues[DF_Macro_TwoPhase + LBM_Symbols::LBM_f_macro] = 1.0;
					}
		}
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			using namespace fluid_lbm_solver_funcs;
			fluid_lbm_solver.init(phaseMesh, DF_Macro_Density);
			fluid_lbm_solver.solver_name = "density";
			Nx = phaseMesh.limit_x;
			Ny = phaseMesh.limit_y;
			Nz = phaseMesh.limit_z;
			PCT_dt = Solvers::get_instance()->parameters.dt;
			dr = phaseMesh.dr;
			double cc = phaseMesh.dr / PCT_dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.max_iterate_steps", max_iterate_steps, infile_debug);
			if (InputFileReader::get_instance()->read_bool_value("Postprocess.FluidDynamics.LatticeBoltzmann.debug_solver", debug_solver, infile_debug))
				InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.debug_output_step", debug_output_step, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.momentum_accuracy", momentum_accuracy, infile_debug);
			lbm_boundary_condition::init(phaseMesh, fluid_lbm_solver);
			lbm_source::init(phaseMesh, fluid_lbm_solver);
			lbm_equilibrium_distribution_function::init(phaseMesh, fluid_lbm_solver);
			lbm_macro_variable::init(phaseMesh, fluid_lbm_solver);
			InputFileReader::get_instance()->read_bool_value("Postprocess.FluidDynamics.LatticeBoltzmann.two_phase_flow", is_two_phase_flow, infile_debug);
			if (is_two_phase_flow) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Preprocess.Microstructure.geometry_layer_x.custom_double = [(-882,liquid_fraction)] \n", InputFileReader::get_instance()->debug_file);
				field_lbm_two_phase_solver.init(phaseMesh, DF_Macro_TwoPhase);
				//init_two_phasae(phaseMesh);
				field_lbm_two_phase_solver.solver_name = "twoPhase";
				lbm_boundary_condition::init_two_phase_solver(phaseMesh, field_lbm_two_phase_solver);
				lbm_source::init_two_phase_solver(phaseMesh, field_lbm_two_phase_solver);
				fluid_two_phase_source_i = lbm_source::fluid_two_phase_source_i;
				lbm_equilibrium_distribution_function::init_two_phase(phaseMesh, field_lbm_two_phase_solver);
				f_eq_two_phase_i = lbm_equilibrium_distribution_function::f_eq_two_phase_i;
				lbm_macro_variable::init_two_phase(phaseMesh, fluid_lbm_solver, field_lbm_two_phase_solver);
				if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
					field_lbm_two_phase_solver._init_distribution_functions = fluid_lbm_solver_funcs::init_two_phase_d2q9;
					field_lbm_two_phase_solver._collision = fluid_lbm_solver_funcs::collision_two_phase_d2q9;
				}
				else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
					field_lbm_two_phase_solver._init_distribution_functions = fluid_lbm_solver_funcs::init_two_phase_d3q19;
					field_lbm_two_phase_solver._collision = fluid_lbm_solver_funcs::collision_two_phase_d3q19;
				}
			}
			fluid_source_i = lbm_source::fluid_source_i;
			f_eq_i = lbm_equilibrium_distribution_function::f_eq_i;

			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				w.resize(9);
				w[0] = 4.0 / 9.0;
				w[1] = 1.0 / 9.0; w[2] = 1.0 / 9.0; w[3] = 1.0 / 9.0; w[4] = 1.0 / 9.0;
				w[5] = 1.0 / 36.0; w[6] = 1.0 / 36.0; w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0;
				fluid_lbm_solver._init_distribution_functions = fluid_lbm_solver_funcs::init_distribution_functions_d2q9;
				fluid_lbm_solver._collision = fluid_lbm_solver_funcs::collision_SRT_d2q9;
				d2q9_w.push_back(Vector3(0.0, 0.0, 0.0));	 // f0  c( 0,  0)
				d2q9_w.push_back(Vector3(1.0, 0.0, 0.0));	 // f1  c( 1,  0)
				d2q9_w.push_back(Vector3(0.0, 1.0, 0.0));	 // f2  c( 0,  1)
				d2q9_w.push_back(Vector3(-1.0, 0.0, 0.0));	 // f3  c(-1,  0)
				d2q9_w.push_back(Vector3(0.0, -1.0, 0.0));	 // f4  c( 0, -1)
				d2q9_w.push_back(Vector3(1.0, 1.0, 0.0));	 // f5  c( 1,  1)
				d2q9_w.push_back(Vector3(-1.0, 1.0, 0.0));	 // f6  c(-1,  1)
				d2q9_w.push_back(Vector3(-1.0, -1.0, 0.0));	 // f7  c(-1, -1)
				d2q9_w.push_back(Vector3(1.0, -1.0, 0.0));	 // f8  c( 1, -1)
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				w.resize(19);
				w[0] = 12.0 / 36.0;
				w[1] = 2.0 / 36.0; w[2] = 2.0 / 36.0; w[3] = 2.0 / 36.0; w[4] = 2.0 / 36.0; w[5] = 2.0 / 36.0; w[6] = 2.0 / 36.0;
				w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0; w[9] = 1.0 / 36.0; w[10] = 1.0 / 36.0; w[11] = 1.0 / 36.0; w[12] = 1.0 / 36.0;
				w[13] = 1.0 / 36.0; w[14] = 1.0 / 36.0; w[15] = 1.0 / 36.0; w[16] = 1.0 / 36.0; w[17] = 1.0 / 36.0; w[18] = 1.0 / 36.0;
				fluid_lbm_solver._init_distribution_functions = fluid_lbm_solver_funcs::init_distribution_functions_d3q19;
				fluid_lbm_solver._collision = fluid_lbm_solver_funcs::collision_SRT_d3q19;
				d3q19_w.push_back(Vector3(0.0, 0.0, 0.0));	 // f0   c( 0,  0,  0)
				d3q19_w.push_back(Vector3(1.0, 0.0, 0.0));	 // f1   c( 1,  0,  0)
				d3q19_w.push_back(Vector3(-1.0, 0.0, 0.0));	 // f2   c(-1,  0,  0)
				d3q19_w.push_back(Vector3(0.0, 1.0, 0.0));	 // f3   c( 0,  1,  0)
				d3q19_w.push_back(Vector3(0.0, -1.0, 0.0));	 // f4   c( 0, -1,  0)
				d3q19_w.push_back(Vector3(0.0, 0.0, 1.0));	 // f5   c( 0,  0,  1)
				d3q19_w.push_back(Vector3(0.0, 0.0, -1.0));	 // f6   c( 0,  0, -1)
				d3q19_w.push_back(Vector3(1.0, 1.0, 0.0));	 // f7   c( 1,  1,  0)
				d3q19_w.push_back(Vector3(-1.0, 1.0, 0.0));	 // f8   c(-1,  1,  0)
				d3q19_w.push_back(Vector3(-1.0, -1.0, 0.0)); // f9   c(-1, -1,  0)
				d3q19_w.push_back(Vector3(1.0, -1.0, 0.0));	 // f10  c( 1, -1,  0)
				d3q19_w.push_back(Vector3(1.0, 0.0, 1.0));	 // f11  c( 1,  0,  1)
				d3q19_w.push_back(Vector3(-1.0, 0.0, 1.0));	 // f12  c(-1,  0,  1)
				d3q19_w.push_back(Vector3(-1.0, 0.0, -1.0)); // f13  c(-1,  0, -1)
				d3q19_w.push_back(Vector3(1.0, 0.0, -1.0));	 // f14  c( 1,  0, -1)
				d3q19_w.push_back(Vector3(0.0, 1.0, 1.0));	 // f15  c( 0,  1,  1)
				d3q19_w.push_back(Vector3(0.0, -1.0, 1.0));	 // f16  c( 0, -1,  1)
				d3q19_w.push_back(Vector3(0.0, -1.0, -1.0)); // f17  c( 0, -1, -1)
				d3q19_w.push_back(Vector3(0.0, 1.0, -1.0));	 // f18  c( 0,  1, -1)
			}

		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			using namespace fluid_lbm_solver_funcs;
			stringstream output;
			fluid_lbm_solver.init_distribution_functions();
			lbm_boundary_condition::cal_fluid_domain(phaseMesh);
			fluid_lbm_solver.boundary_condition();
			fluid_lbm_solver.cal_macro_variables();
			if (is_two_phase_flow) {
				field_lbm_two_phase_solver.init_distribution_functions();
				field_lbm_two_phase_solver.boundary_condition();
				double buff = 0.0;
				field_lbm_two_phase_solver.cal_macro_variables(buff);
			}
			if (debug_solver) {
				output << "> Fluid field solver debug:" << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
			Vector3 MAX_MOMENTUM_CHANGE; double MAX_F_MACRO_CHANGE = 0.0;
			int istep = 0;
			for (istep = 1; istep <= max_iterate_steps; istep++) {
				fluid_lbm_solver.collision();
				fluid_lbm_solver.streaming();
				fluid_lbm_solver.boundary_condition();
				MAX_MOMENTUM_CHANGE = fluid_lbm_solver.cal_macro_variables();
				if (is_two_phase_flow) {
					field_lbm_two_phase_solver.collision();
					field_lbm_two_phase_solver.streaming();
					field_lbm_two_phase_solver.boundary_condition();
					field_lbm_two_phase_solver.cal_macro_variables(MAX_F_MACRO_CHANGE);
				}
				if (debug_solver && (istep % debug_output_step == 0 || istep == 1)) {
					output.str("");
					output << "> Fluid field iterate step:" << istep << endl;
					output << "		MAX_MOMENTUM_CHANGE_X = " << MAX_MOMENTUM_CHANGE[0] << endl;
					output << "		MAX_MOMENTUM_CHANGE_Y = " << MAX_MOMENTUM_CHANGE[1] << endl;
					output << "		MAX_MOMENTUM_CHANGE_Z = " << MAX_MOMENTUM_CHANGE[2] << endl;
					if (is_two_phase_flow)
						output << "		TWO_PHASE_VARIATION   = " << MAX_F_MACRO_CHANGE << endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
				}
				if (MAX_MOMENTUM_CHANGE[0] < momentum_accuracy && MAX_MOMENTUM_CHANGE[1] < momentum_accuracy && MAX_MOMENTUM_CHANGE[2] < momentum_accuracy)
					break;
			}
			istep--;
			output.str("");
			output << "> Fluid field solver:" << endl;
			output << "		MAX_ITERATE_TIMES = " << istep << endl;
			output << "		MAX_MOMENTUM_CHANGE_X = " << MAX_MOMENTUM_CHANGE[0] << endl;
			output << "		MAX_MOMENTUM_CHANGE_Y = " << MAX_MOMENTUM_CHANGE[1] << endl;
			output << "		MAX_MOMENTUM_CHANGE_Z = " << MAX_MOMENTUM_CHANGE[2] << endl;
			if (is_two_phase_flow)
				output << "		TWO_PHASE_VARIATION   = " << MAX_F_MACRO_CHANGE << endl;
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			using namespace fluid_lbm_solver_funcs;
			lbm_properties_automatically_change(phaseMesh);
			stringstream report;
			Vector3 MAX_MOMENTUM_CHANGE; double MAX_F_MACRO_CHANGE = 0.0;
			lbm_boundary_condition::cal_fluid_domain(phaseMesh);
			int istep = 0;
			for (istep = 1; istep <= max_iterate_steps; istep++) {
				fluid_lbm_solver.collision();
				fluid_lbm_solver.streaming();
				fluid_lbm_solver.boundary_condition();
				MAX_MOMENTUM_CHANGE = fluid_lbm_solver.cal_macro_variables();
				if (is_two_phase_flow) {
					field_lbm_two_phase_solver.collision();
					field_lbm_two_phase_solver.streaming();
					field_lbm_two_phase_solver.boundary_condition();
					field_lbm_two_phase_solver.cal_macro_variables(MAX_F_MACRO_CHANGE);
				}
				if (MAX_MOMENTUM_CHANGE[0] < momentum_accuracy && MAX_MOMENTUM_CHANGE[1] < momentum_accuracy && MAX_MOMENTUM_CHANGE[2] < momentum_accuracy)
					break;
			}
			istep--;
			report.str("");
			report << "> Fluid field solver:" << endl;
			report << "		MAX_ITERATE_TIMES = " << istep << endl;
			report << "		MAX_MOMENTUM_CHANGE_X = " << MAX_MOMENTUM_CHANGE[0] << endl;
			report << "		MAX_MOMENTUM_CHANGE_Y = " << MAX_MOMENTUM_CHANGE[1] << endl;
			report << "		MAX_MOMENTUM_CHANGE_Z = " << MAX_MOMENTUM_CHANGE[2] << endl;
			if (is_two_phase_flow)
				report << "		TWO_PHASE_VARIATION   = " << MAX_F_MACRO_CHANGE << endl;
			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			using namespace fluid_lbm_solver_funcs;
			fluid_lbm_solver.free();
			field_lbm_two_phase_solver.free();
			field_variable_init.clear();
			w.clear();
			d2q9_w.clear();
			d3q19_w.clear();
			lbm_boundary_condition::deinit();
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			using namespace fluid_lbm_solver_funcs;
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_speed" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity].abs() << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_pressure" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customValues[ExternalFields::FLUID_pressure] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_density" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customValues[ExternalFields::FLUID_mass] << endl;
					}
			fout << "</DataArray>" << endl;
			if (is_two_phase_flow) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "two_phase" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							fout << phaseMesh(i, j, k).customValues[ExternalFieldsPlus::DF_Macro_TwoPhase + LBM_Symbols::LBM_f_macro] << endl;
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			string name;
			name = "\"fluid_velocity\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity][0] << " "
							<< phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity][1] << " "
							<< phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity][2] << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
}