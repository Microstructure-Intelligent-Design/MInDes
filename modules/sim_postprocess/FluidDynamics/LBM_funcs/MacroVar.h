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
#include "../../../Base.h"
#include "BoundaryCondition.h"
#include "Source.h"

namespace pf {
	namespace lbm_macro_variable {
		static double PCT_dt = 0.0;
		static double Cs2 = 1.0 / 3.0;
		static double Cs4 = 1.0 / 9.0;
		static double dr = 1.0;
		static vector<double> w;
		static vector<Vector3> d2q9_w;
		static vector<Vector3> d3q19_w;
		static int index_h_liang = 0;
		namespace macro_variable_funcs {
			static vector<Vector3(*)(pf::PhaseNode&, int)> fluid_force_list;
			static void load_forces() {
				fluid_force_list.clear();
				string force_key = "Postprocess.FluidDynamics.LatticeBoltzmann.force", force_input = "()";
				InputFileReader::get_instance()->read_string_value(force_key, force_input, false);
				vector<input_value> force_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, force_key, force_input, false);
				for (auto force = force_value.begin(); force < force_value.end(); force++) {
					bool is_already_load = false;
					for (auto old_force = force_value.begin(); old_force < force; old_force++)
						if (force->int_value == old_force->int_value)
							is_already_load = true;
					if (!is_already_load) {
						switch (LBM_Force_Type(force->int_value))
						{
						case LBM_FM_ThermalExpansion:
							fluid_force_list.push_back(lbm_source::force_funcs::Fluid_Force_Thermal_Expansion);
							break;
						case LBM_FM_Gravity:
							fluid_force_list.push_back(lbm_source::force_funcs::Fluid_Force_Gravity);
							break;
						case LBM_FM_H_Liang_SurfaceTension:
							fluid_force_list.push_back(lbm_source::force_funcs::Fluid_Force_H_Liang_Surface_Tension);
							break;
						default:
							break;
						}
					}
				}
			}

			static void cal_macro_variables_d2q9_standard(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double buff = 0.0;
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++) {
						buff = node.customValues[LBM_F_INDEX + index];
						node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] += buff;
						node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += d2q9_w[index] * buff;
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(node, LBM_F_INDEX);
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += force * PCT_dt * 0.5; // force
					node.customVec3s[ExternalFields::FLUID_velocity] = node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] / node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					node.customValues[ExternalFields::FLUID_pressure] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] * Cs2;
					node.customValues[ExternalFields::FLUID_mass] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					node.customValues[ExternalFields::FLUID_pressure] = 0.0;
					node.customValues[ExternalFields::FLUID_mass] = 0.0;
				}
			}
			static void cal_macro_variables_d3q19_standard(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double buff = 0.0;
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++) {
						buff = node.customValues[LBM_F_INDEX + index];
						node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] += buff;
						node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += d3q19_w[index] * buff;
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(node, LBM_F_INDEX);
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += force * PCT_dt * 0.5; // force
					node.customVec3s[ExternalFields::FLUID_velocity] = node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] / node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
					node.customValues[ExternalFields::FLUID_pressure] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] * Cs2;
					node.customValues[ExternalFields::FLUID_mass] = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					node.customValues[ExternalFields::FLUID_pressure] = 0.0;
					node.customValues[ExternalFields::FLUID_mass] = 0.0;
				}
			}

			const double w0_d2q9 = 4.0 / 9.0;
			static double s_0_d2q9_two_phase_flow(Vector3& U) {
				return -w0_d2q9 * (U * U) / Cs2 / 2.0;
			}
			static void cal_macro_variables_d2q9_H_Liang(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double buff = 0.0, density = lbm_boundary_condition::density(node);
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = density;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++) {
						buff += node.customValues[LBM_F_INDEX + index];
						node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += d2q9_w[index] * node.customValues[LBM_F_INDEX + index];
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(node, LBM_F_INDEX);
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += force * PCT_dt * 0.5; // force
					Vector3 Velocity = node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] / density;
					node.customVec3s[ExternalFields::FLUID_velocity] = Velocity;
					node.customValues[ExternalFields::FLUID_mass] = density;
					Vector3 delt_density = node.cal_customValues_gradient(ExternalFields::FLUID_mass, dr);
					buff -= node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_0];
					node.customValues[ExternalFields::FLUID_pressure] = Cs2 / (1.0 - w0_d2q9) 
						* (buff + PCT_dt * 0.5 * (Velocity * delt_density) + density * s_0_d2q9_two_phase_flow(Velocity));
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					node.customValues[ExternalFields::FLUID_pressure] = 0.0;
					node.customValues[ExternalFields::FLUID_mass] = 0.0;
				}
			}
			const double w0_d3q19 = 12.0 / 36.0;
			static double s_0_d3q19_two_phase_flow(Vector3& U) {
				return -w0_d3q19 * (U * U) / Cs2 / 2.0;
			}
			static void cal_macro_variables_d3q19_H_Liang(pf::PhaseNode& node, int LBM_F_INDEX) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					double buff = 0.0, density = lbm_boundary_condition::density(node);
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = density;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++) {
						buff += node.customValues[LBM_F_INDEX + index];
						node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += d3q19_w[index] * node.customValues[LBM_F_INDEX + index];
					}
					Vector3 force(0.0, 0.0, 0.0);
					for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
						force += (*f)(node, LBM_F_INDEX);
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] += force * PCT_dt * 0.5; // force
					Vector3 Velocity = node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro] / density;
					node.customVec3s[ExternalFields::FLUID_velocity] = Velocity;
					node.customValues[ExternalFields::FLUID_mass] = density;
					Vector3 delt_density = node.cal_customValues_gradient(ExternalFields::FLUID_mass, dr);
					buff -= node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_0];
					node.customValues[ExternalFields::FLUID_pressure] = Cs2 / (1.0 - w0_d3q19)
						* (buff + PCT_dt * 0.5 * (Velocity * delt_density) + density * s_0_d3q19_two_phase_flow(Velocity));
				}
				else {
					node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
					node.customVec3s[LBM_F_INDEX + LBM_Symbols::LBM_fv_macro].set_to_zero();
					node.customValues[ExternalFields::FLUID_pressure] = 0.0;
					node.customValues[ExternalFields::FLUID_mass] = 0.0;
				}
			}
			static void cal_macro_variables_d2q9_two_phase(pf::PhaseNode& node, int LBM_F_INDEX) {
				node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_8; index++)
						node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] += node.customValues[LBM_F_INDEX + index];
				}
			}
			static void cal_macro_variables_d3q19_two_phase(pf::PhaseNode& node, int LBM_F_INDEX) {
				node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] = 0.0;
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] >= lbm_boundary_condition::solid_liquid_interface_threshold) {
					for (int index = LBM_Symbols::LBM_f_0; index <= LBM_Symbols::LBM_f_18; index++)
						node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro] += node.customValues[LBM_F_INDEX + index];
				}
			}
		}

		static void (*cal_macro_variables)(pf::PhaseNode& node, int LBM_F_INDEX);

		static void (*cal_macro_variables_two_phase)(pf::PhaseNode& node, int LBM_F_INDEX);

		static void lbm_properties_automatically_change(FieldStorage_forPhaseNode& phaseMesh) {
			PCT_dt = Solvers::get_instance()->parameters.dt;
			double cc = phaseMesh.dr / PCT_dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh, LBM& fluid_lbm_solver) {
			// init parameters
			PCT_dt = Solvers::get_instance()->parameters.dt;
			dr = phaseMesh.dr;
			double cc = dr / PCT_dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d2q9_standard;
				w.resize(9);
				w[0] = 4.0 / 9.0;
				w[1] = 1.0 / 9.0; w[2] = 1.0 / 9.0; w[3] = 1.0 / 9.0; w[4] = 1.0 / 9.0;
				w[5] = 1.0 / 36.0; w[6] = 1.0 / 36.0; w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0;
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
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d3q19_standard;
				w.resize(19);
				w[0] = 12.0 / 36.0;
				w[1] = 2.0 / 36.0; w[2] = 2.0 / 36.0; w[3] = 2.0 / 36.0; w[4] = 2.0 / 36.0; w[5] = 2.0 / 36.0; w[6] = 2.0 / 36.0;
				w[7] = 1.0 / 36.0; w[8] = 1.0 / 36.0; w[9] = 1.0 / 36.0; w[10] = 1.0 / 36.0; w[11] = 1.0 / 36.0; w[12] = 1.0 / 36.0;
				w[13] = 1.0 / 36.0; w[14] = 1.0 / 36.0; w[15] = 1.0 / 36.0; w[16] = 1.0 / 36.0; w[17] = 1.0 / 36.0; w[18] = 1.0 / 36.0;
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
			// init cal_macro_variables function
			macro_variable_funcs::load_forces();
			// init solver
			fluid_lbm_solver._cal_macro_variables = cal_macro_variables;
		}
		static void init_two_phase(FieldStorage_forPhaseNode& phaseMesh, LBM& fluid_lbm_solver, LBM& field_lbm_two_phase_solver) {
			// init parameters
			if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				cal_macro_variables_two_phase = macro_variable_funcs::cal_macro_variables_d2q9_two_phase;
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d2q9_H_Liang;
			}
			else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				cal_macro_variables_two_phase = macro_variable_funcs::cal_macro_variables_d3q19_two_phase;
				cal_macro_variables = macro_variable_funcs::cal_macro_variables_d3q19_H_Liang;
			}
			// init solver
			fluid_lbm_solver._cal_macro_variables = cal_macro_variables;
			field_lbm_two_phase_solver._cal_macro_variables = cal_macro_variables_two_phase;
		}
	}
}