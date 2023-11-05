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
#include "../../../sim_models/InterfaceEnergy.h"
#include "../../../sim_models/BulkEnergy.h"

namespace pf {
	enum LBM_Source_Type { LBM_ST_Forces };
	enum LBM_Force_Term_Model { LBM_FTM_NONE, LBM_FTM_LGA, LBM_FTM_ZS_Guo, LBM_FTM_H_Liang };
	enum LBM_Force_Type { LBM_FM_ThermalExpansion, LBM_FM_Gravity, LBM_FM_H_Liang_SurfaceTension };
	namespace lbm_source {
		static vector<Vector3> d2q9_w;
		static vector<Vector3> d3q19_w;
		static vector<double> w;
		static double Cs2 = 1.0 / 3.0;
		static double Cs4 = 1.0 / 9.0;
		static double dr = 1.0;
		// two phase flow
		static double surface_tension = 0.0;
		static double interface_thickness = 1.0;
		static double beta = 0.0;
		static double kappa = 0.0;
		static double get_interface_thickness() {
			return interface_thickness;
		}
		namespace force_funcs {
			// external force parameters
			static Vector3 gravitational_acceleration(0.0, 0.0, 0.0);
			// fluid force
			static double ref_density = 0.0;
			static Vector3 Fluid_Force_Gravity(pf::PhaseNode& node, int LBM_F_INDEX) {
				return gravitational_acceleration * (node.customValues[ExternalFields::FLUID_mass] - ref_density);
			}
			static Vector3 Fluid_Force_H_Liang_Surface_Tension(pf::PhaseNode& node, int LBM_F_INDEX) {
				double phi = node.customValues[LBM_F_INDEX + LBM_Symbols::LBM_f_macro];
				return node.cal_customValues_gradient(LBM_F_INDEX + LBM_Symbols::LBM_f_macro, dr) * 
					(4.0 * beta * phi * (phi - 1) * (phi - 0.5) - kappa * node.cal_customValues_laplace(LBM_F_INDEX + LBM_Symbols::LBM_f_macro, dr));
			}
			static double ref_temp = 0.0;
			static double thermal_expansion_parameter = 0.0;
			static Vector3 Fluid_Force_Thermal_Expansion(pf::PhaseNode& node, int LBM_F_INDEX) {
				return gravitational_acceleration * thermal_expansion_parameter * (ref_temp - node.temperature.T) * node.customValues[ExternalFields::FLUID_mass];
			}
			static vector<Vector3(*)(pf::PhaseNode&, int)> fluid_force_list;

			// fluid force model
			static double Fluid_Modified_LGA_d2q9_Force_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 force(0.0, 0.0, 0.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(node, LBM_F_INDEX);
				return w[LBM_F_i] * (d2q9_w[LBM_F_i] * force) / Cs2;
			}
			static double Fluid_Modified_LGA_d3q19_Force_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 force(0.0, 0.0, 0.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(node, LBM_F_INDEX);
				return w[LBM_F_i] * (d3q19_w[LBM_F_i] * force) / Cs2;
			}
			static double Fluid_Modified_GZS_d2q9_Force_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = d2q9_w[LBM_F_i], u = node.customVec3s[ExternalFields::FLUID_velocity];
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(node, LBM_F_INDEX);
				return (1.0 - 1.0 / 2.0 / tau) * w[LBM_F_i] * (((ci - u) / Cs2 + ci * (ci * u) / Cs4 / 2.0) * force);
			}
			static double Fluid_Modified_GZS_d3q19_Force_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = d3q19_w[LBM_F_i], u = node.customVec3s[ExternalFields::FLUID_velocity];
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(node, LBM_F_INDEX);
				return (1.0 - 1.0 / 2.0 / tau) * w[LBM_F_i] * (((ci - u) / Cs2 + ci * (ci * u) / Cs4 / 2.0) * force);
			}
			static double Fluid_Modified_HL_d2q9_Force_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = d2q9_w[LBM_F_i], u = node.customVec3s[ExternalFields::FLUID_velocity],
					delt_density = node.cal_customValues_gradient(ExternalFields::FLUID_mass, dr) * (-1.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(node, LBM_F_INDEX);
				return (1.0 - 1.0 / 2.0 / tau) * w[LBM_F_i] * ((ci * force) / Cs2 + (u * ci) * (delt_density * ci) / Cs2);
			}
			static double Fluid_Modified_HL_d3q19_Force_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 force(0.0, 0.0, 0.0), ci = d3q19_w[LBM_F_i], u = node.customVec3s[ExternalFields::FLUID_velocity],
					delt_density = node.cal_customValues_gradient(ExternalFields::FLUID_mass, dr) * (-1.0);
				for (auto f = fluid_force_list.begin(); f < fluid_force_list.end(); f++)
					force += (*f)(node, LBM_F_INDEX);
				return (1.0 - 1.0 / 2.0 / tau) * w[LBM_F_i] * ((ci * force) / Cs2 + (u * ci) * (delt_density * ci) / Cs2);
			}
			static double Fluid_Modified_two_phase_d2q9_Force_i(pf::PhaseNode& node, double tau, Vector3 prefactor, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 ci = d2q9_w[LBM_F_i];
				return (1.0 - 1.0 / 2.0 / tau) * w[LBM_F_i] * (ci * prefactor) / Cs2;
			}
			static double Fluid_Modified_two_phase_d3q19_Force_i(pf::PhaseNode& node, double tau, Vector3 prefactor, int LBM_F_INDEX, int LBM_F_i) {
				Vector3 ci = d3q19_w[LBM_F_i];
				return (1.0 - 1.0 / 2.0 / tau) * w[LBM_F_i] * (ci * prefactor) / Cs2;
			}
			static void load_forces(bool infile_debug){
				fluid_force_list.clear();
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.FluidDynamics.LatticeBoltzmann.force = () \n", InputFileReader::get_instance()->debug_file);
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("#             0 - ThermalExpansion, 1 - Gravity, 2 - H_Liang_SurfaceTension \n", InputFileReader::get_instance()->debug_file);
				string force_key = "Postprocess.FluidDynamics.LatticeBoltzmann.force", force_input = "()";
				InputFileReader::get_instance()->read_string_value(force_key, force_input, infile_debug);
				vector<input_value> force_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, force_key, force_input, infile_debug);
				for (auto force = force_value.begin(); force < force_value.end(); force++) {
					bool is_already_load = false;
					for (auto old_force = force_value.begin(); old_force < force; old_force++)
						if (force->int_value == old_force->int_value)
							is_already_load = true;
					if (!is_already_load) {
						switch (LBM_Force_Type(force->int_value))
						{
						case LBM_FM_ThermalExpansion:
							fluid_force_list.push_back(Fluid_Force_Thermal_Expansion);
							InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.ThermalExpansionForce.reference_temperature", ref_temp, infile_debug);
							InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.ThermalExpansionForce.thermal_expansion_parameter", thermal_expansion_parameter, infile_debug);
							break;
						case LBM_FM_Gravity:
							InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.Gravity.reference_density", ref_density, infile_debug);
							fluid_force_list.push_back(Fluid_Force_Gravity);
							break;
						case LBM_FM_H_Liang_SurfaceTension:
							fluid_force_list.push_back(Fluid_Force_H_Liang_Surface_Tension);
							break;
						default:
							break;
						}
					}
				}
			}
		}

		static vector<double(*)(pf::PhaseNode&, double, int, int)> fluid_source_list;
		static double(*fluid_two_phase_source_list)(pf::PhaseNode&, double, Vector3, int, int) ;

		// source model
		static double fluid_source_i(pf::PhaseNode& node, double tau, int LBM_F_INDEX, int LBM_F_i) {
			double buff = 0.0;
			for (auto func = fluid_source_list.begin(); func < fluid_source_list.end(); func++)
				buff += (*func)(node, tau, LBM_F_INDEX, LBM_F_i);
			return buff;
		}
		static double fluid_two_phase_source_i(pf::PhaseNode& node, double tau, Vector3 prefactor, int LBM_F_INDEX, int LBM_F_i) {
			return fluid_two_phase_source_list(node, tau, prefactor, LBM_F_INDEX, LBM_F_i);
		}
		static void lbm_properties_automatically_change(FieldStorage_forPhaseNode& phaseMesh) {
			double cc = phaseMesh.dr / Solvers::get_instance()->parameters.dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
		}
		static void init(FieldStorage_forPhaseNode& phaseMesh, LBM& fluid_lbm_solver) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			dr = phaseMesh.dr;
			double cc = phaseMesh.dr / Solvers::get_instance()->parameters.dt;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
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
			// load source
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.FluidDynamics.LatticeBoltzmann.source = () \n", InputFileReader::get_instance()->debug_file);
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("#             0 - Forces \n", InputFileReader::get_instance()->debug_file);
			string source_key = "Postprocess.FluidDynamics.LatticeBoltzmann.source", source_input = "()";
			InputFileReader::get_instance()->read_string_value(source_key, source_input, infile_debug);
			vector<input_value> source_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, source_key, source_input, infile_debug);
			int force_model_type = LBM_Force_Term_Model::LBM_FTM_NONE;
			for (auto source = source_value.begin(); source < source_value.end(); source++) {
				bool is_already_load = false;
				for (auto old_source = source_value.begin(); old_source < source; old_source++)
					if (source->int_value == old_source->int_value)
						is_already_load = true;
				if (!is_already_load) {
					int gravity_direction = 0;
					switch (LBM_Source_Type(source->int_value))
					{
					case pf::LBM_ST_Forces:
						if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .gravity_direction : 0 - x_down, 1 - x_up, 2 - y_down, 3 - y_up, 4 - z_dowm, 5 - z_up \n", InputFileReader::get_instance()->debug_file);
						InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.Force.gravity_direction", gravity_direction, infile_debug);
						switch (gravity_direction)
						{
						case 0:  // down_x
							force_funcs::gravitational_acceleration = Vector3(-9.8, 0.0, 0.0);
							break;
						case 1:  // up_x
							force_funcs::gravitational_acceleration = Vector3(9.8, 0.0, 0.0);
							break;
						case 2:  // down_y
							force_funcs::gravitational_acceleration = Vector3(0.0, -9.8, 0.0);
							break;
						case 3:  // up_y
							force_funcs::gravitational_acceleration = Vector3(0.0, 9.8, 0.0);
							break;
						case 4:  // down_z
							force_funcs::gravitational_acceleration = Vector3(0.0, 0.0, -9.8);
							break;
						case 5:  // up_z
							force_funcs::gravitational_acceleration = Vector3(0.0, 0.0, 9.8);
							break;
						default:
							// down_x
							force_funcs::gravitational_acceleration = Vector3(-9.8, 0.0, 0.0);
							break;
						}
						if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .LatticeBoltzmann.force_model : 0 - NONE, 1 - LGA, 2 - ZS_Guo, 3 - H_Liang \n", InputFileReader::get_instance()->debug_file);
						InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.LatticeBoltzmann.force_model", force_model_type, infile_debug);
						if (force_model_type == LBM_Force_Term_Model::LBM_FTM_LGA) {
							if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_LGA_d2q9_Force_i);
							else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_LGA_d3q19_Force_i);
						}
						else if (force_model_type == LBM_Force_Term_Model::LBM_FTM_ZS_Guo) {
							if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_GZS_d2q9_Force_i);
							else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_GZS_d3q19_Force_i);
						}
						else if (force_model_type == LBM_Force_Term_Model::LBM_FTM_H_Liang) {
							if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_HL_d2q9_Force_i);
							else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
								fluid_source_list.push_back(force_funcs::Fluid_Modified_HL_d3q19_Force_i);
						}
						force_funcs::load_forces(infile_debug);
						break;
					default:
						break;
					}
				}
			}
		}
		static void init_two_phase_solver(FieldStorage_forPhaseNode& phaseMesh, LBM& field_lbm_two_phase_solver) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			if (infile_debug) {
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# kappa = 3.0 / 2.0 * surface_tension * interface_thickness \n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# beta  = 12.0 * surface_tension / interface_thickness \n", InputFileReader::get_instance()->debug_file);
			}
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.surface_tension", surface_tension, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.interface_thickness", interface_thickness, infile_debug);
			beta = 12.0 * surface_tension / interface_thickness;
			kappa = 3.0 / 2.0 * surface_tension * interface_thickness;

			if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9)
				fluid_two_phase_source_list = force_funcs::Fluid_Modified_two_phase_d2q9_Force_i;
			else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19)
				fluid_two_phase_source_list = force_funcs::Fluid_Modified_two_phase_d3q19_Force_i;
		}
	}
}