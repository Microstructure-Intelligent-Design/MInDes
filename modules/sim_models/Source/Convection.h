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
	namespace convection {
		static double dr = 1.0;
		static double threshold = 0.5;
		// Phi
		static double convection_a_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		static double convection_a_standard(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double div = (node.get_neighbor_node(Direction::x_up).customVec3s[ExternalFields::FLUID_velocity][0] -
				node.get_neighbor_node(Direction::x_down).customVec3s[ExternalFields::FLUID_velocity][0]) / 2.0 / dr +
				(node.get_neighbor_node(Direction::y_up).customVec3s[ExternalFields::FLUID_velocity][1] -
					node.get_neighbor_node(Direction::y_down).customVec3s[ExternalFields::FLUID_velocity][1]) / 2.0 / dr + 
				(node.get_neighbor_node(Direction::z_up).customVec3s[ExternalFields::FLUID_velocity][2] -
					node.get_neighbor_node(Direction::z_down).customVec3s[ExternalFields::FLUID_velocity][2]) / 2.0 / dr;
			return -(node.customVec3s[ExternalFields::FLUID_velocity] * phase.phi_grad + phase.phi * div);
		}
		static double (*convection_a)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function

		static double convection_ab_none(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return 0.0;
		}
		static double convection_ab_standard(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta) {
			return -(node.customVec3s[ExternalFields::FLUID_velocity] * (alpha.phi_grad - beta.phi_grad) / 2.0);
		}
		static double (*convection_ab)(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta);  // main function

		// Concentration
		static void convection_A_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static void convection_A_standard(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			for (auto x = phase.x.begin(); x < phase.x.end(); x++) {
				Vector3 grad_x;
				PhaseEntry* x_down_phase = &node.get_neighbor_node(Direction::x_down)[phase.index],
					* x_up_phase = &node.get_neighbor_node(Direction::x_up)[phase.index],
					* y_down_phase = &node.get_neighbor_node(Direction::y_down)[phase.index],
					* y_up_phase = &node.get_neighbor_node(Direction::y_up)[phase.index],
					* z_down_phase = &node.get_neighbor_node(Direction::z_down)[phase.index],
					* z_up_phase = &node.get_neighbor_node(Direction::z_up)[phase.index];
				if (x_down_phase->phi < Phi_Num_Cut_Off)
					x_down_phase = &phase;
				if (x_up_phase->phi < Phi_Num_Cut_Off)
					x_up_phase = &phase;
				if (y_down_phase->phi < Phi_Num_Cut_Off)
					y_down_phase = &phase;
				if (y_up_phase->phi < Phi_Num_Cut_Off)
					y_up_phase = &phase;
				if (z_down_phase->phi < Phi_Num_Cut_Off)
					z_down_phase = &phase;
				if (z_up_phase->phi < Phi_Num_Cut_Off)
					z_up_phase = &phase;
				grad_x[0] = (x_up_phase->x[x->index].value - x_down_phase->x[x->index].value) / 2.0 / dr;
				grad_x[1] = (y_up_phase->x[x->index].value - y_down_phase->x[x->index].value) / 2.0 / dr;
				grad_x[2] = (z_up_phase->x[x->index].value - z_down_phase->x[x->index].value) / 2.0 / dr;
				x->ChemicalReactionFlux -= node.customVec3s[ExternalFields::FLUID_velocity] * grad_x;
			}
		}
		static void (*convection_A)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function

		static double convection_i_none(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		static double convection_i_con_standard(pf::PhaseNode& node, int con_i) {
			Vector3 grad_x;
			PhaseNode* x_down_node = &node.get_neighbor_node(Direction::x_down),
				* x_up_node = &node.get_neighbor_node(Direction::x_up),
				* y_down_node = &node.get_neighbor_node(Direction::y_down),
				* y_up_node = &node.get_neighbor_node(Direction::y_up),
				* z_down_node = &node.get_neighbor_node(Direction::z_down),
				* z_up_node = &node.get_neighbor_node(Direction::z_up);
			if (x_down_node->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
				x_down_node = &node;
			if (x_up_node->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
				x_up_node = &node;
			if (y_down_node->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
				y_down_node = &node;
			if (y_up_node->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
				y_up_node = &node;
			if (z_down_node->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
				z_down_node = &node;
			if (z_up_node->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
				z_up_node = &node;
			grad_x[0] = (x_up_node->x[con_i].value - x_down_node->x[con_i].value) / 2.0 / dr;
			grad_x[1] = (y_up_node->x[con_i].value - y_down_node->x[con_i].value) / 2.0 / dr;
			grad_x[2] = (z_up_node->x[con_i].value - z_down_node->x[con_i].value) / 2.0 / dr;
			return node.customVec3s[ExternalFields::FLUID_velocity] * grad_x * (-1.0);
		}
		static double convection_i_con_reverse(pf::PhaseNode& node, int con_i) {
			Vector3 grad_x;
			PhaseNode* x_down_node = &node.get_neighbor_node(Direction::x_down),
				* x_up_node = &node.get_neighbor_node(Direction::x_up),
				* y_down_node = &node.get_neighbor_node(Direction::y_down),
				* y_up_node = &node.get_neighbor_node(Direction::y_up),
				* z_down_node = &node.get_neighbor_node(Direction::z_down),
				* z_up_node = &node.get_neighbor_node(Direction::z_up);
			if (x_down_node->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
				x_down_node = &node;
			if (x_up_node->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
				x_up_node = &node;
			if (y_down_node->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
				y_down_node = &node;
			if (y_up_node->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
				y_up_node = &node;
			if (z_down_node->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
				z_down_node = &node;
			if (z_up_node->customValues[ExternalFields::CON_Smooth_Phi] > threshold)
				z_up_node = &node;
			grad_x[0] = (x_up_node->x[con_i].value - x_down_node->x[con_i].value) / 2.0 / dr;
			grad_x[1] = (y_up_node->x[con_i].value - y_down_node->x[con_i].value) / 2.0 / dr;
			grad_x[2] = (z_up_node->x[con_i].value - z_down_node->x[con_i].value) / 2.0 / dr;
			return node.customVec3s[ExternalFields::FLUID_velocity] * grad_x * (-1.0);
		}
		static double (*convection_i)(pf::PhaseNode& node, int con_i);  // main function

		// Temperature
		static double convection_T_none(pf::PhaseNode& node) {
			return 0.0;
		}
		static double convection_T_standard(pf::PhaseNode& node) {
			Vector3 T_grand = node.cal_temperature_gradient(dr);
			double div = (node.get_neighbor_node(Direction::x_up).customVec3s[ExternalFields::FLUID_velocity][0] -
				node.get_neighbor_node(Direction::x_down).customVec3s[ExternalFields::FLUID_velocity][0]) / 2.0 / dr +
				(node.get_neighbor_node(Direction::y_up).customVec3s[ExternalFields::FLUID_velocity][1] -
					node.get_neighbor_node(Direction::y_down).customVec3s[ExternalFields::FLUID_velocity][1]) / 2.0 / dr +
				(node.get_neighbor_node(Direction::z_up).customVec3s[ExternalFields::FLUID_velocity][2] -
					node.get_neighbor_node(Direction::z_down).customVec3s[ExternalFields::FLUID_velocity][2]) / 2.0 / dr;
			return -(node.customVec3s[ExternalFields::FLUID_velocity] * T_grand + node.temperature.T * div);
		}
		static double (*convection_T)(pf::PhaseNode& node);  // main function

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			dr = phaseMesh.dr;
			threshold = Solvers::get_instance()->C_Solver.threshold;
			convection_a = convection_a_none;
			convection_ab = convection_ab_none;
			convection_A = convection_A_none;
			convection_i = convection_i_none;
			convection_T = convection_T_none;
			bool is_model_load = false;
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				is_model_load = false;
				InputFileReader::get_instance()->read_bool_value("ModelsManager.Phi.convection", is_model_load, infile_debug);
				if (is_model_load)
					convection_ab = convection_ab_standard;
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				is_model_load = false;
				InputFileReader::get_instance()->read_bool_value("ModelsManager.Phi.convection", is_model_load, infile_debug);
				if (is_model_load)
					convection_a = convection_a_standard;
			}

			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				is_model_load = false;
				InputFileReader::get_instance()->read_bool_value("ModelsManager.Con.convection", is_model_load, infile_debug);
				if (is_model_load)
					convection_A = convection_A_standard;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX || Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				is_model_load = false;
				InputFileReader::get_instance()->read_bool_value("ModelsManager.Con.convection", is_model_load, infile_debug);
				if (is_model_load) {
					if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Standard)
						convection_i = convection_i_con_standard;
					else if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse)
						convection_i = convection_i_con_reverse;
				}
			}

			if (Solvers::get_instance()->parameters.TempEType == pf::TemperatureEquationType::TType_Standard) {
				is_model_load = false;
				InputFileReader::get_instance()->read_bool_value("ModelsManager.Temp.convection", is_model_load, infile_debug);
				if (is_model_load)
					convection_T = convection_T_standard;
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Convection !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			dr = phaseMesh.dr;
			threshold = Solvers::get_instance()->C_Solver.threshold;
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}