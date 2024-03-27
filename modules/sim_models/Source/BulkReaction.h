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
#include "../../sim_postprocess/ElectricField.h"

namespace pf {
	namespace bulk_reaction {
		// - 

		//------phi source-----//

		static double (*reaction_a)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function

		static double reaction_a_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}

		// From http://dx.doi.org/10.1016/j.jpowsour.2015.09.055
		// Li dendrite growth, Butler Volmer type source
		static double reaction_a_electrode_reaction(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			bool is_cal = false;
			for (auto index = electric_field::electrode_index.begin(); index < electric_field::electrode_index.end(); index++)
				if (*index == phase.index)
					is_cal = true;  // phase is the solid phase

			double result = 0.0;
			if (is_cal) {
				double charge_trans_coeff{ 0.5 }, & alpha{ charge_trans_coeff };
				double& L_eta{ electric_field::reaction_constant }, & xi{ phase.phi }, & n{ electric_field::electron_num };

				auto h_ = 30.0 * xi * xi * (1 - xi) * (1 - xi);

				double phi_electrode{ electric_field::fix_domain_phi_value(phase.property) };
				double phi_solution{ node.customValues[ElectricFieldIndex::ElectricalPotential] };
				auto eta_a = (phi_electrode - phi_solution - electric_field::E_std);
				eta_a < -0.2 ? eta_a = -0.2 : eta_a>0 ? eta_a = 0 : 0;

				result = -L_eta * h_ * (std::exp(0.5 * n * FaradayConstant * eta_a / (GAS_CONSTANT * ROOM_TEMP)) -
					node.x[phase.index].value * std::exp(-0.5 * n * FaradayConstant * eta_a / (GAS_CONSTANT * ROOM_TEMP)));
			}
			return result;
		}

		//------concentration source-----//

		static void (*reaction_A)(pf::PhaseNode& node, pf::PhaseEntry& phase);   // main function

		static void reaction_A_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}

		static double (*reaction_i)(pf::PhaseNode& node, int con_i);   // main function

		static double reaction_i_none(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}

		static double reaction_i_electrode_reaction(pf::PhaseNode& node, int con_i) {
			double result = 0.0;
			// iterate through every components, pick up needed elements (cation)
			if (con_i != electric_field::active_component_index)
				return result;
			// iterate through every phases, pick up solid phases
			double dxi_dt{};
			for (auto phase = node.begin(); phase < node.end(); phase++)
				for (auto index = electric_field::electrode_index.begin(); index < electric_field::electrode_index.end(); index++)
					if (*index == phase->index) {
						dxi_dt += (phase->phi - phase->old_phi) / electric_field::time_interval;
					}
			Vector3 grad_phi{ node.cal_customValues_gradient(ElectricFieldIndex::ElectricalPotential) };
			double lap_phi{ node.cal_customValues_laplace(ElectricFieldIndex::ElectricalPotential) };

			double con{ node.x[con_i].value };
			Vector3 grad_con{ node.potential[con_i].gradient };

			double D_eff{ node.kinetics_coeff(con_i,con_i).value };
			Vector3 grad_D_eff{ node.kinetics_coeff.get_gradientVec3(con_i, con_i) };

			double& n{ electric_field::electron_num };
			double temp_const{ n * FaradayConstant / (GAS_CONSTANT * ROOM_TEMP) };

			auto source_potential{ temp_const * (D_eff * (grad_con * grad_phi) + con * (grad_phi * grad_D_eff)) + con * D_eff * lap_phi };
			auto source_xi{ -electric_field::c_s / electric_field::c_0 * dxi_dt };

			return source_potential + source_xi;
		}

		//------Temperature source-----//

		static double reaction_T_none(pf::PhaseNode& node) {
			return 0.0;
		}
		static double (*reaction_T)(pf::PhaseNode& node);   // main function


		//------ INTERFACE -----//
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			reaction_a = reaction_a_none;
			reaction_A = reaction_A_none;
			reaction_i = reaction_i_none;
			reaction_T = reaction_T_none;

			string active_comp_name{};
			bool is_electric_field_on{};
			InputFileReader::get_instance()->read_bool_value("Postprocess.PhysicalFields.electric", is_electric_field_on, false);
			if (InputFileReader::get_instance()->read_string_value("ModelsManager.PhiCon.ElectroDeposition.active_component", active_comp_name, infile_debug) && is_electric_field_on) {
				reaction_a = reaction_a_electrode_reaction;
				reaction_i = reaction_i_electrode_reaction;
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Bulk Reaction !\n", LOG_FILE_NAME);
		}

		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++)
				for (auto index = electric_field::electrode_index.begin(); index < electric_field::electrode_index.end(); index++)
					if (*index == phase->index) {
						bool is_defined = false;
						for (auto elecpotential = electric_field::fix_domain_phi_value.begin();
							elecpotential < electric_field::fix_domain_phi_value.end(); elecpotential++)
							if (elecpotential->index == phase->property)
								is_defined = true;
						if (!is_defined) {
							Solvers::get_instance()->writer.add_string_to_txt_and_screen("# ERROR : The electrode electric potential hasn't been defined in Modules.ElectricField.fix_phi ! \n", LOG_FILE_NAME);
							exit(0);
						}
					}
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