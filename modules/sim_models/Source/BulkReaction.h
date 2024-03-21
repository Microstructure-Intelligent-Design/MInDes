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
		static vector<int> electrode_index{};
		static int active_component_index = 0;
		// - 
		static double time_interval = 0.0;
		// - 
		static double reaction_constant{};
		static double electron_num{};
		static double E_std{};
		static double diff_coef_ele{}, diff_coef_sol{};
		static double c_s{}, c_0{};
		//------phi source-----//

		static double (*reaction_a)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function

		static double reaction_a_none(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}

		static inline double interpolate_func(double xi) {
			return 30.0 * xi * xi * (1 - xi) * (1 - xi);
		}
		auto& h_ = interpolate_func;
		// From http://dx.doi.org/10.1016/j.jpowsour.2015.09.055
		// Li dendrite growth, Butler Volmer type source
		static double reaction_a_electrode_reaction(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			bool is_cal = false;
			for (auto index = electrode_index.begin(); index < electrode_index.end(); index++)
				if (*index == phase.index)
					is_cal = true;  // phase is the solid phase
			double result = 0.0;
			if (is_cal) {
				double charge_trans_coeff{ 0.5 }, & alpha{ charge_trans_coeff };
				double& L_eta{ reaction_constant }, & xi{ phase.phi }, & n{ electron_num };

				auto interpolate_func = [&xi]()->double {return 30.0 * xi * xi * (1 - xi) * (1 - xi); }, & h_{ interpolate_func };

				double phi_electrode{ electric_field::fix_domain_phi_value(phase.property) };
				double phi_solution{ node.customValues[ElectricFieldIndex::ElectricalPotential] };
				auto eta_a = [&phi_electrode, &phi_solution]()->double { return phi_electrode - phi_solution - E_std; };

				auto BV_exp = [&eta_a, &n](double _alpha) -> double {return std::exp(_alpha * n * FaradayConstant * eta_a() / (GAS_CONSTANT * ROOM_TEMP)); };

				result = -L_eta * h_() * (BV_exp(1 - alpha) - node.x[phase.index].value * BV_exp(-alpha));
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
			if (con_i != active_component_index)
				return result;
			// iterate through every phases, pick up solid phases
			double dxi_dt{};
			for(auto phase = node.begin(); phase < node.end(); phase++)
				for(auto index = electrode_index.begin(); index < electrode_index.end(); index++)
					if (*index == phase->index) {
						dxi_dt += (phase->phi - phase->old_phi) / time_interval;
					}

			Vector3 grad_phi{ node.cal_customValues_gradient(ElectricFieldIndex::ElectricalPotential) };
			double lap_phi{ node.cal_customValues_laplace(ElectricFieldIndex::ElectricalPotential) };

			double con{ node.potential[con_i].value };
			Vector3 grad_con{ node.potential[con_i].gradient };

			double D_eff{node.kinetics_coeff(con_i,con_i).value};
			Vector3 grad_D_eff{ node.kinetics_coeff.get_gradientVec3(con_i, con_i) };

			//go kinetics.h initialise De and Ds in that module

			double& n{ electron_num };
			double temp_const{ n * FaradayConstant / (GAS_CONSTANT * ROOM_TEMP) };
			auto source_potential = [&]()->double {return temp_const * (D_eff * (grad_con * grad_phi) + con * (grad_phi * grad_D_eff)) + con * D_eff * lap_phi; };

			auto source_xi = [&dxi_dt]()->double {return -c_s / c_0 * dxi_dt; };
			return source_potential() + source_xi();
		}

		//------Temperature source-----//

		static double reaction_T_none(pf::PhaseNode& node) {
			return 0.0;
		}
		static double (*reaction_T)(pf::PhaseNode& node);   // main function

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			reaction_a = reaction_a_none;
			reaction_A = reaction_A_none;
			reaction_i = reaction_i_none;
			reaction_T = reaction_T_none;

			time_interval = Solvers::get_instance()->parameters.dt;
			// active_component_index
			string active_comp_name = "";
			if (InputFileReader::get_instance()->read_string_value("ModelsManager.PhiCon.ElectroDeposition.active_component_index", active_comp_name, infile_debug)) {
				active_component_index = Solvers::get_instance()->parameters.Components[active_comp_name].index;
				
				string electrode_key = "ModelsManager.PhiCon.ElectroDeposition.electrode_index", electrode_input = "()";
				InputFileReader::get_instance()->read_string_value(electrode_key, electrode_input, infile_debug);
				vector<pf::input_value> electrode_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, electrode_key, electrode_input, infile_debug);
				for (int index = 0; index < electrode_value.size(); index++)
					electrode_index.push_back(electrode_value[index].int_value);

				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Butler_Volmer.Reaction_Constant", reaction_constant, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Butler_Volmer.Reaction_Electron_Num", electron_num, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Bulter_Volmer.Standard_Potential", E_std, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.DiffusionCoefficient.Electrode", diff_coef_ele, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.DiffusionCoefficient.Solution", diff_coef_sol, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.Bulter_Volmer.Electrode_Metal_SiteDensity", c_s, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.Bulter_Volmer.Electrolyte_Cation_Con", c_0, infile_debug);
				reaction_a = reaction_a_electrode_reaction;
				reaction_i = reaction_i_electrode_reaction;
			}

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Bulk Reaction !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			
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