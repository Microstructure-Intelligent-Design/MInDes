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
	enum dfdphiType { // phi , phase con , grand potential , temperature
		dfdphi_Const,
		dfdphi_DOUBEL_WELL,
		dfdphi_LQ_CHEN,
		dfdcon_1_CRACK_WELL,
		dfdcon_1_CRACK_OBSTACLE,
		dfdcon_N_CRACK_WELL,
		dfdcon_N_CRACK_OBSTACLE,
		dfdphi_HighOrder,
		dfdphi_DOUBLE_WELL_SIMPLE
	};
	enum dfdconType { // phi , total con , temperature
		dfdcon_Const,
		dfdcon_HighOrder,
	};
	namespace chemical_energy {
		pf::ConEquationDomain _domain = pf::ConEquationDomain::CEDomain_Standard;
		vector<int> smooth_phases;
		//-----------------------------------------------------------------------------------------------
		// fchem for each phase
		static double dfchem_dphi_const(pf::PhaseNode& node, pf::PhaseEntry& phase){
			return 0.0;
		}
		// double well : V * phi_a * phi_a * (1.0 - phi_a) * (1.0 - phi_a) + W * sum_a{ sum_b!=a{ phi_a * phi_a * phi_b * phi_b } }
		static double DOUBLE_WELL_A = 0.0, DOUBLE_WELL_B = 0.0;
		static double dfchem_dphi_double_well(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double phi = phase.phi, sum_phi_j = 0.0;
			for (auto p = node.begin(); p < node.end(); p++)
				if (p->index != phase.index)
					sum_phi_j += p->phi * p->phi;
			return 2.0 * DOUBLE_WELL_A * phi * (1.0 - phi) * (1.0 - 2.0 * phi) + 2.0 * DOUBLE_WELL_B * phi * sum_phi_j;
		}
		// LQ. Chen : sum_a{ - A/2.0 * phi_a * phi_a + B/4.0 * phi_a * phi_a * phi_a * phi_a } + C * sum_a{ sum_b!=a{ phi_a * phi_a * phi_b * phi_b } }
		static double LQ_Chen_A = 0.0, LQ_Chen_B = 0.0, LQ_Chen_C = 0.0;
		static double dfchem_dphi_LQ_Chen(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double phi = phase.phi, sum_phi_j = 0.0;
			for (auto p = node.begin(); p < node.end(); p++)
				if (p->index != phase.index)
					sum_phi_j += p->phi * p->phi;
			return -LQ_Chen_A * phi + LQ_Chen_B * phi * phi * phi + 2.0 * LQ_Chen_C * phi * sum_phi_j;
		}
		//-----------------------------------------------------------------------------------------------
		// const : potential = con
		static void dfchem_dphase_con_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			for (auto p = phase.potential.begin(); p < phase.potential.end(); p++)
				p->value += phase.x[p->index].value;
		}
		//-----------------------------------------------------------------------------------------------
		// fchem for total con
		static double dfchem_dcon_i_const(pf::PhaseNode& node, int con_i) {
			return node.x[con_i].value;
		}
		//-----------------------------------------------------------------------------------------------
		// dx_i^a / du_i
		static double dphase_con_i_du_i_const(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i) {
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						return 1.0;
				return 0.0;
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						return 0.0;
				return 1.0;
			}
			return 0.0;
		}
		//-----------------------------------------------------------------------------------------------
		// x_i^a = func( u_i )
		static void phase_con_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++) {
					if (phase.index == *index) {
						for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
							comp->value = node.potential[comp->index].value;
						return;
					}
				}
				for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
					comp->value = 0.0;
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
					comp->value = node.potential[comp->index].value;
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
							comp->value = 0.0;
			}
		}

		static double (*fchem_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*dfchem_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of this model

		static void (*dfchem_dphase_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of phase con

		static double (*dfchem_dcon_i)(pf::PhaseNode& node, int con_i);  // main function of total con & grand potential

		static double (*dphase_con_i_du_i)(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i);  // grand potential

		static void (*phase_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // grand potential

		static void load_chemical_energy_model(bool infile_debug) {
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				int model_type = 0;
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.PhiCon.BulkEnergy.type : 1 - DoubleWell, 2 - LQ_Chen, 3 - H_Liang , 7 - HighOrder, 8 - SimpleDoubleWell\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.PhiCon.BulkEnergy.type", model_type, infile_debug);
				switch (model_type)
				{
				case pf::dfdphi_DOUBEL_WELL:
					dfchem_dphi = dfchem_dphi_double_well;
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.DoubleWell.A", DOUBLE_WELL_A, infile_debug);
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.DoubleWell.B", DOUBLE_WELL_B, infile_debug);
					break;
				case pf::dfdphi_LQ_CHEN:
					dfchem_dphi = dfchem_dphi_LQ_Chen;
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.LQ_Chen.A", LQ_Chen_A, infile_debug);
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.LQ_Chen.B", LQ_Chen_B, infile_debug);
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.LQ_Chen.C", LQ_Chen_C, infile_debug);
					break;
				case pf::dfdphi_HighOrder:
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# No models for HighOrder with standard AC and CH equations\n", InputFileReader::get_instance()->debug_file);
					std::exit(0);
					break;
				case pf::dfdphi_DOUBLE_WELL_SIMPLE:
					dfchem_dphi = dfchem_dphi_double_well;
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.SimpleDoubleWell.A", DOUBLE_WELL_A, infile_debug);
				default:
					break;
				}
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise || Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				int model_type = 0;
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.PhiCon.BulkEnergy.type : 7 - HighOrder\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.PhiCon.BulkEnergy.type", model_type, infile_debug);
				switch (model_type)
				{
				case pf::PEType_Const:
					dfchem_dphi = dfchem_dphi_const;
					break;
				case pf::dfdcon_1_CRACK_WELL:
					break;
				case pf::dfdcon_1_CRACK_OBSTACLE:
					break;
				case pf::dfdcon_N_CRACK_WELL:
					break;
				case pf::dfdcon_N_CRACK_OBSTACLE:
					break;
				case pf::dfdphi_HighOrder:
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# No models for HighOrder with pair-wise AC and grand-potential equations\n", InputFileReader::get_instance()->debug_file);
					std::exit(0);
					break;
				default:
					break;
				}
			}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			fchem_density = dfchem_dphi_const;
			dfchem_dphi = dfchem_dphi_const;
			dfchem_dphase_con = dfchem_dphase_con_const;
			dfchem_dcon_i = dfchem_dcon_i_const;
			dphase_con_i_du_i = dphase_con_i_du_i_const;
			phase_con = phase_con_const;
			smooth_phases = Solvers::get_instance()->C_Solver.phase_indexes;
			_domain = Solvers::get_instance()->parameters.ConEDomain;
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard)
				dfchem_dphi = dfchem_dphi_double_well;

			load_chemical_energy_model(infile_debug);

		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}