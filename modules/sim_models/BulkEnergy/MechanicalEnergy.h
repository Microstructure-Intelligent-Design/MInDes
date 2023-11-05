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
#include "../../sim_postprocess/Mechanics/ElasticSolver.h"

namespace pf {
	enum MechanicalEnergyType {
		MechType_CONST,   // no stain stress
		MechType_ELASTIC_STANDARD,  // standard elastic energy
		MechType_ELASTIC_PLASTIC_STANDARD,  // standard elastic + plastic energy
	};
	namespace mechanical_energy {
		static double dfmech_dphi_const(pf::PhaseNode& node, pf::PhaseEntry& phase){
			return 0.0;
		}
		static void dfmech_dphase_con_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static double dfmech_dcon_i_const(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		// get chemical energy density for each phase
		static double fmech_elastic(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.5 * (node.customVec6s[ExternalFields::MECH_stress] * (node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain]));
		}
		static double fmech_elastoplastic(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.5 * (node.customVec6s[ExternalFields::MECH_stress] * (node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain] - node.customVec6s[ExternalFields::MECH_plastic_strain])
				+ plastic_solver::get_hardening_modulus(node) * node.customValues[ExternalFields::MECH_ave_plastic_strain] * node.customValues[ExternalFields::MECH_ave_plastic_strain]);
		}
		static double dfmech_dphi_elastic(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			Vector6 strain = node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain];
			return 0.5 * (stiffness_eigenstrain::get_phi_stiffness(phase) * strain * strain);
		}
		static double dfmech_dphi_elastoplastic(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			Vector6 strain = node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain] - node.customVec6s[ExternalFields::MECH_plastic_strain];
			return 0.5 * (stiffness_eigenstrain::get_phi_stiffness(phase) * strain * strain
				+ plastic_solver::get_phase_hardening_modulus(phase) * node.customValues[ExternalFields::MECH_ave_plastic_strain] * node.customValues[ExternalFields::MECH_ave_plastic_strain]);
		}

		static double (*fmech_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*dfmech_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of this model

		static void (*dfmech_dphase_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of phase con

		static double (*dfmech_dcon_i)(pf::PhaseNode& node, int con_i);  // main function of total con & grand potential

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			fmech_density = dfmech_dphi_const;
			dfmech_dphi = dfmech_dphi_const;
			dfmech_dphase_con = dfmech_dphase_con_const;
			dfmech_dcon_i = dfmech_dcon_i_const;

			string fix_boundary_key = "Postprocess.physical_fields", fix_boundary_input = "(false,false,false)";
			InputFileReader::get_instance()->read_string_value(fix_boundary_key, fix_boundary_input, false);
			vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_BOOL, fix_boundary_key, fix_boundary_input, false);
			bool is_mechanical_field_on = fix_boundary_value[0].bool_value;
			if (is_mechanical_field_on) {
				bool is_fmech = false, is_plastic = false;
				InputFileReader::get_instance()->read_bool_value("ModelsManager.Phi.BulkEnergy.Fmech", is_fmech, infile_debug);
				if (is_fmech) {
					InputFileReader::get_instance()->read_bool_value("Postprocess.SolidMechanics.plasticity", is_plastic, false);
					if (is_plastic) {
						fmech_density = fmech_elastoplastic;
						dfmech_dphi = dfmech_dphi_elastoplastic;
					}
					else {
						fmech_density = fmech_elastic;
						dfmech_dphi = dfmech_dphi_elastic;
					}
				}
			}
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}