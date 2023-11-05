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
#include "../Base.h"
#include "BulkEnergy/ChemicalEnergy.h"
#include "BulkEnergy/MechanicalEnergy.h"
#include "BulkEnergy/ElectricEnergy.h"
#include "BulkEnergy/MagneticEnergy.h"

namespace pf {
	namespace bulk_energy {
		// output
		int_box bulk_energy_standard_output;
		PairFlag bulk_energy_pairwise_output;
		static bool is_energy_density_merge_output = false;
		static bool is_energy_density_output = false;
		// model
		static double_box const_fbulk; // (phi.property, value)
		static double dfbulk_dphi_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return const_fbulk[phase.property];
		}
		// energy density
		static double (*fchem_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*fmech_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*felec_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*fmag_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		// phi driving force
		static double (*dfchem_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*dfmech_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*dfelec_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static double (*dfmag_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		// phase con potential
		static void (*dfchem_dcon)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static void (*dfmech_dcon)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static void (*dfelec_dcon)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		static void (*dfmag_dcon)(pf::PhaseNode& node, pf::PhaseEntry& phase);
		// total con & grand potential
		static double (*dfchem_dconi)(pf::PhaseNode& node, int con_i);
		static double (*dfmech_dconi)(pf::PhaseNode& node, int con_i);
		static double (*dfelec_dconi)(pf::PhaseNode& node, int con_i);
		static double (*dfmag_dconi)(pf::PhaseNode& node, int con_i);
		// + grand potential
		static double (*dPhiCon_du)(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i);
		static void (*Phi_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double bulk_energy_density(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return dfbulk_dphi_const(node, phase)
				+ fchem_density(node, phase) + fmech_density(node, phase)
				+ felec_density(node, phase) + fmag_density(node, phase);
		}

		static double dfbulk_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return dfbulk_dphi_const(node, phase)
				+ dfchem_dphi(node, phase) + dfmech_dphi(node, phase)
				+ dfelec_dphi(node, phase) + dfmag_dphi(node, phase);
		}

		static void dfbulk_dcon(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			dfchem_dcon(node, phase);
			dfmech_dcon(node, phase);
			dfelec_dcon(node, phase);
			dfmag_dcon(node, phase);
		}

		static double dfbulk_dconi(pf::PhaseNode& node, int con_i) {
			return dfchem_dconi(node, con_i) +
				dfmech_dconi(node, con_i) +
				dfelec_dconi(node, con_i) +
				dfmag_dconi(node, con_i);
		}

		static double dphase_con_du(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i) {
			return dPhiCon_du(node, phase, con_i);
		}

		static void phase_con(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			Phi_con(node, phase);
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise ||
				Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard)
				Solvers::get_instance()->Phi_Solver_AC.dfbulk_dphi = dfbulk_dphi;
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard)
				Solvers::get_instance()->Phi_Solver_CH.dfbulk_dphi = dfbulk_dphi;

			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				Solvers::get_instance()->C_Solver.dfbulk_dx = dfbulk_dcon;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				Solvers::get_instance()->C_Solver.df_dx = dfbulk_dconi;
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				Solvers::get_instance()->C_Solver.dfbulk_dx = dfbulk_dcon;
				Solvers::get_instance()->C_Solver.dphase_x_du = dphase_con_du;
				Solvers::get_instance()->C_Solver.phase_x = phase_con;
			}
			// bulk const
			for (auto phi_property = Solvers::get_instance()->parameters.Phases.begin();
				phi_property < Solvers::get_instance()->parameters.Phases.end(); phi_property++) {
				const_fbulk.add_double(phi_property->phi_property, 0.0);
			}
			if (Solvers::get_instance()->parameters.PhiEType != PhiEquationType::PEType_Const) {
				string bulk_energy_const_key = "ModelsManager.Phi.BulkEnergy.const", bulk_energy_const_input = "[()]";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.BulkEnergy.const = [(phi_name, bulk_energy), ... ]\n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value(bulk_energy_const_key, bulk_energy_const_input, infile_debug)) {
					// BulkEnergy.const = [(phi_name,value), ... ]
					vector<InputValueType> bulk_energy_const_structure; bulk_energy_const_structure.push_back(InputValueType::IVType_STRING); bulk_energy_const_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> bulk_energy_const_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(bulk_energy_const_structure, bulk_energy_const_key, bulk_energy_const_input, infile_debug);
					for (int index = 0; index < bulk_energy_const_value.size(); index++)
						const_fbulk.add_double(Solvers::get_instance()->parameters.Phases[bulk_energy_const_value[index][0].string_value].phi_property, bulk_energy_const_value[index][1].double_value);
				}
			}

			chemical_energy::init(phaseMesh);
			dfchem_dphi = chemical_energy::dfchem_dphi;
			dfchem_dcon = chemical_energy::dfchem_dphase_con;
			dfchem_dconi = chemical_energy::dfchem_dcon_i;
			dPhiCon_du = chemical_energy::dphase_con_i_du_i;
			Phi_con = chemical_energy::phase_con;
			fchem_density = chemical_energy::fchem_density;

			mechanical_energy::init(phaseMesh);
			dfmech_dphi = mechanical_energy::dfmech_dphi;
			dfmech_dcon = mechanical_energy::dfmech_dphase_con;
			dfmech_dconi = mechanical_energy::dfmech_dcon_i;
			fmech_density = mechanical_energy::fmech_density;

			electric_energy::init(phaseMesh);
			dfelec_dphi = electric_energy::dfelec_dphi;
			dfelec_dcon = electric_energy::dfelec_dphase_con;
			dfelec_dconi = electric_energy::dfelec_dcon_i;
			felec_density = electric_energy::felec_density;

			magnetic_energy::init(phaseMesh);
			dfmag_dphi = magnetic_energy::dfmag_dphi;
			dfmag_dcon = magnetic_energy::dfmag_dphase_con;
			dfmag_dconi = magnetic_energy::dfmag_dcon_i;
			fmag_density = magnetic_energy::fmag_density;

			// output
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				string df_dphi_key = "ModelsManager.Phi.BulkEnergy.vts_output", df_dphi_input = "[()]";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.BulkEnergy.vts_output = [(phi_index_0, phi_index_1), ... ]\n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value(df_dphi_key, df_dphi_input, infile_debug)) {
					vector<vector<input_value>> df_dphi_value;
					df_dphi_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, df_dphi_key, df_dphi_input, infile_debug);
					for (int i = 0; i < df_dphi_value.size(); i++)
						bulk_energy_pairwise_output.set(df_dphi_value[i][0].int_value, df_dphi_value[i][1].int_value, 1);
				}
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard ||
				Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				string df_dphi_key = "ModelsManager.Phi.BulkEnergy.vts_output", df_dphi_input = "()";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.BulkEnergy.vts_output = (phi_index_0, phi_index_1, ...) \n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value(df_dphi_key, df_dphi_input, infile_debug)) {
					vector<input_value> df_dphi_value;
					df_dphi_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, df_dphi_key, df_dphi_input, infile_debug);
					for (int i = 0; i < df_dphi_value.size(); i++)
						bulk_energy_standard_output.add_int(df_dphi_value[i].int_value, 1);
				}
			}
			//InputFileReader::get_instance()->read_bool_value("ModelsManager.Phi.BulkEnergyDensityMerge.vts_output", is_energy_density_merge_output, infile_debug);
			//InputFileReader::get_instance()->read_bool_value("ModelsManager.Phi.BulkEnergyDensity.vts_output", is_energy_density_output, infile_debug);
			
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : BulkEnergy !\n", LOG_FILE_NAME);
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			dfchem_dphi = nullptr;
			dfmech_dphi = nullptr;
			dfelec_dphi = nullptr;
			dfmag_dphi = nullptr;
			chemical_energy::deinit(phaseMesh);
			mechanical_energy::deinit(phaseMesh);
			electric_energy::deinit(phaseMesh);
			magnetic_energy::deinit(phaseMesh);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (bulk_energy_pairwise_output.pairValue.size() != 0 && Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "dfbulkdphi_merge" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double val = 0.0, merge = 0.0;
							for (auto f = bulk_energy_pairwise_output.begin(); f < bulk_energy_pairwise_output.end(); f++) {
								if (node[f->pairIndex_1]._flag == pf_INTERFACE && node[f->pairIndex_2]._flag == pf_INTERFACE) {
									val += (dfbulk_dphi(node, node[f->pairIndex_2]) - dfbulk_dphi(node, node[f->pairIndex_1]))
										* (node[f->pairIndex_1].phi + node[f->pairIndex_2].phi);
									merge += node[f->pairIndex_1].phi + node[f->pairIndex_2].phi;
								}
							}
							if (merge > SYS_EPSILON)
								fout << val / merge << endl;
							else
								fout << 0.0 << endl;
						}
				fout << "</DataArray>" << endl;
			}
			if (bulk_energy_standard_output.size() != 0 && (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard ||
				Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard)) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "dfbulkdphi_merge" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < phaseMesh.limit_z; k++)
					for (int j = 0; j < phaseMesh.limit_y; j++)
						for (int i = 0; i < phaseMesh.limit_x; i++) {
							PhaseNode& node = phaseMesh(i, j, k);
							double val = 0.0, merge = 0.0;
							for (auto f = bulk_energy_standard_output.begin(); f < bulk_energy_standard_output.end(); f++) {
								if (abs(node[f->index].laplacian) > SYS_EPSILON) {
									val += dfbulk_dphi(node, node[f->index]) * node[f->index].phi;
									merge += node[f->index].phi;
								}
							}
							if (merge > SYS_EPSILON)
								fout << val / merge << endl;
							else
								fout << 0.0 << endl;
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}