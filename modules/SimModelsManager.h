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
#include "Base.h"
#include "sim_models/InterfaceEnergy.h"
#include "sim_models/Mobility.h"
#include "sim_models/BulkEnergy.h"
#include "sim_models/Kinetics.h"
#include "sim_models/Source.h"

namespace pf {
	enum PrefabricatedCases {
		SpinodalDecomposition, GrainGrowth, Solidification, Sintering, CellMovement, 
		Precipitation, Rafting, ElectroDeposition, Intercalation,
		CrackPropagation, Bubble, Droplet,
	};
	namespace sim_models_manager {

		static void load_phi_evolution_equation(PhiEquationType _type) {
			Solvers::get_instance()->parameters.PhiEType = _type;
		}

		static void load_con_evolution_equation(ConEquationType _type, ConEquationDomain _domain_type) {
			Solvers::get_instance()->parameters.ConEType = _type;
			Solvers::get_instance()->parameters.ConEDomain = _domain_type;
		}

		static void load_temp_evolution_equation(TemperatureEquationType _type) {
			Solvers::get_instance()->parameters.TempEType = _type;
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			int prefabricated_case = 0;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false); 
			/*if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.prefabricated_case : 0 - SpinodalDecomposition , 1 - GrainGrowth , 2 - Solidification , 3 - Sintering , 4 - CellMovement\n", InputFileReader::get_instance()->debug_file);
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                    5 - Precipitation , 6 - Rafting , 7 - ElectroDeposition , 8 - Intercalation , 9 - CrackPropagation\n", InputFileReader::get_instance()->debug_file);
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                                    10 - BubbleDynamics \n", InputFileReader::get_instance()->debug_file);
			if (InputFileReader::get_instance()->read_int_value("ModelsManager.prefabricated_case", prefabricated_case, infile_debug)) {
				switch (prefabricated_case)
				{
				case pf::SpinodalDecomposition:
					break;
				case pf::GrainGrowth:
					break;
				case pf::Solidification:
					break;
				case pf::Sintering:
					break;
				case pf::CellMovement:
					break;
				case pf::Precipitation:
					break;
				case pf::Rafting:
					break;
				case pf::ElectroDeposition:
					break;
				case pf::Intercalation:
					break;
				case pf::CrackPropagation:
					break;
				case pf::Bubble:
					break;
				case pf::Droplet:
					break;
				default:
					break;
				}
			}*/
			// load evolution equation
			int _PhiEType = Solvers::get_instance()->parameters.PhiEType, _ConEType = Solvers::get_instance()->parameters.ConEType,
				_ConEDomain = Solvers::get_instance()->parameters.ConEDomain, _TEType = Solvers::get_instance()->parameters.TempEType;

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.equation : 0 - Const, 1 - AllenCahn Standard, 2 - AllenCahn Pairwise, 3 - CahnHilliard Standard\n", InputFileReader::get_instance()->debug_file);
			if (InputFileReader::get_instance()->read_int_value("ModelsManager.Phi.equation", _PhiEType, infile_debug))
				Solvers::get_instance()->parameters.PhiEType = PhiEquationType(_PhiEType);
			switch (PhiEquationType(_PhiEType))
			{
			case pf::PEType_AC_Standard:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->Phi_Solver_AC.print_phi_model_normal(), LOG_FILE_NAME);
				break;
			case pf::PEType_AC_Pairwise:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->Phi_Solver_AC.print_phi_model_pairwise(), LOG_FILE_NAME);
				break;
			case pf::PEType_CH_Standard:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->Phi_Solver_CH.print_phi_model_cahnhilliard(), LOG_FILE_NAME);
				break;
			default:
				break;
			}

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Con.equation : 0 - Const, 1 - TotalConcentration, 2 - PhaseConcentration, 3 - GrandPotential\n", InputFileReader::get_instance()->debug_file);
			if (InputFileReader::get_instance()->read_int_value("ModelsManager.Con.equation", _ConEType, infile_debug))
				Solvers::get_instance()->parameters.ConEType = ConEquationType(_ConEType);
			switch (ConEquationType(_ConEType))
			{
			case pf::CEType_TotalX:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->C_Solver.print_total_con_model(), LOG_FILE_NAME);
				break;
			case pf::CEType_PhaseX:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->C_Solver.print_phase_con_model(), LOG_FILE_NAME);
				break;
			case pf::CEType_GrandP:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->C_Solver.print_grand_potential_model(), LOG_FILE_NAME);
				break;
			default:
				break;
			}

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Con.valid_domain : 0 - Standard, 1 - Reverse\n", InputFileReader::get_instance()->debug_file);
			if (InputFileReader::get_instance()->read_int_value("ModelsManager.Con.valid_domain", _ConEDomain, infile_debug))
				Solvers::get_instance()->parameters.ConEDomain = ConEquationDomain(_ConEDomain);

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Temp.equation : 0 - Const, 1 - Standard\n", InputFileReader::get_instance()->debug_file);
			if (InputFileReader::get_instance()->read_int_value("ModelsManager.Temp.equation", _TEType, infile_debug))
				Solvers::get_instance()->parameters.TempEType = TemperatureEquationType(_TEType);
			switch (TemperatureEquationType(_TEType))
			{
			case pf::TType_Standard:
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(Solvers::get_instance()->T_Solver.print_temperature_model(), LOG_FILE_NAME);
				break;
			default:
				break;
			}

			if (_ConEType == ConEquationType::CEType_TotalX || _ConEType == ConEquationType::CEType_GrandP) {
				string phase_indexes_key = "ModelsManager.Con.ValidDomain.phase_indexes", phase_indexes_input = "()";
				InputFileReader::get_instance()->read_string_value(phase_indexes_key, phase_indexes_input, infile_debug);
				vector<input_value> phase_indexes_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, phase_indexes_key, phase_indexes_input, infile_debug);
				for (int index = 0; index < phase_indexes_value.size(); index++)
					Solvers::get_instance()->C_Solver.phase_indexes.push_back(phase_indexes_value[index].int_value);

			}
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				string _solvent = "";
				if (InputFileReader::get_instance()->read_string_value("ModelsManager.Con.solvent", _solvent, infile_debug)) {
					bool is_setting = false;
					for (auto comp = Solvers::get_instance()->parameters.Components.begin();
						comp < Solvers::get_instance()->parameters.Components.end(); comp++)
						if (comp->name.compare(_solvent) == 0) {
							is_setting = true;
							Solvers::get_instance()->C_Solver.solvent = comp->index;
						}
					if (!is_setting) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : This solvent hasn't been defined ! \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
				}
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				string _solvent = "";
				if (InputFileReader::get_instance()->read_string_value("ModelsManager.Con.solvent", _solvent, infile_debug)) {
					bool is_setting = false;
					for (auto comp = Solvers::get_instance()->parameters.Components.begin();
						comp < Solvers::get_instance()->parameters.Components.end(); comp++)
						if (comp->name.compare(_solvent) == 0) {
							is_setting = true;
							Solvers::get_instance()->C_Solver.solvent = comp->index;
						}
					if (!is_setting) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : This solvent hasn't been defined ! \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
				}
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ValidDomain.threshold", Solvers::get_instance()->C_Solver.threshold, infile_debug);
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {

				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ValidDomain.threshold", Solvers::get_instance()->C_Solver.threshold, infile_debug);

				string grand_range_key = "ModelsManager.Con.GrandPotential.range", grand_range_input = "(-100000.0,100000.0)";
				InputFileReader::get_instance()->read_string_value(grand_range_key, grand_range_input, infile_debug);
				vector<input_value> grand_range_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, grand_range_key, grand_range_input, infile_debug);
				Solvers::get_instance()->C_Solver.grand_potential_range[0] = grand_range_value[0].double_value - SYS_EPSILON;
				Solvers::get_instance()->C_Solver.grand_potential_range[1] = grand_range_value[1].double_value + SYS_EPSILON;
			}

			// custom models
			pf::interface_mobility::init(phaseMesh);
			pf::interface_energy::init(phaseMesh);
			pf::bulk_energy::init(phaseMesh);
			pf::kinetics::init(phaseMesh);
			pf::extern_source::init(phaseMesh);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			pf::interface_mobility::exec_pre(phaseMesh);
			pf::interface_energy::exec_pre(phaseMesh);
			pf::extern_source::exec_pre(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			report += pf::interface_mobility::exec_loop(phaseMesh);
			report += pf::interface_energy::exec_loop(phaseMesh);
			report += pf::extern_source::exec_loop(phaseMesh);
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			pf::interface_mobility::deinit(phaseMesh);
			pf::interface_energy::deinit(phaseMesh);
			pf::bulk_energy::deinit(phaseMesh);
			pf::kinetics::deinit(phaseMesh);
			pf::extern_source::deinit(phaseMesh);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			pf::interface_mobility::write_scalar(fout, phaseMesh);
			pf::interface_energy::write_scalar(fout, phaseMesh);
			pf::bulk_energy::write_scalar(fout, phaseMesh);
			pf::extern_source::write_scalar(fout, phaseMesh);
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			pf::interface_mobility::write_vec3(fout, phaseMesh);
			pf::interface_energy::write_vec3(fout, phaseMesh);
			pf::bulk_energy::write_vec3(fout, phaseMesh);
			pf::extern_source::write_vec3(fout, phaseMesh);
		}
		static void load_module() {
			Solvers::get_instance()->create_a_new_module(init, exec_pre, exec_loop, deinit, write_scalar, write_vec3);
		};
	}
}