#pragma once
#include "../../Module.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "DS_Params.h"
#include "DS_Functions.h"
namespace pf {
	namespace dendritic_solidification_model {
		/*
			模型来自文献：
			Kobayashi R (1993) Modeling and numerical simulations of dendritic crystal-growth. Physica D 63:410
			DOI: 10.1016/0167-2789(93)90120-P
		*/
		inline void init_model_modules() {
			// - check
			if (main_field::con_number != 1 || main_field::is_temp_field_on != false) {
				string error_report = "> ERROR, the PCT command line settings do not meet the requirements : PCT = (1,0,true) !\n";
				WriteLog(error_report);
				std::exit(0);
			}
			WriteLog("> \n");
			WriteLog("> Simulation Model - Dendritic Solidification - is Activated ! \n");
			WriteLog(u8"> Reference : Kobayashi R (1993) Modeling and numerical simulations of dendritic crystal-growth. Physica D 63:410 \n");
			WriteLog("> DOI: 10.1016/0167-2789(93)90120-P \n");
			WriteLog("> \n");
			WriteDebugFile("# Funcs.phi  : d eta / d t = - L * (m(c) * (eta_i^3 - eta_i + 2 * eta_i * epsilon * \\sum_{j != i}{eta_j^2}) - 2 * kappa_eta * \\nabla^2 eta_i) \n");
			WriteDebugFile("# Funcs.temp : d T / d t   = \\nabla \\cdot M(eta) \\nabla (c - 5 * c * (1 - c^2) * (3 * c^2 - 1)) * f(eta) + 2 * A * c * (1 - c) * (1 - 2 * c) - 2 * kappa_con * \\nabla^2 c \n");
			WriteDebugFile("#             m(c)          = 1 + 0.5 * c^2 - 2.5 * c^2 * (1 - c^2)^2 \n");
			WriteDebugFile("#             f(eta)        = 0.25 + \\sum_i {eta_i^4 / 4 - eta_i^2 / 2} + epsilon * \\sum_i \\sum_{j>i} {eta_i^2 * eta_j^2} \n");
			WriteDebugFile("#             M(eta)        = Mb + Mg * (\\sum_i \\sum_{j>i} {eta_i^2 * eta_j^2})^0.5 \n");
			// - 
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.L", parameters::L, true);
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.epsilon", parameters::epsilon, true);
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.kappa_eta", parameters::kappa_eta, true);
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.Mb", parameters::Mb, true);
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.Mg", parameters::Mg, true);
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.A", parameters::A, true);
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.kappa_con", parameters::kappa_con, true);
			if (infile_reader::read_int_value("SimulationModels.GrainGrowsSpinodal.noise_seed", parameters::noise_seed, true))
				parameters::is_noise_rand = false;
			infile_reader::read_real_value("SimulationModels.GrainGrowsSpinodal.noise_amplitude", parameters::noise_amplitude, true);
			load_a_new_module(nullptr, nullptr, dendritic_solidification_model::exec_pre_iii,  // exec_pre_i   exec_pre_ii    exec_pre_iii
				dendritic_solidification_model::exec_i, nullptr, nullptr,  // exec_i   exec_ii    exec_iii
				nullptr, nullptr, nullptr,   // exec_pos_i   exec_pos_ii    exec_pos_iii
				dendritic_solidification_model::deinit);  // deinit
		}
	}
}
