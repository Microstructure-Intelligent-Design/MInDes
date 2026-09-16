#pragma once
#include <cmath>
#include <cstdlib>
#include "../../Module.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "DS_Params.h"
#include "DS_Functions.h"

namespace pf {
	namespace dendrite_solidification_model {
		inline void init_model_modules() {
			if (main_field::phi_number != 1 || main_field::con_number != 0 ||
				!main_field::is_temp_field_on || mesh_parameters::MESH_NZ != 1) {
				WriteLog("> ERROR, Dendrite Solidification requires Nx,Ny >= 3, Nz = 1 and PCT = (1,0,true).\n");
				std::exit(EXIT_FAILURE);
			}
			WriteLog("> Simulation Model - Dendrite Solidification - is Activated !\n");
			WriteDebugFile("# DS: phi_t = (d_y(epsilon*epsilon_theta*phi_x)-d_x(epsilon*epsilon_theta*phi_y)+epsilon^2*lap(phi)+phi*(1-phi)*(phi-0.5+m))/tau\n");
			WriteDebugFile("# DS: T_t = lap(T)+kappa*phi_t; m = alpha/pi*atan(gamma*(teq-T)).\n");
			WriteDebugFile("# DS follows the MATLAB epsilon^2*lap(phi), rather than the prose divergence form.\n");
			// - 
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.tau", parameters::tau, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.epsilonb", parameters::epsilonb, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.kappa", parameters::kappa, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.delta", parameters::delta, true);
			infile_reader::read_int_value("SimulationModels.DendriteSolidification.aniso", parameters::aniso, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.alpha", parameters::alpha, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.gamma", parameters::gamma, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.teq", parameters::teq, true);
			infile_reader::read_real_value("SimulationModels.DendriteSolidification.theta0", parameters::theta0, true);
			if (!std::isfinite(parameters::tau) || parameters::tau <= 0 ||
				!std::isfinite(parameters::epsilonb) || parameters::epsilonb <= 0 ||
				!std::isfinite(parameters::kappa) || parameters::kappa < 0 ||
				!std::isfinite(parameters::delta) || std::abs(parameters::delta) >= 1 ||
				parameters::aniso <= 0 ||
				!std::isfinite(parameters::alpha) || parameters::alpha < 0 || parameters::alpha >= 1 ||
				!std::isfinite(parameters::gamma) || parameters::gamma < 0 ||
				!std::isfinite(parameters::teq) || !std::isfinite(parameters::theta0) ||
				!std::isfinite(mesh_parameters::delt_r) || mesh_parameters::delt_r <= 0 ||
				!std::isfinite(time_parameters::delt_t) || time_parameters::delt_t <= 0) {
				WriteLog("> ERROR, invalid Dendrite Solidification parameters, grid spacing or time step.\n");
				std::exit(EXIT_FAILURE);
			}
			load_a_new_module(nullptr, nullptr, exec_pre_iii,
				exec_i, nullptr, nullptr, nullptr, nullptr, nullptr, deinit);
		}
	}
}
