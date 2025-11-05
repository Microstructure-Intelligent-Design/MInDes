#pragma once
#include "../../Module.h"
#include "Model_Params.h"
#include "Model_Functions.h"
namespace pf {
	namespace grain_grows_spinodal_model {
		/*
			模型来自文献：
			郭灿, 赵玉平, 邓英远, 张忠明, 徐春杰, 运动晶界与调幅分解相互作用过程的相场法研究, 物理学报 Acta Phys. Sin. Vol.71, No.7 (2022) 078101
			DOI: 10.7498/aps.71.20211973
		*/
		inline void init_model_modules() {
			WriteLog("> \n");
			WriteLog("> Simulation Model - Grain Grows Spinodal - is Activated ! \n");
			WriteLog(u8"> Reference : 郭灿, 等, 运动晶界与调幅分解相互作用过程的相场法研究, 物理学报 Acta Phys. Sin. Vol.71, No.7 (2022) 078101 \n");
			WriteLog("> DOI: 10.7498/aps.71.20211973 \n");
			WriteLog("> \n");
			WriteDebugFile("# Funcs.phi : d eta_i / d t = - L * (m(c) * (eta_i^3 - eta_i + 2 * eta_i * epsilon * \\sum_{j != i}{eta_j^2}) - 2 * kappa_eta * \\nabla^2 eta_i) \n");
			WriteDebugFile("# Funcs.con : d c / d t     = \\nabla \\cdot M(eta) \\nabla (c - 5 * c * (1 - c) * (1 - 2 * c)) * f(eta) + 2 * A * c * (1 - c) * (1 - 2 * c) - 2 * kappa_con * \\nabla^2 c \n");
			WriteDebugFile("#             m(c)          = 1 + 0.5 * c^2 - 2.5 * c * (1 - c)^2 \n");
			WriteDebugFile("#             f(eta)        = 0.25 + \\sum_i {eta_i^4 / 4 - eta_i^2 / 2} + epsilon * \\sum_i \\sum_{j>i} {eta_i^2 * eta_j^2} \n");
			WriteDebugFile("#             M(eta)        = Mb + Mg * (\\sum_i \\sum_{j>i} {eta_i^2 * eta_j^2})^0.5 \n");
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
			load_a_new_module(default_module_function, default_module_function, grain_grows_spinodal_model::exec_pre_iii,  // exec_pre_i   exec_pre_ii    exec_pre_iii
				grain_grows_spinodal_model::exec_i, default_module_function, default_module_function,  // exec_i   exec_ii    exec_iii
				default_module_function, default_module_function, default_module_function,   // exec_pos_i   exec_pos_ii    exec_pos_iii
				grain_grows_spinodal_model::deinit);  // deinit
		}
	}
}
