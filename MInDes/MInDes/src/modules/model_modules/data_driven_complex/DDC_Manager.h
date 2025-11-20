#pragma once
#include "../../Module.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "DDC_Params.h"
#include "DDC_Functions.h"
namespace pf {
	namespace data_driven_complex_model {
		/*
			模型来自文献：
			  原创
		*/
		inline void init_model_modules() {
			WriteLog("> \n");
			WriteLog("> Simulation Model - Data Driven Complex - is Activated ! \n");
			WriteLog(u8"> Reference : 原创 \n");
			WriteLog(u8"> DOI: 暂无 \n");
			WriteLog("> \n");
			WriteDebugFile("# Funcs.phi : d phi_i / d t = XXX \n");
			WriteDebugFile("# Funcs.con : d c_i / d t     = XXX \n");
			// - 
			parameters::phase_field = &main_field::phase_field;
			parameters::concentration_field = &main_field::concentration_field;
			parameters::MESH_NX = mesh_parameters::MESH_NX;
			parameters::MESH_NY = mesh_parameters::MESH_NY;
			parameters::MESH_NZ = mesh_parameters::MESH_NZ;
			parameters::x_down = mesh_parameters::x_down;
			parameters::y_down = mesh_parameters::y_down;
			parameters::z_down = mesh_parameters::z_down;
			parameters::x_up = mesh_parameters::x_up;
			parameters::y_up = mesh_parameters::y_up;
			parameters::z_up = mesh_parameters::z_up;
			parameters::delt_r = mesh_parameters::delt_r;
			parameters::delt_t = &time_parameters::delt_t;
			parameters::phi_number = main_field::phi_number;
			parameters::PHI_MAX_VARIATION = &main_field::PHI_MAX_VARIATION;
			parameters::con_number = main_field::con_number;
			parameters::CON_MAX_VARIATION = &main_field::CON_MAX_VARIATION;
			// - 
			infile_reader::read_real_value("SimulationModels.DataDrivenComplex.L", parameters::L, true);
			load_a_new_module(nullptr, nullptr, data_driven_complex_model::exec_pre_iii,  // exec_pre_i   exec_pre_ii    exec_pre_iii
				data_driven_complex_model::exec_i, nullptr, nullptr,  // exec_i   exec_ii    exec_iii
				nullptr, nullptr, nullptr,   // exec_pos_i   exec_pos_ii    exec_pos_iii
				data_driven_complex_model::deinit);  // deinit
		}
	}
}
