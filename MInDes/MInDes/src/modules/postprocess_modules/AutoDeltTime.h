#pragma once
#include "../Module.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../model_modules/Model_Params.h"
#include "ShowLoopInfo.h"
#include "../../MainIterator_Params.h"
namespace pf {
	namespace automatic_change_delt_time {
		inline REAL time_interval = 1;
		inline REAL dt_scale = 1;
		inline size_t delt_step = 0;
		inline REAL max_scale = REAL(1e3);
		inline REAL small_scale = REAL(1e-6);
		inline bool is_reduce_output = true;
		inline REAL phi_increment_limit = REAL(1e-3);
		inline REAL con_increment_limit = REAL(1e-3);
		inline REAL temp_increment_limit = REAL(1e-3);
		void exec_post_iii() {
			REAL MAX_dphi_versus_limitPhi = main_field::PHI_MAX_VARIATION / phi_increment_limit;
			REAL MAX_dcon_versus_limitCon = main_field::CON_MAX_VARIATION / con_increment_limit;
			REAL MAX_dT_versus_limitT = main_field::TEMP_MAX_VARIATION / temp_increment_limit;

			if (MAX_dcon_versus_limitCon > 1) {
				dt_scale /= MAX_dcon_versus_limitCon + REAL(0.05);
				MAX_dcon_versus_limitCon = 1;
				if (dt_scale < small_scale)
					dt_scale = small_scale;
				if (!is_reduce_output) {
					string log = "> Adjust time interval for concentration evolving stability, dt_scale = " + to_string(dt_scale) + "\n";
					WriteLog(log);
				}
			}
			else if (MAX_dphi_versus_limitPhi > 1) {
				dt_scale /= MAX_dphi_versus_limitPhi + REAL(0.05);
				MAX_dphi_versus_limitPhi = 1;
				if (dt_scale < small_scale)
					dt_scale = small_scale;
				if (!is_reduce_output) {
					string log = "> Adjust time interval for phi evolving stability, dt_scale = " + to_string(dt_scale) + "\n";
					WriteLog(log);
				}
			}
			else if (MAX_dT_versus_limitT > 1) {
				dt_scale /= MAX_dT_versus_limitT + REAL(0.05);
				MAX_dT_versus_limitT = 1.0;
				if (dt_scale < small_scale)
					dt_scale = small_scale;
				if (!is_reduce_output) {
					string log = "> Adjust time interval for temperature evolving stability, dt_scale = " + to_string(dt_scale) + "\n";
					WriteLog(log);
				}
			}
			else if (main_iterator::Current_ITE_step % delt_step == 0 && (dt_scale * REAL(1.05)) <= max_scale) {
				dt_scale *= REAL(1.05);
				if (!is_reduce_output) {
					string log = "> Adjust time interval for fields evolving quickly, dt_scale = " + to_string(dt_scale) + "\n";
					WriteLog(log);
				}
			}
			else if (main_iterator::Current_ITE_step % delt_step == 0 && (dt_scale * REAL(1.05)) > max_scale) {
				dt_scale = max_scale;
				if (!is_reduce_output) {
					string log = "> Adjust time interval for fields evolving quickly, dt_scale = " + to_string(dt_scale) + "\n";
					WriteLog(log);
				}
			}
			time_parameters::delt_t = time_interval * dt_scale;

			if (main_iterator::Current_ITE_step % show_loop_information::screen_output_step == 0) {
				string log = "# Time interval is automatically adjusting, time interval = " + to_string(time_parameters::delt_t) + "\n";
				WriteLog(log);
			}
		}
		void init_auto_time() {
			infile_reader::read_int_value("Postprocess.PCT.AutoTimeInterval.frequence", delt_step, true);
			if (delt_step > 0) {
				time_interval = time_parameters::delt_t;
				InputFileReader::get_instance()->read_REAL_value("Postprocess.PCT.AutoTimeInterval.large_scale", max_scale, true);
				InputFileReader::get_instance()->read_REAL_value("Postprocess.PCT.AutoTimeInterval.small_scale", small_scale, true);
				InputFileReader::get_instance()->read_bool_value("Postprocess.PCT.AutoTimeInterval.is_reduce_output", is_reduce_output, true);
				if (main_field::is_phi_field_on)
					InputFileReader::get_instance()->read_REAL_value("Postprocess.PCT.AutoTimeInterval.phi_increment_limit", phi_increment_limit, true);
				if (main_field::is_con_field_on)
					InputFileReader::get_instance()->read_REAL_value("Postprocess.PCT.AutoTimeInterval.con_increment_limit", con_increment_limit, true);
				if (main_field::is_temp_field_on)
					InputFileReader::get_instance()->read_REAL_value("Postprocess.PCT.AutoTimeInterval.temp_increment_limit", temp_increment_limit, true);
				load_a_new_module(default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, exec_post_iii, default_module_function);
				WriteLog("> MODULE INIT : Time interval for Phi, Con and Temp will automatically change during simulation !\n");
			}
		}

	}
}