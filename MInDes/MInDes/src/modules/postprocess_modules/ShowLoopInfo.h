#pragma once
#include <climits>
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
#include "../../MainIterator_Params.h"
#include "../model_modules/Model_Params.h"
#include "../base/timer.h"
#include "../postprocess_modules/DataStatistics.h"
namespace pf {
	namespace show_loop_information {
		inline size_t screen_loop_step   = 0;
		inline size_t screen_output_step = 0;
		inline void exec_pre_iii() {
			std::vector<REAL> phi_info, con_info; REAL temp_info = 0;
			if (main_field::is_phi_field_on)
				phi_info = data_statistics_functions::statistical_phi();
			if (main_field::is_con_field_on)
				con_info = data_statistics_functions::statistical_con();
			if (main_field::is_temp_field_on)
				temp_info = data_statistics_functions::statistical_temp();
			stringstream log;
			log << "> PHI CON TEMP information :" << endl;
			if (main_field::is_phi_field_on)
				for (int index = 0; index < main_field::phi_number; index++)
					log << "#  Phi " << index << " = "
					<< setprecision(5) << phi_info[index] << endl;
			if (main_field::is_con_field_on)
				for (int index = 0; index < main_field::con_number; index++)
					log << "#  Con " << index << " = "
					<< setprecision(5) << con_info[index] << endl;
			if (main_field::is_temp_field_on)
				log << "#  Temp " << " = "
				<< setprecision(5) << temp_info << endl;
			WriteLog(log.str());
		}

		inline void exec_pos_i() {
			if (screen_output_step == 0)
				return;
			if (main_iterator::Current_ITE_step % screen_output_step == 0) {
				std::vector<REAL> phi_info, con_info; REAL temp_info = 0;
				if (main_field::is_phi_field_on)
					phi_info = data_statistics_functions::statistical_phi();
				if (main_field::is_con_field_on)
					con_info = data_statistics_functions::statistical_con();
				if (main_field::is_temp_field_on)
					temp_info = data_statistics_functions::statistical_temp();
				stringstream log;
				log << "#====================================================================================================" << endl;
				log << timer::return_cunrrent_time_by_string();
				log.setf(ios::fixed);
				log << "# Simulation step " << main_iterator::Current_ITE_step << " has been finished!" << endl;
				log << "# This " << screen_output_step << " steps used " << setprecision(3) << timer::interval_end(main_iterator::t_interval_begin) << "(secs.), " << endl;
				log << "# Total " << main_iterator::Current_ITE_step << " steps used " << setprecision(3) << timer::total_duration_sec(main_iterator::t_total_begin) << "(secs.)." << endl;
				log << "# Real simulation time is " << setprecision(3) << time_parameters::Real_Time << " (secs.)" << endl;
				log << "#----------------------------------------------------------------------------------------------------" << endl;
				if (main_field::is_phi_field_on)
					for (int index = 0; index < main_field::phi_number; index++)
						log << "#  Phi " << index << " = "
						<< setprecision(5) << phi_info[index] << endl;
				if (main_field::is_con_field_on)
					for (int index = 0; index < main_field::con_number; index++)
						log << "#  Con " << index << " = "
						<< setprecision(5) << con_info[index] << endl;
				if (main_field::is_temp_field_on)
					log << "#  Temp " << " = "
					<< setprecision(5) << temp_info << endl;
				WriteLog(log.str());
				timer::interval_begin(main_iterator::t_interval_begin);
			}
		}

		inline void exec_pos_ii() {
			bool output = false;
			if (screen_loop_step != 0)
				if (main_iterator::Current_ITE_step % screen_loop_step == 0)
					output = true;
			if (screen_output_step != 0)
				if (main_iterator::Current_ITE_step % screen_output_step == 0)
					output = true;
			if (output) {
				stringstream log;
				log << "#------------------------------------------ PCT Field -----------------------------------------------" << endl;
				log << "# CURRENT STEP = " << main_iterator::Current_ITE_step << ", REAL TIME = " << time_parameters::Real_Time << endl;
				if (main_field::is_phi_field_on)
					log << "# MAX PHI  INCREMENT = " << setprecision(5) << main_field::PHI_MAX_VARIATION << endl;
				if (main_field::is_con_field_on)
					log << "# MAX CON  INCREMENT = " << setprecision(5) << main_field::CON_MAX_VARIATION << endl;
				if (main_field::is_temp_field_on)
					log << "# MAX TEMP INCREMENT = " << setprecision(5) << main_field::TEMP_MAX_VARIATION << endl;
				log << "#----------------------------------------------------------------------------------------------------" << endl;
				if (screen_output_step != 0)
					if (main_iterator::Current_ITE_step % screen_output_step != 0)
						log << endl << endl;
				WriteLog(log.str());
			}
		}

		// info end
		inline void exec_pos_iii() {
			if (screen_output_step == 0)
				return;
			if (main_iterator::Current_ITE_step % screen_output_step == 0) {
				stringstream log;
				log << "#====================================================================================================" << endl;
				log << endl << endl;
				WriteLog(log.str());
			}
		}

		inline void init() {
			infile_reader::read_int_value("Solver.Output.LOG.loop_info_step", screen_loop_step, true);
			infile_reader::read_int_value("Solver.Output.LOG.screen_output_step", screen_output_step, true);
			load_a_new_module(default_module_function, default_module_function, exec_pre_iii,
				default_module_function, default_module_function, default_module_function,
				exec_pos_i, exec_pos_ii, exec_pos_iii, default_module_function);
		}
	}
}