#pragma once
#include <climits>
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
#include "../../MainIterator_Params.h"
#include "../Modules_Params.h"
#include "../base/timer.h"
#include "DataStatistics.h"
namespace pf {
	namespace show_loop_information {
		inline size_t screen_loop_step   = 0;
		inline size_t screen_output_step = 0;
		inline std::string statistic_file_name = "Statistics.txt";

		std::vector<REAL> statistical_phi();

		std::vector<REAL> statistical_con();

		REAL statistical_temp();

		void exec_pre_iii();

		void exec_pos_i();

		void exec_pos_ii();

		// info end
		void exec_pos_iii();

		void init_show_loop_information();
	}
}