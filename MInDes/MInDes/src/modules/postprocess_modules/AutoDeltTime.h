#pragma once
#include "../base/MACRO_DEF.h"
#include "../Module.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Modules_Params.h"
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
		void exec_post_iii();
		void init_auto_time();
	}
}