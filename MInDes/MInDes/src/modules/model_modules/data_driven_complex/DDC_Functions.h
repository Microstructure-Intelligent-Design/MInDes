#pragma once
#include "DDC_Phi_Functions.h"
#include "DDC_Con_Functions.h"
#include "DDC_Temp_Functions.h"
#include "../../postprocess_modules/ShowLoopInfo.h"
namespace pf {
	namespace data_driven_complex_model {
		// - main functions
		void exec_pre_ii();
		void exec_pre_iii();
		void exec_i();
		void exec_pos_i();
		void deinit();
	}
}