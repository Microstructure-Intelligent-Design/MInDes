#pragma once
#include "DDCCPAI_Phi_Functions.h"
#include "DDCCPAI_Con_Functions.h"
#include "DDCCPAI_Temp_Functions.h"
#include "../../postprocess_modules/ShowLoopInfo.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// - main functions
		void exec_pre_ii();
		void exec_pre_iii();
		void exec_i();
		void exec_pos_i();
		void deinit();
	}
}