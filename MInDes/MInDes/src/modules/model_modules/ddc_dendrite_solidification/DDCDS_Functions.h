#pragma once
#include "DDCDS_Phi_Functions.h"
#include "DDCDS_Con_Functions.h"
#include "DDCDS_Temp_Functions.h"
#include "DDCDS_BulkEnergy.h"
#include "../../postprocess_modules/ShowLoopInfo.h"
namespace pf {
	namespace ddc_dendrite_solidification {
		// - main functions
		void exec_pre_ii();
		void exec_pre_iii();
		void exec_i();
		void exec_pos_i();
		void exec_pos_iii();
		void deinit();
	}
}