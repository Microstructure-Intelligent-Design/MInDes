#pragma once
#include "Model_Params.h"
#include <random>
namespace pf {
	namespace grain_grows_spinodal_model {
		// - main functions
		// - init noise con field
		void exec_pre_iii();
		void exec_i();
		void deinit();
	}
}