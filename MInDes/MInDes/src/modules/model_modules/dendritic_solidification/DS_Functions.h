#pragma once
#include "DS_Params.h"
#include <random>
namespace pf {
	namespace dendritic_solidification_model {
		// - main functions
		// - init noise con field
		void exec_pre_iii();
		void exec_i();
		void deinit();
	}
}