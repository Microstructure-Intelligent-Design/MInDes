#pragma once
#include "../base/MACRO_DEF.h"
#include "../input_modules/ioFiles_Params.h"
#include "FluidDynamics/LatticeBoltzmann.h"
namespace pf {
	namespace fluid_field {
		// statement functons
		void init();
		void exec_pre();
		std::string exec_loop();
		void deinit();
	}
}