#pragma once
#include "MicrostructureInit.h"
namespace pf {
	inline void init_preprocess_modules() {
		microstructure_init::init();
		WriteLog("> MODULE INIT : Preprocess Ready !\n");
	}
}
