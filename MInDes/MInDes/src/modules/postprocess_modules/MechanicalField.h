#pragma once
#include "../base/MACRO_DEF.h"
namespace pf {
	namespace mechanical_field {
		void init();
		void exec_pre_iii();
		std::string exec_loop_iii();
		void deinit();
	}
}