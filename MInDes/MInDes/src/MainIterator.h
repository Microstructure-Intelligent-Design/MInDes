#pragma once
#include "modules/base/MACRO_DEF.h"
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <climits>
#include "MainIterator_Params.h"
#include "modules/Module.h"
#include "modules/base/timer.h"
#include "modules/input_modules/InputModulesManager.h"
namespace pf {
	namespace iterator_times {
		std::string print_time_interval();
	}

	namespace main_iterator {

		void init_modules(int argc, char* argv[]);

		void run();

	};

}
