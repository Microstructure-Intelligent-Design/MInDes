#include "Module.h"
#include "preprocess_modules/MicrostructureInit.h"
#include "postprocess_modules/AutoDeltTime.h"
#include "postprocess_modules/ShowLoopInfo.h"
#include "postprocess_modules/WriteVTS.h"
#include "postprocess_modules/CpuMemoryUsage.h"
namespace pf {
	void register_all_modules() {
		// - 
		microstructure_init::init_microstructure();
		// - 
		automatic_change_delt_time::init_auto_time();
		// - 
		show_loop_information::init_show_loop_information();
		// - 
		write_vts::init_write_vts();
		// - 
		cpu_memory_usage::init_cpu_memory_usage();
	}
}