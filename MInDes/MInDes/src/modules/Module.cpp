#include "Module.h"
namespace pf {
	void register_all_modules() {
		// - models
		microstructure_init::init_microstructure();
		// - preprocess
		automatic_change_delt_time::init_auto_time();
		// - posprocess
		show_loop_information::init_show_loop_information();
		write_vts::init_write_vts();
		cpu_memory_usage::init_cpu_memory_usage();
	}
}