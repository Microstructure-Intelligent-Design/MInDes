#pragma once
#include "ShowLoopInfo.h"
#include "WriteVTS.h"
#include "AutoDeltTime.h"
#include "CpuMemoryUsage.h"
namespace pf {
	inline void init_postprocess_modules() {
		automatic_change_delt_time::init_auto_time();
		// output
		show_loop_information::init();
		write_vts::init();
		cpu_memory_usage::init();
		WriteLog("> MODULE INIT : Postprocess Ready !\n");
	}
}
