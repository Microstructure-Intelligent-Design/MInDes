#include "Module.h"
#include "preprocess_modules/MicrostructureInit.h"
#include "postprocess_modules/AutoDeltTime.h"
#include "postprocess_modules/ShowLoopInfo.h"
#include "postprocess_modules/WriteVTS.h"
#include "postprocess_modules/CpuMemoryUsage.h"
// - simlulation models

namespace pf {
	enum SimulationModels { SM_None };
	void register_all_modules() {
		// - models
		// - init models
		WriteDebugFile("# SimulationModels.model =  0 - None \n");
		int sm_model = SimulationModels::SM_None;
		infile_reader::read_int_value("SimulationModels.model", sm_model, true);
		switch (SimulationModels(sm_model)) {
		case SimulationModels::SM_None: {
			// - model settings
			break;
		}
		}
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