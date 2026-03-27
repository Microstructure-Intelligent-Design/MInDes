#include "Module.h"
// - preprocess
#include "preprocess_modules/MicrostructureInit.h"
// - postprocess
#include "postprocess_modules/AutoDeltTime.h"
#include "postprocess_modules/ShowLoopInfo.h"
#include "postprocess_modules/WriteVTS.h"
#include "postprocess_modules/CpuMemoryUsage.h"
#include "postprocess_modules/MachineLearning.h"
// - simlulation models
#include "model_modules/data_driven_complex/DDC_Manager.h"
#include "model_modules/grain_grows_spinodal/GGS_Manager.h"
namespace pf {
	enum SimulationModels { SM_None, SM_GGS_MODEL, SM_DDC };
	void register_all_modules() {
		// - models
		// - init models
		WriteDebugFile("# SimulationModels.model =  0 - None \n");
		WriteDebugFile("#                           1 - Grain Grows Spinodal , PCT = (N,1,false) \n");
		WriteDebugFile("#                           2 - Data Driven Complex Model , PCT = (N, K, true/false) \n");
		int sm_model = SimulationModels::SM_None;
		infile_reader::read_int_value("SimulationModels.model", sm_model, true);
		switch (SimulationModels(sm_model)) {
		case SimulationModels::SM_None: {
			// - model settings
			break;
		}
		case SimulationModels::SM_GGS_MODEL: {
			// - model settings
			grain_grows_spinodal_model::init_model_modules();
			break;
		}
		case SimulationModels::SM_DDC: {
			// - model settings
			data_driven_complex_model::init_model_modules();
			break;
		}
		// - preprocess
		microstructure_init::init_microstructure();
		// - posprocess
		automatic_change_delt_time::init_auto_time();
		show_loop_information::init_show_loop_information();
		write_vts::init_write_vts();
		cpu_memory_usage::init_cpu_memory_usage();
		machine_learning::init_machine_learning();
		}
	}
}