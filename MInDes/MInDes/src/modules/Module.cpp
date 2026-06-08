#include "Module.h"
// - preprocess
#include "preprocess_modules/MicrostructureInit.h"
#include "preprocess_modules/Pretreatment.h"
// - postprocess
#include "postprocess_modules/AutoDeltTime.h"
#include "postprocess_modules/ShowLoopInfo.h"
#include "postprocess_modules/WriteVTS.h"
#include "postprocess_modules/CpuMemoryUsage.h"
#include "postprocess_modules/MachineLearning.h"
// - simlulation models template
#include "model_modules/Data_Driven_Complex/DDC_Manager.h"
// - simlulation models 
#include "model_modules/Grain_Grows_Spinodal/GGS_Manager.h"
#include "model_modules/ddc_calphad_ai/DDCCPAI_Manager.h"
namespace pf {
	enum SimulationModels { SM_None, SM_GGS_MODEL, SM_DDC, SM_DDC_CPAI };
	void register_all_modules() {
		// - basic functions
		automatic_change_delt_time::init_auto_time();
		microstructure_init::init_microstructure();
		pretreatment::init_pretreatment();
		// - models
		// - init models
		WriteDebugFile("========================================================================================= \n");
		WriteDebugFile("# SimulationModels.model =  0 - None \n");
		WriteDebugFile("#                           1 - Grain Grows Spinodal , PCT = (N,1,false) \n");
		WriteDebugFile("#                           2 - Data Driven Complex Model , PCT = (N > 0, K, true/false) \n");
		WriteDebugFile("#                           3 - Data Driven Complex Model - coupled with CALPHAD & AI , PCT = (N > 0, K, true/false) \n");
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
		case SimulationModels::SM_DDC_CPAI: {
			// - model settings
			ddc_calphad_ai_model::init_model_modules();
			break;
		}
		}
		WriteDebugFile("========================================================================================= \n");
		// - other method
		machine_learning::init_machine_learning();
		// - tail
		WriteDebugFile("========================================================================================= \n");
		show_loop_information::init_show_loop_information();
		write_vts::init_write_vts();
		cpu_memory_usage::init_cpu_memory_usage();
		WriteDebugFile("========================================================================================= \n");
		WriteDebugFile("=============================== Parameters Definition End ===============================\n");
		WriteDebugFile("=========================================================================================\n");
	}
}