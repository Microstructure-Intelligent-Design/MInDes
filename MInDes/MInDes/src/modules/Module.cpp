#include "Module.h"
// - preprocess
#include "preprocess_modules/MicrostructureInit.h"
#include "preprocess_modules/Pretreatment.h"
// - postprocess
#include "postprocess_modules/AutoDeltTime.h"
#include "postprocess_modules/ShowLoopInfo.h"
#include "postprocess_modules/WriteVTS.h"
#include "postprocess_modules/CpuMemoryUsage.h"
// - simlulation models template
#include "model_modules/data_driven_complex/DDC_Manager.h"
// - simlulation models
#include "model_modules/grain_grows_spinodal/GGS_Manager.h"
#include "model_modules/dendrite_solidification/DS_Manager.h"
#include "model_modules/solid_state_sintering/SSS_Manager.h"
// - external field
#include "postprocess_modules/FluidDynamics/LatticeBoltzmann.h"
#include "postprocess_modules/Mechanics/ElasticSolver.h"
namespace pf {
	enum SimulationModels { SM_None, SM_GGS, SM_DS, SM_SSS, SM_DDC };
	void register_all_modules() {
		// - basic functions
		microstructure_init::init_microstructure();
		pretreatment::init_pretreatment();
		// - models
		// - init models
		WriteDebugFile("========================================================================================= \n");
		WriteDebugFile("# SimulationModels.model =  0 - None \n");
		WriteDebugFile("#                           1 - Grain Grows Spinodal , PCT = (N,1,false) \n");
		WriteDebugFile("#                           2 - Dendrite Solidification , PCT = (1,0,true), 2D \n");
		WriteDebugFile("#                           3 - Solid State Sintering , PCT = (N >= 1,1,false) \n");
		WriteDebugFile("#                           4 - Data Driven Complex Model , PCT = (N > 0, K, true/false) \n");
		int sm_model = SimulationModels::SM_None;
		infile_reader::read_int_value("SimulationModels.model", sm_model, true);
		switch (SimulationModels(sm_model)) {
		case SimulationModels::SM_None: {
			// - model settings
			break;
		}
		case SimulationModels::SM_GGS: {
			// - model settings
			grain_grows_spinodal_model::init_model_modules();
			break;
		}
		case SimulationModels::SM_DS: {
			dendrite_solidification_model::init_model_modules();
			break;
		}
		case SimulationModels::SM_SSS: {
			solid_state_sintering_model::init_model_modules();
			break;
		}
		case SimulationModels::SM_DDC: {
			// - model settings
			data_driven_complex_model::init_model_modules();
			break;
		}
		}
		WriteDebugFile("========================================================================================= \n");
		show_loop_information::init_show_loop_information();
		// - other method
		automatic_change_delt_time::init_auto_time();
		if (external_physical_field::is_fluid_field_on)
			lattice_boltzmann::init();
		if (external_physical_field::is_mech_field_on)
			elastic_solver::init();
		// - tail
		WriteDebugFile("========================================================================================= \n");
		cpu_memory_usage::init_cpu_memory_usage();
		write_vts::init_write_vts();
		WriteDebugFile("========================================================================================= \n");
		WriteDebugFile("=============================== Parameters Definition End ===============================\n");
		WriteDebugFile("=========================================================================================\n");
	}
}
