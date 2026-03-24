#pragma once
#include "../base/MACRO_DEF.h"
#include "Mechanics/StiffnessEigenStrain.h"
#include "Mechanics/ElasticSolver.h"
#include "Mechanics/PlasticSolver.h"
#include "../model_modules/GrainsOrientations.h"
#include "../model_modules/PhiProperties.h"
#include "../input_modules/ioFiles_Params.h"
#include "../Modules_Params.h"
#include "WriteVTS.h"
namespace pf {
	namespace mechanical_field {
		void init();
		void exec_pre_iii();
		std::string exec_loop_iii();
		void deinit();
	}
}