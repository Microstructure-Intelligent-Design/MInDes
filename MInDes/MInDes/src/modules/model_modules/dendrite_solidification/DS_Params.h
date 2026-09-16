#pragma once
#include "../../Modules_Params.h"

namespace pf {
	namespace dendrite_solidification_model {
		enum FIELD { EPSILON, FLUX_X, FLUX_Y, NUM };
		namespace parameters {
			inline Mesh_Boundry<Matrix1D<REAL>> anisotropy;
			inline Mesh_Boundry<REAL> phi_increment;
			inline Mesh_Boundry<REAL> temperature_increment;
			inline REAL tau = 0.0003;
			inline REAL epsilonb = 0.01;
			inline REAL kappa = 1.8;
			inline REAL delta = 0.02;
			inline int aniso = 4;
			inline REAL alpha = 0.9;
			inline REAL gamma = 10.0;
			inline REAL teq = 1.0;
			inline REAL theta0 = 0.2;
		}
	}
}
