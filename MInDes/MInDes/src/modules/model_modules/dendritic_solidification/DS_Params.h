#pragma once
#include "../../base/Mesh_0.h"
#include "../../Modules_Params.h"
namespace pf {
	namespace dendritic_solidification_model {
		enum FIELD { MOB, DFDCON, NUM };
		namespace parameters {
			// - field
			inline Mesh_Boundry<Matrix1D<REAL>>  con_field_variavles;
			inline Mesh_Boundry<Matrix1D<REAL>>  phi_increment;
			// - parameters
			inline REAL A = 1;
			inline REAL epsilon = 20;
			inline REAL kappa_eta = 1;
			inline REAL kappa_con = 1;
			inline REAL L = 1;
			inline REAL Mb = 1;
			inline REAL Mg = 0;
			inline bool is_noise_rand = true;
			inline int noise_seed = 0;
			inline REAL noise_amplitude = 0;
		}
	}
}