#pragma once
#include "../../base/Mesh_0.h"
namespace pf {
	namespace grain_grows_spinodal_model {
		enum FIELD { MOB, DFDCON, NUM };
		namespace parameters {
			// - field
			inline Mesh_Boundry<std::vector<REAL>>* phase_field;
			inline Mesh_Boundry<std::vector<REAL>>* concentration_field;
			inline Mesh_Boundry<std::vector<REAL>>  con_field_variavles;
			inline Mesh_Boundry<std::vector<REAL>>  phi_increment;
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
			// - 
			inline size_t MESH_NX = 1;
			inline size_t MESH_NY = 1;
			inline size_t MESH_NZ = 1;
			inline BoundaryCondition x_down = BoundaryCondition::PERIODIC;
			inline BoundaryCondition y_down = BoundaryCondition::PERIODIC;
			inline BoundaryCondition z_down = BoundaryCondition::PERIODIC;
			inline BoundaryCondition x_up = BoundaryCondition::PERIODIC;
			inline BoundaryCondition y_up = BoundaryCondition::PERIODIC;
			inline BoundaryCondition z_up = BoundaryCondition::PERIODIC;
			inline REAL delt_r = 1.0;
			inline REAL* delt_t;
			inline size_t phi_number = 0;
			inline REAL* PHI_MAX_VARIATION;
			inline size_t con_number = 0;
			inline REAL* CON_MAX_VARIATION;
		}
	}
}