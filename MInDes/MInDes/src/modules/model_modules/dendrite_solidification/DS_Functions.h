#pragma once
#include "DS_Params.h"

namespace pf {
	namespace dendrite_solidification_model {
		namespace field_functions {
			// Default Lij callback with cubic anisotropy in the canonical beta crystal frame.
			REAL interface_mobility_const(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi,
				const Vector3& alpha_grad, const Vector3& beta_grad, REAL temperature);
			REAL concentration_mobility_const(long long x, long long y, long long z, size_t con_index);
			REAL temperature_mobility_mixture(long long x, long long y, long long z);
			// Unified entry: concentration plus temperature for K>0, temperature only for K=0.
			REAL linearized_driving_force(long long x, long long y, long long z);
			REAL solidification_driving_pair_rate(long long x, long long y, long long z, size_t solid_index, REAL mobility);

			void init_dendrite_field();
			void prepare_step_snapshot();
			void prepare_thermodynamic_fields();
			void calculate_right_hand_sides();
			void solve_fields();
		}

		void exec_pre_iii();
		void exec_i();
		void deinit();
	}
}
