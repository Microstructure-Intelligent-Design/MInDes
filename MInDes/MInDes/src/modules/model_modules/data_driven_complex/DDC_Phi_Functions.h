#pragma once
#include "DDC_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
namespace pf {
	namespace data_driven_complex_model {
		// - default functions
		namespace phase_field_functions {
			// - phi
			inline REAL interpolation_func(REAL phi) {
				return phi * phi * (3 - 2 * phi);
			}

			inline REAL dinterpolation_func(REAL phi) {
				return 6 * phi * (1 - phi);
			}
			void interphase_gradient_lapace_calculation_7P(size_t x, size_t y, size_t z, size_t phi_index);
			void interphase_gradient_lapace_calculation_19P(size_t x, size_t y, size_t z, size_t phi_index);
			inline void (*interphase_gradient_lapace_calculation)(size_t x, size_t y, size_t z, size_t phi_index);
			InterfaceFlag currentFlag_7P(size_t x, size_t y, size_t z, size_t phi_index);
			InterfaceFlag currentFlag_19P(size_t x, size_t y, size_t z, size_t phi_index);
			inline InterfaceFlag(*currentFlag)(size_t x, size_t y, size_t z, size_t phi_index);
			void upgradeFlag_7P(size_t x, size_t y, size_t z, size_t phi_index);
			void upgradeFlag_19P(size_t x, size_t y, size_t z, size_t phi_index);
			inline void(*upgradeFlag)(size_t x, size_t y, size_t z, size_t phi_index);
			// - 
			Vector3 normals(REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			// - mobility
			REAL Lij(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_temp(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_temp_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_temp_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_temp_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			REAL Lij_temp_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			// - interface energy
			REAL xi_ab(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			REAL xi_ab_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			REAL xi_ab_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			REAL xi_ab_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			REAL xi_ab_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			REAL xi_abc(size_t alpha_index, size_t beta_index, size_t gamma_index);
			// - 
			inline REAL(*_Lij)(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature);
			inline REAL(*_xi_ab)(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad);
			inline REAL(*_xi_abc)(size_t alpha_index, size_t beta_index, size_t gamma_index);
			// - interface energy model
			REAL dfint_dphi_grad_S1996_acc(FIELD_PhiTemp& point, size_t phi_index);
			REAL dfint_dphi_grad_S1999_acc(FIELD_PhiTemp& point, size_t phi_index);
			REAL dfint_dphi_pot_Nobstacle_acc(FIELD_PhiTemp& point, size_t phi_index);
			REAL dfint_dphi_pot_Nwell_acc(FIELD_PhiTemp& point, size_t phi_index);
			REAL dfint_dphi_pairwise_S2009(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index);
			REAL dfint_dphi_pairwise_S1996_NO(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index);
			REAL dfint_dphi_pairwise_S1996_NW(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index);
			REAL dfint_dphi_pairwise_S1999_NO(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index);
			REAL dfint_dphi_pairwise_S1999_NW(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index);
			// - 
			inline REAL(*dfint_dphi_pairwise_acc)(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index);
			// - 
			void pairwise_normalize_acc(std::vector<REAL>& active_index, std::vector<bool>& interphase, const Matrix1D<REAL>& old_phi, std::vector<REAL>& phi_increment);
			REAL dfbulk_dphi_pairwise_acc(long long x, long long y, long long z, size_t phi_index);
			REAL source_pairwise_acc(long long x, long long y, long long z, size_t alpha_index, size_t beta_index);
			// - default functions
			REAL dfbulk_dphi_0(size_t x, size_t y, size_t z, size_t phi_index);
			//=========================================================================================================================================
			// - phi pair-wise evolution equation
			void init_phi_pair_wise();
			void pre_calculation_phi_pair_wise();
			REAL solve_phi_pair_wise();
			//=========================================================================================================================================
			// - output
			void write_scalar_active_phi_number(std::ofstream& fout);
		}
	}
}
