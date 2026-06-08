#pragma once
#include "DDCCPAI_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace concentration_field_functions {

			void init_con_in_moving_region_default(size_t x, size_t y, size_t z, size_t region_index);

			void deinit_con_in_moving_region_default(size_t x, size_t y, size_t z, size_t region_index);

			void cal_mob_miu_grad_lap_7P(size_t x, size_t y, size_t z, size_t region_index,
				Vector3& grad_region, std::vector<Vector3>& grad_miu, std::vector<Vector3>& grad_mob, std::vector<REAL>& lap_miu);

			void cal_mob_miu_grad_lap_19P(size_t x, size_t y, size_t z, size_t region_index,
				Vector3& grad_region, std::vector<Vector3>& grad_miu, std::vector<Vector3>& grad_mob, std::vector<REAL>& lap_miu);

			REAL dfdcon_con(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL dfdcon_lap_con_7P(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL dfdcon_lap_con_19P(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL surface_reaction_flux(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL bulk_reaction_flux(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL bulk_diffusion_mobility(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL bulk_diffusion_mobility_with_temperature(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL surface_diffusion_mobility(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			REAL interphase_diffusion_mobility(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);

			//========================================================================================================================
			// - evolution equations
			void init_total_concentration();

			void pre_calculation_total_concentration();

			REAL solve_total_concentration();

			//========================================================================================================================
			// - output
			void write_scalar_con_all(std::ofstream& fout);

		}
	}
}
