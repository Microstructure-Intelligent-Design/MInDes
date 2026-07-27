#pragma once
#include "DDCCPAI_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace temperature_functions {

			void cal_mob_temp_grad_lap_7P(size_t x, size_t y, size_t z, Vector3& grad_temp, Vector3& grad_mob, REAL& lap_temp);

			void cal_mob_temp_grad_lap_19P(size_t x, size_t y, size_t z, Vector3& grad_temp, Vector3& grad_mob, REAL& lap_temp);

			REAL temperature_mobility_0(size_t x, size_t y, size_t z);

			REAL temperature_source_phi(size_t x, size_t y, size_t z);

			//========================================================================================================================
			// - evolution equations
			void init_temperature_field();

			void pre_calculation_temperature();

			REAL solve_temperature();

			//========================================================================================================================
			// - output

		}
	}
}