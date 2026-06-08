#pragma once
#include "DDCCPAI_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace chemical_energy_functions {
			void init_phase_concentration();
			// - polynomial
			REAL fchem_polynomial(std::vector<REAL> con, size_t phi_property);
			REAL miu_polynomial(std::vector<REAL> con, size_t phi_property, size_t con_index);
			// - 
			REAL delt_Fchem_delt_phi_polynomial(size_t x, size_t y, size_t z, size_t phi_index);
			REAL delt_Fchem_delt_con_polynomial(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);
			// - AI

			// - Energy minimization
			std::pair<REAL, size_t> energy_minimazation(FIELD_PhaseCon& phase_con, size_t region_index, std::vector<size_t> active_phase,
				size_t active_phase_number, REAL active_phi);
			void calculation_energy_minimazation_pre();
		}
	}
}