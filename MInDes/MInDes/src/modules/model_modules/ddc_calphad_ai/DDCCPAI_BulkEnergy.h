#pragma once
#include "DDCCPAI_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
#include "../../postprocess_modules/ShowLoopInfo.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace chemical_energy_functions {
			// - polynomial
			REAL fchem_polynomial(std::vector<REAL> con, REAL temperature, size_t phi_property);
			REAL miu_polynomial(std::vector<REAL> con, REAL temperature, size_t phi_property, size_t con_index);
			// - 
			REAL delt_Fchem_delt_phi(size_t x, size_t y, size_t z, size_t phi_index);
			REAL delt_Fchem_delt_con(size_t x, size_t y, size_t z, size_t region_index, size_t con_index);
			// - AI

			// - Energy minimization
			void init_phase_concentration();
			// { VARIATION , ITERATION_STEP }
			std::pair<REAL, size_t> energy_minimazation(std::vector<REAL>& phase_phi, std::vector<std::vector<REAL>>& phase_con, std::vector<std::vector<REAL>>& phase_miu, std::vector<size_t> active_phase,
				size_t active_phase_number, REAL active_phi, REAL temperature, size_t region_index);
			// { MAX_VARIATION , MAX_ITERATION_STEP }
			std::pair<REAL, size_t> local_concentration_redistribution(size_t x, size_t y, size_t z);
			void calculation_energy_minimazation_pre();
			// - output
			void write_scalar_phasecon(std::ofstream& fout);
			void write_scalar_con(std::ofstream& fout);
			void write_scalar_phasemiu(std::ofstream& fout);
			void write_scalar_miu(std::ofstream& fout);
			void write_scalar_fchem(std::ofstream& fout);
			// - 
			void write_csv_energy_results();
			void write_csv_energy_minimization_results();
		}
	}
}
