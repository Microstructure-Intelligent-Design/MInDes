#pragma once
#include "../../base/MACRO_DEF.h"
#include "GoverningEquationsF.h"
#include "../../Modules_Params.h"
namespace pf {
	namespace lbm_equilibrium_distribution_function {
		inline double Cs2 = LBM_Cs2;
		inline double Cs4 = LBM_Cs4;

		inline double (*f_eq_i)(size_t INDEX_i, double p_macro, double f_macro, Vector3& U);

		inline double (*f_eq_two_phase_i)(size_t INDEX_i, double p_macro, double f_macro, Vector3& U);

		void lbm_properties_automatically_change();

		void init(LBM& fluid_lbm_solver);

	}
}