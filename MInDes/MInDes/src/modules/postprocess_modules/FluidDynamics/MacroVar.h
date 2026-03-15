#pragma once
#include "../../base/MACRO_DEF.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "BoundaryCondition.h"
#include "SourceF.h"
#include "GoverningEquationsF.h"
namespace pf {
	namespace lbm_macro_variable {
		inline double Cs2 = 1.0 / 3.0;
		inline double Cs4 = 1.0 / 9.0;
		inline int index_h_liang = 0;
		namespace macro_variable_funcs {
			inline vector<Vector3(*)(long long x, long long y, long long z)> fluid_force_list;
			void load_forces();
			const double w0_d2q9 = 4.0 / 9.0;
			const double w0_d3q19 = 12.0 / 36.0;
		}

		inline void (*cal_macro_variables)(long long x, long long y, long long z);

		inline void (*cal_macro_variables_two_phase)(long long x, long long y, long long z);

		void lbm_properties_automatically_change();

		void init(LBM& fluid_lbm_solver);

		void init_two_phase(LBM& fluid_lbm_solver, LBM& field_lbm_two_phase_solver);
	}
}