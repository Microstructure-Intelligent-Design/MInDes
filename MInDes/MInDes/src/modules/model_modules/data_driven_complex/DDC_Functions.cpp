#pragma once
#include "DDC_Functions.h"
namespace pf {
	namespace data_driven_complex_model {
		void exec_pre_iii() {
			parameters::field_variable.init();
			phase_field_functions::init_phi_pair_wise();
		}
		void exec_i() {
			phase_field_functions::pre_calculation_phi_pair_wise();
			phase_field_functions::solve_phi_pair_wise();
		}
		void deinit() {

		}
	}
}