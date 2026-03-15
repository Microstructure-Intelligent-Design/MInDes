#include "EquiDistFunc.h"
namespace pf {
	namespace lbm_equilibrium_distribution_function {
		namespace default_functions {
			static double f_eq_a_d2q9_standard(size_t INDEX_A, double p_macro, double f_macro, Vector3& U) {
				double CU = LBM_d2q9_w_vec[INDEX_A] * U;
				return LBM_d2q9_w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double f_eq_a_d3q19_standard(size_t INDEX_A, double p_macro, double f_macro, Vector3& U) {
				double CU = LBM_d3q19_w_vec[INDEX_A] * U;
				return LBM_d3q19_w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
		}

		void lbm_properties_automatically_change() {
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
		}

		void init(LBM& fluid_lbm_solver) {
			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				f_eq_i = default_functions::f_eq_a_d2q9_standard;
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				f_eq_i = default_functions::f_eq_a_d3q19_standard;
			}
		}
	}
}