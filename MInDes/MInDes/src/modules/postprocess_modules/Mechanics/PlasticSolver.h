#pragma once
#include "StiffnessEigenStrain.h"
#include "../../postprocess_modules/WriteVTS.h"
#include "../../Modules_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "../../model_modules/PhiProperties.h"
namespace pf {
	namespace plastic_solver {
		// Plastic Parameters
		inline size_t phi_number = 0;
		inline size_t phi_property_number = 0;
		inline std::vector<double> phase_yield_stress;
		inline std::vector<double> phi_yield_stress;
		inline std::vector<double> phase_hardening_modulus;
		inline std::vector<double> phi_hardening_modulus;
		inline std::vector<double> phase_shear_modulus;
		inline std::vector<double> phi_shear_modulus;
		inline double get_phase_hardening_modulus(size_t phi_property) {
			return phase_hardening_modulus[phi_property];
		}
		inline double get_phi_hardening_modulus(size_t phi_index) {
			return phi_hardening_modulus[phi_index];
		}
		inline double get_phase_yield_stress(size_t phi_property) {
			return phase_yield_stress[phi_property];
		}
		inline double get_phi_yield_stress(size_t phi_index) {
			return phi_yield_stress[phi_index];
		}
		inline double get_phase_shear_modulus(size_t phi_property) {
			return phase_shear_modulus[phi_property];
		}
		inline double get_phi_shear_modulus(size_t phi_index) {
			return phi_shear_modulus[phi_index];
		}
		// stiffness and eigenstain
		inline void(*cal_plastic_parameters)(long long x, long long y, long long z,
			double& yield_stress, double& hardening_modulus, double& shear_modulus);
		// 
		double Prandtl_Reuss_formula(double deviatoric_part_of_mises_stress_norm, double yield_stress
			, double hardening_modulus, double ave_plastic_strain);

		void solve_plastic_flow(double& max_dplastic_strain, double& max_dave_plastic_strain);

		void init();

		void deinit();

		void write_scalar(std::ofstream& fout);
	}
}