#pragma once
#include "../../base/MACRO_DEF.h"
#include "../WriteVTS.h"
#include "GoverningEquationsF.h"
#include "BoundaryCondition.h"
#include "SourceF.h"
#include "EquiDistFunc.h"
#include "MacroVar.h"
namespace pf {
	namespace lattice_boltzmann {
		// FST_LBM_Difference
		inline LBM fluid_lbm_solver;
		inline double momentum_accuracy = 1e-4;
		inline bool debug_solver = false;
		inline int debug_output_step = 1000;
		inline int max_iterate_steps = 0;
		// normal parameters
		// two phase flow
		inline bool is_two_phase_flow = false;
		// inline LBM fluid_lbm_two_phase_solver;
		// Source
		inline double (*fluid_source_i)(long long x, long long y, long long z, double tau, size_t LBM_F_i);
		inline double (*fluid_two_phase_source_i)(long long x, long long y, long long z, double tau, Vector3 prefactor, size_t LBM_F_i);
		inline double (*f_eq_i)(size_t INDEX_i, double p_macro, double f_macro, Vector3& U);
		inline double (*f_eq_two_phase_i)(size_t INDEX_i, double p_macro, double f_macro, Vector3& U);
		void init_distribution_functions_d2q9();
		void init_distribution_functions_d3q19();
		void lbm_properties_automatically_change();
		void init();
		void exec_pre();
		std::string exec_loop();
		void deinit();
		void write_velocity(std::ofstream& fout);
		void write_abs_velocity(std::ofstream& fout);
		void write_pressure(std::ofstream& fout);
		void write_density(std::ofstream& fout);
	}
}