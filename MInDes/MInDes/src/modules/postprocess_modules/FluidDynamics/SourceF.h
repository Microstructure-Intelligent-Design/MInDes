#pragma once
#include "../../base/MACRO_DEF.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "GoverningEquationsF.h"
namespace pf {
	enum LBM_Source_Type { LBM_ST_Forces };
	enum LBM_Force_Term_Model { LBM_FTM_NONE, LBM_FTM_LGA, LBM_FTM_ZS_Guo, LBM_FTM_H_Liang };
	enum LBM_Force_Type { LBM_FM_ThermalExpansion, LBM_FM_Gravity, LBM_FM_H_Liang_SurfaceTension };
	namespace lbm_source {
		inline double Cs2 = 1.0 / 3.0;
		inline double Cs4 = 1.0 / 9.0;
		// two phase flow
		inline double surface_tension = 0.0;
		inline double interface_thickness = 1.0;
		inline double beta = 0.0;
		inline double kappa = 0.0;
		inline double get_interface_thickness() {
			return interface_thickness;
		}
		inline std::vector<double(*)(long long x, long long y, long long z, double tau, size_t LBM_i)> fluid_source_list;
		inline double(*fluid_two_phase_source_list)(long long x, long long y, long long z, double tau, Vector3 prefactor, size_t LBM_i);

		namespace force_funcs {
			// external force parameters
			inline Vector3 gravitational_acceleration(0.0, 0.0, 0.0);
			// fluid force
			inline double ref_density = 0.0;
			Vector3 Fluid_Force_Gravity(long long x, long long y, long long z);
			// 
			inline double ref_temp = 0.0;
			inline double thermal_expansion_parameter = 0.0;
			Vector3 Fluid_Force_Thermal_Expansion(long long x, long long y, long long z);
			Vector3 Fluid_Force_H_Liang_Surface_Tension(long long x, long long y, long long z);
			//
			inline std::vector<Vector3(*)(long long x, long long y, long long z)> fluid_force_list;
			// fluid force model
			void load_forces();
		}
		// source model
		double fluid_source_i(long long x, long long y, long long z, double tau, size_t LBM_i);
		double fluid_two_phase_source_i(long long x, long long y, long long z, double tau, Vector3 prefactor, size_t LBM_i);
		void lbm_properties_automatically_change();
		void init(LBM& fluid_lbm_solver);
		void init_two_phase_solver(LBM& field_lbm_two_phase_solver);
	}
}