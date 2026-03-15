#pragma once
#include "../../base/MACRO_DEF.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../model_modules/PhiProperties.h"
#include "GoverningEquationsF.h"
namespace pf {
	namespace lbm_boundary_condition {
		enum Fluid_Domain_Boundary_Condition { FDBC_Wall, FDBC_Period, FDBC_Free, FDBC_Pressure, FDBC_Normal_Flow, FDBC_SIZE };
		enum Fluid_Boundary_Condition { FBC_UP_X, FBC_UP_Y, FBC_UP_Z, FBC_DOWN_X, FBC_DOWN_Y, FBC_DOWN_Z, FBC_SIZE };
		enum Fluid_Boundary_Property { FBP_WallRoughness, FBP_WallSpeed, FBP_DensityValue, FBP_NormalFlowSpeed, FBP_SIZE };
		// boundary condition
		inline std::vector<std::vector<REAL>> fluid_boundary_condition;
		// smooth boundary
		const double solid_liquid_interface_threshold = 0.5;
		const double q_c = 0.75;
		inline std::vector<bool> is_solid_phases;   // phase index
		inline double Cs2 = 1.0 / 3.0;
		inline double Cs4 = 1.0 / 9.0;
		// viscosity
		inline double viscosity_liquid = 0.0;
		inline double viscosity_gas = 0.0;
		inline double viscosity_one_phase(long long x, long long y, long long z) {
			return viscosity_liquid;
		}
		inline double (*viscosity)(long long x, long long y, long long z);
		// density
		inline double density_liquid = 1.0;
		inline double density_gas = 1.0;
		inline double density_one_phase(long long x, long long y, long long z) {
			return density_liquid;
		}
		inline double (*density)(long long x, long long y, long long z);
		// tau
		inline double mobility_two_phase = 0.0;
		inline double _tau_const = DBL_MAX;
		inline double _tau_two_phase = DBL_MAX;
		inline double tau_const(long long x, long long y, long long z) {
			return _tau_const;
		}
		inline double tau_standard(long long x, long long y, long long z) {
			return viscosity(x, y, z) / time_parameters::delt_t / Cs2 + 0.5;
		}
		inline double tau_two_phase_const(long long x, long long y, long long z) {
			return _tau_two_phase;
		}
		inline double (*tau)(long long x, long long y, long long z);
		inline double (*tau_two_phase)(long long x, long long y, long long z);

		inline void (*d2q9_domain_boundary_x_down)(long long x, long long y, long long z);
		inline void (*d2q9_domain_boundary_x_up)(long long x, long long y, long long z);
		inline void (*d2q9_domain_boundary_y_down)(long long x, long long y, long long z);
		inline void (*d2q9_domain_boundary_y_up)(long long x, long long y, long long z);
		inline void(*d2q9_fluid_solid_boundary)(long long x, long long y, long long z);
		void boundary_condition_d2q9(long long x, long long y, long long z);

		inline void (*d3q19_domain_boundary_x_down)(long long x, long long y, long long z);
		inline void (*d3q19_domain_boundary_x_up)(long long x, long long y, long long z);
		inline void (*d3q19_domain_boundary_y_down)(long long x, long long y, long long z);
		inline void (*d3q19_domain_boundary_y_up)(long long x, long long y, long long z);
		inline void (*d3q19_domain_boundary_z_down)(long long x, long long y, long long z);
		inline void (*d3q19_domain_boundary_z_up)(long long x, long long y, long long z);
		inline void(*d3q19_fluid_solid_boundary)(long long x, long long y, long long z);
		void boundary_condition_d3q19(long long x, long long y, long long z);

		void cal_fluid_domain();

		void init(LBM& fluid_lbm_solver);

		void init_two_phase_solver(LBM& field_lbm_two_phase_solver);

		void lbm_properties_automatically_change();

		void deinit();
	}
}