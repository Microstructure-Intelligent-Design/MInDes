#include "BoundaryCondition.h"
namespace pf {
	namespace lbm_boundary_condition {
		namespace bc_funcs {
			static void default_domain_boundary_condition(long long x, long long y, long long z) {};
			static double f_eq_a_virtual_d2q9(int INDEX_A, double f_macro, Vector3& U) {
				double CU = LBM_d2q9_w_vec[INDEX_A] * U;
				return LBM_d2q9_w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static double f_eq_a_virtual_d3q19(int INDEX_A, double f_macro, Vector3& U) {
				double CU = LBM_d3q19_w_vec[INDEX_A] * U;
				return LBM_d3q19_w[INDEX_A] * f_macro * (1.0 + CU / Cs2 + CU * CU / Cs4 / 2.0 - U * U / Cs2 / 2.0);
			}
			static void d2q9_fluid_solid_boundary_Guo2002(long long x, long long y, long long z) {
				double f_eq = 0.0, f_neq = 0.0, _tau = tau(x, y, z);
				Vector3 U;
				// LBM_f_0
				LBMPoint& point = external_physical_field::lbm_field(x, y, z);
				point.F[LBM_0] = 0.0;
				LBMPoint* near_point;
				// LBM_f_1
				near_point = &external_physical_field::lbm_field(x + 1, y, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_1] - f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_1] - f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_1] - f_eq_a_virtual_d2q9(LBM_1, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_1] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_1] = 0.0;
				}
				// LBM_f_2
				near_point = &external_physical_field::lbm_field(x, y + 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_2] - f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y + 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_2] - f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_2] - f_eq_a_virtual_d2q9(LBM_2, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_2] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_2] = 0.0;
				}
				// LBM_f_3
				near_point = &external_physical_field::lbm_field(x - 1, y, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_3] - f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_3] - f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_3] - f_eq_a_virtual_d2q9(LBM_3, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_3] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_3] = 0.0;
				}
				// LBM_f_4
				near_point = &external_physical_field::lbm_field(x, y - 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_4] - f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y - 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_4] - f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_4] - f_eq_a_virtual_d2q9(LBM_4, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_4] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_4] = 0.0;
				}
				// LBM_f_5
				near_point = &external_physical_field::lbm_field(x + 1, y + 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_5] - f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y + 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_5] - f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_5] - f_eq_a_virtual_d2q9(LBM_5, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_5] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_5] = 0.0;
				}
				// LBM_f_6
				near_point = &external_physical_field::lbm_field(x - 1, y + 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_6] - f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y + 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_6] - f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_6] - f_eq_a_virtual_d2q9(LBM_6, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_6] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_6] = 0.0;
				}
				// LBM_f_7
				near_point = &external_physical_field::lbm_field(x - 1, y - 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_7] - f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y - 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_7] - f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_7] - f_eq_a_virtual_d2q9(LBM_7, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_7] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_7] = 0.0;
				}
				// LBM_f_8
				near_point = &external_physical_field::lbm_field(x + 1, y - 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_8] - f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y - 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_8] - f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_8] - f_eq_a_virtual_d2q9(LBM_8, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_8] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_8] = 0.0;
				}
			}
			namespace lbm_bc_d2q9 {
				// FDBC_Wall_No_Slip
				static void wall_no_slip_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_1] = point.F[LBM_3];
						point.F[LBM_5] = point.F[LBM_7] * roughness + point.F[LBM_6] * (1.0 - roughness);
						point.F[LBM_8] = point.F[LBM_6] * roughness + point.F[LBM_7] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_3] = point.F[LBM_1];
						point.F[LBM_6] = point.F[LBM_8] * roughness + point.F[LBM_5] * (1.0 - roughness);
						point.F[LBM_7] = point.F[LBM_5] * roughness + point.F[LBM_8] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_2] = point.F[LBM_4];
						point.F[LBM_5] = point.F[LBM_7] * roughness + point.F[LBM_8] * (1.0 - roughness);
						point.F[LBM_6] = point.F[LBM_8] * roughness + point.F[LBM_7] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_4] = point.F[LBM_2];
						point.F[LBM_7] = point.F[LBM_5] * roughness + point.F[LBM_6] * (1.0 - roughness);
						point.F[LBM_8] = point.F[LBM_6] * roughness + point.F[LBM_5] * (1.0 - roughness);
					}
				}
				// FDBC_Wall_Slip
				static void wall_slip_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double loc_density = point.F[LBM_0] + point.F[LBM_2] + point.F[LBM_4]
							+ 2 * (point.F[LBM_3] + point.F[LBM_6] + point.F[LBM_7]),
							tail = loc_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_1] = point.F[LBM_3];
						point.F[LBM_5] = point.F[LBM_7] * roughness + point.F[LBM_6] * (1.0 - roughness) + tail;
						point.F[LBM_8] = point.F[LBM_6] * roughness + point.F[LBM_7] * (1.0 - roughness) - tail;
					}
				}
				static void wall_slip_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double loc_density = point.F[LBM_0] + point.F[LBM_2] + point.F[LBM_4]
							+ 2 * (point.F[LBM_1] + point.F[LBM_5] + point.F[LBM_8]),
							tail = loc_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_3] = point.F[LBM_1];
						point.F[LBM_6] = point.F[LBM_8] * roughness + point.F[LBM_5] * (1.0 - roughness) + tail;
						point.F[LBM_7] = point.F[LBM_5] * roughness + point.F[LBM_8] * (1.0 - roughness) - tail;
					}
				}
				static void wall_slip_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double loc_density = point.F[LBM_0] + point.F[LBM_1] + point.F[LBM_3]
							+ 2 * (point.F[LBM_4] + point.F[LBM_7] + point.F[LBM_8]),
							tail = loc_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_2] = point.F[LBM_4];
						point.F[LBM_5] = point.F[LBM_7] * roughness + point.F[LBM_8] * (1.0 - roughness) + tail;
						point.F[LBM_6] = point.F[LBM_8] * roughness + point.F[LBM_7] * (1.0 - roughness) - tail;
					}
				}
				static void wall_slip_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double loc_density = point.F[LBM_0] + point.F[LBM_1] + point.F[LBM_3]
							+ 2 * (point.F[LBM_2] + point.F[LBM_6] + point.F[LBM_5]),
							tail = loc_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallSpeed] / 6.0,
							roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_4] = point.F[LBM_2];
						point.F[LBM_7] = point.F[LBM_5] * roughness + point.F[LBM_6] * (1.0 - roughness) - tail;
						point.F[LBM_8] = point.F[LBM_6] * roughness + point.F[LBM_5] * (1.0 - roughness) + tail;
					}
				}
				// FDBC_Period
				static void period_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_x_down = external_physical_field::lbm_field(x - 1, y, z);
						point.F[LBM_1] = point_x_down.F[LBM_1];
						point.F[LBM_5] = point_x_down.F[LBM_5];
						point.F[LBM_8] = point_x_down.F[LBM_8];
					}
				}
				static void period_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_x_up = external_physical_field::lbm_field(x + 1, y, z);
						point.F[LBM_3] = point_x_up.F[LBM_3];
						point.F[LBM_7] = point_x_up.F[LBM_7];
						point.F[LBM_6] = point_x_up.F[LBM_6];
					}
				}
				static void period_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_y_down = external_physical_field::lbm_field(x, y - 1, z);
						point.F[LBM_2] = point_y_down.F[LBM_2];
						point.F[LBM_5] = point_y_down.F[LBM_5];
						point.F[LBM_6] = point_y_down.F[LBM_6];
					}
				}
				static void period_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_y_up = external_physical_field::lbm_field(x, y - 1, z);
						point.F[LBM_4] = point_y_up.F[LBM_4];
						point.F[LBM_7] = point_y_up.F[LBM_7];
						point.F[LBM_8] = point_y_up.F[LBM_8];
					}
				}
				// FDBC_Free
				static void free_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_x_up = external_physical_field::lbm_field(x + 1, y, z);
						point.F[LBM_1] = point_x_up.F[LBM_1];
						point.F[LBM_5] = point_x_up.F[LBM_5];
						point.F[LBM_8] = point_x_up.F[LBM_8];
					}
				}
				static void free_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_x_down = external_physical_field::lbm_field(x - 1, y, z);
						point.F[LBM_3] = point_x_down.F[LBM_3];
						point.F[LBM_7] = point_x_down.F[LBM_7];
						point.F[LBM_6] = point_x_down.F[LBM_6];
					}
				}
				static void free_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_y_up = external_physical_field::lbm_field(x, y - 1, z);
						point.F[LBM_2] = point_y_up.F[LBM_2];
						point.F[LBM_5] = point_y_up.F[LBM_5];
						point.F[LBM_6] = point_y_up.F[LBM_6];
					}
				}
				static void free_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& point_y_down = external_physical_field::lbm_field(x, y - 1, z);
						point.F[LBM_4] = point_y_down.F[LBM_4];
						point.F[LBM_7] = point_y_down.F[LBM_7];
						point.F[LBM_8] = point_y_down.F[LBM_8];
					}
				}
				// FDBC_Pressure
				static void pressure_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double pu = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] - (point.F[LBM_0] + point.F[LBM_2] +
							point.F[LBM_4] + 2.0 * point.F[LBM_6] + 2.0 * point.F[LBM_7]),
							diff = (point.F[LBM_4] - point.F[LBM_2]) / 2.0;
						point.F[LBM_1] = point.F[LBM_3] + 2.0 / 3.0 * pu;
						point.F[LBM_5] = point.F[LBM_7] + diff + pu / 6.0;
						point.F[LBM_8] = point.F[LBM_6] - diff + pu / 6.0;
					}
				}
				static void pressure_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double pu = point.F[LBM_0] + point.F[LBM_2] +
							point.F[LBM_4] + 2.0 * point.F[LBM_5] + 2.0 * point.F[LBM_8] - fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue],
							diff = (point.F[LBM_4] - point.F[LBM_2]) / 2.0;
						point.F[LBM_3] = point.F[LBM_1] - 2.0 / 3.0 * pu;
						point.F[LBM_6] = point.F[LBM_8] + diff - pu / 6.0;
						point.F[LBM_7] = point.F[LBM_5] - diff - pu / 6.0;
					}
				}
				static void pressure_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double pu = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] - (point.F[LBM_0] + point.F[LBM_1] +
							point.F[LBM_3] + 2.0 * point.F[LBM_7] + 2.0 * point.F[LBM_8]),
							diff = (point.F[LBM_3] - point.F[LBM_1]) / 2.0;
						point.F[LBM_2] = point.F[LBM_4] + 2.0 / 3.0 * pu;
						point.F[LBM_5] = point.F[LBM_7] + diff + pu / 6.0;
						point.F[LBM_6] = point.F[LBM_8] - diff + pu / 6.0;
					}
				}
				static void pressure_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double pu = (point.F[LBM_0] + point.F[LBM_1] +
							point.F[LBM_3] + 2.0 * point.F[LBM_5] + 2.0 * point.F[LBM_6]) - fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue],
							diff = (point.F[LBM_3] - point.F[LBM_1]) / 2.0;
						point.F[LBM_4] = point.F[LBM_2] - 2.0 / 3.0 * pu;
						point.F[LBM_7] = point.F[LBM_5] - diff - pu / 6.0;
						point.F[LBM_8] = point.F[LBM_6] + diff - pu / 6.0;
					}
				}
				// FDBC_Normal_Micro_Flow
				static void normal_micro_flow_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double local_density = (point.F[LBM_0] + point.F[LBM_2] + point.F[LBM_4]
							+ 2 * (point.F[LBM_3] + point.F[LBM_6] + point.F[LBM_7])) / (1.0 - fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						point.F[LBM_1] = point.F[LBM_3] + 2.0 / 3.0 * tail;
						point.F[LBM_5] = point.F[LBM_7] + tail / 6.0;
						point.F[LBM_8] = point.F[LBM_6] + tail / 6.0;
					}
				}
				static void normal_micro_flow_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double local_density = (point.F[LBM_0] + point.F[LBM_2] + point.F[LBM_4]
							+ 2 * (point.F[LBM_1] + point.F[LBM_5] + point.F[LBM_8])) / (1.0 - fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						point.F[LBM_3] = point.F[LBM_1] - 2.0 / 3.0 * tail;
						point.F[LBM_6] = point.F[LBM_8] - tail / 6.0;
						point.F[LBM_7] = point.F[LBM_5] - tail / 6.0;
					}
				}
				static void normal_micro_flow_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double local_density = (point.F[LBM_0] + point.F[LBM_1] + point.F[LBM_3]
							+ 2 * (point.F[LBM_4] + point.F[LBM_7] + point.F[LBM_8])) / (1.0 - fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						point.F[LBM_2] = point.F[LBM_4] + 2.0 / 3.0 * tail;
						point.F[LBM_5] = point.F[LBM_7] + tail / 6.0;
						point.F[LBM_6] = point.F[LBM_8] + tail / 6.0;
					}
				}
				static void normal_micro_flow_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double local_density = (point.F[LBM_0] + point.F[LBM_1] + point.F[LBM_3]
							+ 2 * (point.F[LBM_2] + point.F[LBM_5] + point.F[LBM_6])) / (1.0 - fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed]),
							tail = local_density * fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed];
						point.F[LBM_4] = point.F[LBM_2] - 2.0 / 3.0 * tail;
						point.F[LBM_7] = point.F[LBM_5] - tail / 6.0;
						point.F[LBM_8] = point.F[LBM_6] - tail / 6.0;
					}
				}
			}

			static void d3q19_fluid_solid_boundary_Guo2002(long long x, long long y, long long z) {
				double f_eq = 0.0, f_neq = 0.0, _tau = tau(x, y, z);
				Vector3 U;
				// LBM_f_0
				LBMPoint& point = external_physical_field::lbm_field(x, y, z);
				point.F[LBM_0] = 0.0;
				LBMPoint* near_point;
				// LBM_f_1
				near_point = &external_physical_field::lbm_field(x + 1, y, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_1] - f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_1] - f_eq_a_virtual_d2q9(LBM_1, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_1] - f_eq_a_virtual_d2q9(LBM_1, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_1] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_1] = 0.0;
				}
				// LBM_f_2
				near_point = &external_physical_field::lbm_field(x, y + 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_2] - f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y + 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_2] - f_eq_a_virtual_d2q9(LBM_2, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_2] - f_eq_a_virtual_d2q9(LBM_2, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_2] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_2] = 0.0;
				}
				// LBM_f_3
				near_point = &external_physical_field::lbm_field(x - 1, y, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_3] - f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_3] - f_eq_a_virtual_d2q9(LBM_3, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_3] - f_eq_a_virtual_d2q9(LBM_3, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_3] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_3] = 0.0;
				}
				// LBM_f_4
				near_point = &external_physical_field::lbm_field(x, y - 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_4] - f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y - 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_4] - f_eq_a_virtual_d2q9(LBM_4, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_4] - f_eq_a_virtual_d2q9(LBM_4, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_4] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_4] = 0.0;
				}
				// LBM_f_5
				near_point = &external_physical_field::lbm_field(x, y, z + 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_5] - f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y, z + 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_5] - f_eq_a_virtual_d2q9(LBM_5, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_5] - f_eq_a_virtual_d2q9(LBM_5, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_5] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_5] = 0.0;
				}
				// LBM_f_6
				near_point = &external_physical_field::lbm_field(x, y, z - 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_6] - f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y, z - 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_6] - f_eq_a_virtual_d2q9(LBM_6, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_6] - f_eq_a_virtual_d2q9(LBM_6, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_6] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_6] = 0.0;
				}
				// LBM_f_7
				near_point = &external_physical_field::lbm_field(x + 1, y + 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_7] - f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y + 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_7] - f_eq_a_virtual_d2q9(LBM_7, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_7] - f_eq_a_virtual_d2q9(LBM_7, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_7] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_7] = 0.0;
				}
				// LBM_f_8
				near_point = &external_physical_field::lbm_field(x - 1, y + 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_8] - f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y + 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_8] - f_eq_a_virtual_d2q9(LBM_8, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_8] - f_eq_a_virtual_d2q9(LBM_8, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_8] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_8] = 0.0;
				}
				// LBM_f_9
				near_point = &external_physical_field::lbm_field(x - 1, y - 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_9, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_9] - f_eq_a_virtual_d2q9(LBM_9, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y - 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_9, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_9] - f_eq_a_virtual_d2q9(LBM_9, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_9] - f_eq_a_virtual_d2q9(LBM_9, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_9] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_9] = 0.0;
				}
				// LBM_f_10
				near_point = &external_physical_field::lbm_field(x + 1, y - 1, z);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_10, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_10] - f_eq_a_virtual_d2q9(LBM_10, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y - 2, z);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_10, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_10] - f_eq_a_virtual_d2q9(LBM_10, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_10] - f_eq_a_virtual_d2q9(LBM_10, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_10] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_10] = 0.0;
				}
				// LBM_f_11
				near_point = &external_physical_field::lbm_field(x + 1, y, z + 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_11, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_11] - f_eq_a_virtual_d2q9(LBM_11, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y, z + 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_11, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_11] - f_eq_a_virtual_d2q9(LBM_11, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_11] - f_eq_a_virtual_d2q9(LBM_11, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_11] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_11] = 0.0;
				}
				// LBM_f_12
				near_point = &external_physical_field::lbm_field(x - 1, y, z + 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_12, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_12] - f_eq_a_virtual_d2q9(LBM_12, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y, z + 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_12, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_12] - f_eq_a_virtual_d2q9(LBM_12, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_12] - f_eq_a_virtual_d2q9(LBM_12, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_12] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_12] = 0.0;
				}
				// LBM_f_13
				near_point = &external_physical_field::lbm_field(x - 1, y, z - 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_13, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_13] - f_eq_a_virtual_d2q9(LBM_13, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x - 2, y, z - 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_13, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_13] - f_eq_a_virtual_d2q9(LBM_13, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_13] - f_eq_a_virtual_d2q9(LBM_13, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_13] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_13] = 0.0;
				}
				// LBM_f_14
				near_point = &external_physical_field::lbm_field(x + 1, y, z - 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_14, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_14] - f_eq_a_virtual_d2q9(LBM_14, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x + 2, y, z - 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_14, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_14] - f_eq_a_virtual_d2q9(LBM_14, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_14] - f_eq_a_virtual_d2q9(LBM_14, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_14] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_14] = 0.0;
				}
				// LBM_f_15
				near_point = &external_physical_field::lbm_field(x, y + 1, z + 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_15, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_15] - f_eq_a_virtual_d2q9(LBM_15, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y + 2, z + 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_15, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_15] - f_eq_a_virtual_d2q9(LBM_15, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_15] - f_eq_a_virtual_d2q9(LBM_15, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_15] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_15] = 0.0;
				}
				// LBM_f_16
				near_point = &external_physical_field::lbm_field(x, y - 1, z + 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_16, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_16] - f_eq_a_virtual_d2q9(LBM_16, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y - 2, z + 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_16, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_16] - f_eq_a_virtual_d2q9(LBM_16, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_16] - f_eq_a_virtual_d2q9(LBM_16, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_16] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_16] = 0.0;
				}
				// LBM_f_17
				near_point = &external_physical_field::lbm_field(x, y - 1, z - 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_17, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_17] - f_eq_a_virtual_d2q9(LBM_17, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y - 2, z - 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_17, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_17] - f_eq_a_virtual_d2q9(LBM_17, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_17] - f_eq_a_virtual_d2q9(LBM_17, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_17] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_17] = 0.0;
				}
				// LBM_f_18
				near_point = &external_physical_field::lbm_field(x, y + 1, z - 1);
				if (near_point->fluid_region >= solid_liquid_interface_threshold) {
					//double qqq = (solid_liquid_interface_threshold - near_point->solid_region) / (point.solid_region - near_point->solid_region);
					double qqq = 1.0 - (solid_liquid_interface_threshold - point.fluid_region) / (near_point->fluid_region - point.fluid_region);
					if (qqq >= q_c) {
						U = (point.velocity + near_point->velocity * (qqq - 1.0)) / qqq;
						f_eq = f_eq_a_virtual_d2q9(LBM_18, near_point->F_MACRO, U);
						f_neq = near_point->F[LBM_18] - f_eq_a_virtual_d2q9(LBM_18, near_point->F_MACRO, near_point->velocity);
					}
					else {
						LBMPoint& near_near_point = external_physical_field::lbm_field.at2(x, y + 2, z - 2);
						U = point.velocity + near_point->velocity * (qqq - 1.0)
							+ (point.velocity * 2.0 + near_near_point.velocity * (qqq - 1.0)) * (1.0 - qqq) / (1.0 + qqq);
						f_eq = f_eq_a_virtual_d2q9(LBM_18, near_point->F_MACRO, U);
						f_neq = (near_point->F[LBM_18] - f_eq_a_virtual_d2q9(LBM_18, near_point->F_MACRO, near_point->velocity)) * qqq
							+ (near_near_point.F[LBM_18] - f_eq_a_virtual_d2q9(LBM_18, near_near_point.F_MACRO, near_near_point.velocity)) * (1.0 - qqq);
					}
					point.F[LBM_18] = f_eq + (1.0 - 1.0 / _tau) * f_neq;
				}
				else {
					point.F[LBM_18] = 0.0;
				}
			}
			namespace lbm_bc_d3q19 {
				// FDBC_Wall_No_Slip
				static void wall_no_slip_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_1] = point.F[LBM_2];
						point.F[LBM_7] = point.F[LBM_9] * roughness + point.F[LBM_8] * (1.0 - roughness);
						point.F[LBM_10] = point.F[LBM_8] * roughness + point.F[LBM_9] * (1.0 - roughness);
						point.F[LBM_11] = point.F[LBM_13] * roughness + point.F[LBM_12] * (1.0 - roughness);
						point.F[LBM_14] = point.F[LBM_12] * roughness + point.F[LBM_13] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_2] = point.F[LBM_1];
						point.F[LBM_8] = point.F[LBM_10] * roughness + point.F[LBM_7] * (1.0 - roughness);
						point.F[LBM_9] = point.F[LBM_7] * roughness + point.F[LBM_10] * (1.0 - roughness);
						point.F[LBM_12] = point.F[LBM_14] * roughness + point.F[LBM_11] * (1.0 - roughness);
						point.F[LBM_13] = point.F[LBM_11] * roughness + point.F[LBM_14] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_3] = point.F[LBM_4];
						point.F[LBM_7] = point.F[LBM_9] * roughness + point.F[LBM_10] * (1.0 - roughness);
						point.F[LBM_8] = point.F[LBM_10] * roughness + point.F[LBM_9] * (1.0 - roughness);
						point.F[LBM_15] = point.F[LBM_17] * roughness + point.F[LBM_16] * (1.0 - roughness);
						point.F[LBM_18] = point.F[LBM_16] * roughness + point.F[LBM_17] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_4] = point.F[LBM_3];
						point.F[LBM_9] = point.F[LBM_7] * roughness + point.F[LBM_8] * (1.0 - roughness);
						point.F[LBM_10] = point.F[LBM_8] * roughness + point.F[LBM_7] * (1.0 - roughness);
						point.F[LBM_16] = point.F[LBM_18] * roughness + point.F[LBM_15] * (1.0 - roughness);
						point.F[LBM_17] = point.F[LBM_15] * roughness + point.F[LBM_18] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_z_down(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_5] = point.F[LBM_6];
						point.F[LBM_11] = point.F[LBM_13] * roughness + point.F[LBM_14] * (1.0 - roughness);
						point.F[LBM_12] = point.F[LBM_14] * roughness + point.F[LBM_13] * (1.0 - roughness);
						point.F[LBM_15] = point.F[LBM_17] * roughness + point.F[LBM_18] * (1.0 - roughness);
						point.F[LBM_16] = point.F[LBM_18] * roughness + point.F[LBM_17] * (1.0 - roughness);
					}
				}
				static void wall_no_slip_z_up(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double roughness = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_WallRoughness];
						point.F[LBM_6] = point.F[LBM_5];
						point.F[LBM_13] = point.F[LBM_11] * roughness + point.F[LBM_12] * (1.0 - roughness);
						point.F[LBM_14] = point.F[LBM_12] * roughness + point.F[LBM_11] * (1.0 - roughness);
						point.F[LBM_17] = point.F[LBM_15] * roughness + point.F[LBM_16] * (1.0 - roughness);
						point.F[LBM_18] = point.F[LBM_16] * roughness + point.F[LBM_15] * (1.0 - roughness);
					}
				}
				// FDBC_Period
				static void period_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x - 1, y, z);
						point.F[LBM_1] = near_point.F[LBM_1];
						point.F[LBM_7] = near_point.F[LBM_7];
						point.F[LBM_10] = near_point.F[LBM_10];
						point.F[LBM_11] = near_point.F[LBM_11];
						point.F[LBM_14] = near_point.F[LBM_14];
					}
				}
				static void period_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x + 1, y, z);
						point.F[LBM_2] = near_point.F[LBM_2];
						point.F[LBM_8] = near_point.F[LBM_8];
						point.F[LBM_9] = near_point.F[LBM_9];
						point.F[LBM_12] = near_point.F[LBM_12];
						point.F[LBM_13] = near_point.F[LBM_13];
					}
				}
				static void period_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y - 1, z);
						point.F[LBM_3] = near_point.F[LBM_3];
						point.F[LBM_7] = near_point.F[LBM_7];
						point.F[LBM_8] = near_point.F[LBM_8];
						point.F[LBM_15] = near_point.F[LBM_15];
						point.F[LBM_18] = near_point.F[LBM_18];
					}
				}
				static void period_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y + 1, z);
						point.F[LBM_4] = near_point.F[LBM_4];
						point.F[LBM_9] = near_point.F[LBM_9];
						point.F[LBM_10] = near_point.F[LBM_10];
						point.F[LBM_16] = near_point.F[LBM_16];
						point.F[LBM_17] = near_point.F[LBM_17];
					}
				}
				static void period_z_down(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y, z - 1);
						point.F[LBM_5] = near_point.F[LBM_5];
						point.F[LBM_11] = near_point.F[LBM_11];
						point.F[LBM_12] = near_point.F[LBM_12];
						point.F[LBM_15] = near_point.F[LBM_15];
						point.F[LBM_16] = near_point.F[LBM_16];
					}
				}
				static void period_z_up(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y, z + 1);
						point.F[LBM_6] = near_point.F[LBM_6];
						point.F[LBM_13] = near_point.F[LBM_13];
						point.F[LBM_14] = near_point.F[LBM_14];
						point.F[LBM_17] = near_point.F[LBM_17];
						point.F[LBM_18] = near_point.F[LBM_18];
					}
				}
				// FDBC_Free
				static void free_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x + 1, y, z);
						point.F[LBM_1] = near_point.F[LBM_1];
						point.F[LBM_7] = near_point.F[LBM_7];
						point.F[LBM_10] = near_point.F[LBM_10];
						point.F[LBM_11] = near_point.F[LBM_11];
						point.F[LBM_14] = near_point.F[LBM_14];
					}
				}
				static void free_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x - 1, y, z);
						point.F[LBM_2] = near_point.F[LBM_2];
						point.F[LBM_8] = near_point.F[LBM_8];
						point.F[LBM_9] = near_point.F[LBM_9];
						point.F[LBM_12] = near_point.F[LBM_12];
						point.F[LBM_13] = near_point.F[LBM_13];
					}
				}
				static void free_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y + 1, z);
						point.F[LBM_3] = near_point.F[LBM_3];
						point.F[LBM_7] = near_point.F[LBM_7];
						point.F[LBM_8] = near_point.F[LBM_8];
						point.F[LBM_15] = near_point.F[LBM_15];
						point.F[LBM_18] = near_point.F[LBM_18];
					}
				}
				static void free_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y - 1, z);
						point.F[LBM_4] = near_point.F[LBM_4];
						point.F[LBM_9] = near_point.F[LBM_9];
						point.F[LBM_10] = near_point.F[LBM_10];
						point.F[LBM_16] = near_point.F[LBM_16];
						point.F[LBM_17] = near_point.F[LBM_17];
					}
				}
				static void free_z_down(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y, z + 1);
						point.F[LBM_5] = near_point.F[LBM_5];
						point.F[LBM_11] = near_point.F[LBM_11];
						point.F[LBM_12] = near_point.F[LBM_12];
						point.F[LBM_15] = near_point.F[LBM_15];
						point.F[LBM_16] = near_point.F[LBM_16];
					}
				}
				static void free_z_up(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						LBMPoint& near_point = external_physical_field::lbm_field(x, y, z - 1);
						point.F[LBM_6] = near_point.F[LBM_6];
						point.F[LBM_13] = near_point.F[LBM_13];
						point.F[LBM_14] = near_point.F[LBM_14];
						point.F[LBM_17] = near_point.F[LBM_17];
						point.F[LBM_18] = near_point.F[LBM_18];
					}
				}
				// FDBC_Pressure
				static void pressure_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locDensity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (point.F[LBM_0]
								+ point.F[LBM_3] + point.F[LBM_4] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_15] + point.F[LBM_16] + point.F[LBM_17] + point.F[LBM_18]
								+ 2.0 * (point.F[LBM_2] + point.F[LBM_8] + point.F[LBM_9] + point.F[LBM_12]
									+ point.F[LBM_13])) / locDensity;
						point.F[LBM_1] = point.F[LBM_2] + locDensity * u / 3.0;
						point.F[LBM_7] = point.F[LBM_9] + locDensity * u / 6.0;
						point.F[LBM_10] = point.F[LBM_8] + locDensity * u / 6.0;
						point.F[LBM_11] = point.F[LBM_13] + locDensity * u / 6.0;
						point.F[LBM_14] = point.F[LBM_12] + locDensity * u / 6.0;
					}
				}
				static void pressure_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locDensity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (point.F[LBM_0]
								+ point.F[LBM_3] + point.F[LBM_4] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_15] + point.F[LBM_16] + point.F[LBM_17] + point.F[LBM_18]
								+ 2.0 * (point.F[LBM_1] + point.F[LBM_7] + point.F[LBM_10] + point.F[LBM_11]
									+ point.F[LBM_14])) / locDensity;
						point.F[LBM_2] = point.F[LBM_1] + locDensity * u / 3.0;
						point.F[LBM_8] = point.F[LBM_10] + locDensity * u / 6.0;
						point.F[LBM_9] = point.F[LBM_7] + locDensity * u / 6.0;
						point.F[LBM_12] = point.F[LBM_14] + locDensity * u / 6.0;
						point.F[LBM_13] = point.F[LBM_11] + locDensity * u / 6.0;
					}
				}
				static void pressure_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locDensity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_11] + point.F[LBM_12] + point.F[LBM_13] + point.F[LBM_14]
								+ 2.0 * (point.F[LBM_4] + point.F[LBM_9] + point.F[LBM_10] + point.F[LBM_16]
									+ point.F[LBM_17])) / locDensity;
						point.F[LBM_3] = point.F[LBM_4] + locDensity * u / 3.0;
						point.F[LBM_7] = point.F[LBM_9] + locDensity * u / 6.0;
						point.F[LBM_8] = point.F[LBM_10] + locDensity * u / 6.0;
						point.F[LBM_15] = point.F[LBM_17] + locDensity * u / 6.0;
						point.F[LBM_18] = point.F[LBM_16] + locDensity * u / 6.0;
					}
				}
				static void pressure_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locDensity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_11] + point.F[LBM_12] + point.F[LBM_13] + point.F[LBM_14]
								+ 2.0 * (point.F[LBM_3] + point.F[LBM_7] + point.F[LBM_8] + point.F[LBM_15]
									+ point.F[LBM_18])) / locDensity;
						point.F[LBM_4] = point.F[LBM_3] + locDensity * u / 3.0;
						point.F[LBM_9] = point.F[LBM_7] + locDensity * u / 6.0;
						point.F[LBM_10] = point.F[LBM_8] + locDensity * u / 6.0;
						point.F[LBM_16] = point.F[LBM_18] + locDensity * u / 6.0;
						point.F[LBM_17] = point.F[LBM_15] + locDensity * u / 6.0;
					}
				}
				static void pressure_z_down(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locDensity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_3] + point.F[LBM_4]
								+ point.F[LBM_7] + point.F[LBM_8] + point.F[LBM_9] + point.F[LBM_10]
								+ 2.0 * (point.F[LBM_6] + point.F[LBM_13] + point.F[LBM_14] + point.F[LBM_17]
									+ point.F[LBM_18])) / locDensity;
						point.F[LBM_5] = point.F[LBM_6] + locDensity * u / 3.0;
						point.F[LBM_11] = point.F[LBM_13] + locDensity * u / 6.0;
						point.F[LBM_12] = point.F[LBM_14] + locDensity * u / 6.0;
						point.F[LBM_15] = point.F[LBM_17] + locDensity * u / 6.0;
						point.F[LBM_16] = point.F[LBM_18] + locDensity * u / 6.0;
					}
				}
				static void pressure_z_up(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locDensity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_DensityValue],
							u = 1.0 - (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_3] + point.F[LBM_4]
								+ point.F[LBM_7] + point.F[LBM_8] + point.F[LBM_9] + point.F[LBM_10]
								+ 2.0 * (point.F[LBM_5] + point.F[LBM_11] + point.F[LBM_12] + point.F[LBM_15]
									+ point.F[LBM_16])) / locDensity;
						point.F[LBM_6] = point.F[LBM_5] + locDensity * u / 3.0;
						point.F[LBM_13] = point.F[LBM_11] + locDensity * u / 6.0;
						point.F[LBM_14] = point.F[LBM_12] + locDensity * u / 6.0;
						point.F[LBM_17] = point.F[LBM_15] + locDensity * u / 6.0;
						point.F[LBM_18] = point.F[LBM_16] + locDensity * u / 6.0;
					}
				}
				// FDBC_Normal_Micro_Flow
				static void normal_micro_flow_x_down(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locVelocity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (point.F[LBM_0]
								+ point.F[LBM_3] + point.F[LBM_4] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_15] + point.F[LBM_16] + point.F[LBM_17] + point.F[LBM_18]
								+ 2.0 * (point.F[LBM_2] + point.F[LBM_8] + point.F[LBM_9] + point.F[LBM_12]
									+ point.F[LBM_13])) / (1.0 - locVelocity);
						point.F[LBM_1] = point.F[LBM_2] + locVelocity * locDensity / 3.0;
						point.F[LBM_7] = point.F[LBM_9] + locVelocity * locDensity / 6.0;
						point.F[LBM_10] = point.F[LBM_8] + locVelocity * locDensity / 6.0;
						point.F[LBM_11] = point.F[LBM_13] + locVelocity * locDensity / 6.0;
						point.F[LBM_14] = point.F[LBM_12] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_x_up(long long x, long long y, long long z) {
					if (x == external_physical_field::lbm_field.COMP_X_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locVelocity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (point.F[LBM_0]
								+ point.F[LBM_3] + point.F[LBM_4] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_15] + point.F[LBM_16] + point.F[LBM_17] + point.F[LBM_18]
								+ 2.0 * (point.F[LBM_1] + point.F[LBM_7] + point.F[LBM_10] + point.F[LBM_11]
									+ point.F[LBM_14])) / (1.0 - locVelocity);
						point.F[LBM_2] = point.F[LBM_1] + locVelocity * locDensity / 3.0;
						point.F[LBM_8] = point.F[LBM_10] + locVelocity * locDensity / 6.0;
						point.F[LBM_9] = point.F[LBM_7] + locVelocity * locDensity / 6.0;
						point.F[LBM_12] = point.F[LBM_14] + locVelocity * locDensity / 6.0;
						point.F[LBM_13] = point.F[LBM_11] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_y_down(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locVelocity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_11] + point.F[LBM_12] + point.F[LBM_13] + point.F[LBM_14]
								+ 2.0 * (point.F[LBM_4] + point.F[LBM_9] + point.F[LBM_10] + point.F[LBM_16]
									+ point.F[LBM_17])) / (1.0 - locVelocity);
						point.F[LBM_3] = point.F[LBM_4] + locVelocity * locDensity / 3.0;
						point.F[LBM_7] = point.F[LBM_9] + locVelocity * locDensity / 6.0;
						point.F[LBM_8] = point.F[LBM_10] + locVelocity * locDensity / 6.0;
						point.F[LBM_15] = point.F[LBM_17] + locVelocity * locDensity / 6.0;
						point.F[LBM_18] = point.F[LBM_16] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_y_up(long long x, long long y, long long z) {
					if (y == external_physical_field::lbm_field.COMP_Y_END()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locVelocity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_5] + point.F[LBM_6]
								+ point.F[LBM_11] + point.F[LBM_12] + point.F[LBM_13] + point.F[LBM_14]
								+ 2.0 * (point.F[LBM_3] + point.F[LBM_7] + point.F[LBM_8] + point.F[LBM_15]
									+ point.F[LBM_18])) / (1.0 - locVelocity);
						point.F[LBM_4] = point.F[LBM_3] + locVelocity * locDensity / 3.0;
						point.F[LBM_9] = point.F[LBM_7] + locVelocity * locDensity / 6.0;
						point.F[LBM_10] = point.F[LBM_8] + locVelocity * locDensity / 6.0;
						point.F[LBM_16] = point.F[LBM_18] + locVelocity * locDensity / 6.0;
						point.F[LBM_17] = point.F[LBM_15] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_z_down(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locVelocity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_3] + point.F[LBM_4]
								+ point.F[LBM_7] + point.F[LBM_8] + point.F[LBM_9] + point.F[LBM_10]
								+ 2.0 * (point.F[LBM_6] + point.F[LBM_13] + point.F[LBM_14] + point.F[LBM_17]
									+ point.F[LBM_18])) / (1.0 - locVelocity);
						point.F[LBM_5] = point.F[LBM_6] + locVelocity * locDensity / 3.0;
						point.F[LBM_11] = point.F[LBM_13] + locVelocity * locDensity / 6.0;
						point.F[LBM_12] = point.F[LBM_14] + locVelocity * locDensity / 6.0;
						point.F[LBM_15] = point.F[LBM_17] + locVelocity * locDensity / 6.0;
						point.F[LBM_16] = point.F[LBM_18] + locVelocity * locDensity / 6.0;
					}
				}
				static void normal_micro_flow_z_up(long long x, long long y, long long z) {
					if (z == external_physical_field::lbm_field.COMP_Z_BGN()) {
						LBMPoint& point = external_physical_field::lbm_field(x, y, z);
						double locVelocity = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed],
							locDensity = (point.F[LBM_0]
								+ point.F[LBM_1] + point.F[LBM_2] + point.F[LBM_3] + point.F[LBM_4]
								+ point.F[LBM_7] + point.F[LBM_8] + point.F[LBM_9] + point.F[LBM_10]
								+ 2.0 * (point.F[LBM_5] + point.F[LBM_11] + point.F[LBM_12] + point.F[LBM_15]
									+ point.F[LBM_16])) / (1.0 - locVelocity);
						point.F[LBM_6] = point.F[LBM_5] + locVelocity * locDensity / 3.0;
						point.F[LBM_13] = point.F[LBM_11] + locVelocity * locDensity / 6.0;
						point.F[LBM_14] = point.F[LBM_12] + locVelocity * locDensity / 6.0;
						point.F[LBM_17] = point.F[LBM_15] + locVelocity * locDensity / 6.0;
						point.F[LBM_18] = point.F[LBM_16] + locVelocity * locDensity / 6.0;
					}
				}
			}
		}

		void boundary_condition_d2q9(long long x, long long y, long long z) {
			d2q9_domain_boundary_x_down(x, y, z);
			d2q9_domain_boundary_x_up(x, y, z);
			d2q9_domain_boundary_y_down(x, y, z);
			d2q9_domain_boundary_y_up(x, y, z);
			if (external_physical_field::lbm_field(x, y, z).fluid_region < solid_liquid_interface_threshold)
				d2q9_fluid_solid_boundary(x, y, z);
		}

		void boundary_condition_d3q19(long long x, long long y, long long z) {
			d3q19_domain_boundary_x_down(x, y, z);
			d3q19_domain_boundary_x_up(x, y, z);
			d3q19_domain_boundary_y_down(x, y, z);
			d3q19_domain_boundary_y_up(x, y, z);
			d3q19_domain_boundary_z_down(x, y, z);
			d3q19_domain_boundary_z_up(x, y, z);
			if (external_physical_field::lbm_field(x, y, z).fluid_region < solid_liquid_interface_threshold)
				d3q19_fluid_solid_boundary(x, y, z);
		}

		void cal_fluid_domain() {
#pragma omp parallel for
			for (long long x = 0; x < external_physical_field::lbm_field.Nx(); x++)
				for (long long y = 0; y < external_physical_field::lbm_field.Ny(); y++)
					for (long long z = 0; z < external_physical_field::lbm_field.Nz(); z++) {
						LBMPoint& fluid_point = external_physical_field::lbm_field(x, y, z);
						Matrix1D<REAL>& phi_point = main_field::phase_field(x, y, z);
						fluid_point.fluid_region = 0.0;
						for (size_t index = 0; index < main_field::phi_number; index++)
							if (is_solid_phases[index])
								fluid_point.fluid_region += phi_point[index];
						if (fluid_point.fluid_region > 1.0)
							fluid_point.fluid_region = 1.0;
						else if (fluid_point.fluid_region < 0.0)
							fluid_point.fluid_region = 0.0;
						fluid_point.fluid_region = 1.0 - fluid_point.fluid_region;
					}
		}

		void init(LBM& fluid_lbm_solver) {
			PhiProperties::instance().init();
			tau = tau_const;
			viscosity = viscosity_one_phase;
			density = density_one_phase;
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			bool is_solid_phase_in_simulation = false;
			fluid_boundary_condition.resize(Fluid_Boundary_Condition::FBC_SIZE, std::vector<REAL>(Fluid_Boundary_Property::FBP_SIZE, 0));
			WriteDebugFile("# Postprocess.FluidDynamics.LatticeBoltzmann.solid_phases = (phase_name, ... ) \n");
			std::string fluid_phase_key = "Postprocess.FluidDynamics.LatticeBoltzmann.solid_phases", fluid_phase_input = "()";
			infile_reader::read_string_value(fluid_phase_key, fluid_phase_input, true);
			std::vector<input_value> fluid_phase_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, fluid_phase_key, fluid_phase_input, true);
			is_solid_phases.resize(main_field::phi_number, false);
			for (auto fluid_name = fluid_phase_value.begin(); fluid_name < fluid_phase_value.end(); fluid_name++) {
				size_t property = 0;
				if (PhiProperties::instance().is_phi_property(fluid_name->string_value))
					property = PhiProperties::instance().phi_property(fluid_name->string_value);
				else {
					WriteDebugFile("# ERROR , phase name in Postprocess.FluidDynamics.LatticeBoltzmann.solid_phases is not defined in phi property ! \n");
					SYS_PROGRAM_STOP;
				}
				is_solid_phase_in_simulation = true;
				for (size_t index = 0; index < main_field::phi_number; index++)
					if (PhiProperties::instance().phi_property(index) == property)
						is_solid_phases[index] = true;
			}
			WriteDebugFile("# tau = viscosity / fluid_dt / Cs2 + 0.5 \n");
			infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.liquid_viscosity", viscosity_liquid, true);
			infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.liquid_density", density_liquid, true);
			_tau_const = viscosity_liquid / time_parameters::delt_t / Cs2 + 0.5;

			if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				fluid_lbm_solver._boundary_condition = boundary_condition_d2q9;
				// load boundary condition
				if (is_solid_phase_in_simulation)
					d2q9_fluid_solid_boundary = bc_funcs::d2q9_fluid_solid_boundary_Guo2002;
				else
					d2q9_fluid_solid_boundary = bc_funcs::default_domain_boundary_condition;
				WriteDebugFile("# .LatticeBoltzmann.boundary_condition = (down_x,up_x,down_y,up_y) \n");
				WriteDebugFile("#                                        0 - Wall, 1 - Period, 2 - Free, 3 - Pressure, 4 - Normal_Flow \n");
				WriteDebugFile("#                            .pressure = p0 , density0 = p0 / Cs^2 , Cs = 1 / sqrt(3) \n");
				std::string bc_key = "Postprocess.FluidDynamics.LatticeBoltzmann.boundary_condition", bc_input = "(0,0,0,0)";
				infile_reader::read_string_value(bc_key, bc_input, true);
				std::vector<input_value> bc_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, bc_key, bc_input, true);
				switch (Fluid_Domain_Boundary_Condition(bc_value[0].int_value)) // down_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed] = 0.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness], true);
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.wall_speed", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed], true);
					if (isTwoREALEquality(fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::wall_no_slip_x_down;
					else
						d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::wall_slip_x_down;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::period_x_down;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::free_x_down;
					break;
				case FDBC_Pressure:
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::pressure_x_down;
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d2q9_domain_boundary_x_down = bc_funcs::lbm_bc_d2q9::normal_micro_flow_x_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[1].int_value)) // up_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallSpeed] = 0.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness], true);
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.wall_speed", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallSpeed], true);
					if (isTwoREALEquality(fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::wall_no_slip_x_up;
					else
						d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::wall_slip_x_up;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::period_x_up;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::free_x_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::pressure_x_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d2q9_domain_boundary_x_up = bc_funcs::lbm_bc_d2q9::normal_micro_flow_x_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[2].int_value)) // down_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed] = 0.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness], true);
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.wall_speed", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed], true);
					if (isTwoREALEquality(fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::wall_no_slip_y_down;
					else
						d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::wall_slip_y_down;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::period_y_down;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::free_y_down;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::pressure_y_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d2q9_domain_boundary_y_down = bc_funcs::lbm_bc_d2q9::normal_micro_flow_y_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[3].int_value)) // up_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallSpeed] = 0.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness], true);
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.wall_speed", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallSpeed], true);
					if (isTwoREALEquality(fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallSpeed], 0.0))
						d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::wall_no_slip_y_up;
					else
						d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::wall_slip_y_up;
					break;
				case FDBC_Period:
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::period_y_up;
					break;
				case FDBC_Free:
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::free_y_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::pressure_y_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d2q9_domain_boundary_y_up = bc_funcs::lbm_bc_d2q9::normal_micro_flow_y_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
			}
			else if (fluid_lbm_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				fluid_lbm_solver._boundary_condition = boundary_condition_d3q19;
				if (is_solid_phase_in_simulation)
					d3q19_fluid_solid_boundary = bc_funcs::d3q19_fluid_solid_boundary_Guo2002;
				else
					d3q19_fluid_solid_boundary = bc_funcs::default_domain_boundary_condition;
				// load boundary condition
				WriteDebugFile("# .LatticeBoltzmann.boundary_condition = (down_x,up_x,down_y,up_y,down_z,up_z) \n");
				WriteDebugFile("#                                        0 - Wall, 1 - Period, 2 - Free, 3 - Pressure, 4 - Normal_Flow \n");
				WriteDebugFile("#                            .pressure = p0 , density0 = p0 / Cs^2 , Cs = 1 / sqrt(3) \n");
				std::string bc_key = "Postprocess.FluidDynamics.LatticeBoltzmann.boundary_condition", bc_input = "(0,0,0,0,0,0)";
				infile_reader::read_string_value(bc_key, bc_input, true);
				std::vector<input_value> bc_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, bc_key, bc_input, true);
				switch (Fluid_Domain_Boundary_Condition(bc_value[0].int_value)) // down_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_WallRoughness], true);
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::wall_no_slip_x_down;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::period_x_down;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::free_x_down;
					break;
				case FDBC_Pressure:
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::pressure_x_down;
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d3q19_domain_boundary_x_down = bc_funcs::lbm_bc_d3q19::normal_micro_flow_x_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_X.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[1].int_value)) // up_x
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_WallRoughness], true);
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::wall_no_slip_x_up;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::period_x_up;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::free_x_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::pressure_x_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d3q19_domain_boundary_x_up = bc_funcs::lbm_bc_d3q19::normal_micro_flow_x_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_X.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_X][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[2].int_value)) // down_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_WallRoughness], true);
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::wall_no_slip_y_down;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::period_y_down;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::free_y_down;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::pressure_y_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d3q19_domain_boundary_y_down = bc_funcs::lbm_bc_d3q19::normal_micro_flow_y_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Y.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[3].int_value)) // up_y
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_WallRoughness], true);
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::wall_no_slip_y_up;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::period_y_up;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::free_y_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::pressure_y_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d3q19_domain_boundary_y_up = bc_funcs::lbm_bc_d3q19::normal_micro_flow_y_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Y.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Y][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[4].int_value)) // down_z
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Z.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_WallRoughness], true);
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::wall_no_slip_z_down;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::period_z_down;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::free_z_down;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::pressure_z_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Z.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d3q19_domain_boundary_z_down = bc_funcs::lbm_bc_d3q19::normal_micro_flow_z_down;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Down_Z.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_DOWN_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
				switch (Fluid_Domain_Boundary_Condition(bc_value[5].int_value)) // up_z
				{
				case FDBC_Wall:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_WallRoughness] = 1.0;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Z.wall_roughness", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_WallRoughness], true);
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::wall_no_slip_z_up;
					break;
				case FDBC_Period:
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::period_z_up;
					break;
				case FDBC_Free:
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::free_z_up;
					break;
				case FDBC_Pressure:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_DensityValue] = 1.0;
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::pressure_z_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Z.pressure", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_DensityValue], true);
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_DensityValue] = fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_DensityValue] / Cs2;
					break;
				case FDBC_Normal_Flow:
					fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed] = 0.0;
					d3q19_domain_boundary_z_up = bc_funcs::lbm_bc_d3q19::normal_micro_flow_z_up;
					infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.BC_Up_Z.normal_velocity", fluid_boundary_condition[Fluid_Boundary_Condition::FBC_UP_Z][Fluid_Boundary_Property::FBP_NormalFlowSpeed], true);
					break;
				default:
					break;
				}
			}
		}

		void init_two_phase_solver(LBM& field_lbm_two_phase_solver) {
			// tau = tau_standard;
			// tau_two_phase = tau_two_phase_const;
			// viscosity = viscosity_two_phase;
			// density = density_two_phase;
			// WriteDebugFile("# tau_two_phase = Mobility / fluid_dt / Cs2 + 0.5 \n");
			// infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.Mobility", mobility_two_phase, true);
			// _tau_two_phase = mobility_two_phase / PCT_dt / Cs2 + 0.5;
			// infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.gas_viscosity", viscosity_gas, true);
			// infile_reader::read_real_value("Postprocess.FluidDynamics.LatticeBoltzmann.TwoPhaseFLow.gas_density", density_gas, true);
			// 
			// if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
			// 	field_lbm_two_phase_solver._boundary_condition = boundary_condition_d2q9;
			// }
			// else if (field_lbm_two_phase_solver.lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
			// 	field_lbm_two_phase_solver._boundary_condition = boundary_condition_d3q19;
			// }
		}

		void lbm_properties_automatically_change() {
			double cc = mesh_parameters::delt_r / time_parameters::delt_t;
			Cs2 = cc * cc / 3.0;
			Cs4 = cc * cc / 9.0;
			_tau_const = viscosity_liquid / time_parameters::delt_t / Cs2 + 0.5;
			_tau_two_phase = mobility_two_phase / time_parameters::delt_t / Cs2 + 0.5;
		}

		void deinit() {
			fluid_boundary_condition.clear();
			tau = nullptr;
			viscosity = nullptr;
			density = nullptr;
			d2q9_domain_boundary_x_down = nullptr;
			d2q9_domain_boundary_x_up = nullptr;
			d2q9_domain_boundary_y_down = nullptr;
			d2q9_domain_boundary_y_up = nullptr;
			d3q19_domain_boundary_x_down = nullptr;
			d3q19_domain_boundary_x_up = nullptr;
			d3q19_domain_boundary_y_down = nullptr;
			d3q19_domain_boundary_y_up = nullptr;
			d3q19_domain_boundary_z_down = nullptr;
			d3q19_domain_boundary_z_up = nullptr;
		}
	}
}
