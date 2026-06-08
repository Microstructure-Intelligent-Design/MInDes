#include "DDC_Temp_Functions.h"
namespace pf {
	namespace data_driven_complex_model {
		// - default functions
		namespace temperature_functions {

			void cal_mob_temp_grad_lap_7P(size_t x, size_t y, size_t z, Vector3& grad_temp, Vector3& grad_mob, REAL& lap_temp) {
				FIELD_PhiTemp& point = parameters::PhiTemp_field(x, y, z);
				FIELD_PhiTemp& point_upx = parameters::PhiTemp_field(x + 1, y, z);
				FIELD_PhiTemp& point_downx = parameters::PhiTemp_field(x - 1, y, z);
				FIELD_PhiTemp& point_upy = parameters::PhiTemp_field(x, y + 1, z);
				FIELD_PhiTemp& point_downy = parameters::PhiTemp_field(x, y - 1, z);
				FIELD_PhiTemp& point_upz = parameters::PhiTemp_field(x, y, z + 1);
				FIELD_PhiTemp& point_downz = parameters::PhiTemp_field(x, y, z - 1);
				grad_temp[0] = (point_upx.old_temp - point_downx.old_temp) / 2 / mesh_parameters::delt_r;
				grad_temp[1] = (point_upy.old_temp - point_downy.old_temp) / 2 / mesh_parameters::delt_r;
				grad_temp[2] = (point_upz.old_temp - point_downz.old_temp) / 2 / mesh_parameters::delt_r;
				lap_temp =
					(point_upx.old_temp + point_downx.old_temp + point_upy.old_temp + point_downy.old_temp
						+ point_upz.old_temp + point_downz.old_temp - 6 * point.old_temp)
					/ mesh_parameters::delt_r / mesh_parameters::delt_r;
				grad_mob[0] = (point_upx.mob_temp - point_downx.mob_temp) / 2 / mesh_parameters::delt_r;
				grad_mob[1] = (point_upy.mob_temp - point_downy.mob_temp) / 2 / mesh_parameters::delt_r;
				grad_mob[2] = (point_upz.mob_temp - point_downz.mob_temp) / 2 / mesh_parameters::delt_r;
			}

			void cal_mob_temp_grad_lap_19P(size_t x, size_t y, size_t z, Vector3& grad_temp, Vector3& grad_mob, REAL& lap_temp) {
				FIELD_PhiTemp& point = parameters::PhiTemp_field(x, y, z);
				FIELD_PhiTemp& point_upx = parameters::PhiTemp_field(x + 1, y, z);
				FIELD_PhiTemp& point_downx = parameters::PhiTemp_field(x - 1, y, z);
				FIELD_PhiTemp& point_upy = parameters::PhiTemp_field(x, y + 1, z);
				FIELD_PhiTemp& point_downy = parameters::PhiTemp_field(x, y - 1, z);
				FIELD_PhiTemp& point_upz = parameters::PhiTemp_field(x, y, z + 1);
				FIELD_PhiTemp& point_downz = parameters::PhiTemp_field(x, y, z - 1);
				FIELD_PhiTemp& point_upxupy = parameters::PhiTemp_field(x + 1, y + 1, z);
				FIELD_PhiTemp& point_downxdowny = parameters::PhiTemp_field(x - 1, y - 1, z);
				FIELD_PhiTemp& point_upydownx = parameters::PhiTemp_field(x - 1, y + 1, z);
				FIELD_PhiTemp& point_downyupx = parameters::PhiTemp_field(x + 1, y - 1, z);
				FIELD_PhiTemp& point_upxupz = parameters::PhiTemp_field(x + 1, y, z + 1);
				FIELD_PhiTemp& point_downxdownz = parameters::PhiTemp_field(x - 1, y, z - 1);
				FIELD_PhiTemp& point_upzdownx = parameters::PhiTemp_field(x - 1, y, z + 1);
				FIELD_PhiTemp& point_downzupx = parameters::PhiTemp_field(x + 1, y, z - 1);
				FIELD_PhiTemp& point_upzupy = parameters::PhiTemp_field(x, y + 1, z + 1);
				FIELD_PhiTemp& point_downzdowny = parameters::PhiTemp_field(x, y - 1, z - 1);
				FIELD_PhiTemp& point_upydownz = parameters::PhiTemp_field(x, y + 1, z - 1);
				FIELD_PhiTemp& point_downyupz = parameters::PhiTemp_field(x, y - 1, z + 1);
				grad_temp[0] = (point_upx.old_temp - point_downx.old_temp) / 2 / mesh_parameters::delt_r;
				grad_temp[1] = (point_upy.old_temp - point_downy.old_temp) / 2 / mesh_parameters::delt_r;
				grad_temp[2] = (point_upz.old_temp - point_downz.old_temp) / 2 / mesh_parameters::delt_r;
				lap_temp =
					(4 * (point_upx.old_temp + point_downx.old_temp + point_upy.old_temp + point_downy.old_temp
						+ point_upz.old_temp + point_downz.old_temp)
						+ point_upxupy.old_temp + point_downxdowny.old_temp + point_upydownx.old_temp + point_downyupx.old_temp
						+ point_upxupz.old_temp + point_downxdownz.old_temp + point_upzdownx.old_temp + point_downzupx.old_temp
						+ point_upzupy.old_temp + point_downzdowny.old_temp + point_upydownz.old_temp + point_downyupz.old_temp
						- 36 * point.old_temp) / 6 / mesh_parameters::delt_r / mesh_parameters::delt_r;
				grad_mob[0] = (point_upx.mob_temp - point_downx.mob_temp) / 2 / mesh_parameters::delt_r;
				grad_mob[1] = (point_upy.mob_temp - point_downy.mob_temp) / 2 / mesh_parameters::delt_r;
				grad_mob[2] = (point_upz.mob_temp - point_downz.mob_temp) / 2 / mesh_parameters::delt_r;
			}

			REAL temperature_mobility_0(size_t x, size_t y, size_t z) {
				FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
				REAL Mtemp = 0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_var.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_index = field_var.active_index[index];
					Mtemp += field_var.new_phi[phi_index] * parameters::Mtemp[PhiProperties::instance()[phi_index]];
				}
				return Mtemp;
			}

			//========================================================================================================================
			// - evolution equations
			void init_temperature_field() {
#pragma omp parallel for
				for (long long x = main_field::temperature_field.COMP_X_BGN(); x <= main_field::temperature_field.COMP_X_END(); x++)
					for (long long y = main_field::temperature_field.COMP_Y_BGN(); y <= main_field::temperature_field.COMP_Y_END(); y++)
						for (long long z = main_field::temperature_field.COMP_Z_BGN(); z <= main_field::temperature_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_temp = parameters::PhiTemp_field(x, y, z);
							REAL& temp = main_field::temperature_field(x, y, z);
							field_temp.old_temp = temp;
							field_temp.new_temp = temp;
							field_temp.mob_temp = 0;
							for (auto mobi : parameters::mobility_temp)
								field_temp.mob_temp += mobi(x, y, z);
						}
				parameters::PhiTemp_field.init_boundary_condition();
			}

			void pre_calculation_temperature() {
#pragma omp parallel for
				for (long long x = main_field::temperature_field.COMP_X_BGN(); x <= main_field::temperature_field.COMP_X_END(); x++)
					for (long long y = main_field::temperature_field.COMP_Y_BGN(); y <= main_field::temperature_field.COMP_Y_END(); y++)
						for (long long z = main_field::temperature_field.COMP_Z_BGN(); z <= main_field::temperature_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_temp = parameters::PhiTemp_field(x, y, z);
							field_temp.old_temp = field_temp.new_temp;
							field_temp.new_temp = 0;
							field_temp.mob_temp = 0;
							// - mobility calculation
							for (auto mobi : parameters::mobility_temp)
								field_temp.mob_temp += mobi(x, y, z);
						}
				parameters::PhiTemp_field.do_boundary_condition();
#pragma omp parallel for
				for (long long x = main_field::temperature_field.COMP_X_BGN(); x <= main_field::temperature_field.COMP_X_END(); x++)
					for (long long y = main_field::temperature_field.COMP_Y_BGN(); y <= main_field::temperature_field.COMP_Y_END(); y++)
						for (long long z = main_field::temperature_field.COMP_Z_BGN(); z <= main_field::temperature_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_temp = parameters::PhiTemp_field(x, y, z);
							Vector3 grad_temp, grad_mob;
							REAL lap_temp = 0, source_temp = 0;
							parameters::cal_mob_temp_grad_lap(x, y, z, grad_temp, grad_mob, lap_temp);
							for (auto source : parameters::source_temp)
								source_temp += source(x, y, z);
							field_temp.new_temp = grad_temp * grad_mob + field_temp.mob_temp * lap_temp + source_temp;
						}
			}

			REAL solve_temperature() {
				REAL MAX_T_INCREMENT = 0.0;
#pragma omp parallel for
				for (long long x = main_field::temperature_field.COMP_X_BGN(); x <= main_field::temperature_field.COMP_X_END(); x++)
					for (long long y = main_field::temperature_field.COMP_Y_BGN(); y <= main_field::temperature_field.COMP_Y_END(); y++)
						for (long long z = main_field::temperature_field.COMP_Z_BGN(); z <= main_field::temperature_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_temp = parameters::PhiTemp_field(x, y, z);
							field_temp.new_temp = field_temp.old_temp + time_parameters::delt_t * field_temp.new_temp;
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (abs(field_temp.new_temp - field_temp.old_temp) > MAX_T_INCREMENT)
									MAX_T_INCREMENT = abs(field_temp.new_temp - field_temp.old_temp);
							}
#ifdef _DEBUG
							if (_isnan(field_temp.new_temp)) {
								std::cout << "DEBUG: temperature NaN on position x = " << x << " ,y = " << y << " , z = " << z << std::endl;
								SYS_PROGRAM_STOP;
							}
#endif
							main_field::temperature_field(x, y, z) = field_temp.new_temp;
						}
				main_field::temperature_field.do_boundary_condition();
				return MAX_T_INCREMENT;
			}

			//========================================================================================================================
			// - output

		}
	}
}