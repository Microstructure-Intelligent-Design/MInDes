#pragma once
#include "DDCCPAI_Functions.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		void exec_pre_ii() {
			if (parameters::is_phase_volume_statistic) {
				for (size_t phase_index = 0; phase_index < PhiProperties::instance().phi_property_number(); phase_index++)
					data_statistics_functions::statistics_data.push_back
					({ parameters::statistic_phase_volume.first + "_" + PhiProperties::instance().phi_property_name(phase_index), REAL(0) });
			}
			if (parameters::is_region_volume_statistic) {
				for (size_t reg_index = 0; reg_index < ConRegions::instance().region_number(); reg_index++)
					data_statistics_functions::statistics_data.push_back
					({ parameters::statistic_region_volume.first + "_" + ConRegions::instance().region_name(reg_index), REAL(0) });
			}
			if (parameters::is_average_con_statistic) {
				for (size_t con_index = 0; con_index < main_field::con_number; con_index++)
					data_statistics_functions::statistics_data.push_back
					({ parameters::statistic_average_con.first + "_" + ConRegions::instance().con_name(con_index), REAL(0) });
			}
			if (parameters::is_total_con_statistic) {
				for (size_t con_index = 0; con_index < main_field::con_number; con_index++)
					data_statistics_functions::statistics_data.push_back
					({ parameters::statistic_total_con.first + "_" + ConRegions::instance().con_name(con_index), REAL(0) });
			}
		}
		void exec_pre_iii() {
			// - Memory Allocation & Init each field
			parameters::PhiTemp_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, mesh_parameters::delt_r, mesh_parameters::x_down,
				mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			phase_field_functions::init_phi_pair_wise();
			if (main_field::is_con_field_on) {
				parameters::Con_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, mesh_parameters::delt_r, mesh_parameters::x_down,
					mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
				concentration_field_functions::init_total_concentration();
			}
			if (main_field::is_temp_field_on) {
				temperature_functions::init_temperature_field();
			}
		}
		void exec_i() {
			// - phase-field evolution
			phase_field_functions::pre_calculation_phi_pair_wise();
			main_field::PHI_MAX_VARIATION = phase_field_functions::solve_phi_pair_wise();
			// - concentration-field evolution
			if (main_field::is_con_field_on) {
				concentration_field_functions::pre_calculation_total_concentration();
				main_field::CON_MAX_VARIATION = concentration_field_functions::solve_total_concentration();
			}
			// - temperature-field evolution
			if (main_field::is_temp_field_on) {
				temperature_functions::pre_calculation_temperature();
				main_field::TEMP_MAX_VARIATION = temperature_functions::solve_temperature();
			}
		}
		void exec_pos_i() {
			if (show_loop_information::screen_output_step == 0)
				return;
			if (main_iterator::Current_ITE_step % show_loop_information::screen_output_step == 0) {
				std::vector<REAL> phase_volume_info(PhiProperties::instance().phi_property_number(), 0);
				std::vector<REAL> region_volume_info(ConRegions::instance().region_number(), 0);
				std::vector<REAL> average_con_info(main_field::con_number, 0);
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& ptpoint = parameters::PhiTemp_field(x, y, z);
							FIELD_Con& cpoint = parameters::Con_field(x, y, z);
							if (parameters::is_phase_volume_statistic) {
								for (size_t index = 0; index < main_field::phi_number; index++)
									phase_volume_info[PhiProperties::instance()[index]] += ptpoint.new_phi[index];
							}
							if (parameters::is_region_volume_statistic) {
								for (size_t index = 0; index < ConRegions::instance().region_number(); index++)
									region_volume_info[index] += cpoint.new_region[index];
							}
							if (parameters::is_average_con_statistic) {
								for (size_t index = 0; index < main_field::con_number; index++)
									average_con_info[index] += cpoint.new_con[index] * cpoint.new_region[ConRegions::instance().con_region(index)];
							}
						}
				size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
				if (parameters::is_phase_volume_statistic) {
					for (size_t index = 0; index < PhiProperties::instance().phi_property_number(); index++) {
						data_statistics_functions::update_statistics(parameters::statistic_phase_volume.first + "_" + 
							PhiProperties::instance().phi_property_name(index), phase_volume_info[index] / SIZE);
					}
				}
				if (parameters::is_region_volume_statistic) {
					for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
						data_statistics_functions::update_statistics(parameters::statistic_region_volume.first + "_" +
							ConRegions::instance().region_name(index), region_volume_info[index] / SIZE);
					}
				}
				if (parameters::is_average_con_statistic) {
					for (size_t index = 0; index < main_field::con_number; index++) {
						data_statistics_functions::update_statistics(parameters::statistic_average_con.first + "_" + 
							ConRegions::instance().con_name(index), average_con_info[index] / region_volume_info[ConRegions::instance().con_region(index)]);
					}
				}
				if (parameters::is_total_con_statistic) {
					for (size_t index = 0; index < main_field::con_number; index++) {
						data_statistics_functions::update_statistics(parameters::statistic_total_con.first + "_" +
							ConRegions::instance().con_name(index), average_con_info[index] / SIZE);
					}
				}
			}
		}
		void deinit() {
			// - Deallocate Memory
			main_field::phase_field.clear();
			main_field::concentration_field.clear();
			main_field::temperature_field.clear();
			parameters::PhiTemp_field.clear();
			parameters::Con_field.clear();
		}
	}
}