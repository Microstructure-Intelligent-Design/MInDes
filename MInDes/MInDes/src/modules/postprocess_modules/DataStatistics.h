#pragma once
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
#include "../../MainIterator_Params.h"
#include "../model_modules/Model_Params.h"
namespace pf {
	namespace data_statistics_functions {

		inline std::vector<REAL> statistical_phi() {
			std::vector<REAL> phi_info(main_field::phi_number, 0);
			size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
				for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
					for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
						std::vector<REAL>& phi = main_field::phase_field(x, y, z);
						for (size_t index = 0; index < main_field::phi_number; index++)
							phi_info[index] += phi[index];
					}
			for (size_t index = 0; index < main_field::phi_number; index++)
				phi_info[index] /= SIZE;
			return phi_info;
		}

		inline std::vector<REAL> statistical_con() {
			std::vector<REAL> con_info(main_field::con_number, 0);
			size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
				for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
					for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
						std::vector<REAL>& con = main_field::concentration_field(x, y, z);
						for (size_t index = 0; index < main_field::con_number; index++)
							con_info[index] += con[index];
					}
			for (size_t index = 0; index < main_field::con_number; index++)
				con_info[index] /= SIZE;
			return con_info;
		}

		inline REAL statistical_temp() {
			REAL temp = 0;
			size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (long long x = main_field::temperature_field.COMP_X_BGN(); x <= main_field::temperature_field.COMP_X_END(); x++)
				for (long long y = main_field::temperature_field.COMP_Y_BGN(); y <= main_field::temperature_field.COMP_Y_END(); y++)
					for (long long z = main_field::temperature_field.COMP_Z_BGN(); z <= main_field::temperature_field.COMP_Z_END(); z++) {
						temp += main_field::temperature_field(x, y, z);
					}
			return temp / SIZE;
		}



	}

	inline void init() {

	}

}