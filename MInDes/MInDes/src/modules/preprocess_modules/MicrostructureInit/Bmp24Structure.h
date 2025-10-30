#pragma once
#include "GeometryStructure.h"
#include <string>
#include <algorithm>
#include <iostream>
namespace pf {
	namespace bmp24_structure {
		inline bool is_read_bmp24file = false;
		// -
		inline std::string bmp24file_path = "fig.bmp";
		inline int bmp24_layer = 0;
		inline std::vector<std::vector<float>> bmp24_threshold;
		inline std::vector<int> bmp24_phi_index;
		inline std::vector<float> bmp24_phi_value;
		inline std::vector<bool> bmp24_phi_normalized;
		inline std::vector<std::vector<float>> bmp24_con;
		inline std::vector<float> bmp24_temperature;
		// -
		extern "C" MINDES_INIT_API void generate_structure_from_BMP_pic(size_t MESH_NX, size_t MESH_NY, size_t MESH_NZ, 
			bool is_phi_field_on, bool is_con_field_on, bool is_temp_field_on, std::string bmp24file_path, int bmp24_layer, std::vector<std::vector<float>> bmp24_threshold, 
			std::vector<int> bmp24_phi_index, std::vector<float> bmp24_phi_value, std::vector<bool> bmp24_phi_normalized, 
			std::vector<std::vector<float>> bmp24_con, std::vector<float> bmp24_temperature, 
			geometry_structure::NucleationBox& nucleation_box);
	}
}