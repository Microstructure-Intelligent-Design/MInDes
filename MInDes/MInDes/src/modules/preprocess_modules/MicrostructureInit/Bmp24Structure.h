#pragma once
#include <string>
#include <algorithm>
#include <iostream>
#include "../../base/VectorMatrix.h"
#include "GeometryStructure.h"
namespace pf {
	namespace bmp24_structure {
		inline bool is_read_bmp24file = false;
		// -
		inline std::string bmp24file_path = "fig.bmp";
		inline int bmp24_layer = 0;
		inline std::vector<std::vector<double>> bmp24_threshold;
		inline std::vector<int> bmp24_phi_index;
		inline std::vector<double> bmp24_phi_value;
		inline std::vector<bool> bmp24_phi_normalized;
		inline std::vector<bool> is_phi;
		inline std::vector<Matrix1D<REAL>> bmp24_con;
		inline std::vector<bool> is_con;
		inline std::vector<double> bmp24_temperature;
		inline std::vector<bool> is_temp;
		// -
		void generate_structure_from_BMP_pic(size_t MESH_NX, size_t MESH_NY, size_t MESH_NZ, 
			bool is_phi_field_on, bool is_con_field_on, bool is_temp_field_on);
	}
}