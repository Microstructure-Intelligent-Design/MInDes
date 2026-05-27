#pragma once
#include <random>
#include <vector>
#include <string>
#include "GeometryStructure.h"
namespace pf {
	namespace porous_structure {
		// - 
		inline bool is_porous = false;
		// -
		inline bool is_porous_rand = true;
		inline int porous_rand_seed = 0;
		inline double porosity = 2;
		inline double porous_init_noise = 0;
		inline size_t porous_first_phi_index = 0;
		inline size_t porous_second_phi_index = 0;
		inline Matrix1D<REAL> porous_first_con;
		inline Matrix1D<REAL> porous_second_con;
		inline double porous_first_temperature = 0;
		inline double porous_second_temperature = 0;
		inline bool is_porous_normalized = true;
		inline bool is_phi = false;
		inline bool is_con = false;
		inline bool is_temp = false;
		// -
		inline double porous_TwoD_d1 = double(0.05), porous_TwoD_d5 = double(0.0125);
		inline double porous_ThreeD_d1 = double(0.02), porous_ThreeD_d7 = porous_ThreeD_d1 / 2, porous_ThreeD_d19 = porous_ThreeD_d1 / 8;
		// - 
		void quartet_structure_generation(size_t NX, size_t NY, size_t NZ);

		// - aim_phi storage aim phis' fraction , which need to be set zero first.
		void quartet_structure_generation_in_phis(size_t NX, size_t NY, size_t NZ, std::vector<std::vector<std::vector<double>>>& aim_phi);
	}
}