#pragma once
#include "GeometryStructure.h"
#include <random>
#include <vector>
#include <string>
namespace pf {
	namespace porous_structure {
		// - 
		inline bool is_porous = false;
		// -
		inline bool is_porous_rand = true;
		inline int porous_rand_seed = 0;
		inline float porosity = 2;
		inline float porous_init_noise = 0;
		inline size_t porous_first_phi_index = 0;
		inline size_t porous_second_phi_index = 0;
		inline std::vector<float> porous_first_con;
		inline std::vector<float> porous_second_con;
		inline float porous_first_temperature = 0;
		inline float porous_second_temperature = 0;
		inline bool is_porous_normalized = true;
		// -
		inline float porous_TwoD_d1 = float(0.05), porous_TwoD_d5 = float(0.0125);
		inline float porous_ThreeD_d1 = float(0.02), porous_ThreeD_d7 = porous_ThreeD_d1 / 2, porous_ThreeD_d19 = porous_ThreeD_d1 / 8;
		// - 
		extern "C" MINDES_INIT_API void quartet_structure_generation(size_t NX, size_t NY, size_t NZ, bool is_porous_rand, 
			int porous_rand_seed, float porosity, float porous_init_noise, size_t porous_first_phi_index, size_t porous_second_phi_index,
			std::vector<float> porous_first_con, std::vector<float> porous_second_con, float porous_first_temperature, 
			float porous_second_temperature, bool is_porous_normalized, geometry_structure::NucleationBox& nucleation_box);

		// - aim_phi storage aim phis' fraction , which need to be set zero first.
		extern "C" MINDES_INIT_API void quartet_structure_generation_in_phis(size_t NX, size_t NY, size_t NZ, std::vector<std::vector<std::vector<float>>>& aim_phi,
			bool is_porous_rand, int porous_rand_seed, float porosity, float porous_init_noise, size_t porous_first_phi_index, size_t porous_second_phi_index,
			std::vector<float> porous_first_con, std::vector<float> porous_second_con, float porous_first_temperature,
			float porous_second_temperature, bool is_porous_normalized, geometry_structure::NucleationBox& nucleation_box);
	}
}