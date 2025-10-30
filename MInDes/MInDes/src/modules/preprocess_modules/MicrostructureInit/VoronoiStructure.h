#pragma once
#include "GeometryStructure.h"
#include <string>
#include <iostream>
#include <limits>
namespace pf {
	namespace voronoi_structure {
		// -
		enum VoronoiType { VDT_CONST, VDT_REF_DOT, VDT_REF_SURFACE, VDT_DOTS_MATRIX };
		inline bool is_voronoi = false;
		// -
		inline bool is_voronoi_mirror_generation = true;
		inline bool is_voronoi_rand = true;
		inline int voronoi_rand_seed = 0;
		inline std::vector<int> voronoi_box_position = { 0, 0, 0 };
		inline std::vector<int> voronoi_box_size = { 0, 0, 0 };
		inline std::vector<size_t> voronoi_phi_index_range = { 0, 0 };
		inline std::vector<float> voronoi_con;
		inline float voronoi_temperature = 0.0;
		inline VoronoiType voronoi_type = VoronoiType::VDT_CONST;
		inline std::vector<geometry_structure::Point> voronoi_points;
		inline float voronoi_const_pointsDistance = -1.0;
		inline geometry_structure::Point voronoi_reference_dot;
		inline float voronoi_reference_dot_distance = -1.0;
		inline float voronoi_reference_dot_min_pointsDistance = -1.0;
		inline float voronoi_reference_dot_max_pointsDistance = -1.0;
		inline geometry_structure::surf_func_3D voronoi_reference_surface;
		inline float voronoi_reference_surface_distance = -1.0;
		inline float voronoi_reference_surface_min_pointsDistance = -1.0;
		inline float voronoi_reference_surface_max_pointsDistance = -1.0;
		inline std::vector<geometry_structure::Point> voronoi_matrix_dots;
		inline std::vector<float> voronoi_matrix_dots_pointsDistance;
		// -
		void generate_voronoi_structure();
		void generate_voronoi_structure_new_method();
		// -
		void generate_voronoi_structure_in_phis(std::vector<std::vector<std::vector<float>>>& aim_phi);
	}
}