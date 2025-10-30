#include "VoronoiStructure.h"
namespace pf {
	namespace voronoi_structure {

		inline float voronoi_points_distance(geometry_structure::Point point) {
			float points_distance = -1;
			if (voronoi_type == VoronoiType::VDT_CONST) {
				points_distance = voronoi_const_pointsDistance;
			}
			else if (voronoi_type == VoronoiType::VDT_REF_DOT) {
				float pos = voronoi_reference_dot.to_length(point.x, point.y, point.z);
				if (pos > voronoi_reference_dot_distance) {
					points_distance = voronoi_reference_dot_max_pointsDistance;
				}
				else {
					points_distance = pos / voronoi_reference_dot_distance * voronoi_reference_dot_max_pointsDistance + (1 - pos / voronoi_reference_dot_distance) * voronoi_reference_dot_min_pointsDistance;
				}
			}
			else if (voronoi_type == VoronoiType::VDT_REF_SURFACE) {
				float pos = voronoi_reference_surface.distance(point.x, point.y, point.z);
				if (pos > voronoi_reference_surface_distance) {
					points_distance = voronoi_reference_surface_max_pointsDistance;
				}
				else {
					points_distance = pos / voronoi_reference_surface_distance * voronoi_reference_surface_max_pointsDistance + (1 - pos / voronoi_reference_surface_distance) * voronoi_reference_surface_min_pointsDistance;
				}
			}
			else if (voronoi_type == VoronoiType::VDT_DOTS_MATRIX) {
				std::vector<float> dots_lenghts;
				points_distance = 0;
				float total_lenght = 0;
				for (auto dot = voronoi_matrix_dots.begin(); dot < voronoi_matrix_dots.end(); dot++) {
					dots_lenghts.push_back(dot->to_length(point.x, point.y, point.z));
					total_lenght += dot->to_length(point.x, point.y, point.z);
				}
				total_lenght = total_lenght * float(voronoi_matrix_dots.size() - 1);
				if (total_lenght < 1e-6f) {
					std::cout << "ERROR: VDT_DOTS_MATRIX method need more than one dot ! \n";
					exit(0);
				}
				for (auto index = 0; index < voronoi_matrix_dots.size(); index++) {
					float potential = 0;
					for (auto index2 = 0; index2 < voronoi_matrix_dots.size(); index2++)
						if (index != index2)
							potential += dots_lenghts[index2];
					potential /= total_lenght;
					points_distance += potential * voronoi_matrix_dots_pointsDistance[index];
				}
			}
			return points_distance;
		}

		void generate_voronoi_structure() {
			int dimention = 3;
			if (voronoi_box_size[0] == 0 && voronoi_box_size[1] == 0 && voronoi_box_size[2] == 0)
				return;
			else if (voronoi_box_size[0] == 0 || voronoi_box_size[1] == 0 || voronoi_box_size[2] == 0)
				dimention = 2;
			size_t grain_number = voronoi_phi_index_range[1] - voronoi_phi_index_range[0] + 1;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!is_voronoi_rand) {
				gen.seed(voronoi_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			size_t grain_index = 0;
			while (grain_index < grain_number) {
				bool is_point_add = true;
				float rand_x = float(real_dist(gen)), rand_y = float(real_dist(gen)), rand_z = float(real_dist(gen));
				geometry_structure::Point p(rand_x * voronoi_box_size[0] + voronoi_box_position[0],
					rand_y * voronoi_box_size[1] + voronoi_box_position[1],
					rand_z * voronoi_box_size[2] + voronoi_box_position[2]);
				for (auto ip = voronoi_points.begin(); ip < voronoi_points.end(); ip++) {
					float D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
					float distance = voronoi_points_distance(geometry_structure::Point((p.x + ip->x) / 2, (p.y + ip->y) / 2, (p.z + ip->z) / 2));
					if (distance < 0)
						distance = 0;
					if (D2 < (distance * distance)) {
						is_point_add = false;
					}
				}
				if (is_point_add) {
					voronoi_points.push_back(p);
					grain_index++;
					std::string str_report = "> Voronoi: generate point index : " + std::to_string(grain_index) + ", at : (x, y, z) (" 
						+ std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")\n";
					std::cout << str_report;
				}
			}
			// > periodic boundary condition
			std::vector<std::vector<geometry_structure::Point>> mirror_points; // = 27 * points.size()
			size_t region_number = 0;
			if (dimention == 3)
				region_number = 27;
			else
				region_number = 9;
			if (not is_voronoi_mirror_generation)
				region_number = 1;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, 0);
				if (is_voronoi_mirror_generation) {
					mirror_points[1][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, 0);
					mirror_points[2][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, 0);
					mirror_points[3][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), 0);
					mirror_points[4][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), 0);
					mirror_points[5][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[6][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[7][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					mirror_points[8][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					if (dimention == 3) {
						mirror_points[9][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, int(voronoi_box_size[2]));
						mirror_points[10][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, int(-voronoi_box_size[2]));
						mirror_points[11][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[12][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[13][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[14][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[15][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[16][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[17][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[18][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[19][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[20][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[23][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[24][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[26][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[21][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[22][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[25][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
					}
				}
			}
			for (size_t region = 0; region < region_number; region++)
				for (size_t grain = 0; grain < grain_number; grain++) {
					geometry_structure::Polyhedron poly(mirror_points[region][grain]);
					std::vector<geometry_structure::point_in_region_index> record_points;
					record_points.push_back(geometry_structure::point_in_region_index(region, grain));
					for (size_t region_index = 0; region_index < mirror_points.size(); region_index++)
						for (size_t grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
							// Avoid inclusion points
							bool is_point_contained = false;
							for (auto re = record_points.begin(); re < record_points.end(); re++)
								if (re->region == region_index && re->grain_index == grain_index)
									is_point_contained = true;
							if (is_point_contained)
								continue;
							// prepare
							geometry_structure::Vector3 norm = get_vector(mirror_points[region_index][grain_index], poly.point_inside_polyhedron);
							// Vector3 norm = get_vector(poly.point_inside_polyhedron, mirror_points[region_index][grain_index]);
							geometry_structure::Point mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
							// Judged ipsilateral to poly center
							if (poly.check_point_inside_polyhedron(mid_point) == false)
								continue;
							// < add point
							poly.add_surf(norm, mid_point);
							// < Eliminate meaningless points in poly
							for (auto re = record_points.begin(); re < record_points.end();) {
								geometry_structure::Point check = (voronoi_points[grain] + mirror_points[re->region][re->grain_index]) / 2;
								if (poly.check_point_inside_polyhedron(check) == false) {
									re = record_points.erase(re);
								}
								else {
									++re;
								}
							}
						}
					std::string str_report = "> Voronoi: One polyhedron in region : " + std::to_string(region) + ", grain : "
						+ std::to_string(grain) + " has been generated ! \n";
					std::cout << str_report;
					geometry_structure::GeometricRegion geo;
					geo.geometryProperty = geometry_structure::Geometry::Geo_Polyhedron;
					geo.generate_step = 0;
					geo.polyhedron = poly;
					geo.phaseIndex = voronoi_phi_index_range[0] + grain;
					geo.temperature = voronoi_temperature;
					geo.con = voronoi_con;
					geo.phi = 1;
					geo.isNormalized = true;
					geometry_structure::nucleation_box.geometry_box.push_back(geo);
				}

		}

		void generate_voronoi_structure_new_method() {
			int dimention = 3;
			if (voronoi_box_size[0] == 0 && voronoi_box_size[1] == 0 && voronoi_box_size[2] == 0)
				return;
			else if (voronoi_box_size[0] == 0 || voronoi_box_size[1] == 0 || voronoi_box_size[2] == 0)
				dimention = 2;
			size_t grain_number = voronoi_phi_index_range[1] - voronoi_phi_index_range[0] + 1;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!is_voronoi_rand) {
				gen.seed(voronoi_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			size_t grain_index = 0;
			while (grain_index < grain_number) {
				bool is_point_add = true;
				float rand_x = float(real_dist(gen)), rand_y = float(real_dist(gen)), rand_z = float(real_dist(gen));
				geometry_structure::Point p(rand_x * voronoi_box_size[0] + voronoi_box_position[0],
					rand_y * voronoi_box_size[1] + voronoi_box_position[1],
					rand_z * voronoi_box_size[2] + voronoi_box_position[2]);
				for (auto ip = voronoi_points.begin(); ip < voronoi_points.end(); ip++) {
					float D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
					float distance = voronoi_points_distance(geometry_structure::Point((p.x + ip->x) / 2, (p.y + ip->y) / 2, (p.z + ip->z) / 2));
					if (distance < 0)
						distance = 0;
					if (D2 < (distance * distance)) {
						is_point_add = false;
					}
				}
				if (is_point_add) {
					voronoi_points.push_back(p);
					grain_index++;
					std::string str_report = "> Voronoi: generate point index : " + std::to_string(grain_index) + ", at : (x, y, z) (" 
						+ std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")\n";
					std::cout << str_report;
				}
			}
			// > periodic boundary condition
			std::vector<std::vector<geometry_structure::Point>> mirror_points; // = 27 * points.size()
			size_t region_number = 0;
			if (dimention == 3)
				region_number = 27;
			else
				region_number = 9;
			if (not is_voronoi_mirror_generation)
				region_number = 1;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, 0);
				if (is_voronoi_mirror_generation) {
					mirror_points[1][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, 0);
					mirror_points[2][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, 0);
					mirror_points[3][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), 0);
					mirror_points[4][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), 0);
					mirror_points[5][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[6][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[7][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					mirror_points[8][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					if (dimention == 3) {
						mirror_points[9][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, int(voronoi_box_size[2]));
						mirror_points[10][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, int(-voronoi_box_size[2]));
						mirror_points[11][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[12][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[13][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[14][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[15][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[16][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[17][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[18][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[19][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[20][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[23][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[24][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[26][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[21][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[22][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[25][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
					}
				}
			}
			for (size_t region = 0; region < region_number; region++)
				for (size_t grain = 0; grain < grain_number; grain++) {
					geometry_structure::Polyhedron poly(mirror_points[region][grain]);
					std::vector<geometry_structure::point_in_region_index> record_points;
					record_points.push_back(geometry_structure::point_in_region_index(region, grain));
					bool find_new_point = false;
					do
					{
						find_new_point = false;
						float point_distance2 = std::numeric_limits<float>::max();;
						size_t buff_region = 0, buff_grain = 0;
						for (size_t region_index = 0; region_index < mirror_points.size(); region_index++)
							for (size_t grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
								// Avoid inclusion points
								bool is_point_contained = false;
								for (auto re = record_points.begin(); re < record_points.end(); re++)
									if (re->region == region_index && re->grain_index == grain_index)
										is_point_contained = true;
								if (is_point_contained)
									continue;
								// check mid point
								geometry_structure::Point mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
								// Judged ipsilateral to poly center
								if (poly.check_point_inside_polyhedron(mid_point) == false)
									continue;
								find_new_point = true;
								float distance2 = (poly.point_inside_polyhedron.x - mirror_points[region_index][grain_index].x) * (poly.point_inside_polyhedron.x - mirror_points[region_index][grain_index].x)
									+ (poly.point_inside_polyhedron.y - mirror_points[region_index][grain_index].y) * (poly.point_inside_polyhedron.y - mirror_points[region_index][grain_index].y)
									+ (poly.point_inside_polyhedron.z - mirror_points[region_index][grain_index].z) * (poly.point_inside_polyhedron.z - mirror_points[region_index][grain_index].z);
								if (distance2 < point_distance2) {
									point_distance2 = distance2;
									buff_region = region_index;
									buff_grain = grain_index;
								}
							}
						if (find_new_point) {
							geometry_structure::Point mid_point = (poly.point_inside_polyhedron + mirror_points[buff_region][buff_grain]) / 2;
							geometry_structure::Vector3 norm = get_vector(mirror_points[buff_region][buff_grain], poly.point_inside_polyhedron);
							poly.add_surf(norm, mid_point);
							record_points.push_back(geometry_structure::point_in_region_index(buff_region, buff_grain));
						}
					} while (find_new_point);
					std::string str_report = "> Voronoi: One polyhedron in region : " + std::to_string(region) + ", grain : "
						+ std::to_string(grain) + " has been generated ! \n";
					std::cout << str_report;
					geometry_structure::GeometricRegion geo;
					geo.geometryProperty = geometry_structure::Geometry::Geo_Polyhedron;
					geo.generate_step = 0;
					geo.polyhedron = poly;
					geo.phaseIndex = voronoi_phi_index_range[0] + grain;
					geo.temperature = voronoi_temperature;
					geo.con = voronoi_con;
					geo.phi = 1;
					geo.isNormalized = true;
					geometry_structure::nucleation_box.geometry_box.push_back(geo);
				}

		}

		void generate_voronoi_structure_in_phis(std::vector<std::vector<std::vector<float>>>& aim_phi) {
			int dimention = 3;
			if (voronoi_box_size[0] == 0 && voronoi_box_size[1] == 0 && voronoi_box_size[2] == 0)
				return;
			else if (voronoi_box_size[0] == 0 || voronoi_box_size[1] == 0 || voronoi_box_size[2] == 0)
				dimention = 2;
			size_t grain_number = voronoi_phi_index_range[1] - voronoi_phi_index_range[0] + 1;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!is_voronoi_rand) {
				gen.seed(voronoi_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			size_t grain = 0;
			while (grain < grain_number) {
				bool is_point_add = true;
				float rand_x = float(real_dist(gen)), rand_y = float(real_dist(gen)), rand_z = float(real_dist(gen));
				geometry_structure::Point p(rand_x * voronoi_box_size[0] + voronoi_box_position[0],
					rand_y * voronoi_box_size[1] + voronoi_box_position[1],
					rand_z * voronoi_box_size[2] + voronoi_box_position[2]);
				for (auto ip = voronoi_points.begin(); ip < voronoi_points.end(); ip++) {
					float D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
					float distance = voronoi_points_distance(geometry_structure::Point((p.x + ip->x) / 2, (p.y + ip->y) / 2, (p.z + ip->z) / 2));
					if (distance < 0)
						distance = 0;
					if (D2 < (distance * distance)) {
						is_point_add = false;
					}
				}
				if (is_point_add && aim_phi[geometry_structure::REAL_to_int(p.x)]
					[geometry_structure::REAL_to_int(p.y)][geometry_structure::REAL_to_int(p.z)] > 0.5) {
					grain++;
					voronoi_points.push_back(p);
					std::string str_report = "> Voronoi: generate point index : " + std::to_string(grain) + ",  at : (x, y, z) (" 
						+ std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")\n";
					std::cout << str_report;
				}
			}
			// > periodic boundary condition
			std::vector<std::vector<geometry_structure::Point>> mirror_points; // = 27 * points.size()
			size_t region_number = 0;
			if (dimention == 3)
				region_number = 27;
			else
				region_number = 9;
			if (not is_voronoi_mirror_generation)
				region_number = 1;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, 0);
				if (is_voronoi_mirror_generation) {
					mirror_points[1][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, 0);
					mirror_points[2][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, 0);
					mirror_points[3][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), 0);
					mirror_points[4][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), 0);
					mirror_points[5][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[6][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[7][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					mirror_points[8][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					if (dimention == 3) {
						mirror_points[9][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, int(voronoi_box_size[2]));
						mirror_points[10][grain] = voronoi_points[grain] + geometry_structure::Point(0, 0, int(-voronoi_box_size[2]));
						mirror_points[11][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[12][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[13][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[14][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[15][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[16][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[17][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[18][grain] = voronoi_points[grain] + geometry_structure::Point(0, int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[19][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[20][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[23][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[24][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[26][grain] = voronoi_points[grain] + geometry_structure::Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[21][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[22][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[25][grain] = voronoi_points[grain] + geometry_structure::Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
					}
				}
			}
			for (size_t region = 0; region < region_number; region++)
				for (size_t grain = 0; grain < grain_number; grain++) {
					geometry_structure::Polyhedron poly(mirror_points[region][grain]);
					std::vector<geometry_structure::point_in_region_index> record_points;
					geometry_structure::point_in_region_index rp(region, grain);
					record_points.push_back(rp);
					for (size_t region_index = 0; region_index < mirror_points.size(); region_index++)
						for (size_t grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
							// Avoid inclusion points
							bool is_point_contained = false;
							for (auto re = record_points.begin(); re < record_points.end(); re++)
								if (re->region == region_index && re->grain_index == grain_index)
									is_point_contained = true;
							if (is_point_contained)
								continue;
							// prepare
							// Vector3 norm = get_vector(mirror_points[region_index][grain_index], poly.point_inside_polyhedron);
							geometry_structure::Vector3 norm = get_vector(poly.point_inside_polyhedron, mirror_points[region_index][grain_index]);
							geometry_structure::Point mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
							// Judged ipsilateral to poly center
							if (poly.check_point_inside_polyhedron(mid_point) == false)
								continue;
							// < add point
							poly.add_surf(norm, mid_point);
							// < Eliminate meaningless points in poly
							for (auto re = record_points.begin(); re < record_points.end();) {
								geometry_structure::Point check = (voronoi_points[grain] + mirror_points[re->region][re->grain_index]) / 2;
								if (poly.check_point_inside_polyhedron(check) == false) {
									re = record_points.erase(re);
								}
								else {
									++re;
								}
							}
						}
					geometry_structure::PointSet set;
					float sum_sum_phi = 0;
					for (size_t z = 0; z <= voronoi_box_size[2]; z++)
						for (size_t y = 0; y <= voronoi_box_size[1]; y++)
							for (size_t x = 0; x <= voronoi_box_size[0]; x++) {
								if (poly.check_point_inside_polyhedron(geometry_structure::Point(x, y, z)) == true) {
									if (aim_phi[x][y][z] > 1e-6) {
										set.add_point(x, y, z, aim_phi[x][y][z]);
										sum_sum_phi += aim_phi[x][y][z];
									}
								}
							}
					std::string str_report = "> Voronoi: One polyhedron in region : " + std::to_string(region) + ", grain : " 
						+ std::to_string(grain) + " could been generated ! \n";
					std::cout << str_report;
					set.generate_step = 0;
					set.phaseIndex = voronoi_phi_index_range[0] + grain;
					set.temperature = voronoi_temperature;
					set.con = voronoi_con;
					set.is_normalized = false;
					if (sum_sum_phi > 1e-6) {
						geometry_structure::nucleation_box.point_set_box.push_back(set);
					}
				}
		}

	}
}