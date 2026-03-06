#include "PorousStructure.h"
namespace pf {
	namespace porous_structure {

		static void quartet_structure_grow_2D(std::vector<std::vector<size_t>>& Solid, std::mt19937& gen) {
			std::vector<std::vector<size_t>> buff_solid;
			const size_t Nx = Solid.size(), Ny = Solid[0].size();
			buff_solid.resize(Nx);
			for (int i = 0; i < Nx; i++) {
				buff_solid[i].resize(Ny);
				for (int j = 0; j < Ny; j++) {
					buff_solid[i][j] = Solid[i][j];
				}
			}
			/// \brief Grow in eight directions.
			/// Grow directions
			///*****    6    2    5   *****
			///*****                  *****
			///*****    3    C    1   *****
			///*****                  *****
			///*****    7    4    8   *****
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					if (i == 0) {
						if (j == 0) {
							if (Solid[i][j] == 1)
							{
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
							}
						}
						else if (j == Ny - 1) {
							if (Solid[i][j] == 1)
							{
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
						else {
							if (Solid[i][j] == 1)
							{
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
					}
					else if (i == Nx - 1) {
						if (j == 0) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
							}
						}
						else if (j == Ny - 1) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
							}
						}
						else {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
							}
						}
					}
					else {
						if (j == 0) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
							}
						}
						else if (j == Ny - 1) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
						else {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
					}
				}
			}
#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					Solid[i][j] = buff_solid[i][j];
				}
			}
		}

		static void quartet_structure_grow_3D(std::vector<std::vector<std::vector<size_t>>>& arrgrid, std::vector<std::vector<size_t>>& soild, size_t& Tnumsoild, std::mt19937& gen) {
			const size_t NX = arrgrid.size(), NY = arrgrid[0].size(), NZ = arrgrid[0][0].size();
			size_t numsoild = Tnumsoild;
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			for (size_t index_soild = 0; index_soild < Tnumsoild; index_soild++) {
				size_t index_i = soild[index_soild][0];
				size_t index_j = soild[index_soild][1];
				size_t index_k = soild[index_soild][2];
				//1向右方向生长
				if (index_j < NY - 1) {
					size_t i = index_i;
					size_t j = index_j + 1;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//2向后方向生长
				if (index_i > 0) {
					size_t i = index_i - 1;
					size_t j = index_j;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//3向左方向生长
				if (index_j > 0) {
					size_t i = index_i;
					size_t j = index_j - 1;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//4向前方向生长
				if (index_i < NX - 1) {
					size_t i = index_i + 1;
					size_t j = index_j;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//5向上方向生长
				if (index_k < NZ - 1) {
					size_t i = index_i;
					size_t j = index_j;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//6向下方向生长		
				if (index_k > 0) {
					size_t i = index_i;
					size_t j = index_j;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//7向水平方向右前生长
				if (index_i < NX - 1 && index_j < NY - 1) {
					size_t i = index_i + 1;
					size_t j = index_j + 1;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//8向水平方向左前生长
				if (index_i < NX - 1 && index_j > 0) {
					size_t i = index_i + 1;
					size_t j = index_j - 1;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//9向水平方向右后生长
				if (index_i > 0 && index_j < NY - 1) {
					size_t i = index_i - 1;
					size_t j = index_j + 1;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//10向水平方向左后生长
				if (index_i > 0 && index_j > 0) {
					size_t i = index_i - 1;
					size_t j = index_j - 1;
					size_t k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//11向右上方向生长
				if (index_j < NY - 1 && index_k < NZ - 1) {
					size_t i = index_i;
					size_t j = index_j + 1;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//12向右下方向生长
				if (index_j < NY - 1 && index_k >0) {
					size_t i = index_i;
					size_t j = index_j + 1;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//13向左上方向生长
				if (index_j > 0 && index_k < NZ - 1) {
					size_t i = index_i;
					size_t j = index_j - 1;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//14向左下方向生长
				if (index_j > 0 && index_k > 0) {
					size_t i = index_i;
					size_t j = index_j - 1;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//15向前上方向生长
				if (index_i < NX - 1 && index_k < NZ - 1) {
					size_t i = index_i + 1;
					size_t j = index_j;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//16向前下方向生长
				if (index_i < NX - 1 && index_k >0) {
					size_t i = index_i + 1;
					size_t j = index_j;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//17向后上方向生长
				if (index_i > 0 && index_k < NZ - 1) {
					size_t i = index_i - 1;
					size_t j = index_j;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//18向后下方向生长
				if (index_i > 0 && index_k > 0) {
					size_t i = index_i - 1;
					size_t j = index_j;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//19向右前上对角线方向生长
				if (index_i < NX - 1 && index_j < NY - 1 && index_k < NZ - 1) {
					size_t i = index_i + 1;
					size_t j = index_j + 1;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//20向右后上对角线方向生长
				if (index_i > 0 && index_j < NY - 1 && index_k < NZ - 1) {
					size_t i = index_i - 1;
					size_t j = index_j + 1;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//21向左后上对角线方向生长
				if (index_i > 0 && index_j > 0 && index_k < NZ - 1) {
					size_t i = index_i - 1;
					size_t j = index_j - 1;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//22向左前上对角线方向生长
				if (index_i < NX - 1 && index_j>0 && index_k < NZ - 1) {
					size_t i = index_i + 1;
					size_t j = index_j - 1;
					size_t k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//23向右前下对角线方向生长
				if (index_i < NX - 1 && index_j < NY - 1 && index_k>0) {
					size_t i = index_i + 1;
					size_t j = index_j + 1;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//24向右后下对角线方向生长
				if (index_i > 0 && index_j < NY - 1 && index_k>0) {
					size_t i = index_i - 1;
					size_t j = index_j + 1;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//25向左后下对角线方向生长
				if (index_i > 0 && index_j > 0 && index_k > 0) {
					size_t i = index_i - 1;
					size_t j = index_j - 1;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//26向左前下对角线方向生长
				if (index_i < NX - 1 && index_j>0 && index_k > 0) {
					size_t i = index_i + 1;
					size_t j = index_j - 1;
					size_t k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						std::vector<size_t> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
			}
			Tnumsoild = numsoild;
		}

		void quartet_structure_generation(size_t MESH_NX, size_t MESH_NY, size_t MESH_NZ) {
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!is_porous_rand) {
				gen.seed(porous_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			std::cout << "> Quartet structure generation:" << std::endl;
			if (MESH_NZ == 1) {
				std::vector<std::vector<size_t>> Solid;
				Solid.resize(MESH_NX);
				/// No solid in the beginning
				for (size_t i = 0; i < MESH_NX; i++) {
					Solid[i].resize(MESH_NY);
					for (size_t j = 0; j < MESH_NY; j++) {
						Solid[i][j] = 0;
					}
				}
				/// Produce growth core
				double core_p = 0.0;
				for (size_t i = 0; i < MESH_NX; i++) {
					for (size_t j = 0; j < MESH_NY; j++) {
						if (real_dist(gen) < porous_init_noise) {
							Solid[i][j] = 1;
							core_p += 1.0;
						}
					}
				}
				core_p /= MESH_NX * MESH_NY;
				std::string out_ptr = "  percentage of core in 2D = " + std::to_string(core_p);
				std::cout << out_ptr << std::endl;
				///  Produce process
				size_t step = 0;
				do {
					step++;
					quartet_structure_grow_2D(Solid, gen);
					/// Calculate porosity
					core_p = 0.0;
					for (size_t i = 0; i < MESH_NX; i++) {
						for (size_t j = 0; j < MESH_NY; j++) {
							core_p += Solid[i][j];
						}
					}
					core_p = core_p / (MESH_NX * MESH_NY);
					out_ptr = "  grow step: " + std::to_string(step) + ", percentage of particle/pore = " + std::to_string(core_p);
					std::cout << out_ptr << std::endl;
				} while (core_p < porosity);
				out_ptr = "  End of growth ! \n";
				std::cout << out_ptr;
				geometry_structure::PointSet set_phi_1, set_phi_2;
				for (size_t i = 0; i < MESH_NX; i++) {
					for (size_t j = 0; j < MESH_NY; j++) {
						if (Solid[i][j]) {
							set_phi_1.add_point(float(i + 1), float(j + 1), float(1), float(1.0));
						}
						else {
							set_phi_2.add_point(float(i + 1), float(j + 1), float(1), float(1.0));
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = porous_first_phi_index;
				set_phi_1.temperature = porous_first_temperature;
				set_phi_1.con = porous_first_con;
				set_phi_1.is_normalized = is_porous_normalized;
				
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = porous_second_phi_index;
				set_phi_2.temperature = porous_second_temperature;
				set_phi_2.con = porous_second_con;
				set_phi_2.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else if (MESH_NZ > 1) {
				std::vector<std::vector<std::vector<size_t>>> arrgrid; arrgrid.resize(MESH_NX);
				std::vector<std::vector<size_t>> soild;
				size_t numsoild = 0, Nxyz = MESH_NX * MESH_NY * MESH_NZ, numtotal_need = size_t(porosity * Nxyz);
				for (size_t i = 0; i < MESH_NX; i++) {
					arrgrid[i].resize(MESH_NY);
					for (size_t j = 0; j < MESH_NY; j++) {
						arrgrid[i][j].resize(MESH_NZ);
						for (size_t k = 0; k < MESH_NZ; k++) {
							arrgrid[i][j][k] = 0;
							if (real_dist(gen) < porous_init_noise) {
								arrgrid[i][j][k] = 1;
								std::vector<size_t> buff = { i, j, k };
								soild.push_back(buff);
								numsoild = numsoild + 1;
							}
						}
					}
				}
				double core_p = double(numsoild) / double(Nxyz);
				std::string out_ptr = "  percentage of core in 3D = " + std::to_string(core_p) + "\n";
				std::cout << out_ptr;
				size_t Tnumsoild = numsoild, step = 0;
				while (Tnumsoild < numtotal_need) {
					step++;
					quartet_structure_grow_3D(arrgrid, soild, Tnumsoild, gen);
					core_p = double(Tnumsoild) / double(Nxyz);
					out_ptr = "  grow step: " + std::to_string(step) + ", percentage of particle/pore = " + std::to_string(core_p) + "\n";
					std::cout << out_ptr;
				}
				out_ptr = "  End of growth ! \n";
				std::cout << out_ptr;
				geometry_structure::PointSet set_phi_1, set_phi_2;
				for (size_t i = 0; i < MESH_NX; i++) {
					for (size_t j = 0; j < MESH_NY; j++) {
						for (size_t k = 0; k < MESH_NZ; k++) {
							if (arrgrid[i][j][k]) {
								set_phi_1.add_point(float(i + 1), float(j + 1), float(k + 1), float(1.0));
							}
							else {
								set_phi_2.add_point(float(i + 1), float(j + 1), float(k + 1), float(1.0));
							}
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = porous_first_phi_index;
				set_phi_1.temperature = porous_first_temperature;
				set_phi_1.con = porous_first_con;
				set_phi_1.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = porous_second_phi_index;
				set_phi_2.temperature = porous_second_temperature;
				set_phi_2.con = porous_second_con;
				set_phi_2.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else {
				std::cout << "> Quartet structure generation Failed, can't generate at one dimension ! \n";
			}
		}

		void quartet_structure_generation_in_phis(size_t MESH_NX, size_t MESH_NY, size_t MESH_NZ, std::vector<std::vector<std::vector<float>>>& aim_phi) {
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!is_porous_rand) {
				gen.seed(porous_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			std::cout << "> Quartet structure generation:";
			if (MESH_NZ == 1) {
				std::vector<std::vector<size_t>> Solid;
				Solid.resize(MESH_NX);
				/// No solid in the beginning
				for (size_t i = 0; i < MESH_NX; i++) {
					Solid[i].resize(MESH_NY);
					for (size_t j = 0; j < MESH_NY; j++) {
						Solid[i][j] = 0;
					}
				}
				/// Produce growth core
				float core_p = 0.0;
				for (size_t i = 0; i < MESH_NX; i++) {
					for (size_t j = 0; j < MESH_NY; j++) {
						if (real_dist(gen) < porous_init_noise) {
							Solid[i][j] = 1;
							core_p += 1.0;
						}
					}
				}
				core_p /= MESH_NX * MESH_NY;
				std::string out_ptr = "  percentage of core in 2D = " + std::to_string(core_p) + "\n";
				std::cout << out_ptr;
				///  Produce process
				size_t step = 0;
				do {
					step++;
					quartet_structure_grow_2D(Solid, gen);
					/// Calculate porosity
					core_p = 0.0;
					for (size_t i = 0; i < MESH_NX; i++) {
						for (size_t j = 0; j < MESH_NY; j++) {
							core_p += Solid[i][j];
						}
					}
					core_p = core_p / (MESH_NX * MESH_NY);
					out_ptr = "  grow step: " + std::to_string(step) + ", percentage of particle/pore = " + std::to_string(core_p) + "\n";
					std::cout << out_ptr;
				} while (core_p < porosity);
				out_ptr = "  End of growth ! \n";
				std::cout << out_ptr;
				geometry_structure::PointSet set_phi_1, set_phi_2;
				for (int i = 0; i < MESH_NX; i++) {
					for (int j = 0; j < MESH_NY; j++) {
						if (aim_phi[i][j][0] > 1e-6) {
							if (Solid[i][j]) {
								set_phi_1.add_point(i + 1, j + 1, 0 + 1, aim_phi[i][j][0]);
							}
							else {
								set_phi_2.add_point(i + 1, j + 1, 0 + 1, aim_phi[i][j][0]);
							}
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = porous_first_phi_index;
				set_phi_1.temperature = porous_first_temperature;
				set_phi_1.con = porous_first_con;
				set_phi_1.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = porous_second_phi_index;
				set_phi_2.temperature = porous_second_temperature;
				set_phi_2.con = porous_second_con;
				set_phi_2.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else if (MESH_NZ > 1) {
				std::vector<std::vector<std::vector<size_t>>> arrgrid; arrgrid.resize(MESH_NX);
				std::vector<std::vector<size_t>> soild;
				size_t numsoild = 0, Nxyz = MESH_NX * MESH_NY * MESH_NZ, numtotal_need = size_t(porosity * Nxyz);
				for (size_t i = 0; i < MESH_NX; i++) {
					arrgrid[i].resize(MESH_NY);
					for (size_t j = 0; j < MESH_NY; j++) {
						arrgrid[i][j].resize(MESH_NZ);
						for (size_t k = 0; k < MESH_NZ; k++) {
							arrgrid[i][j][k] = 0;
							if (real_dist(gen) < porous_init_noise) {
								arrgrid[i][j][k] = 1;
								std::vector<size_t> buff = { i, j, k };
								soild.push_back(buff);
								numsoild = numsoild + 1;
							}
						}
					}
				}
				float core_p = float(numsoild) / float(Nxyz);
				std::string out_ptr = "  percentage of core in 3D = " + std::to_string(core_p) + "\n";
				std::cout << out_ptr;
				size_t Tnumsoild = numsoild, step = 0;
				while (Tnumsoild < numtotal_need) {
					step++;
					quartet_structure_grow_3D(arrgrid, soild, Tnumsoild, gen);
					core_p = float(Tnumsoild) / float(Nxyz);
					out_ptr = "  grow step: " + std::to_string(step) + ", percentage of particle/pore = " + std::to_string(core_p) + "\n";
					std::cout << out_ptr;
				}
				out_ptr = "  End of growth ! \n";
				std::cout << out_ptr;
				geometry_structure::PointSet set_phi_1, set_phi_2;
				for (size_t i = 0; i < MESH_NX; i++) {
					for (size_t j = 0; j < MESH_NY; j++) {
						for (size_t k = 0; k < MESH_NZ; k++) {
							if (aim_phi[i][j][k] > 1e-6) {
								if (arrgrid[i][j][k]) {
									set_phi_1.add_point(i + 1, j + 1, k + 1, aim_phi[i][j][k]);
								}
								else {
									set_phi_2.add_point(i + 1, j + 1, k + 1, aim_phi[i][j][k]);
								}
							}
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = porous_first_phi_index;
				set_phi_1.temperature = porous_first_temperature;
				set_phi_1.con = porous_first_con;
				set_phi_1.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = porous_second_phi_index;
				set_phi_2.temperature = porous_second_temperature;
				set_phi_2.con = porous_second_con;
				set_phi_2.is_normalized = is_porous_normalized;
				geometry_structure::nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else {
				std::cout << "> Quartet structure generation Failed, can't generate at one dimension ! \n";
			}
		}

	}
}