#pragma once
#include "Model_Functions.h"
namespace pf {
	namespace grain_grows_spinodal_model {
		void exec_pre_iii() {
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!parameters::is_noise_rand) {
				gen.seed(parameters::noise_seed);
			}
			std::uniform_real_distribution<> real_dist(-1.0, 1.0); // [-1.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			for (long long x = 0; x < parameters::concentration_field->Nx(); x++)
				for (long long y = 0; y < parameters::concentration_field->Ny(); y++)
					for (long long z = 0; z < parameters::concentration_field->Nz(); z++) {
						parameters::concentration_field->at(x, y, z)[0] += parameters::noise_amplitude * REAL(real_dist(gen));
					}
			parameters::con_field_variavles.init(parameters::MESH_NX, parameters::MESH_NY, parameters::MESH_NZ,
				parameters::delt_r, parameters::x_down, parameters::x_up, parameters::y_down, parameters::y_up,
				parameters::z_down, parameters::z_up);
			parameters::phi_increment.init(parameters::MESH_NX, parameters::MESH_NY, parameters::MESH_NZ,
				parameters::delt_r, parameters::x_down, parameters::x_up, parameters::y_down, parameters::y_up,
				parameters::z_down, parameters::z_up);
			for (long long x = 0; x < parameters::phi_increment.Nx(); x++)
				for (long long y = 0; y < parameters::phi_increment.Ny(); y++)
					for (long long z = 0; z < parameters::phi_increment.Nz(); z++) {
						parameters::con_field_variavles(x, y, z).resize(FIELD::NUM, 0);
						parameters::phi_increment(x, y, z).resize(parameters::phi_number, 0);
					}
			parameters::phase_field->do_boundary_condition();
			parameters::concentration_field->do_boundary_condition();
			// - init boundary
#pragma omp parallel for
			for (long long x = parameters::phase_field->COMP_X_BGN(); x <= parameters::phase_field->COMP_X_END(); x++)
				for (long long y = parameters::phase_field->COMP_Y_BGN(); y <= parameters::phase_field->COMP_Y_END(); y++)
					for (long long z = parameters::phase_field->COMP_Z_BGN(); z <= parameters::phase_field->COMP_Z_END(); z++) {
						std::vector<REAL>& phi = parameters::phase_field->at(x, y, z);
						REAL& con = parameters::concentration_field->at(x, y, z)[0];
						std::vector<REAL>& con_var = parameters::con_field_variavles(x, y, z);
						std::vector<REAL>& phi_incre = parameters::phi_increment(x, y, z);
						REAL etai2j2 = 0, etai42 = 0,
							m = REAL(1 + 0.5 * con * con - 2.5 * con * con * (1 - con) * (1 - con));
						for (size_t i = 0; i < parameters::phi_number; i++) {
							etai42 += phi[i] * phi[i] * phi[i] * phi[i] / 4 - phi[i] * phi[i] / 2;
							for (size_t j = i + 1; j < parameters::phi_number; j++)
								etai2j2 += phi[i] * phi[i] * phi[j] * phi[j];
						}
						for (size_t index = 0; index < parameters::phi_number; index++) {
							phi_incre[index] = (parameters::phase_field->at(x + 1, y, z)[index] + parameters::phase_field->at(x - 1, y, z)[index]
								+ parameters::phase_field->at(x, y + 1, z)[index] + parameters::phase_field->at(x, y - 1, z)[index]
								+ parameters::phase_field->at(x, y, z + 1)[index] + parameters::phase_field->at(x, y, z - 1)[index]
								- 6 * phi[index]) / parameters::delt_r / parameters::delt_r;
							REAL etai1j2 = 0;
							for (size_t j = 0; j < parameters::phi_number; j++)
								if (j != index)
									etai1j2 += 2 * phi[index] * phi[j] * phi[j];
							phi_incre[index] = -parameters::L * (m * (phi[index] * phi[index] * phi[index] - phi[index] + parameters::epsilon * etai1j2) - 2 * parameters::kappa_eta * phi_incre[index]);
						}
						con_var[FIELD::MOB] = parameters::Mb + std::sqrt(etai2j2) * parameters::Mg;
						con_var[FIELD::DFDCON] = (parameters::concentration_field->at(x + 1, y, z)[0] + parameters::concentration_field->at(x - 1, y, z)[0]
							+ parameters::concentration_field->at(x, y + 1, z)[0] + parameters::concentration_field->at(x, y - 1, z)[0]
							+ parameters::concentration_field->at(x, y, z + 1)[0] + parameters::concentration_field->at(x, y, z - 1)[0]
							- 6 * con) / parameters::delt_r / parameters::delt_r;
						con_var[FIELD::DFDCON] = (con - 5 * con * (1 - con) * (1 - 2 * con)) * (REAL(0.25) + etai42 + parameters::epsilon * etai2j2)
							+ 2 * parameters::A * con * (1 - con) * (1 - 2 * con) - 2 * parameters::kappa_con * con_var[FIELD::DFDCON];
					}
			parameters::phase_field->init_boundary_condition();
			parameters::concentration_field->init_boundary_condition();
			parameters::con_field_variavles.init_boundary_condition();
		}
		void exec_i() {
			// - boundary condition
			parameters::phase_field->do_boundary_condition();
			parameters::concentration_field->do_boundary_condition();
			// - phi increment + dFdcon + Mob
#pragma omp parallel for
			for (long long x = parameters::phase_field->COMP_X_BGN(); x <= parameters::phase_field->COMP_X_END(); x++)
				for (long long y = parameters::phase_field->COMP_Y_BGN(); y <= parameters::phase_field->COMP_Y_END(); y++)
					for (long long z = parameters::phase_field->COMP_Z_BGN(); z <= parameters::phase_field->COMP_Z_END(); z++) {
						std::vector<REAL>& phi = parameters::phase_field->at(x, y, z);
						REAL& con = parameters::concentration_field->at(x, y, z)[0];
						std::vector<REAL>& con_var = parameters::con_field_variavles(x, y, z);
						std::vector<REAL>& phi_incre = parameters::phi_increment(x, y, z);
						REAL etai2j2 = 0, etai42 = 0, 
							m = REAL(1 + 0.5 * con * con - 2.5 * con * con * (1 - con) * (1 - con));
						for (size_t i = 0; i < parameters::phi_number; i++) {
							etai42 += phi[i] * phi[i] * phi[i] * phi[i] / 4 - phi[i] * phi[i] / 2;
							for (size_t j = i + 1; j < parameters::phi_number; j++)
								etai2j2 += phi[i] * phi[i] * phi[j] * phi[j];
						}
						for (size_t index = 0; index < parameters::phi_number; index++) {
							phi_incre[index] = (parameters::phase_field->at(x + 1, y, z)[index] + parameters::phase_field->at(x - 1, y, z)[index]
								+ parameters::phase_field->at(x, y + 1, z)[index] + parameters::phase_field->at(x, y - 1, z)[index]
								+ parameters::phase_field->at(x, y, z + 1)[index] + parameters::phase_field->at(x, y, z - 1)[index]
								- 6 * phi[index]) / parameters::delt_r / parameters::delt_r;
							REAL etai1j2 = 0;
							for (size_t j = 0; j < parameters::phi_number; j++)
								if (j != index)
									etai1j2 += 2 * phi[index] * phi[j] * phi[j];
							phi_incre[index] = -parameters::L * (m * (phi[index] * phi[index] * phi[index] - phi[index] + parameters::epsilon * etai1j2) - 2 * parameters::kappa_eta * phi_incre[index]);
						}
						con_var[FIELD::MOB] = parameters::Mb + std::sqrt(etai2j2) * parameters::Mg;
						con_var[FIELD::DFDCON] = (parameters::concentration_field->at(x + 1, y, z)[0] + parameters::concentration_field->at(x - 1, y, z)[0]
							+ parameters::concentration_field->at(x, y + 1, z)[0] + parameters::concentration_field->at(x, y - 1, z)[0]
							+ parameters::concentration_field->at(x, y, z + 1)[0] + parameters::concentration_field->at(x, y, z - 1)[0]
							- 6 * con) / parameters::delt_r / parameters::delt_r;
						con_var[FIELD::DFDCON] = (con - 5 * con * (1 - con) * (1 - 2 * con)) * (REAL(0.25) + etai42 + parameters::epsilon * etai2j2)
							+ 2 * parameters::A * con * (1 - con) * (1 - 2 * con) - 2 * parameters::kappa_con * con_var[FIELD::DFDCON];
					}
			parameters::con_field_variavles.do_boundary_condition();
			// - phi + con
#pragma omp parallel for
			for (long long x = parameters::phase_field->COMP_X_BGN(); x <= parameters::phase_field->COMP_X_END(); x++)
				for (long long y = parameters::phase_field->COMP_Y_BGN(); y <= parameters::phase_field->COMP_Y_END(); y++)
					for (long long z = parameters::phase_field->COMP_Z_BGN(); z <= parameters::phase_field->COMP_Z_END(); z++) {
						std::vector<REAL>& phi = parameters::phase_field->at(x, y, z);
						REAL& con = parameters::concentration_field->at(x, y, z)[0];
						std::vector<REAL>& con_var = parameters::con_field_variavles(x, y, z);
						std::vector<REAL>& phi_incre = parameters::phi_increment(x, y, z);
						for (size_t index = 0; index < parameters::phi_number; index++) {
							REAL old_phi = phi[index];
							phi[index] += *parameters::delt_t * phi_incre[index];
							if (std::abs(phi[index] - old_phi) > *parameters::PHI_MAX_VARIATION)
								*parameters::PHI_MAX_VARIATION = std::abs(phi[index] - old_phi);
						}
						REAL old_con = con;
						std::vector<REAL> grad_M = { 
							 (parameters::con_field_variavles.at(x + 1, y, z)[FIELD::MOB]
							- parameters::con_field_variavles.at(x - 1, y, z)[FIELD::MOB]) / 2 / parameters::delt_r,
							 (parameters::con_field_variavles.at(x, y + 1, z)[FIELD::MOB]
							- parameters::con_field_variavles.at(x, y - 1, z)[FIELD::MOB]) / 2 / parameters::delt_r,
							 (parameters::con_field_variavles.at(x, y, z + 1)[FIELD::MOB]
							- parameters::con_field_variavles.at(x, y, z - 1)[FIELD::MOB]) / 2 / parameters::delt_r
						}, grad_DFDCON = {
							 (parameters::con_field_variavles.at(x + 1, y, z)[FIELD::DFDCON]
							- parameters::con_field_variavles.at(x - 1, y, z)[FIELD::DFDCON]) / 2 / parameters::delt_r,
							 (parameters::con_field_variavles.at(x, y + 1, z)[FIELD::DFDCON]
							- parameters::con_field_variavles.at(x, y - 1, z)[FIELD::DFDCON]) / 2 / parameters::delt_r,
							 (parameters::con_field_variavles.at(x, y, z + 1)[FIELD::DFDCON]
							- parameters::con_field_variavles.at(x, y, z - 1)[FIELD::DFDCON]) / 2 / parameters::delt_r 
						};
						REAL lap_DFDCON = (parameters::con_field_variavles.at(x + 1, y, z)[FIELD::DFDCON]
							+ parameters::con_field_variavles.at(x - 1, y, z)[FIELD::DFDCON] + parameters::con_field_variavles.at(x, y + 1, z)[FIELD::DFDCON]
							+ parameters::con_field_variavles.at(x, y - 1, z)[FIELD::DFDCON] + parameters::con_field_variavles.at(x, y, z + 1)[FIELD::DFDCON]
							+ parameters::con_field_variavles.at(x, y, z - 1)[FIELD::DFDCON] - 6 * con_var[FIELD::DFDCON]) / parameters::delt_r / parameters::delt_r;
						con += *parameters::delt_t * (grad_M[0] * grad_DFDCON[0] + grad_M[1] * grad_DFDCON[1] 
							+ grad_M[2] * grad_DFDCON[2] + con_var[FIELD::MOB] * lap_DFDCON);
						if (std::abs(con - old_con) > *parameters::CON_MAX_VARIATION)
							*parameters::CON_MAX_VARIATION = std::abs(con - old_con);
					}
		}
		void deinit() {
			parameters::con_field_variavles.clear();
			parameters::phi_increment.clear();
			parameters::phase_field->clear();
			parameters::phase_field = nullptr;
			parameters::concentration_field->clear();
			parameters::concentration_field = nullptr;
		}
	}
}