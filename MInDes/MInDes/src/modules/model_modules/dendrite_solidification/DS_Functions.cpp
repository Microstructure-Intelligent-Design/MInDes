#include "DS_Functions.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <new>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace pf {
	namespace dendrite_solidification_model {
		namespace {
			constexpr REAL obstacle_factor = REAL(16.0 / PI2);

			inline size_t checked_mul(size_t a, size_t b) {
				if (a != 0 && b > (std::numeric_limits<size_t>::max)() / a) throw std::bad_alloc();
				return a * b;
			}
			template<typename T> inline size_t vector_bytes(const std::vector<T>& value) { return value.capacity() * sizeof(T); }
			inline void mix_signature(std::uintptr_t& result, const void* data, size_t capacity) {
				result ^= reinterpret_cast<std::uintptr_t>(data) + std::uintptr_t(0x9e3779b9U) + (result << 6) + (result >> 2);
				result ^= static_cast<std::uintptr_t>(capacity) + std::uintptr_t(0x85ebca6bU) + (result << 6) + (result >> 2);
			}
			inline size_t thread_index() {
#ifdef _OPENMP
				return static_cast<size_t>(omp_get_thread_num());
#else
				return 0;
#endif
			}
			inline bool callbacks_unchanged() {
				using namespace parameters;
				return phi_bulk_data == phi_bulk_terms.data() && phi_bulk_size == phi_bulk_terms.size()
					&& phi_pair_data == phi_pair_sources.data() && phi_pair_size == phi_pair_sources.size()
					&& con_mu_data == con_mu_terms.data() && con_mu_size == con_mu_terms.size()
					&& con_source_data == con_sources.data() && con_source_size == con_sources.size()
					&& temp_source_data == temp_sources.data() && temp_source_size == temp_sources.size();
			}
			inline void freeze_callbacks() {
				using namespace parameters;
				phi_bulk_data = phi_bulk_terms.data(); phi_bulk_size = phi_bulk_terms.size();
				phi_pair_data = phi_pair_sources.data(); phi_pair_size = phi_pair_sources.size();
				con_mu_data = con_mu_terms.data(); con_mu_size = con_mu_terms.size();
				con_source_data = con_sources.data(); con_source_size = con_sources.size();
				temp_source_data = temp_sources.data(); temp_source_size = temp_sources.size();
				callbacks_frozen = true;
			}

			inline REAL divergence_of_phase_flux(long long x, long long y, long long z, size_t i) {
				const auto& flux = parameters::Dendrite_workspace.interface_flux;
				const REAL inv2 = REAL(1) / (REAL(2) * mesh_parameters::delt_r);
				return ((flux(x + 1, y, z, i)[0] - flux(x - 1, y, z, i)[0])
					+ (flux(x, y + 1, z, i)[1] - flux(x, y - 1, z, i)[1])
					+ (flux(x, y, z + 1, i)[2] - flux(x, y, z - 1, i)[2])) * inv2;
			}

			inline int current_flag(long long x, long long y, long long z, size_t i) {
				const REAL phi = parameters::old_phi(x, y, z, i);
				if (phi >= parameters::Phi_Cut_Off && phi <= parameters::Phi_Cut_Off_R) return InterfaceFlag::IF_INTERFACE;
				const REAL neighbours[6] = {
					parameters::old_phi(x + 1, y, z, i), parameters::old_phi(x - 1, y, z, i),
					parameters::old_phi(x, y + 1, z, i), parameters::old_phi(x, y - 1, z, i),
					parameters::old_phi(x, y, z + 1, i), parameters::old_phi(x, y, z - 1, i)
				};
				for (REAL neighbour : neighbours)
					if ((phi < parameters::Phi_Cut_Off && neighbour >= parameters::Phi_Cut_Off)
						|| (phi > parameters::Phi_Cut_Off_R && neighbour <= parameters::Phi_Cut_Off_R))
						return InterfaceFlag::IF_NEAR_INTERFACE;
				return InterfaceFlag::IF_BULK;
			}

			inline void scalar_gradient_laplacian(long long x, long long y, long long z, size_t i,
				Mesh_Boundry<Matrix1D<REAL>>& field, Vector3& grad, REAL& lap) {
				const REAL dr = mesh_parameters::delt_r;
				grad[0] = (field(x + 1, y, z)[i] - field(x - 1, y, z)[i]) / (REAL(2) * dr);
				grad[1] = (field(x, y + 1, z)[i] - field(x, y - 1, z)[i]) / (REAL(2) * dr);
				grad[2] = (field(x, y, z + 1)[i] - field(x, y, z - 1)[i]) / (REAL(2) * dr);
				lap = (field(x + 1, y, z)[i] + field(x - 1, y, z)[i] + field(x, y + 1, z)[i]
					+ field(x, y - 1, z)[i] + field(x, y, z + 1)[i] + field(x, y, z - 1)[i]
					- REAL(6) * field(x, y, z)[i]) / (dr * dr);
			}
			inline void temperature_gradient_laplacian(long long x, long long y, long long z, Vector3& grad, REAL& lap) {
				auto& field = main_field::temperature_field;
				const REAL dr = mesh_parameters::delt_r;
				grad[0] = (field(x + 1, y, z) - field(x - 1, y, z)) / (REAL(2) * dr);
				grad[1] = (field(x, y + 1, z) - field(x, y - 1, z)) / (REAL(2) * dr);
				grad[2] = (field(x, y, z + 1) - field(x, y, z - 1)) / (REAL(2) * dr);
				lap = (field(x + 1, y, z) + field(x - 1, y, z) + field(x, y + 1, z) + field(x, y - 1, z)
					+ field(x, y, z + 1) + field(x, y, z - 1) - REAL(6) * field(x, y, z)) / (dr * dr);
			}

			inline void accumulate_interface_pair(long long x, long long y, long long z, size_t cell, size_t alpha, size_t beta) {
				auto& workspace = parameters::Dendrite_workspace;
				const size_t oa = workspace.phi_offset(cell, alpha), ob = workspace.phi_offset(cell, beta);
				const REAL phi_a = parameters::old_phi(x, y, z, alpha), phi_b = parameters::old_phi(x, y, z, beta);
				const Vector3& grad_a = workspace.grad_phi[oa]; const Vector3& grad_b = workspace.grad_phi[ob];
				Vector3 p = grad_a * phi_b - grad_b * phi_a;
				const REAL length = std::sqrt(p * p), sigma0 = parameters::sigma_ab(alpha, beta);
				REAL sigma = sigma0; Vector3 h(0, 0, 0);
				if (length >= SYS_EPSILON) {
					Vector3 normal = p / length;
					Vector3 crystal_normal = parameters::grain_rotation_matrix[beta] * normal;
					const REAL fourth_sum = std::pow(crystal_normal[0], 4.0) + std::pow(crystal_normal[1], 4.0) + std::pow(crystal_normal[2], 4.0);
					sigma *= REAL(1) + parameters::anisotropy_delta * (REAL(1.5) - REAL(2.5) * fourth_sum);
					Vector3 dsigma_dcrystal(-REAL(10) * sigma0 * parameters::anisotropy_delta * std::pow(crystal_normal[0], 3.0),
						-REAL(10) * sigma0 * parameters::anisotropy_delta * std::pow(crystal_normal[1], 3.0),
						-REAL(10) * sigma0 * parameters::anisotropy_delta * std::pow(crystal_normal[2], 3.0));
					Vector3 dsigma_dnormal = parameters::grain_rotation_matrix[beta].get_transposed() * dsigma_dcrystal;
					h = (dsigma_dnormal - normal * (normal * dsigma_dnormal)) / length;
				}
				const REAL eta = parameters::interface_width, s = grad_a * grad_b;
				const REAL h_dot_a = h * grad_a, h_dot_b = h * grad_b;
				workspace.mu_phi[oa] += eta * (sigma0 * workspace.lap_phi_or_rhs[ob] + s * h_dot_b)
					+ obstacle_factor / eta * (sigma * phi_b - phi_a * phi_b * h_dot_b);
				workspace.interface_flux(x, y, z, alpha) += (grad_b * (sigma - sigma0) + h * (s * phi_b)) * eta
					- h * (obstacle_factor / eta * phi_a * phi_b * phi_b);
				workspace.mu_phi[ob] += eta * (sigma0 * workspace.lap_phi_or_rhs[oa] - s * h_dot_a)
					+ obstacle_factor / eta * (sigma * phi_a + phi_a * phi_b * h_dot_a);
				workspace.interface_flux(x, y, z, beta) += (grad_a * (sigma - sigma0) - h * (s * phi_a)) * eta
					+ h * (obstacle_factor / eta * phi_a * phi_a * phi_b);
			}

			inline void normalize_rhs(long long x, long long y, long long z, size_t cell, REAL* candidate) {
				auto& workspace = parameters::Dendrite_workspace;
				REAL sum = 0;
				for (size_t i = 0; i < workspace.phi_number; ++i) {
					const REAL rhs = workspace.interface_flux(x, y, z, i)[0];
					candidate[i] = std::max(REAL(0), std::min(REAL(1), parameters::old_phi(x, y, z, i) + time_parameters::delt_t * rhs));
					sum += candidate[i];
				}
				if (sum > SYS_EPSILON) for (size_t i = 0; i < workspace.phi_number; ++i) candidate[i] /= sum;
				else for (size_t i = 0; i < workspace.phi_number; ++i) candidate[i] = parameters::old_phi(x, y, z, i);
				for (size_t i = 0; i < workspace.phi_number; ++i)
					workspace.interface_flux(x, y, z, i)[0] = (candidate[i] - parameters::old_phi(x, y, z, i)) / time_parameters::delt_t;
			}
		}

		void DendriteWorkspace::init(size_t nx_value, size_t ny_value, size_t nz_value, size_t phi_count, size_t con_count,
			size_t active_count_value, size_t threads) {
			clear();
			nx = nx_value; ny = ny_value; nz = nz_value; phi_number = phi_count; con_number = con_count;
			active_capacity = active_count_value; thread_number = std::max(size_t(1), threads);
			cell_number = checked_mul(checked_mul(nx, ny), nz);
			active_indices.assign(checked_mul(cell_number, active_capacity), phi_number);
			active_flags.assign(checked_mul(cell_number, active_capacity), uint8_t(InterfaceFlag::IF_BULK));
			active_count.assign(cell_number, 0);
			grad_phi.assign(checked_mul(cell_number, phi_number), Vector3(0, 0, 0));
			lap_phi_or_rhs.assign(checked_mul(cell_number, phi_number), 0);
			mu_phi.assign(checked_mul(cell_number, phi_number), 0);
			grad_con.assign(checked_mul(cell_number, con_number), Vector3(0, 0, 0));
			lap_con_or_rhs.assign(checked_mul(cell_number, con_number), 0);
			grad_temp.assign(cell_number, Vector3(0, 0, 0)); lap_temp.assign(cell_number, 0);
			delta_g.assign(cell_number, 0); temp_rhs.assign(cell_number, 0);
			interface_flux.init(nx, ny, nz, phi_number, mesh_parameters::x_down, mesh_parameters::x_up,
				mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			mu_con.init(nx, ny, nz, con_number, mesh_parameters::x_down, mesh_parameters::x_up,
				mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			mob_con.init(nx, ny, nz, con_number, mesh_parameters::x_down, mesh_parameters::x_up,
				mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			mob_temp.init(nx, ny, nz, 1, mesh_parameters::x_down, mesh_parameters::x_up,
				mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			thread_phi_scratch.assign(checked_mul(thread_number, phi_number), 0);
			initialized = true; allocation_signature = current_allocation_signature();
		}

		void DendriteWorkspace::clear() {
			std::vector<size_t>().swap(active_indices); std::vector<size_t>().swap(active_count); std::vector<uint8_t>().swap(active_flags);
			std::vector<Vector3>().swap(grad_phi); std::vector<Vector3>().swap(grad_con); std::vector<Vector3>().swap(grad_temp);
			std::vector<REAL>().swap(lap_phi_or_rhs); std::vector<REAL>().swap(mu_phi); std::vector<REAL>().swap(lap_con_or_rhs);
			std::vector<REAL>().swap(lap_temp); std::vector<REAL>().swap(delta_g); std::vector<REAL>().swap(temp_rhs);
			std::vector<REAL>().swap(thread_phi_scratch);
			interface_flux.clear(); mu_con.clear(); mob_con.clear(); mob_temp.clear();
			nx = ny = nz = phi_number = con_number = active_capacity = cell_number = 0; thread_number = 1;
			allocation_signature = 0; initialized = false;
		}

		size_t DendriteWorkspace::cell(long long x, long long y, long long z) const {
			return static_cast<size_t>(x - 1) + static_cast<size_t>(y - 1) * nx + static_cast<size_t>(z - 1) * nx * ny;
		}

		std::uintptr_t DendriteWorkspace::current_allocation_signature() const {
			std::uintptr_t result = 0;
			mix_signature(result, active_indices.data(), active_indices.capacity()); mix_signature(result, active_flags.data(), active_flags.capacity());
			mix_signature(result, active_count.data(), active_count.capacity()); mix_signature(result, grad_phi.data(), grad_phi.capacity());
			mix_signature(result, lap_phi_or_rhs.data(), lap_phi_or_rhs.capacity()); mix_signature(result, mu_phi.data(), mu_phi.capacity());
			mix_signature(result, grad_con.data(), grad_con.capacity()); mix_signature(result, lap_con_or_rhs.data(), lap_con_or_rhs.capacity());
			mix_signature(result, grad_temp.data(), grad_temp.capacity()); mix_signature(result, lap_temp.data(), lap_temp.capacity());
			mix_signature(result, delta_g.data(), delta_g.capacity()); mix_signature(result, temp_rhs.data(), temp_rhs.capacity());
			mix_signature(result, interface_flux.data(), interface_flux.capacity()); mix_signature(result, mu_con.data(), mu_con.capacity());
			mix_signature(result, mob_con.data(), mob_con.capacity()); mix_signature(result, mob_temp.data(), mob_temp.capacity());
			mix_signature(result, thread_phi_scratch.data(), thread_phi_scratch.capacity());
			return result;
		}

		size_t DendriteWorkspace::allocated_bytes() const {
			return vector_bytes(active_indices) + vector_bytes(active_flags) + vector_bytes(active_count)
				+ vector_bytes(grad_phi) + vector_bytes(lap_phi_or_rhs) + vector_bytes(mu_phi)
				+ vector_bytes(grad_con) + vector_bytes(lap_con_or_rhs) + vector_bytes(grad_temp)
				+ vector_bytes(lap_temp) + vector_bytes(delta_g) + vector_bytes(temp_rhs) + vector_bytes(thread_phi_scratch)
				+ interface_flux.bytes() + mu_con.bytes() + mob_con.bytes() + mob_temp.bytes();
		}

		namespace field_functions {
			REAL interface_mobility_const(size_t a, size_t b, REAL phi_a, REAL phi_b,
				const Vector3& grad_a, const Vector3& grad_b, REAL) {
				const REAL base_mobility = parameters::Lij(a, b);
				const Vector3 p = grad_a * phi_b - grad_b * phi_a;
				const REAL length = std::sqrt(p * p);
				if (length < SYS_EPSILON) return base_mobility;
				const Vector3 crystal_normal = parameters::grain_rotation_matrix[b] * (p / length);
				const REAL fourth_sum = std::pow(crystal_normal[0], 4.0)
					+ std::pow(crystal_normal[1], 4.0) + std::pow(crystal_normal[2], 4.0);
				return base_mobility * (REAL(1) - parameters::interface_mobility_anisotropy_delta
					* (REAL(1.5) - REAL(2.5) * fourth_sum));
			}
			REAL concentration_mobility_const(long long, long long, long long, size_t i) {
				return parameters::con_mobility_field.empty() ? parameters::con_mobility_const : parameters::con_mobility_field[i];
			}
			REAL linearized_driving_force(long long x, long long y, long long z) {
				REAL result = parameters::driving_force_dtemp * (parameters::old_temp(x, y, z) - parameters::equilibrium_temperature);
				if (parameters::has_components())
					for (size_t i = 0; i < parameters::component_number; ++i)
						result += parameters::driving_force_dcon[i] * (parameters::old_con(x, y, z, parameters::liquid_con_index(i)) - parameters::liquid_con_equilibrium[i]);
				return result;
			}
			REAL solidification_driving_pair_rate(long long x, long long y, long long z, size_t solid, REAL mobility) {
				if (solid == 0 || solid >= main_field::phi_number) return 0;
				return -mobility * REAL(PI) / parameters::interface_width
					* std::sqrt(std::max(REAL(0), parameters::old_phi(x, y, z, 0) * parameters::old_phi(x, y, z, solid)))
					* parameters::driving_force(x, y, z);
			}
			REAL temperature_mobility_mixture(long long x, long long y, long long z) {
				if (parameters::temp_mobility_phase.empty()) return parameters::temp_mobility_const;
				REAL mobility = 0;
				for (size_t i = 0; i < main_field::phi_number; ++i) mobility += parameters::old_phi(x, y, z, i) * parameters::temp_mobility_phase[i];
				return mobility;
			}

			void init_dendrite_field() {
				size_t threads = 1;
#ifdef _OPENMP
				threads = static_cast<size_t>(omp_get_max_threads());
#endif
				try {
					parameters::Dendrite_workspace.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ,
						main_field::phi_number, main_field::con_number, parameters::PHI_ACC_NUMBER, threads);
				}
				catch (const std::bad_alloc&) {
					WriteLog("> ERROR, Dendrite Solidification: workspace size overflow or allocation failure.\n"); SYS_PROGRAM_STOP;
				}
				if (parameters::is_thermal_only()) {
					const auto& workspace = parameters::Dendrite_workspace;
					if (workspace.con_number != 0 || workspace.grad_con.capacity() != 0
						|| workspace.lap_con_or_rhs.capacity() != 0 || workspace.mu_con.bytes() != 0 || workspace.mob_con.bytes() != 0) {
						WriteLog("> ERROR, Dendrite Solidification: thermal-only workspace allocated concentration storage.\n"); SYS_PROGRAM_STOP;
					}
					WriteLog("> DS thermal-only workspace: concentration cache capacity = 0 bytes.\n");
				}
				freeze_callbacks();
				const size_t bytes = parameters::Dendrite_workspace.allocated_bytes();
				WriteLog("> DS fixed workspace allocated: " + std::to_string(bytes/1000) + " Kbytes (" + std::to_string(bytes / REAL(1024 * 1024)) + " MiB).\n");
			}

			void prepare_step_snapshot() {
				auto& workspace = parameters::Dendrite_workspace;
				main_field::phase_field.do_boundary_condition();
				if (parameters::has_components()) main_field::concentration_field.do_boundary_condition();
				main_field::temperature_field.do_boundary_condition();
				std::fill(workspace.active_indices.begin(), workspace.active_indices.end(), parameters::PAIRWISE_ACC_STOP);
				std::fill(workspace.active_flags.begin(), workspace.active_flags.end(), uint8_t(InterfaceFlag::IF_BULK));
				std::fill(workspace.active_count.begin(), workspace.active_count.end(), size_t(0));
				std::fill(workspace.mu_phi.begin(), workspace.mu_phi.end(), REAL(0));
				workspace.interface_flux.fill(Vector3(0, 0, 0));
				if (parameters::has_components()) { workspace.mu_con.fill(0); workspace.mob_con.fill(0); }
				workspace.mob_temp.fill(0);
				bool overflow = false; long long overflow_x = 0, overflow_y = 0, overflow_z = 0; size_t overflow_count = 0;
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z); size_t active = 0;
							for (size_t i = 0; i < workspace.phi_number; ++i) {
								const size_t offset = workspace.phi_offset(cell, i);
								scalar_gradient_laplacian(x, y, z, i, main_field::phase_field, workspace.grad_phi[offset], workspace.lap_phi_or_rhs[offset]);
								const int flag = current_flag(x, y, z, i);
								if (flag || parameters::old_phi(x, y, z, i) > parameters::PhiCon_Cut_Off) {
									if (active < workspace.active_capacity) {
										const size_t ao = workspace.active_offset(cell, active);
										workspace.active_indices[ao] = i; workspace.active_flags[ao] = static_cast<uint8_t>(flag);
									}
									++active;
								}
							}
							workspace.active_count[cell] = std::min(active, workspace.active_capacity);
							if (active > workspace.active_capacity) {
#ifdef _OPENMP
#pragma omp critical(ds_active_overflow)
#endif
								if (!overflow) { overflow = true; overflow_x = x; overflow_y = y; overflow_z = z; overflow_count = active; }
							}
							for (size_t i = 0; i < workspace.con_number; ++i) {
								const size_t offset = workspace.con_offset(cell, i);
								scalar_gradient_laplacian(x, y, z, i, main_field::concentration_field, workspace.grad_con[offset], workspace.lap_con_or_rhs[offset]);
							}
							temperature_gradient_laplacian(x, y, z, workspace.grad_temp[cell], workspace.lap_temp[cell]);
							workspace.delta_g[cell] = linearized_driving_force(x, y, z); workspace.temp_rhs[cell] = 0;
						}
				if (overflow) {
					WriteLog("> ERROR, Dendrite Solidification: active phase capacity exceeded at (" + std::to_string(overflow_x) + ","
						+ std::to_string(overflow_y) + "," + std::to_string(overflow_z) + "): required " + std::to_string(overflow_count)
						+ ", configured " + std::to_string(workspace.active_capacity) + ".\n"); SYS_PROGRAM_STOP;
				}
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z), count = workspace.active_count[cell];
							for (size_t ai = 0; ai < count; ++ai) for (size_t bi = ai + 1; bi < count; ++bi) {
								const size_t ao = workspace.active_offset(cell, ai), bo = workspace.active_offset(cell, bi);
								if (workspace.active_flags[ao] && workspace.active_flags[bo])
									accumulate_interface_pair(x, y, z, cell, workspace.active_indices[ao], workspace.active_indices[bo]);
							}
							if (parameters::triple_junction_energy != 0) {
								const REAL coeff = parameters::triple_junction_energy / parameters::interface_width;
								for (size_t ai = 0; ai < count; ++ai) for (size_t bi = ai + 1; bi < count; ++bi) for (size_t ci = bi + 1; ci < count; ++ci) {
									const size_t ao = workspace.active_offset(cell, ai), bo = workspace.active_offset(cell, bi), co = workspace.active_offset(cell, ci);
									if (workspace.active_flags[ao] && workspace.active_flags[bo] && workspace.active_flags[co]) {
										const size_t a = workspace.active_indices[ao], b = workspace.active_indices[bo], c = workspace.active_indices[co];
										workspace.mu_phi[workspace.phi_offset(cell, a)] += coeff * parameters::old_phi(x, y, z, b) * parameters::old_phi(x, y, z, c);
										workspace.mu_phi[workspace.phi_offset(cell, b)] += coeff * parameters::old_phi(x, y, z, a) * parameters::old_phi(x, y, z, c);
										workspace.mu_phi[workspace.phi_offset(cell, c)] += coeff * parameters::old_phi(x, y, z, a) * parameters::old_phi(x, y, z, b);
									}
								}
							}
						}
				workspace.interface_flux.do_boundary_condition();
			}

			void prepare_thermodynamic_fields() {
				auto& workspace = parameters::Dendrite_workspace;
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z);
							for (size_t i = 0; i < workspace.phi_number; ++i) {
								REAL& mu = workspace.mu_phi[workspace.phi_offset(cell, i)]; mu += divergence_of_phase_flux(x, y, z, i);
								for (auto term : parameters::phi_bulk_terms) mu += term(x, y, z, i);
							}
							for (size_t i = 0; i < workspace.con_number; ++i) {
								REAL mu = 0; for (auto term : parameters::con_mu_terms) mu += term(x, y, z, i);
								workspace.mu_con(x, y, z, i) = mu;
								workspace.mob_con(x, y, z, i) = parameters::con_mobility ? parameters::con_mobility(x, y, z, i) : 0;
							}
							workspace.mob_temp(x, y, z) = parameters::temp_mobility ? parameters::temp_mobility(x, y, z) : 0;
						}
				if (parameters::has_components()) {
					workspace.mu_con.do_boundary_condition(); workspace.mob_con.do_boundary_condition();
				}
				workspace.mob_temp.do_boundary_condition();
			}

			void calculate_right_hand_sides() {
				auto& workspace = parameters::Dendrite_workspace;
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z), count = workspace.active_count[cell];
							for (size_t i = 0; i < workspace.phi_number; ++i) workspace.interface_flux(x, y, z, i)[0] = 0;
							for (size_t ai = 0; ai < count; ++ai) for (size_t bi = ai + 1; bi < count; ++bi) {
								const size_t ao = workspace.active_offset(cell, ai), bo = workspace.active_offset(cell, bi);
								if (!workspace.active_flags[ao] || !workspace.active_flags[bo]) continue;
								const size_t a = workspace.active_indices[ao], b = workspace.active_indices[bo];
								const REAL mobility = parameters::interface_mobility ? parameters::interface_mobility(a, b,
									parameters::old_phi(x, y, z, a), parameters::old_phi(x, y, z, b),
									workspace.grad_phi[workspace.phi_offset(cell, a)], workspace.grad_phi[workspace.phi_offset(cell, b)], parameters::old_temp(x, y, z)) : 0;
								REAL rate = mobility / parameters::interface_width * (workspace.mu_phi[workspace.phi_offset(cell, b)] - workspace.mu_phi[workspace.phi_offset(cell, a)]);
								if (a == 0) rate += solidification_driving_pair_rate(x, y, z, b, mobility);
								for (auto source : parameters::phi_pair_sources) rate += source(x, y, z, a, b);
								const REAL phi_a = parameters::old_phi(x, y, z, a), phi_b = parameters::old_phi(x, y, z, b);
								if ((rate > SYS_EPSILON && (phi_a >= parameters::Phi_Cut_Off_R || phi_b <= parameters::Phi_Cut_Off))
									|| (rate < -SYS_EPSILON && (phi_a <= parameters::Phi_Cut_Off || phi_b >= parameters::Phi_Cut_Off_R))) rate = 0;
								workspace.interface_flux(x, y, z, a)[0] += rate; workspace.interface_flux(x, y, z, b)[0] -= rate;
							}
							if (parameters::is_phi_normalized) normalize_rhs(x, y, z, cell, workspace.phi_scratch(thread_index()));
						}
				// The phase Laplace buffer is now dead and becomes the accepted phase RHS.
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z);
							for (size_t i = 0; i < workspace.phi_number; ++i)
								workspace.lap_phi_or_rhs[workspace.phi_offset(cell, i)] = workspace.interface_flux(x, y, z, i)[0];
						}
				const REAL dr = mesh_parameters::delt_r, inv2 = REAL(1) / (REAL(2) * dr), invsq = REAL(1) / (dr * dr);
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z);
							for (size_t i = 0; i < workspace.con_number; ++i) {
								Vector3 grad_mu, grad_mob;
								grad_mu[0] = (workspace.mu_con(x + 1, y, z, i) - workspace.mu_con(x - 1, y, z, i)) * inv2;
								grad_mu[1] = (workspace.mu_con(x, y + 1, z, i) - workspace.mu_con(x, y - 1, z, i)) * inv2;
								grad_mu[2] = (workspace.mu_con(x, y, z + 1, i) - workspace.mu_con(x, y, z - 1, i)) * inv2;
								grad_mob[0] = (workspace.mob_con(x + 1, y, z, i) - workspace.mob_con(x - 1, y, z, i)) * inv2;
								grad_mob[1] = (workspace.mob_con(x, y + 1, z, i) - workspace.mob_con(x, y - 1, z, i)) * inv2;
								grad_mob[2] = (workspace.mob_con(x, y, z + 1, i) - workspace.mob_con(x, y, z - 1, i)) * inv2;
								const REAL lap_mu = (workspace.mu_con(x + 1, y, z, i) + workspace.mu_con(x - 1, y, z, i)
									+ workspace.mu_con(x, y + 1, z, i) + workspace.mu_con(x, y - 1, z, i)
									+ workspace.mu_con(x, y, z + 1, i) + workspace.mu_con(x, y, z - 1, i) - REAL(6) * workspace.mu_con(x, y, z, i)) * invsq;
								REAL rhs = grad_mob * grad_mu + workspace.mob_con(x, y, z, i) * lap_mu;
								for (auto source : parameters::con_sources) rhs += source(x, y, z, i);
								workspace.lap_con_or_rhs[workspace.con_offset(cell, i)] = rhs;
							}
							Vector3 grad_mt;
							grad_mt[0] = (workspace.mob_temp(x + 1, y, z) - workspace.mob_temp(x - 1, y, z)) * inv2;
							grad_mt[1] = (workspace.mob_temp(x, y + 1, z) - workspace.mob_temp(x, y - 1, z)) * inv2;
							grad_mt[2] = (workspace.mob_temp(x, y, z + 1) - workspace.mob_temp(x, y, z - 1)) * inv2;
							REAL rhs = grad_mt * workspace.grad_temp[cell] + workspace.mob_temp(x, y, z) * workspace.lap_temp[cell];
							for (size_t i = 0; i < workspace.phi_number; ++i) rhs += parameters::latent_heat_phase[i] * workspace.lap_phi_or_rhs[workspace.phi_offset(cell, i)];
							for (auto source : parameters::temp_sources) rhs += source(x, y, z);
							workspace.temp_rhs[cell] = rhs;
						}
			}

			void solve_fields() {
				auto& workspace = parameters::Dendrite_workspace; REAL max_phi = 0, max_con = 0, max_temp = 0;
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
							const size_t cell = workspace.cell(x, y, z); REAL local_phi = 0, local_con = 0;
							for (size_t i = 0; i < workspace.phi_number; ++i) {
								const REAL old_value = parameters::old_phi(x, y, z, i);
								const REAL new_value = old_value + time_parameters::delt_t * workspace.lap_phi_or_rhs[workspace.phi_offset(cell, i)];
								local_phi = std::max(local_phi, std::abs(new_value - old_value)); main_field::phase_field(x, y, z)[i] = new_value;
							}
							for (size_t i = 0; i < workspace.con_number; ++i) {
								const REAL old_value = parameters::old_con(x, y, z, i);
								const REAL new_value = old_value + time_parameters::delt_t * workspace.lap_con_or_rhs[workspace.con_offset(cell, i)];
								local_con = std::max(local_con, std::abs(new_value - old_value)); main_field::concentration_field(x, y, z)[i] = new_value;
							}
							const REAL old_t = parameters::old_temp(x, y, z), new_t = old_t + time_parameters::delt_t * workspace.temp_rhs[cell];
							main_field::temperature_field(x, y, z) = new_t;
#ifdef _OPENMP
#pragma omp critical(ds_max_variation)
#endif
							{ max_phi = std::max(max_phi, local_phi); max_con = std::max(max_con, local_con); max_temp = std::max(max_temp, std::abs(new_t - old_t)); }
#ifdef _DEBUG
							bool finite = std::isfinite(new_t);
							for (size_t i = 0; i < workspace.con_number; ++i) finite = finite && std::isfinite(main_field::concentration_field(x, y, z)[i]);
							if (!finite) { std::cout << "DEBUG: DS non-finite value at (" << x << "," << y << "," << z << ")" << std::endl; SYS_PROGRAM_STOP; }
#endif
						}
				main_field::PHI_MAX_VARIATION = max_phi; main_field::CON_MAX_VARIATION = max_con; main_field::TEMP_MAX_VARIATION = max_temp;
				main_field::phase_field.do_boundary_condition();
				if (parameters::has_components()) main_field::concentration_field.do_boundary_condition();
				main_field::temperature_field.do_boundary_condition();
			}
		}

		void exec_pre_iii() { field_functions::init_dendrite_field(); }
		void exec_i() {
			if (!parameters::Dendrite_workspace.initialized
				|| parameters::Dendrite_workspace.current_allocation_signature() != parameters::Dendrite_workspace.allocation_signature) {
				WriteLog("> ERROR, Dendrite Solidification: workspace allocation changed after initialization.\n"); SYS_PROGRAM_STOP;
			}
			if (!parameters::callbacks_frozen || !callbacks_unchanged()) {
				WriteLog("> ERROR, Dendrite Solidification: callback tables changed after workspace initialization.\n"); SYS_PROGRAM_STOP;
			}
			field_functions::prepare_step_snapshot(); field_functions::prepare_thermodynamic_fields();
			field_functions::calculate_right_hand_sides(); field_functions::solve_fields();
		}
		void deinit() {
			parameters::Dendrite_workspace.clear(); parameters::callbacks_frozen = false;
			main_field::phase_field.clear();
			if (parameters::has_components()) main_field::concentration_field.clear();
			main_field::temperature_field.clear();
		}
	}
}
