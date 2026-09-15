#pragma once
#include "../../base/Mesh_0.h"
#include "../../base/RotationMatrix.h"
#include "../../Modules_Params.h"
#include "../GrainsOrientations.h"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <new>
#include <string>
#include <vector>

namespace pf {
	namespace dendrite_solidification_model {
		enum InterfaceFlag { IF_BULK, IF_NEAR_INTERFACE, IF_INTERFACE };

		// One allocation for all cells/components: no per-cell dynamic storage.
		template<typename T>
		class DS_HaloField {
		public:
			void init(size_t nx, size_t ny, size_t nz, size_t components,
				BoundaryCondition x_down, BoundaryCondition x_up,
				BoundaryCondition y_down, BoundaryCondition y_up,
				BoundaryCondition z_down, BoundaryCondition z_up) {
				nx_ = checked_add(nx, 2); ny_ = checked_add(ny, 2); nz_ = checked_add(nz, 2); components_ = components;
				bc_x_down_ = x_down; bc_x_up_ = x_up; bc_y_down_ = y_down;
				bc_y_up_ = y_up; bc_z_down_ = z_down; bc_z_up_ = z_up;
				data_.assign(checked_mul(checked_mul(checked_mul(nx_, ny_), nz_), components_), T());
			}
			void clear() { std::vector<T>().swap(data_); nx_ = ny_ = nz_ = components_ = 0; }
			T& operator()(long long x, long long y, long long z, size_t c = 0) { return data_[index(x, y, z, c)]; }
			const T& operator()(long long x, long long y, long long z, size_t c = 0) const { return data_[index(x, y, z, c)]; }
			void fill(const T& value) { std::fill(data_.begin(), data_.end(), value); }
			const T* data() const { return data_.data(); }
			size_t capacity() const { return data_.capacity(); }
			size_t bytes() const { return data_.capacity() * sizeof(T); }

			void do_boundary_condition() {
				for (size_t c = 0; c < components_; ++c) {
					for (long long y = 0; y < static_cast<long long>(ny_); ++y)
						for (long long z = 0; z < static_cast<long long>(nz_); ++z) {
							apply((*this)(0, y, z, c), (*this)(1, y, z, c), (*this)(2, y, z, c), (*this)(static_cast<long long>(nx_) - 2, y, z, c), bc_x_down_);
							apply((*this)(static_cast<long long>(nx_) - 1, y, z, c), (*this)(static_cast<long long>(nx_) - 2, y, z, c), (*this)(static_cast<long long>(nx_) - 3, y, z, c), (*this)(1, y, z, c), bc_x_up_);
						}
					for (long long x = 0; x < static_cast<long long>(nx_); ++x)
						for (long long z = 0; z < static_cast<long long>(nz_); ++z) {
							apply((*this)(x, 0, z, c), (*this)(x, 1, z, c), (*this)(x, 2, z, c), (*this)(x, static_cast<long long>(ny_) - 2, z, c), bc_y_down_);
							apply((*this)(x, static_cast<long long>(ny_) - 1, z, c), (*this)(x, static_cast<long long>(ny_) - 2, z, c), (*this)(x, static_cast<long long>(ny_) - 3, z, c), (*this)(x, 1, z, c), bc_y_up_);
						}
					for (long long x = 0; x < static_cast<long long>(nx_); ++x)
						for (long long y = 0; y < static_cast<long long>(ny_); ++y) {
							apply((*this)(x, y, 0, c), (*this)(x, y, 1, c), (*this)(x, y, 2, c), (*this)(x, y, static_cast<long long>(nz_) - 2, c), bc_z_down_);
							apply((*this)(x, y, static_cast<long long>(nz_) - 1, c), (*this)(x, y, static_cast<long long>(nz_) - 2, c), (*this)(x, y, static_cast<long long>(nz_) - 3, c), (*this)(x, y, 1, c), bc_z_up_);
						}
				}
			}

		private:
			static size_t checked_add(size_t a, size_t b) { if (a > (std::numeric_limits<size_t>::max)() - b) throw std::bad_alloc(); return a + b; }
			static size_t checked_mul(size_t a, size_t b) { if (a && b > (std::numeric_limits<size_t>::max)() / a) throw std::bad_alloc(); return a * b; }
			static void apply(T& ghost, const T& adjacent, const T& next, const T& periodic, BoundaryCondition bc) {
				if (bc == BoundaryCondition::PERIODIC) ghost = periodic;
				else if (bc == BoundaryCondition::ZEROFLUX) ghost = adjacent;
				else if (bc == BoundaryCondition::OPENFLUX) ghost = adjacent * REAL(2) - next;
				else if (bc == BoundaryCondition::FIXED) ghost = adjacent;
			}
			size_t index(long long x, long long y, long long z, size_t c) const {
				return (((c * nz_ + static_cast<size_t>(z)) * ny_ + static_cast<size_t>(y)) * nx_ + static_cast<size_t>(x));
			}
			size_t nx_ = 0, ny_ = 0, nz_ = 0, components_ = 0;
			BoundaryCondition bc_x_down_ = BoundaryCondition::PERIODIC, bc_x_up_ = BoundaryCondition::PERIODIC;
			BoundaryCondition bc_y_down_ = BoundaryCondition::PERIODIC, bc_y_up_ = BoundaryCondition::PERIODIC;
			BoundaryCondition bc_z_down_ = BoundaryCondition::PERIODIC, bc_z_up_ = BoundaryCondition::PERIODIC;
			std::vector<T> data_;
		};

		struct DendriteWorkspace {
			size_t nx = 0, ny = 0, nz = 0, phi_number = 0, con_number = 0, active_capacity = 0;
			size_t cell_number = 0, thread_number = 1;
			std::vector<size_t> active_indices, active_count;
			std::vector<uint8_t> active_flags;
			std::vector<Vector3> grad_phi, grad_con, grad_temp;
			std::vector<REAL> lap_phi_or_rhs, mu_phi, lap_con_or_rhs, lap_temp, delta_g, temp_rhs;
			DS_HaloField<Vector3> interface_flux;
			DS_HaloField<REAL> mu_con, mob_con, mob_temp;
			std::vector<REAL> thread_phi_scratch;
			std::uintptr_t allocation_signature = 0;
			bool initialized = false;

			void init(size_t nx_value, size_t ny_value, size_t nz_value, size_t phi_count, size_t con_count, size_t active_count_value, size_t threads);
			void clear();
			size_t cell(long long x, long long y, long long z) const;
			size_t phi_offset(size_t c, size_t i) const { return c * phi_number + i; }
			size_t con_offset(size_t c, size_t i) const { return c * con_number + i; }
			size_t active_offset(size_t c, size_t slot) const { return c * active_capacity + slot; }
			REAL* phi_scratch(size_t thread) { return thread_phi_scratch.data() + thread * phi_number; }
			std::uintptr_t current_allocation_signature() const;
			size_t allocated_bytes() const;
		};

		namespace parameters {
			using PhiBulkTerm = REAL(*)(long long, long long, long long, size_t);
			using PhiPairSource = REAL(*)(long long, long long, long long, size_t, size_t);
			using ConTerm = REAL(*)(long long, long long, long long, size_t);
			using ScalarTerm = REAL(*)(long long, long long, long long);
			using InterfaceMobility = REAL(*)(size_t, size_t, REAL, REAL, const Vector3&, const Vector3&, REAL);

			inline DendriteWorkspace Dendrite_workspace;
			inline size_t PAIRWISE_ACC_STOP = 0, PHI_ACC_NUMBER = 0;
			inline REAL Phi_Cut_Off = REAL(0.001), Phi_Cut_Off_R = REAL(0.999), PhiCon_Cut_Off = REAL(0.1);
			inline bool is_phi_normalized = true;
			inline REAL interface_width = REAL(4.0), anisotropy_delta = 0, interface_mobility_anisotropy_delta = 0, triple_junction_energy = 0;
			inline Matrix2D<REAL> sigma_ab, Lij;
			inline std::vector<Matrix3x3> grain_rotation_matrix;
			inline REAL con_mobility_const = 0, temp_mobility_const = 0;
			inline size_t component_number = 0;
			inline std::vector<std::string> component_names;
			inline std::vector<REAL> liquid_con_equilibrium, driving_force_dcon, con_mobility_field;
			inline REAL equilibrium_temperature = 0, driving_force_dtemp = 0;
			inline std::vector<REAL> temp_mobility_phase, latent_heat_phase;
			inline std::vector<PhiBulkTerm> phi_bulk_terms;
			inline std::vector<PhiPairSource> phi_pair_sources;
			inline std::vector<ConTerm> con_mu_terms, con_sources;
			inline std::vector<ScalarTerm> temp_sources;
			inline InterfaceMobility interface_mobility = nullptr;
			inline ConTerm con_mobility = nullptr;
			inline ScalarTerm temp_mobility = nullptr;
			inline bool callbacks_frozen = false;
			inline const void* phi_bulk_data = nullptr, * phi_pair_data = nullptr, * con_mu_data = nullptr;
			inline const void* con_source_data = nullptr, * temp_source_data = nullptr;
			inline size_t phi_bulk_size = 0, phi_pair_size = 0, con_mu_size = 0, con_source_size = 0, temp_source_size = 0;

			inline bool has_components() { return component_number != 0; }
			inline bool is_thermal_only() { return !has_components(); }
			inline size_t liquid_con_index(size_t i) { return i; }
			inline size_t solid_con_index(size_t i) { return component_number + i; }
			inline size_t component_index(size_t i) {
				return has_components() ? i % component_number : (std::numeric_limits<size_t>::max)();
			}
			inline bool is_liquid_con(size_t i) { return has_components() && i < component_number; }
			inline REAL old_phi(long long x, long long y, long long z, size_t i) { return main_field::phase_field(x, y, z)[i]; }
			inline REAL old_con(long long x, long long y, long long z, size_t i) {
				return has_components() ? main_field::concentration_field(x, y, z)[i] : REAL(0);
			}
			inline REAL old_temp(long long x, long long y, long long z) { return main_field::temperature_field(x, y, z); }
			inline const Vector3& phi_gradient(long long x, long long y, long long z, size_t i) { auto c = Dendrite_workspace.cell(x, y, z); return Dendrite_workspace.grad_phi[Dendrite_workspace.phi_offset(c, i)]; }
			inline REAL phi_laplacian(long long x, long long y, long long z, size_t i) { auto c = Dendrite_workspace.cell(x, y, z); return Dendrite_workspace.lap_phi_or_rhs[Dendrite_workspace.phi_offset(c, i)]; }
			inline const Vector3& con_gradient(long long x, long long y, long long z, size_t i) { auto c = Dendrite_workspace.cell(x, y, z); return Dendrite_workspace.grad_con[Dendrite_workspace.con_offset(c, i)]; }
			inline REAL con_laplacian(long long x, long long y, long long z, size_t i) { auto c = Dendrite_workspace.cell(x, y, z); return Dendrite_workspace.lap_con_or_rhs[Dendrite_workspace.con_offset(c, i)]; }
			inline const Vector3& temp_gradient(long long x, long long y, long long z) { return Dendrite_workspace.grad_temp[Dendrite_workspace.cell(x, y, z)]; }
			inline REAL temp_laplacian(long long x, long long y, long long z) { return Dendrite_workspace.lap_temp[Dendrite_workspace.cell(x, y, z)]; }
			inline REAL phase_potential(long long x, long long y, long long z, size_t i) { auto c = Dendrite_workspace.cell(x, y, z); return Dendrite_workspace.mu_phi[Dendrite_workspace.phi_offset(c, i)]; }
			inline REAL concentration_potential(long long x, long long y, long long z, size_t i) { return Dendrite_workspace.mu_con(x, y, z, i); }
			inline REAL concentration_mobility(long long x, long long y, long long z, size_t i) { return Dendrite_workspace.mob_con(x, y, z, i); }
			inline REAL temperature_mobility(long long x, long long y, long long z) { return Dendrite_workspace.mob_temp(x, y, z); }
			inline REAL driving_force(long long x, long long y, long long z) { return Dendrite_workspace.delta_g[Dendrite_workspace.cell(x, y, z)]; }
			inline size_t active_phase_number(long long x, long long y, long long z) { return Dendrite_workspace.active_count[Dendrite_workspace.cell(x, y, z)]; }
			inline size_t active_phase_index(long long x, long long y, long long z, size_t slot) { auto c = Dendrite_workspace.cell(x, y, z); return Dendrite_workspace.active_indices[Dendrite_workspace.active_offset(c, slot)]; }
			inline InterfaceFlag active_phase_flag(long long x, long long y, long long z, size_t slot) { auto c = Dendrite_workspace.cell(x, y, z); return static_cast<InterfaceFlag>(Dendrite_workspace.active_flags[Dendrite_workspace.active_offset(c, slot)]); }
		}
	}
}
