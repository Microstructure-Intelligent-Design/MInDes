#include "DS_Functions.h"
#include <algorithm>
#include <cmath>

namespace pf {
	namespace dendrite_solidification_model {
		using namespace parameters;
		namespace {

			void init_scalar_mesh(Mesh_Boundry<REAL>& mesh) {
				mesh.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ,
					mesh_parameters::delt_r, mesh_parameters::x_down, mesh_parameters::x_up,
					mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down,
					mesh_parameters::z_up);
			}

			// On a nonperiodic ghost layer use a one-sided normal derivative of the
			// actual phase ghost value. Tangential derivatives remain centered.
			void evaluate_anisotropy(long long x, long long y) {
				auto& phi = main_field::phase_field;
				const long long nx = static_cast<long long>(mesh_parameters::MESH_NX);
				const long long ny = static_cast<long long>(mesh_parameters::MESH_NY);
				const REAL h = mesh_parameters::delt_r;
				const REAL gx = x == 0 ? (phi(x + 1, y, 1LL)[0] - phi(x, y, 1LL)[0]) / h
					: x == nx + 1 ? (phi(x, y, 1LL)[0] - phi(x - 1, y, 1LL)[0]) / h
					: (phi(x + 1, y, 1LL)[0] - phi(x - 1, y, 1LL)[0]) / (2 * h);
				const REAL gy = y == 0 ? (phi(x, y + 1, 1LL)[0] - phi(x, y, 1LL)[0]) / h
					: y == ny + 1 ? (phi(x, y, 1LL)[0] - phi(x, y - 1, 1LL)[0]) / h
					: (phi(x, y + 1, 1LL)[0] - phi(x, y - 1, 1LL)[0]) / (2 * h);
				const REAL angle = std::atan2(gy, gx) - theta0;
				const REAL epsilon = epsilonb * (1 + delta * std::cos(aniso * angle));
				const REAL epsilon_theta = -epsilonb * aniso * delta * std::sin(aniso * angle);
				auto& values = anisotropy(x, y, 1LL);
				values[EPSILON] = epsilon;
				values[FLUX_X] = epsilon * epsilon_theta * gx;
				values[FLUX_Y] = epsilon * epsilon_theta * gy;
			}

			void compute_anisotropy() {
				const long long nx = static_cast<long long>(mesh_parameters::MESH_NX);
				const long long ny = static_cast<long long>(mesh_parameters::MESH_NY);
#pragma omp parallel for
				for (long long x = 1; x <= nx; ++x)
					for (long long y = 1; y <= ny; ++y)
						evaluate_anisotropy(x, y);

				// Only side ghosts are needed by centered derivatives of the interior flux.
				for (long long y = 1; y <= ny; ++y) {
					if (main_field::phase_field.BC_X_DOWN() == BoundaryCondition::PERIODIC)
						anisotropy(0LL, y, 1LL) = anisotropy(nx, y, 1LL);
					else
						evaluate_anisotropy(0, y);
					if (main_field::phase_field.BC_X_UP() == BoundaryCondition::PERIODIC)
						anisotropy(nx + 1, y, 1LL) = anisotropy(1LL, y, 1LL);
					else
						evaluate_anisotropy(nx + 1, y);
				}
				for (long long x = 1; x <= nx; ++x) {
					if (main_field::phase_field.BC_Y_DOWN() == BoundaryCondition::PERIODIC)
						anisotropy(x, 0LL, 1LL) = anisotropy(x, ny, 1LL);
					else
						evaluate_anisotropy(x, 0);
					if (main_field::phase_field.BC_Y_UP() == BoundaryCondition::PERIODIC)
						anisotropy(x, ny + 1, 1LL) = anisotropy(x, 1LL, 1LL);
					else
						evaluate_anisotropy(x, ny + 1);
				}
			}
		}

		void exec_pre_iii() {
			anisotropy.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ,
				mesh_parameters::delt_r, mesh_parameters::x_down, mesh_parameters::x_up,
				mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down,
				mesh_parameters::z_up);
			init_scalar_mesh(phi_increment);
			init_scalar_mesh(temperature_increment);
			for (long long x = 0; x < anisotropy.Nx(); ++x)
				for (long long y = 0; y < anisotropy.Ny(); ++y)
					for (long long z = 0; z < anisotropy.Nz(); ++z)
						anisotropy(x, y, z).resize(NUM, 0);
			// FIXED boundaries capture the initialized values exactly once.
			main_field::phase_field.init_boundary_condition();
			main_field::temperature_field.init_boundary_condition();
			main_field::phase_field.do_boundary_condition();
			main_field::temperature_field.do_boundary_condition();
		}

		void exec_i() {
			auto& phi = main_field::phase_field;
			auto& temperature = main_field::temperature_field;
			phi.do_boundary_condition();
			temperature.do_boundary_condition();
			compute_anisotropy();
			const REAL h = mesh_parameters::delt_r;
			const REAL inv_h2 = 1 / (h * h);
			const REAL inv_2h = 1 / (2 * h);
			const REAL dt = time_parameters::delt_t;
			const REAL pi = std::acos(REAL(-1));
			const long long nx = static_cast<long long>(mesh_parameters::MESH_NX);
			const long long ny = static_cast<long long>(mesh_parameters::MESH_NY);
#pragma omp parallel for
			for (long long x = 1; x <= nx; ++x)
				for (long long y = 1; y <= ny; ++y) {
					const REAL old_phi = phi(x, y, 1LL)[0];
					const REAL old_temp = temperature(x, y, 1LL);
					const REAL lap_phi = (phi(x + 1, y, 1LL)[0] + phi(x - 1, y, 1LL)[0] +
						phi(x, y + 1, 1LL)[0] + phi(x, y - 1, 1LL)[0] - 4 * old_phi) * inv_h2;
					const REAL lap_temp = (temperature(x + 1, y, 1LL) + temperature(x - 1, y, 1LL) +
						temperature(x, y + 1, 1LL) + temperature(x, y - 1, 1LL) - 4 * old_temp) * inv_h2;
					const REAL directional = (anisotropy(x, y + 1, 1LL)[FLUX_X] -
						anisotropy(x, y - 1, 1LL)[FLUX_X] -
						anisotropy(x + 1, y, 1LL)[FLUX_Y] +
						anisotropy(x - 1, y, 1LL)[FLUX_Y]) * inv_2h;
					const REAL epsilon = anisotropy(x, y, 1LL)[EPSILON];
					const REAL m = alpha / pi * std::atan(gamma * (teq - old_temp));
					// Match the reference MATLAB code: epsilon^2 * lap(phi), not
					// the divergence form printed in the accompanying text.
					const REAL dphi = dt / tau * (directional + epsilon * epsilon * lap_phi +
						old_phi * (1 - old_phi) * (old_phi - REAL(0.5) + m));
					phi_increment(x, y, 1LL) = dphi;
					temperature_increment(x, y, 1LL) = dt * lap_temp + kappa * dphi;
				}

			REAL max_phi = 0;
			REAL max_temp = 0;
#pragma omp parallel
			{
				REAL local_phi = 0;
				REAL local_temp = 0;
#pragma omp for nowait
				for (long long x = 1; x <= nx; ++x)
					for (long long y = 1; y <= ny; ++y) {
						const REAL dphi = phi_increment(x, y, 1LL);
						const REAL dtemp = temperature_increment(x, y, 1LL);
						phi(x, y, 1LL)[0] += dphi;
						temperature(x, y, 1LL) += dtemp;
						local_phi = std::max(local_phi, std::abs(dphi));
						local_temp = std::max(local_temp, std::abs(dtemp));
					}
#pragma omp critical
				{
					max_phi = std::max(max_phi, local_phi);
					max_temp = std::max(max_temp, local_temp);
				}
			}
			main_field::PHI_MAX_VARIATION = std::max(main_field::PHI_MAX_VARIATION, max_phi);
			main_field::TEMP_MAX_VARIATION = std::max(main_field::TEMP_MAX_VARIATION, max_temp);
		}

		void deinit() {
			anisotropy.clear();
			phi_increment.clear();
			temperature_increment.clear();
			main_field::phase_field.clear();
			main_field::temperature_field.clear();
		}
	}
}



