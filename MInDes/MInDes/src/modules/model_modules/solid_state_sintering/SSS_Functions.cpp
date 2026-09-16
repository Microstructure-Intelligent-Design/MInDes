#include "SSS_Functions.h"
#include <algorithm>
#include <cmath>

namespace pf {
    namespace solid_state_sintering_model {
        namespace {
            REAL laplacian_con(long long x, long long y, long long z) {
                return (main_field::concentration_field.at(x + 1, y, z)[0]
                    + main_field::concentration_field.at(x - 1, y, z)[0]
                    + main_field::concentration_field.at(x, y + 1, z)[0]
                    + main_field::concentration_field.at(x, y - 1, z)[0]
                    + main_field::concentration_field.at(x, y, z + 1)[0]
                    + main_field::concentration_field.at(x, y, z - 1)[0]
                    - 6 * main_field::concentration_field.at(x, y, z)[0]) 
                    / (mesh_parameters::delt_r * mesh_parameters::delt_r);
            }

            REAL laplacian_eta(long long x, long long y, long long z, size_t i) {
                return (main_field::phase_field.at(x + 1, y, z)[i]
                    + main_field::phase_field.at(x - 1, y, z)[i]
                    + main_field::phase_field.at(x, y + 1, z)[i]
                    + main_field::phase_field.at(x, y - 1, z)[i]
                    + main_field::phase_field.at(x, y, z + 1)[i]
                    + main_field::phase_field.at(x, y, z - 1)[i]
                    - 6 * main_field::phase_field.at(x, y, z)[i]) 
                    / (mesh_parameters::delt_r * mesh_parameters::delt_r);
            }

            REAL laplacian_mu(long long x, long long y, long long z) {
                return (parameters::concentration_variables.at(x + 1, y, z)[FIELD::CHEMICAL_POTENTIAL]
                    + parameters::concentration_variables.at(x - 1, y, z)[FIELD::CHEMICAL_POTENTIAL]
                    + parameters::concentration_variables.at(x, y + 1, z)[FIELD::CHEMICAL_POTENTIAL]
                    + parameters::concentration_variables.at(x, y - 1, z)[FIELD::CHEMICAL_POTENTIAL]
                    + parameters::concentration_variables.at(x, y, z + 1)[FIELD::CHEMICAL_POTENTIAL]
                    + parameters::concentration_variables.at(x, y, z - 1)[FIELD::CHEMICAL_POTENTIAL]
                    - 6 * parameters::concentration_variables.at(x, y, z)[FIELD::CHEMICAL_POTENTIAL]) 
                    / (mesh_parameters::delt_r * mesh_parameters::delt_r);
            }

            REAL clamped_rho(REAL value) {
                return std::max(REAL(0.00001), std::min(REAL(0.9999), value));
            }

            REAL clamped_eta(REAL value) {
                return std::max(REAL(0.0001), std::min(REAL(0.9999), value));
            }
        }

        void exec_pre_iii() {
            parameters::concentration_variables.init(
                mesh_parameters::MESH_NX, mesh_parameters::MESH_NY,
                mesh_parameters::MESH_NZ, mesh_parameters::delt_r,
                mesh_parameters::x_down, mesh_parameters::x_up,
                mesh_parameters::y_down, mesh_parameters::y_up,
                mesh_parameters::z_down, mesh_parameters::z_up);
            parameters::eta_increment.init(
                mesh_parameters::MESH_NX, mesh_parameters::MESH_NY,
                mesh_parameters::MESH_NZ, mesh_parameters::delt_r,
                mesh_parameters::x_down, mesh_parameters::x_up,
                mesh_parameters::y_down, mesh_parameters::y_up,
                mesh_parameters::z_down, mesh_parameters::z_up);
            for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
                for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
                    for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
                        parameters::concentration_variables.at(x, y, z).resize(FIELD::NUM, 0);
                        parameters::eta_increment.at(x, y, z).resize(main_field::phi_number, 0);
                    }
            main_field::phase_field.init_boundary_condition();
            main_field::concentration_field.init_boundary_condition();
        }

        void exec_i() {
            main_field::phase_field.do_boundary_condition();
            main_field::concentration_field.do_boundary_condition();
            // MATLAB: compute mu and the pointwise mobility from the old rho and eta.
#pragma omp parallel for
            for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
                for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
                    for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
                        const Matrix1D<REAL>& eta = main_field::phase_field.at(x, y, z);
                        const REAL rho = main_field::concentration_field.at(x, y, z)[0];
                        REAL eta_sum = 0, eta_square_sum = 0, eta_cube_sum = 0;
                        for (size_t i = 0; i < main_field::phi_number; ++i) {
                            eta_sum += eta[i];
                            eta_square_sum += eta[i] * eta[i];
                            eta_cube_sum += eta[i] * eta[i] * eta[i];
                        }
                        const REAL phi = rho * rho * rho * (10 - 15 * rho + 6 * rho * rho);
                        Matrix1D<REAL>& variables = parameters::concentration_variables.at(x, y, z);
                        variables[FIELD::MOB] = parameters::Dvol * phi
                            + parameters::Dvap * (1 - phi)
                            + parameters::Dsurf * rho * (1 - rho)
                            + parameters::Dgb * (eta_sum * eta_sum - eta_square_sum);
                        const REAL dfd_rho = parameters::B
                            * (2 * rho + 4 * eta_cube_sum - 6 * eta_square_sum)
                            + 2 * parameters::A * rho * (1 - rho) * (1 - 2 * rho);
                        variables[FIELD::CHEMICAL_POTENTIAL] = dfd_rho
                            - REAL(0.5) * parameters::kappa_rho * laplacian_con(x, y, z);
                    }
            parameters::concentration_variables.init_boundary_condition();
            parameters::concentration_variables.do_boundary_condition();

            REAL max_con_delta = 0;
#pragma omp parallel
            {
                REAL thread_max_con_delta = 0;
#pragma omp for
                for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); ++x)
                    for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); ++y)
                        for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); ++z) {
                            const Matrix1D<REAL>& variables = parameters::concentration_variables.at(x, y, z);
                            REAL& rho = main_field::concentration_field.at(x, y, z)[0];
                            const REAL old_rho = rho;
                            rho = clamped_rho(rho + time_parameters::delt_t
                                * variables[FIELD::MOB] * laplacian_mu(x, y, z));
                            thread_max_con_delta = std::max(thread_max_con_delta, std::abs(rho - old_rho));
                        }
#pragma omp critical(sss_con_max)
                max_con_delta = std::max(max_con_delta, thread_max_con_delta);
            }
            main_field::CON_MAX_VARIATION = std::max(main_field::CON_MAX_VARIATION, max_con_delta);
            main_field::concentration_field.do_boundary_condition();

            // MATLAB updates eta components in order, using the new rho and already updated eta components.
            REAL max_eta_delta = 0;
            for (size_t i = 0; i < main_field::phi_number; ++i) {
                main_field::phase_field.do_boundary_condition();
#pragma omp parallel for
                for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
                    for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
                        for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
                            const Matrix1D<REAL>& eta = main_field::phase_field.at(x, y, z);
                            const REAL rho = main_field::concentration_field.at(x, y, z)[0];
                            REAL eta_square_sum = 0;
                            for (size_t j = 0; j < main_field::phi_number; ++j)
                                eta_square_sum += eta[j] * eta[j];
                            const REAL dfd_eta = 12 * parameters::B * eta[i]
                                * ((1 - rho) - (2 - rho) * eta[i] + eta_square_sum);
                            parameters::eta_increment.at(x, y, z)[i] = -parameters::L
                                * (dfd_eta - REAL(0.5) * parameters::kappa_eta
                                    * laplacian_eta(x, y, z, i));
                        }
#pragma omp parallel
                {
                    REAL thread_max_eta_delta = 0;
#pragma omp for
                    for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); ++x)
                        for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); ++y)
                            for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); ++z) {
                                REAL& eta = main_field::phase_field.at(x, y, z)[i];
                                const REAL old_eta = eta;
                                eta = clamped_eta(eta + time_parameters::delt_t
                                    * parameters::eta_increment.at(x, y, z)[i]);
                                thread_max_eta_delta = std::max(thread_max_eta_delta, std::abs(eta - old_eta));
                            }
#pragma omp critical(sss_eta_max)
                    max_eta_delta = std::max(max_eta_delta, thread_max_eta_delta);
                }
            }
            main_field::PHI_MAX_VARIATION = std::max(main_field::PHI_MAX_VARIATION, max_eta_delta);
        }

        void deinit() {
            parameters::concentration_variables.clear();
            parameters::eta_increment.clear();
            main_field::phase_field.clear();
            main_field::concentration_field.clear();
        }
    }
}
