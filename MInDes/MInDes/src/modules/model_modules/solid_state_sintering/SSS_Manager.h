#pragma once
#include <cmath>
#include <cstdlib>
#include "../../Module.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "SSS_Params.h"
#include "SSS_Functions.h"

namespace pf {
    namespace solid_state_sintering_model {
        inline void init_model_modules() {
            if (main_field::phi_number < 1 || main_field::con_number != 1 ||
                main_field::is_temp_field_on) {
                WriteLog("> ERROR, Solid State Sintering requires PCT = (N >= 1,1,false).\n");
                std::exit(EXIT_FAILURE);
            }

            WriteLog("> Simulation Model - Solid State Sintering - is Activated !\n");
            WriteDebugFile("# SSS MATLAB update: rho_t = D(rho,eta) lap(mu); mu = df/drho - 0.5*kappa_rho*lap(rho).\n");
            WriteDebugFile("# eta_i,t = -L*(df/deta_i - 0.5*kappa_eta*lap(eta_i)); eta_i are updated in index order.\n");
            WriteDebugFile("# D = Dvol*phi(rho) + Dvap*(1-phi(rho)) + Dsurf*rho*(1-rho) + Dgb*sum(i!=j)eta_i*eta_j.\n");
            WriteDebugFile("# phi(rho) = rho^3*(10-15*rho+6*rho^2); MATLAB bounds are applied after each field update.\n");

            infile_reader::read_real_value("SimulationModels.SolidStateSintering.A", parameters::A, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.B", parameters::B, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.kappa_rho", parameters::kappa_rho, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.kappa_eta", parameters::kappa_eta, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.L", parameters::L, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.Dvol", parameters::Dvol, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.Dvap", parameters::Dvap, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.Dsurf", parameters::Dsurf, true);
            infile_reader::read_real_value("SimulationModels.SolidStateSintering.Dgb", parameters::Dgb, true);

            if (!std::isfinite(parameters::A) || parameters::A < 0 ||
                !std::isfinite(parameters::B) || parameters::B < 0 ||
                !std::isfinite(parameters::kappa_rho) || parameters::kappa_rho < 0 ||
                !std::isfinite(parameters::kappa_eta) || parameters::kappa_eta < 0 ||
                !std::isfinite(parameters::L) || parameters::L < 0 ||
                !std::isfinite(parameters::Dvol) || parameters::Dvol < 0 ||
                !std::isfinite(parameters::Dvap) || parameters::Dvap < 0 ||
                !std::isfinite(parameters::Dsurf) || parameters::Dsurf < 0 ||
                !std::isfinite(parameters::Dgb) || parameters::Dgb < 0 ||
                !std::isfinite(mesh_parameters::delt_r) || mesh_parameters::delt_r <= 0 ||
                !std::isfinite(time_parameters::delt_t) || time_parameters::delt_t <= 0) {
                WriteLog("> ERROR, invalid Solid State Sintering parameters, grid spacing or time step.\n");
                std::exit(EXIT_FAILURE);
            }

            load_a_new_module(nullptr, nullptr, exec_pre_iii,
                exec_i, nullptr, nullptr, nullptr, nullptr, nullptr, deinit);
        }
    }
}
