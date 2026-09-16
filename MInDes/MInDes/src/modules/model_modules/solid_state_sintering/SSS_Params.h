#pragma once
#include "../../Modules_Params.h"

namespace pf {
    namespace solid_state_sintering_model {
        enum FIELD { MOB, CHEMICAL_POTENTIAL, NUM };
        namespace parameters {
            inline Mesh_Boundry<Matrix1D<REAL>> concentration_variables;
            inline Mesh_Boundry<Matrix1D<REAL>> eta_increment;

            // Defaults from Table 4.1 and the MATLAB listing in SSS模型.pdf.
            inline REAL A = 16.0;
            inline REAL B = 1.0;
            inline REAL kappa_rho = 5.0;
            inline REAL kappa_eta = 2.0;
            inline REAL L = 10.0;
            inline REAL Dvol = 0.040;
            inline REAL Dvap = 0.002;
            inline REAL Dsurf = 16.0;
            inline REAL Dgb = 1.6;
        }
    }
}
