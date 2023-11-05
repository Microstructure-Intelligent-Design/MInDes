/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include "../../Base.h"

namespace pf {
	enum MagneticEnergyType {
		MagType_CONST,   // no magnetic potential
		MagType_STANDARD,  // standard electric energy
	};
	namespace magnetic_energy {
		static double dfmag_dphi_const(pf::PhaseNode& node, pf::PhaseEntry& phase){
			return 0.0;
		}
		static void dfmag_dphase_con_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static double dfmag_dcon_i_const(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}

		static double (*fmag_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*dfmag_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of this model

		static void (*dfmag_dphase_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of phase con

		static double (*dfmag_dcon_i)(pf::PhaseNode& node, int con_i);  // main function of total con & grand potential

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			fmag_density = dfmag_dphi_const;
			dfmag_dphi = dfmag_dphi_const;
			dfmag_dphase_con = dfmag_dphase_con_const;
			dfmag_dcon_i = dfmag_dcon_i_const;

			
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}