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
// include base
#include "../solvers/Solvers.h"
#ifdef _WIN32
#include "base/BMP24Reader.h"
#endif
#include "base/NucleationTools.h"
#include "base/InputFileMath.h"
#include "base/InputFileReader.h"
#include "base/PresetFunctions.h"

namespace pf {
	enum ExternalFieldsPlus {
		DF_Macro_Density = ExternalFields::LBM_Symbols_INDEX_0
		, DF_Macro_TwoPhase = ExternalFields::LBM_Symbols_INDEX_0 + LBM_Symbols::LBM_SIZE * 2
		, DF_Macro_TwoPhase_f_macro_old = ExternalFields::LBM_Symbols_INDEX_0 + LBM_Symbols::LBM_SIZE * 3, DF_Macro_TwoPhase_velocity_old
		, EFP_Crack = 0, EFP_Crack_Incre = 1000, EFP_GrainEigenStrain_Region = 2000 };

	static inline double interpolation_func(double phi) {
		return phi * phi * phi * (6.0 * phi * phi - 15.0 * phi + 10.0);
	}

	static inline double dinterpolation_func_dphi(double phi) {
		return 30.0 * phi * phi * (1.0 - phi) * (1.0 - phi);
	}

}