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
// include solvers
#include "Base.h"
// include modules
#include "InputManager.h"
#include "SimPreprocess.h"
#include "SimPostprocess.h"
#include "SimModelsManager.h"

#ifdef _WIN32
//windowsƽ̨ x86 or x68
#ifdef _WIN64
 //x64
#ifdef _DEBUG
//#pragma comment(lib,"lib/x64/debug/MID.lib")
#else
//#pragma comment(lib,"lib/x64/release/MID.lib")
#endif
#else
 //x86
#ifdef _DEBUG
//#pragma comment(lib,"lib/x86/debug/MID.lib")
#else
//#pragma comment(lib,"lib/x86/release/MID.lib")
#endif
#endif //_WIN64
#endif //_WIN32

namespace pf {
	static void init_modules(string input_file_path, int input_file_lines_buff = 1000, char split = ' ') {
		// load modules
		pf::input_manager::load_module(input_file_path, input_file_lines_buff, split);  // work in init() to read material/physical functions from input file  ( need to be reconstructed from InputFileMath.h and InputFileReader.h )

		pf::sim_models_manager::load_module();  // work in init(), to choose a phase-field model and load material/physical functions

		pf::sim_preprocess::load_module();  // mainly work in init() and exec_pre() functions, prepare for simulation

		pf::sim_postprocess::load_module();  // mainly work in exec_pre() and exec_loop() functions , work with simulation and prepare for special output

	}
}