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
#include "../Base.h"
#include "FluidDynamics/PressureCorrection.h"
#include "FluidDynamics/LatticeBoltzmann.h"

namespace pf {
	namespace fluid_field {
		static Fluid_Solver_Type fluid_solver_type = Fluid_Solver_Type::FST_NONE;
		// statement functons
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			int _fluid_solver_type = Fluid_Solver_Type::FST_NONE;
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.FluidDynamics.solver = 0 - None , 1 - Pressure_Correction , 2 - Lattice_Boltzmann \n", InputFileReader::get_instance()->debug_file);
			InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.solver", _fluid_solver_type, infile_debug);
			fluid_solver_type = Fluid_Solver_Type(_fluid_solver_type);
			if (fluid_solver_type == Fluid_Solver_Type::FST_Pressure_Correction) {
				pressure_correction::init(phaseMesh);
			}
			else if (fluid_solver_type == Fluid_Solver_Type::FST_Lattice_Boltzmann) {
				lattice_boltzmann::init(phaseMesh);
			}
			
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : FluidField !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (fluid_solver_type == Fluid_Solver_Type::FST_Pressure_Correction)
				pressure_correction::exec_pre(phaseMesh);
			else if (fluid_solver_type == Fluid_Solver_Type::FST_Lattice_Boltzmann)
				lattice_boltzmann::exec_pre(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			if (fluid_solver_type == Fluid_Solver_Type::FST_Pressure_Correction)
				report = pressure_correction::exec_loop(phaseMesh);
			else if (fluid_solver_type == Fluid_Solver_Type::FST_Lattice_Boltzmann)
				report = lattice_boltzmann::exec_loop(phaseMesh);
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			if (fluid_solver_type == Fluid_Solver_Type::FST_Pressure_Correction) {
				pressure_correction::deinit(phaseMesh);
			}
			else if (fluid_solver_type == Fluid_Solver_Type::FST_Lattice_Boltzmann) {
				lattice_boltzmann::deinit(phaseMesh);
			}
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (fluid_solver_type == Fluid_Solver_Type::FST_Pressure_Correction) {
				pressure_correction::write_scalar(fout, phaseMesh);
			}
			else if (fluid_solver_type == Fluid_Solver_Type::FST_Lattice_Boltzmann) {
				lattice_boltzmann::write_scalar(fout, phaseMesh);
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (fluid_solver_type == Fluid_Solver_Type::FST_Pressure_Correction) {
				pressure_correction::write_vec3(fout, phaseMesh);
			}
			else if (fluid_solver_type == Fluid_Solver_Type::FST_Lattice_Boltzmann) {
				lattice_boltzmann::write_vec3(fout, phaseMesh);
			}
		}
	}
}