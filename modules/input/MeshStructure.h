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

namespace pf {
	namespace mesh_structure {
		static int Nx = 1;
		static int Ny = 1;
		static int Nz = 1;
		static double dr = 1.0;
		static int x_up_bc = BoundaryCondition::PERIODIC;
		static int y_up_bc = BoundaryCondition::PERIODIC;
		static int z_up_bc = BoundaryCondition::PERIODIC;
		static int x_down_bc = BoundaryCondition::PERIODIC;
		static int y_down_bc = BoundaryCondition::PERIODIC;
		static int z_down_bc = BoundaryCondition::PERIODIC;
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			InputFileReader::get_instance()->read_int_value("Solver.Mesh.Nx", Nx, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.Ny", Ny, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.Nz", Nz, infile_debug);
			InputFileReader::get_instance()->read_double_value("Solver.Mesh.dr", dr, infile_debug);

			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Solver.Mesh.BoundaryCondition : 0 - FIXED , 1 - PERIODIC , 2 - ADIABATIC\n", InputFileReader::get_instance()->debug_file);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_up", x_up_bc, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_down", x_down_bc, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_up", y_up_bc, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_down", y_down_bc, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_up", z_up_bc, infile_debug);
			InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_dowm", z_down_bc, infile_debug);

			// some dimension rules
			if (Nx == 0 || Ny == 0 || Nz == 0) {
				string error_report = "> error, one edge length (Nx or Ny or Nz) of domain is zero, domain does not exist ! Please set Nx > 0 and Ny > 0 and Nz > 0 !\n";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				std::exit(0);
			}
			else if (Nx == 1 && (Ny > 1 || Nz > 1)) {
				string error_report = "> error, edge length (Nx, Ny, Nz) of domain should be set (Nx > 1, Ny = 1, Nz = 1) in 1 dimension or (Nx > 1, Ny > 1, Nz = 1) in 2 dimension !\n";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				std::exit(0);
			}
			else if (Nx > 1 && Ny == 1 && Nz > 1) {
				string error_report = "> error, edge length (Nx, Ny, Nz) of domain should be set as (Nx > 1, Ny > 1, Nz = 1) in 2 dimension !\n";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				std::exit(0);
			}
			// init main mesh
			phaseMesh.init(Nx, Ny, Nz, dr, BoundaryCondition(x_up_bc), BoundaryCondition(y_up_bc), BoundaryCondition(z_up_bc), BoundaryCondition(x_down_bc), BoundaryCondition(y_down_bc), BoundaryCondition(z_down_bc));
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			phaseMesh.free();
		}
	}
}