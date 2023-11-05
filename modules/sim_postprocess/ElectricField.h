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
	enum ElectricFieldIndex { ElectricalConductivity = 100, ElectricalPotential, ChargeDensity};
	namespace electric_field_funcs {
		static void rhs_cal(pf::PhaseNode& node, int rhs_index) {
			node.customValues[rhs_index] = 0.0;
		}
		static void lhs_cal(pf::PhaseNode& node, int lhs_index) {
			node.customValues[lhs_index] = 0.0;
		}
		static void boundary(pf::PhaseNode& node, int r_index) {
			// node.customValues[r_index] = custom;
		}
	}
	namespace electric_field {
		static PoissonEquationSolver electric_field_solver;
		static vector<bool> fix_domain_boundary;
		static vector<double> fix_domain_boundary_value;
		static tensor1_double fix_domain_phi_value;
		static tensor1_double conductivity_phi;
		static double solver_accuracy = 1e-3;
		static int solver_max_iterate_times = 100;
		static bool solver_debug = false;
		static int solver_debug_output_steps = 100;
		// - 
		static int Nx = 1;
		static int Ny = 1;
		static int Nz = 1;
		static void charge_density_phix(pf::PhaseNode& node, int rhs_index) {
			node.customValues[rhs_index] = 0.0;
		}
		static void charge_density_x(pf::PhaseNode& node, int rhs_index) {
			node.customValues[rhs_index] = 0.0;
		}
		static void conductivity(pf::PhaseNode& node, int lhs_index) {
			double conductivity = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++)
				conductivity += phase->phi * conductivity_phi(phase->property);
			node.customValues[lhs_index] = conductivity;
		}
		static void boundary(pf::PhaseNode& node, int r_index) {
			// node.customValues[r_index] = custom;
			if(fix_domain_boundary[Boundary::DOWN_X] && node._x == 0)
				node.customValues[r_index] = fix_domain_boundary_value[Boundary::DOWN_X];
			if (fix_domain_boundary[Boundary::UP_X] && node._x == Nx - 1)
				node.customValues[r_index] = fix_domain_boundary_value[Boundary::UP_X];
			if (fix_domain_boundary[Boundary::DOWN_Y] && node._y == 0)
				node.customValues[r_index] = fix_domain_boundary_value[Boundary::DOWN_Y];
			if (fix_domain_boundary[Boundary::UP_Y] && node._y == Ny - 1)
				node.customValues[r_index] = fix_domain_boundary_value[Boundary::UP_Y];
			if (fix_domain_boundary[Boundary::DOWN_Z] && node._z == 0)
				node.customValues[r_index] = fix_domain_boundary_value[Boundary::DOWN_Z];
			if (fix_domain_boundary[Boundary::UP_Z] && node._z == Nz - 1)
				node.customValues[r_index] = fix_domain_boundary_value[Boundary::UP_Z];
			// 1
			/*double phi_sum = 0.0, fix_conduc = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++)
				for (auto c_phi = fix_domain_phi_value.begin(); c_phi < fix_domain_phi_value.end(); c_phi++)
					if (phase->property == c_phi->index) {
						fix_conduc += phase->phi * c_phi->val;
						phi_sum += phase->phi;
					}
			if(phi_sum > SYS_EPSILON)
				fix_conduc /= phi_sum;
			node.customValues[r_index] = node.customValues[r_index] * (1.0 - phi_sum) + fix_conduc * phi_sum;*/
			// 2
			for (auto phase = node.begin(); phase < node.end(); phase++)
				for (auto c_phi = fix_domain_phi_value.begin(); c_phi < fix_domain_phi_value.end(); c_phi++)
					if (phase->property == c_phi->index && phase->phi > (1.0 - Phi_Num_Cut_Off)) {
						node.customValues[r_index] = c_phi->val;
					}
		}
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			electric_field_solver.init_field(Solvers::get_instance()->phaseMesh, ElectricalConductivity, ElectricalPotential, ChargeDensity, "ElectricFieldSolver");
			fix_domain_boundary.resize(6); 
			fix_domain_boundary[Boundary::DOWN_X] = false; 
			fix_domain_boundary[Boundary::DOWN_Y] = false; 
			fix_domain_boundary[Boundary::DOWN_Z] = false;
			fix_domain_boundary[Boundary::UP_X] = false; 
			fix_domain_boundary[Boundary::UP_Y] = false; 
			fix_domain_boundary[Boundary::UP_Z] = false;
			fix_domain_boundary_value.resize(6); 
			fix_domain_boundary_value[Boundary::DOWN_X] = 0.0; 
			fix_domain_boundary_value[Boundary::DOWN_Y] = 0.0; 
			fix_domain_boundary_value[Boundary::DOWN_Z] = 0.0;
			fix_domain_boundary_value[Boundary::UP_X] = 0.0; 
			fix_domain_boundary_value[Boundary::UP_Y] = 0.0; 
			fix_domain_boundary_value[Boundary::UP_Z] = 0.0;
			electric_field_solver.set_LHS_calfunc(conductivity);
			electric_field_solver.set_BoundaryCondition_calfunc(boundary);
			electric_field_solver.set_field_variable(1.0, 0.0, 0.0);
			Nx = phaseMesh.limit_x;
			Ny = phaseMesh.limit_y;
			Nz = phaseMesh.limit_z;
			InputFileReader::get_instance()->read_double_value("Modules.ElectricField.accuracy", solver_accuracy, infile_debug);
			InputFileReader::get_instance()->read_int_value("Modules.ElectricField.max_iteration_steps", solver_max_iterate_times, infile_debug);
			if(InputFileReader::get_instance()->read_bool_value("Modules.ElectricField.debug", solver_debug, infile_debug))
				InputFileReader::get_instance()->read_int_value("Modules.ElectricField.Debug.output_steps", solver_debug_output_steps, infile_debug);
			// (DOWN_X, UP_X, DOWN_Y, UP_Y, DOWN_Z, UP_Z)
			string fix_boundary_key = "Modules.ElectricField.fix_boundary.type", fix_boundary_input = "(false,false,false,false,false,false)";
			InputFileReader::get_instance()->read_string_value(fix_boundary_key, fix_boundary_input, infile_debug);
			vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_BOOL, fix_boundary_key, fix_boundary_input, infile_debug);
			fix_domain_boundary[Boundary::DOWN_X]	= fix_boundary_value[0].bool_value;
			fix_domain_boundary[Boundary::UP_X]		= fix_boundary_value[1].bool_value;
			fix_domain_boundary[Boundary::DOWN_Y]	= fix_boundary_value[2].bool_value;
			fix_domain_boundary[Boundary::UP_Y]		= fix_boundary_value[3].bool_value;
			fix_domain_boundary[Boundary::DOWN_Z]	= fix_boundary_value[4].bool_value;
			fix_domain_boundary[Boundary::UP_Z]		= fix_boundary_value[5].bool_value;
			string fix_boundary_val_key = "Modules.ElectricField.fix_boundary.value", fix_boundary_val_input = "(0,0,0,0,0,0)";
			InputFileReader::get_instance()->read_string_value(fix_boundary_val_key, fix_boundary_val_input, infile_debug);
			vector<input_value> fix_boundary_val_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, fix_boundary_val_key, fix_boundary_val_input, infile_debug);
			fix_domain_boundary_value[Boundary::DOWN_X] = fix_boundary_val_value[0].double_value;
			fix_domain_boundary_value[Boundary::UP_X] = fix_boundary_val_value[1].double_value;
			fix_domain_boundary_value[Boundary::DOWN_Y] = fix_boundary_val_value[2].double_value;
			fix_domain_boundary_value[Boundary::UP_Y] = fix_boundary_val_value[3].double_value;
			fix_domain_boundary_value[Boundary::DOWN_Z] = fix_boundary_val_value[4].double_value;
			fix_domain_boundary_value[Boundary::UP_Z] = fix_boundary_val_value[5].double_value;
			string conductivity_phi_key = "Modules.ElectricField.conductivity", conductivity_phi_input = "()";
			InputFileReader::get_instance()->read_string_value(conductivity_phi_key, conductivity_phi_input, infile_debug);
			vector<input_value> conductivity_phi_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, conductivity_phi_key, conductivity_phi_input, infile_debug);
			int index = 0;
			for(auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++){
				conductivity_phi.add_double(phi->phi_property, conductivity_phi_value[index].double_value);
				index++;
			}
			string fix_phi_key = "Modules.ElectricField.fix_phi", fix_phi_input = "[()]";
			if (InputFileReader::get_instance()->read_string_value(fix_phi_key, fix_phi_input, infile_debug)) {
				vector<InputValueType> fix_phi_structure; fix_phi_structure.push_back(InputValueType::IVType_STRING); fix_phi_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> fix_phi_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(fix_phi_structure, fix_phi_key, fix_phi_input, infile_debug);
					for (auto fix_phi = fix_phi_value.begin(); fix_phi < fix_phi_value.end(); fix_phi++)
						fix_domain_phi_value.add_double(Solvers::get_instance()->parameters.Phases[(*fix_phi)[0].string_value].phi_property, (*fix_phi)[1].double_value);
			}
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE: ElectricField has been initialized !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			int iterate_steps = electric_field_solver.solve_whole_domain(solver_accuracy, solver_max_iterate_times, solver_debug, solver_debug_output_steps);
			output << "> Solver :" << electric_field_solver.solver_name << ", iterate " << iterate_steps << " times." << endl;
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream report;
			int iterate_steps = electric_field_solver.solve_whole_domain(solver_accuracy, solver_max_iterate_times, solver_debug, solver_debug_output_steps);
			report << "> Solver :" << electric_field_solver.solver_name << ", iterate " << iterate_steps << " times." << endl;
			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			fout << "<DataArray type = \"Float64\" Name = \"" << "elec_conductivity" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
			for(int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						fout << node.customValues[ElectricalConductivity] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "elec_potential" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						fout << node.customValues[ElectricalPotential] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "charge_density" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						fout << node.customValues[ChargeDensity] << endl;
					}
			fout << "</DataArray>" << endl;
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}

		static void init_module() {
			Solvers::get_instance()->create_a_new_module(init, exec_pre, exec_loop, deinit, write_scalar, write_vec3);
		};
	}
}