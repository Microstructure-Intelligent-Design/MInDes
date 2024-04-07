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
	enum ElectricFieldSolver { EFS_NONE, EFS_EXPLICITE_DIFFERENCE, EFS_EXPLICITE_FOURIER_SPECTRAL, EFS_IMPLICIT_FOURIER_SPECTRAL };
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
		static int electric_field_solver = ElectricFieldSolver::EFS_NONE;
		// - 
		static vector<int> electrode_index{};
		static int active_component_index = 0;
		// - 
		static double time_interval = 0.0;
		// - 
		static double reaction_constant{};
		static double electron_num{};
		static double E_std{};
		static double c_s{}, c_0{};
		// - 
		static PoissonEquationSolver_Explicit electric_field_solver_diff;
		static FourierTransformSolver electric_field_solver_fourier;
		static double im_dt = 0.0;
		static double im_varient = 0.0;
		static double im_max_varient = 0.0;
		static string im_solver_name = "ElectricFieldSolver_fourier";
		// -
		static vector<bool> fix_domain_boundary;
		static vector<double> fix_domain_boundary_value;
		static tensor1_double fix_domain_phi_value;
		static tensor1_double conductivity_phi;
		static double conductivity_background = 0.0;
		static double solver_accuracy = 1e-3;
		static int solver_max_iterate_times = 100;
		static bool solver_debug = false;
		static int solver_debug_output_steps = 100;
		static double threshold{1.0-1e-3};
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
			double conductivity = 0.0, phi = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++) {
				conductivity += phase->phi * conductivity_phi(phase->property);
				phi += phase->phi;
			}
			if (phi > SYS_EPSILON)
				conductivity = conductivity / phi;
			double h_func = interpolation_func(phi);
			node.customValues[lhs_index] = conductivity * h_func + conductivity_background * (1.0 - h_func);
		}
		
		// From http://dx.doi.org/10.1016/j.jpowsour.2015.09.055
		static void reaction_term(pf::PhaseNode& node, int rhs_index) {
			double dxi_dt{};
			for (auto phase = node.begin(); phase < node.end(); phase++)
				for (auto index = electrode_index.begin(); index < electrode_index.end(); index++)
					if (*index == phase->index) {
						dxi_dt += (phase->phi - phase->old_phi)/ time_interval;
					}
			double result{ c_s * FaradayConstant * electron_num * dxi_dt };

			node.customValues[rhs_index] = result;
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
					if (phase->property == c_phi->index && phase->phi > threshold) {
						node.customValues[r_index] = c_phi->val;
					}
		}

		static void init_real_space(std::complex<double>& basic_real_space, PhaseNode& node, int real_x, int real_y, int real_z) {
			node.customValues[ElectricalPotential] = 0.0;
			// node.customValues[r_index] = custom;
			if (fix_domain_boundary[Boundary::DOWN_X] && node._x == 0)
				node.customValues[ElectricalPotential] = fix_domain_boundary_value[Boundary::DOWN_X];
			if (fix_domain_boundary[Boundary::UP_X] && node._x == Nx - 1)
				node.customValues[ElectricalPotential] = fix_domain_boundary_value[Boundary::UP_X];
			if (fix_domain_boundary[Boundary::DOWN_Y] && node._y == 0)
				node.customValues[ElectricalPotential] = fix_domain_boundary_value[Boundary::DOWN_Y];
			if (fix_domain_boundary[Boundary::UP_Y] && node._y == Ny - 1)
				node.customValues[ElectricalPotential] = fix_domain_boundary_value[Boundary::UP_Y];
			if (fix_domain_boundary[Boundary::DOWN_Z] && node._z == 0)
				node.customValues[ElectricalPotential] = fix_domain_boundary_value[Boundary::DOWN_Z];
			if (fix_domain_boundary[Boundary::UP_Z] && node._z == Nz - 1)
				node.customValues[ElectricalPotential] = fix_domain_boundary_value[Boundary::UP_Z];
			// 2
			for (auto phase = node.begin(); phase < node.end(); phase++)
				for (auto c_phi = fix_domain_phi_value.begin(); c_phi < fix_domain_phi_value.end(); c_phi++)
					if (phase->property == c_phi->index && phase->phi > threshold) {
						node.customValues[ElectricalPotential] = c_phi->val;
					}
			// - 
			basic_real_space._Val[FFTW_REAL] = node.customValues[ElectricalPotential];
			// basic_real_space._Val[FFTW_IMAG] = 0.0;
		}

		static void fill_node_real_space(vector<std::complex<double>>& real_space, PhaseNode& node, int real_x, int real_y, int real_z) {
			// - calculate real space

			// - calculate node
			if (real_x != node._x || real_y != node._y || real_z != node._z)
				return;
			double conductivity = 0.0, phi = 0.0;
			for (auto phase = node.begin(); phase < node.end(); phase++) {
				conductivity += phase->phi * conductivity_phi(phase->property);
				phi += phase->phi;
			}
			if (phi > SYS_EPSILON)
				conductivity = conductivity / phi;
			double h_func = interpolation_func(phi);
			node.customValues[ElectricalConductivity] = conductivity * h_func + conductivity_background * (1.0 - h_func);
		}

		static std::complex<double> dynamic_equation_fourier_space_ex(vector<std::complex<double>> fourier_space, pf::PhaseNode& node, double Q2, double Q4) {
			double parameter = 1.0 - Q2 * node.customValues[ElectricalConductivity] * im_dt;
			return fourier_space[0] * parameter;
		}

		static std::complex<double> dynamic_equation_fourier_space_im(vector<std::complex<double>> fourier_space, pf::PhaseNode& node, double Q2, double Q4) {
			double parameter = 1.0 + Q2 * node.customValues[ElectricalConductivity] * im_dt;
			return fourier_space[0] / parameter;
		}

		static void boundary_condition_real_space(std::complex<double>& basic_real_space, PhaseNode& node, int real_x, int real_y, int real_z) {
			// node.customValues[r_index] = custom;
			if (fix_domain_boundary[Boundary::DOWN_X] && node._x == 0)
				basic_real_space._Val[FFTW_REAL] = fix_domain_boundary_value[Boundary::DOWN_X];
			if (fix_domain_boundary[Boundary::UP_X] && node._x == Nx - 1)
				basic_real_space._Val[FFTW_REAL] = fix_domain_boundary_value[Boundary::UP_X];
			if (fix_domain_boundary[Boundary::DOWN_Y] && node._y == 0)
				basic_real_space._Val[FFTW_REAL] = fix_domain_boundary_value[Boundary::DOWN_Y];
			if (fix_domain_boundary[Boundary::UP_Y] && node._y == Ny - 1)
				basic_real_space._Val[FFTW_REAL] = fix_domain_boundary_value[Boundary::UP_Y];
			if (fix_domain_boundary[Boundary::DOWN_Z] && node._z == 0)
				basic_real_space._Val[FFTW_REAL] = fix_domain_boundary_value[Boundary::DOWN_Z];
			if (fix_domain_boundary[Boundary::UP_Z] && node._z == Nz - 1)
				basic_real_space._Val[FFTW_REAL] = fix_domain_boundary_value[Boundary::UP_Z];
			// 2
			for (auto phase = node.begin(); phase < node.end(); phase++)
				for (auto c_phi = fix_domain_phi_value.begin(); c_phi < fix_domain_phi_value.end(); c_phi++)
					if (phase->property == c_phi->index && phase->phi > threshold) {
						basic_real_space._Val[FFTW_REAL] = c_phi->val;
					}
			/*if (basic_real_space > 1.0)
				basic_real_space = 1.0;
			else if(basic_real_space < 0.0)
				basic_real_space = 0.0;*/
			// varient
			double varient = abs(node.customValues[ElectricalPotential] - basic_real_space._Val[FFTW_REAL]);
			if (varient > im_varient)
				im_varient = varient;
			// set
			node.customValues[ElectricalPotential] = basic_real_space._Val[FFTW_REAL];
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			time_interval = Solvers::get_instance()->parameters.dt;
			InputFileReader::get_instance()->read_int_value("Postprocess.PhysicalFields.electric", electric_field_solver, false);
			if (electric_field_solver == ElectricFieldSolver::EFS_NONE)
				return;
			std::string active_comp_name{};
			if (InputFileReader::get_instance()->read_string_value("ModelsManager.PhiCon.ElectroDeposition.active_component", active_comp_name, infile_debug) && electric_field_solver) {
				active_component_index = Solvers::get_instance()->parameters.Components[active_comp_name].index;

				string electrode_key = "ModelsManager.PhiCon.ElectroDeposition.electrode_index", electrode_input = "()";
				InputFileReader::get_instance()->read_string_value(electrode_key, electrode_input, infile_debug);
				vector<pf::input_value> electrode_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, electrode_key, electrode_input, infile_debug);
				for (int index = 0; index < electrode_value.size(); index++)
					electrode_index.push_back(electrode_value[index].int_value);

				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Butler_Volmer.Reaction_Constant", reaction_constant, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Butler_Volmer.Reaction_Electron_Num", electron_num, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Phi.Bulter_Volmer.Standard_Potential", E_std, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.Bulter_Volmer.Electrode_Metal_SiteDensity", c_s, infile_debug);
				InputFileReader::get_instance()->read_double_value("ModelsManager.Con.Bulter_Volmer.Electrolyte_Cation_Con", c_0, infile_debug);

				//electric_field_solver_ex.set_RHS_calfunc(reaction_term);
			}

			electric_field_solver_diff.init_field(Solvers::get_instance()->phaseMesh, ElectricalConductivity, ElectricalPotential, ChargeDensity, "ElectricFieldSolver_DIFF");
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

			if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_DIFFERENCE) {

				electric_field_solver_diff.set_LHS_calfunc(conductivity);

				electric_field_solver_diff.set_BoundaryCondition_calfunc(boundary);
				electric_field_solver_diff.set_field_variable(1.0, 0.0, 0.0);
			}
			else if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_FOURIER_SPECTRAL || electric_field_solver == ElectricFieldSolver::EFS_IMPLICIT_FOURIER_SPECTRAL) {
				InputFileReader::get_instance()->read_double_value("Modules.ElectricField.dt", im_dt, infile_debug);
				BoundaryCondition x_bc = BoundaryCondition::PERIODIC, y_bc = BoundaryCondition::PERIODIC, z_bc = BoundaryCondition::PERIODIC;
				if (phaseMesh._bc_x_up != BoundaryCondition::PERIODIC || phaseMesh._bc_x_down != BoundaryCondition::PERIODIC)
					x_bc = BoundaryCondition::ADIABATIC;
				if (phaseMesh._bc_y_up != BoundaryCondition::PERIODIC || phaseMesh._bc_y_down != BoundaryCondition::PERIODIC)
					y_bc = BoundaryCondition::ADIABATIC;
				if (phaseMesh._bc_z_up != BoundaryCondition::PERIODIC || phaseMesh._bc_z_down != BoundaryCondition::PERIODIC)
					z_bc = BoundaryCondition::ADIABATIC;
				electric_field_solver_fourier.init(phaseMesh, x_bc, y_bc, z_bc);
				electric_field_solver_fourier.init_real_space = init_real_space;
				electric_field_solver_fourier.fill_node_real_space = fill_node_real_space;
				if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_FOURIER_SPECTRAL)
					electric_field_solver_fourier.dynamic_equation_fourier_space = dynamic_equation_fourier_space_ex;
				else if(electric_field_solver == ElectricFieldSolver::EFS_IMPLICIT_FOURIER_SPECTRAL)
					electric_field_solver_fourier.dynamic_equation_fourier_space = dynamic_equation_fourier_space_im;
				electric_field_solver_fourier.boundary_condition_real_space = boundary_condition_real_space;
			}
			Nx = phaseMesh.limit_x;
			Ny = phaseMesh.limit_y;
			Nz = phaseMesh.limit_z;
			InputFileReader::get_instance()->read_double_value("Modules.ElectricField.accuracy", solver_accuracy, infile_debug);
			InputFileReader::get_instance()->read_double_value("Modules.ElectricField.threshold", threshold, infile_debug);
			InputFileReader::get_instance()->read_int_value("Modules.ElectricField.max_iteration_steps", solver_max_iterate_times, infile_debug);
			if(InputFileReader::get_instance()->read_bool_value("Modules.ElectricField.debug", solver_debug, infile_debug))
				InputFileReader::get_instance()->read_int_value("Modules.ElectricField.Debug.output_steps", solver_debug_output_steps, infile_debug);
			// (DOWN_X, UP_X, DOWN_Y, UP_Y, DOWN_Z, UP_Z)
			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Modules.ElectricField.fix_boundary.value = (x_down,x_up,y_down,y_up,z_down,z_up) \n", InputFileReader::get_instance()->debug_file);
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
			if (conductivity_phi_value.size() != Solvers::get_instance()->parameters.Phases.size()) {
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ERROR : size of Modules.ElectricField.conductivity isn't equal to the number of phases \n", InputFileReader::get_instance()->debug_file);
				exit(0);
			}
			for(auto phi = Solvers::get_instance()->parameters.Phases.begin(); phi < Solvers::get_instance()->parameters.Phases.end(); phi++){
				conductivity_phi.add_double(phi->phi_property, conductivity_phi_value[index].double_value);
				index++;
			}

			InputFileReader::get_instance()->read_double_value("Modules.ElectricField.BackGround.conductivity", conductivity_background, infile_debug);

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Modules.ElectricField.fix_phi = [(phi_name, elec_potential), ... ] \n", InputFileReader::get_instance()->debug_file);
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
			if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_DIFFERENCE) {
				int iterate_steps = electric_field_solver_diff.solve_whole_domain(solver_accuracy, solver_max_iterate_times, solver_debug, solver_debug_output_steps);
				output << "> Solver :" << electric_field_solver_diff.solver_name << ", iterate " << iterate_steps << " times." << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
			else if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_FOURIER_SPECTRAL || electric_field_solver == ElectricFieldSolver::EFS_IMPLICIT_FOURIER_SPECTRAL) {
				electric_field_solver_fourier.init_basic_real_space();
				im_max_varient = 0.0;
				int iterate_steps = 1;
				for (iterate_steps = 1; iterate_steps <= solver_max_iterate_times; iterate_steps++) {
					im_varient = 0.0;
					electric_field_solver_fourier.solve_one_step();
					if (solver_debug && (iterate_steps % solver_debug_output_steps == 0)) {
						output.str("");
						output << im_solver_name << " iterate " << iterate_steps << " times !" << endl;
						output << im_solver_name << " variation = " << im_varient << " !" << endl;
						Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
					}
					if (im_max_varient < im_varient)
						im_max_varient = im_varient;
					if (im_varient < solver_accuracy)
						break;
				}
				if (solver_debug) {
					output.str("");
					output << im_solver_name << " iterate " << iterate_steps << " times !" << endl;
					output << im_solver_name << " max_variation = " << im_max_varient << " !" << endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
				}
				output.str("");
				output << "> Solver :" << im_solver_name << ", iterate " << iterate_steps << " times," << " max_variation = " << im_max_varient << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream report;
			if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_DIFFERENCE) {
				int iterate_steps = electric_field_solver_diff.solve_whole_domain(solver_accuracy, solver_max_iterate_times, false, solver_debug_output_steps);
				report << "> Solver :" << electric_field_solver_diff.solver_name << ", iterate " << iterate_steps << " times." << endl;
			}
			else if (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_FOURIER_SPECTRAL || electric_field_solver == ElectricFieldSolver::EFS_IMPLICIT_FOURIER_SPECTRAL) {
				electric_field_solver_fourier.init_basic_real_space();
				im_max_varient = 0.0;
				int iterate_steps = 1;
				for (iterate_steps = 1; iterate_steps <= solver_max_iterate_times; iterate_steps++) {
					im_varient = 0.0;
					electric_field_solver_fourier.solve_one_step();
					if (im_max_varient < im_varient)
						im_max_varient = im_varient;
					if (im_varient < solver_accuracy)
						break;
				}
				report << "> Solver :" << im_solver_name << ", iterate " << iterate_steps << " times," << " max_variation = " << im_max_varient << endl;
			}
			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			if (electric_field_solver == ElectricFieldSolver::EFS_NONE)
				return;
			fout << "<DataArray type = \"Float64\" Name = \"" << "elec_conductivity" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
			for (int z = 0; z < phaseMesh.limit_z; z++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int x = 0; x < phaseMesh.limit_x; x++) {
						PhaseNode& node = phaseMesh(x, y, z);
						fout << node.customValues[ElectricalConductivity] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "elec_potential" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
			for (int z = 0; z < phaseMesh.limit_z; z++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int x = 0; x < phaseMesh.limit_x; x++) {
						PhaseNode& node = phaseMesh(x, y, z);
						fout << node.customValues[ElectricalPotential] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "charge_density" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
			for (int z = 0; z < phaseMesh.limit_z; z++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int x = 0; x < phaseMesh.limit_x; x++) {
						PhaseNode& node = phaseMesh(x, y, z);
						fout << node.customValues[ChargeDensity] << endl;
					}
			fout << "</DataArray>" << endl;
			// debug Q2 Q4
			if (solver_debug && (electric_field_solver == ElectricFieldSolver::EFS_EXPLICITE_FOURIER_SPECTRAL || electric_field_solver == ElectricFieldSolver::EFS_IMPLICIT_FOURIER_SPECTRAL) && false) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "Q2" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
				for (int z = 0; z < phaseMesh.limit_z; z++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int x = 0; x < phaseMesh.limit_x; x++) {
							fout << electric_field_solver_fourier.Q2(x, y, z) << endl;
						}
				fout << "</DataArray>" << endl;
				fout << "<DataArray type = \"Float64\" Name = \"" << "Q4" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
#pragma omp parallel for
				for (int z = 0; z < phaseMesh.limit_z; z++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int x = 0; x < phaseMesh.limit_x; x++) {
							fout << electric_field_solver_fourier.Q4(x, y, z) << endl;
						}
				fout << "</DataArray>" << endl;
			}
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}

		static void init_module() {
			Solvers::get_instance()->create_a_new_module(init, exec_pre, exec_loop, deinit, write_scalar, write_vec3);
		};
	}
}