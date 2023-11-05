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
	namespace pressure_correction {
		// FST_Traditional_Difference
		FluidField fluid_solver;
		static double viscosity = 1.0;
		static double mass = 1.0;
		static double fluid_dt = 1.0;
		static bool debug_solver = false;
		static int debug_output_step = 1000;
		static int max_solver_iterate_steps = 0;
		static int max_pressure_correction_steps = 1000;
		static double pressure_accuracy = 1e-4;
		static int correct_steps = 1;
		// boundary condition
		static vector<int> solid_phase_properties;
		const double fluid_solid_interface_threshold = 0.5;
		static vector<Vector3> flow_boundary_condition; // domain_boundary_values[Direction]
		static vector<double> pressure_boundary_condition; // domain_boundary_values[Direction]
		static void (*domain_boundary_condition_u_x_down)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_u_y_down)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_u_z_down)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_u_x_up)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_u_y_up)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_u_z_up)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_v_x_down)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_v_y_down)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_v_z_down)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_v_x_up)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_v_y_up)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_v_z_up)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_w_x_down)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_w_y_down)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_w_z_down)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_w_x_up)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_w_y_up)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_w_z_up)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_pressure_x_down)(pf::PhaseNode& node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_pressure_y_down)(pf::PhaseNode& node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_pressure_z_down)(pf::PhaseNode& node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_pressure_x_up)(pf::PhaseNode& node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_pressure_y_up)(pf::PhaseNode& node, int Nx, int Ny, int Nz);
		static void (*domain_boundary_condition_pressure_z_up)(pf::PhaseNode& node, int Nx, int Ny, int Nz);
		static void (*smooth_boundary_condition_u)(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*smooth_boundary_condition_v)(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*smooth_boundary_condition_w)(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz);
		static void (*smooth_boundary_condition_pressure)(pf::PhaseNode& node, int Nx, int Ny, int Nz);

		namespace inse_boundary_condition {
			static void default_boundary_condition_node(pf::PhaseNode& node, int Nx, int Ny, int Nz) {};
			static void default_boundary_condition_stagered(VectorNode& vn, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {};
			// speed x_down
			static void fix_u_x_down(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._x == 0) {
					u.vals[0] = flow_boundary_condition[Boundary::DOWN_X][0];
				}
			}
			static void fix_v_x_down(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._x == 0) {
					v.vals[0] = flow_boundary_condition[Boundary::DOWN_X][1];
				}
			}
			static void fix_w_x_down(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._x == 0) {
					w.vals[0] = flow_boundary_condition[Boundary::DOWN_X][2];
				}
			}
			// speed y_down
			static void fix_u_y_down(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._y == 0) {
					u.vals[0] = flow_boundary_condition[Boundary::DOWN_Y][0];
				}
			}
			static void fix_v_y_down(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._y == 0) {
					v.vals[0] = flow_boundary_condition[Boundary::DOWN_Y][1];
				}
			}
			static void fix_w_y_down(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._y == 0) {
					w.vals[0] = flow_boundary_condition[Boundary::DOWN_Y][2];
				}
			}
			// speed z_down
			static void fix_u_z_down(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._z == 0) {
					u.vals[0] = flow_boundary_condition[Boundary::DOWN_Z][0];
				}
			}
			static void fix_v_z_down(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._z == 0) {
					v.vals[0] = flow_boundary_condition[Boundary::DOWN_Z][1];
				}
			}
			static void fix_w_z_down(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._z == 0) {
					w.vals[0] = flow_boundary_condition[Boundary::DOWN_Z][2];
				}
			}
			// speed x_up
			static void fix_u_x_up(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._x == Nx - 1) {
					u.vals[0] = flow_boundary_condition[Boundary::UP_X][0];
				}
			}
			static void fix_v_x_up(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._x == Nx - 1) {
					v.vals[0] = flow_boundary_condition[Boundary::UP_X][1];
				}
			}
			static void fix_w_x_up(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._x == Nx - 1) {
					w.vals[0] = flow_boundary_condition[Boundary::UP_X][2];
				}
			}
			// speed y_up
			static void fix_u_y_up(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._y == Ny - 1) {
					u.vals[0] = flow_boundary_condition[Boundary::UP_Y][0];
				}
			}
			static void fix_v_y_up(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._y == Ny - 1) {
					v.vals[0] = flow_boundary_condition[Boundary::UP_Y][1];
				}
			}
			static void fix_w_y_up(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._y == Ny - 1) {
					w.vals[0] = flow_boundary_condition[Boundary::UP_Y][2];
				}
			}
			// speed z_up
			static void fix_u_z_up(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._z == Nz - 1) {
					u.vals[0] = flow_boundary_condition[Boundary::UP_Z][0];
				}
			}
			static void fix_v_z_up(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._z == Nz - 1) {
					v.vals[0] = flow_boundary_condition[Boundary::UP_Z][1];
				}
			}
			static void fix_w_z_up(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._z == Nz - 1) {
					w.vals[0] = flow_boundary_condition[Boundary::UP_Z][2];
				}
			}
			// pressure
			static void fix_pressure_x_down(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node._x == 0) {
					node.customValues[ExternalFields::FLUID_pressure] = pressure_boundary_condition[Boundary::DOWN_X];
				}
			}
			static void fix_pressure_x_up(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node._x == Nx - 1) {
					node.customValues[ExternalFields::FLUID_pressure] = pressure_boundary_condition[Boundary::UP_X];
				}
			}
			static void fix_pressure_y_down(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node._y == 0) {
					node.customValues[ExternalFields::FLUID_pressure] = pressure_boundary_condition[Boundary::DOWN_Y];
				}
			}
			static void fix_pressure_y_up(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node._y == Ny - 1) {
					node.customValues[ExternalFields::FLUID_pressure] = pressure_boundary_condition[Boundary::UP_Y];
				}
			}
			static void fix_pressure_z_down(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node._z == 0) {
					node.customValues[ExternalFields::FLUID_pressure] = pressure_boundary_condition[Boundary::DOWN_Z];
				}
			}
			static void fix_pressure_z_up(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node._z == Nz - 1) {
					node.customValues[ExternalFields::FLUID_pressure] = pressure_boundary_condition[Boundary::UP_Z];
				}
			}
			// free boundary
			// x_down
			static void free_u_x_down(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._x == 0) {
					u.vals[0] = u.get_neighbor_node(Direction::x_up).vals[0];
				}
			}
			static void free_v_x_down(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._x == 0) {
					v.vals[0] = v.get_neighbor_node(Direction::x_up).vals[0];
				}
			}
			static void free_w_x_down(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._x == 0) {
					w.vals[0] = w.get_neighbor_node(Direction::x_up).vals[0];
				}
			}
			// x_up
			static void free_u_x_up(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._x == Nx - 1) {
					u.vals[0] = u.get_neighbor_node(Direction::x_down).vals[0];
				}
			}
			static void free_v_x_up(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._x == Nx - 1) {
					v.vals[0] = v.get_neighbor_node(Direction::x_down).vals[0];
				}
			}
			static void free_w_x_up(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._x == Nx - 1) {
					w.vals[0] = w.get_neighbor_node(Direction::x_down).vals[0];
				}
			}
			// y_down
			static void free_u_y_down(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._y == 0) {
					u.vals[0] = u.get_neighbor_node(Direction::y_up).vals[0];
				}
			}
			static void free_v_y_down(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._y == 0) {
					v.vals[0] = v.get_neighbor_node(Direction::y_up).vals[0];
				}
			}
			static void free_w_y_down(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._y == 0) {
					w.vals[0] = w.get_neighbor_node(Direction::y_up).vals[0];
				}
			}
			// y_up
			static void free_u_y_up(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._y == Ny - 1) {
					u.vals[0] = u.get_neighbor_node(Direction::y_down).vals[0];
				}
			}
			static void free_v_y_up(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._y == Ny - 1) {
					v.vals[0] = v.get_neighbor_node(Direction::y_down).vals[0];
				}
			}
			static void free_w_y_up(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._y == Ny - 1) {
					w.vals[0] = w.get_neighbor_node(Direction::y_down).vals[0];
				}
			}
			// z_down
			static void free_u_z_down(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._z == 0) {
					u.vals[0] = u.get_neighbor_node(Direction::z_up).vals[0];
				}
			}
			static void free_v_z_down(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._z == 0) {
					v.vals[0] = v.get_neighbor_node(Direction::z_up).vals[0];
				}
			}
			static void free_w_z_down(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._z == 0) {
					w.vals[0] = w.get_neighbor_node(Direction::z_up).vals[0];
				}
			}
			// z_up
			static void free_u_z_up(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (u._z == Nz - 1) {
					u.vals[0] = u.get_neighbor_node(Direction::z_down).vals[0];
				}
			}
			static void free_v_z_up(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (v._z == Nz - 1) {
					v.vals[0] = v.get_neighbor_node(Direction::z_down).vals[0];
				}
			}
			static void free_w_z_up(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (w._z == Nz - 1) {
					w.vals[0] = w.get_neighbor_node(Direction::z_down).vals[0];
				}
			}
			// smooth condition
			static void smooth_u(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (down_node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold &&
					up_node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold) {
					u.vals[0] = 0.0;
				}
			}
			static void smooth_v(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (down_node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold &&
					up_node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold) {
					v.vals[0] = 0.0;
				}
			}
			static void smooth_w(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
				if (down_node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold &&
					up_node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold) {
					w.vals[0] = 0.0;
				}
			}
			static void smooth_pressure(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
				if (node.customValues[ExternalFields::FLUID_Fluid_Domain] < fluid_solid_interface_threshold) {
					node.customValues[ExternalFields::FLUID_pressure] = 0.0;
				}
			}
		}

		// smooth boundary
		static void boundary_condition_for_U(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			domain_boundary_condition_u_x_down(u, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_u_y_down(u, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_u_z_down(u, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_u_x_up(u, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_u_y_up(u, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_u_z_up(u, down_node, up_node, Nx, Ny, Nz);
			smooth_boundary_condition_u(u, down_node, up_node, Nx, Ny, Nz);
		}
		static void boundary_condition_for_V(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			domain_boundary_condition_v_x_down(v, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_v_y_down(v, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_v_z_down(v, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_v_x_up(v, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_v_y_up(v, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_v_z_up(v, down_node, up_node, Nx, Ny, Nz);
			smooth_boundary_condition_v(v, down_node, up_node, Nx, Ny, Nz);
		}
		static void boundary_condition_for_W(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			domain_boundary_condition_w_x_down(w, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_w_y_down(w, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_w_z_down(w, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_w_x_up(w, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_w_y_up(w, down_node, up_node, Nx, Ny, Nz);
			domain_boundary_condition_w_z_up(w, down_node, up_node, Nx, Ny, Nz);
			smooth_boundary_condition_w(w, down_node, up_node, Nx, Ny, Nz);
		}
		static void boundary_condition_for_pressure(PhaseNode& node, int Nx, int Ny, int Nz) {
			domain_boundary_condition_pressure_x_down(node, Nx, Ny, Nz);
			domain_boundary_condition_pressure_y_down(node, Nx, Ny, Nz);
			domain_boundary_condition_pressure_z_down(node, Nx, Ny, Nz);
			domain_boundary_condition_pressure_x_up(node, Nx, Ny, Nz);
			domain_boundary_condition_pressure_y_up(node, Nx, Ny, Nz);
			domain_boundary_condition_pressure_z_up(node, Nx, Ny, Nz);
			smooth_boundary_condition_pressure(node, Nx, Ny, Nz);
		}
		static void cal_parameters_for_main_domain(PhaseNode& node, int Nx, int Ny, int Nz) {
			double solid_phases_volume_fraction = 0.0;
			for (auto phi = node.begin(); phi < node.end(); phi++)
				for (auto phase = solid_phase_properties.begin(); phase < solid_phase_properties.end(); phase++)
					if (phi->property == *phase) {
						solid_phases_volume_fraction += phi->phi;
					}
			if (solid_phases_volume_fraction > 1.0)
				solid_phases_volume_fraction = 1.0;
			else if (solid_phases_volume_fraction < 0.0)
				solid_phases_volume_fraction = 0.0;
			node.customValues[ExternalFields::FLUID_Fluid_Domain] = 1.0 - solid_phases_volume_fraction;
			node.customValues[ExternalFields::FLUID_viscosity] = viscosity;
			node.customValues[ExternalFields::FLUID_mass] = mass;
			node.customVec3s[ExternalFields::FLUID_volume_force].set_to_zero();
		}
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			fluid_solver.init(phaseMesh, phaseMesh._bc_x_up, phaseMesh._bc_x_down, phaseMesh._bc_y_up, phaseMesh._bc_y_down, phaseMesh._bc_z_up, phaseMesh._bc_z_down);
			fluid_solver.boundary_condition_for_domain_U = boundary_condition_for_U;
			fluid_solver.boundary_condition_for_domain_V = boundary_condition_for_V;
			fluid_solver.boundary_condition_for_domain_W = boundary_condition_for_W;
			fluid_solver.boundary_condition_for_pressure = boundary_condition_for_pressure;
			fluid_solver.cal_parameters_for_main_domain = cal_parameters_for_main_domain;
			flow_boundary_condition.resize(6);
			pressure_boundary_condition.resize(6);

			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.viscosity", viscosity, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.density", mass, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.fluid_dt", fluid_dt, infile_debug);
			if (InputFileReader::get_instance()->read_bool_value("Postprocess.FluidDynamics.PressureCorrection.debug_solver", debug_solver, infile_debug))
				InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.PressureCorrection.debug_output_step", debug_output_step, infile_debug);
			InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.PressureCorrection.Pressure.max_iterate_steps", max_pressure_correction_steps, infile_debug);
			InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.Pressure.accuracy", pressure_accuracy, infile_debug);
			InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.PressureCorrection.Velocity.max_iterate_steps", max_solver_iterate_steps, infile_debug);
			InputFileReader::get_instance()->read_int_value("Postprocess.FluidDynamics.PressureCorrection.Velocity.correct_steps", correct_steps, infile_debug);
			smooth_boundary_condition_pressure = inse_boundary_condition::default_boundary_condition_node;
			smooth_boundary_condition_u = inse_boundary_condition::default_boundary_condition_stagered;
			smooth_boundary_condition_v = inse_boundary_condition::default_boundary_condition_stagered;
			smooth_boundary_condition_w = inse_boundary_condition::default_boundary_condition_stagered;
			if (infile_debug)
			InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.FluidDynamics.PressureCorrection.solid_phases = (phase_name, ... ) \n", InputFileReader::get_instance()->debug_file);
			string solid_phase_key = "Postprocess.FluidDynamics.PressureCorrection.solid_phases", solid_phase_input = "()";
			if (InputFileReader::get_instance()->read_string_value(solid_phase_key, solid_phase_input, infile_debug)) {
				vector<input_value> solid_phase_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, solid_phase_key, solid_phase_input, infile_debug);
				for (auto solid_name = solid_phase_value.begin(); solid_name < solid_phase_value.end(); solid_name++)
					solid_phase_properties.push_back(Solvers::get_instance()->parameters.Phases[solid_name->string_value].phi_property);
				smooth_boundary_condition_pressure = inse_boundary_condition::smooth_pressure;
				smooth_boundary_condition_u = inse_boundary_condition::smooth_u;
				smooth_boundary_condition_v = inse_boundary_condition::smooth_v;
				smooth_boundary_condition_w = inse_boundary_condition::smooth_w;
			}

			domain_boundary_condition_pressure_x_down = inse_boundary_condition::default_boundary_condition_node;
			domain_boundary_condition_u_x_down = inse_boundary_condition::free_u_x_down;
			domain_boundary_condition_v_x_down = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_w_x_down = inse_boundary_condition::default_boundary_condition_stagered;
			if (phaseMesh._bc_x_down != BoundaryCondition::PERIODIC) {
				string domain_boundary_key = "Postprocess.FluidDynamics.PressureCorrection.BC_Down_X.fix_velocity", domain_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(domain_boundary_key, domain_boundary_input, infile_debug)) {
					vector<input_value> domain_boundary_val = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, domain_boundary_key, domain_boundary_input, infile_debug);
					flow_boundary_condition[Boundary::DOWN_X] = Vector3(domain_boundary_val[0].double_value, domain_boundary_val[1].double_value, domain_boundary_val[2].double_value);
					domain_boundary_condition_u_x_down = inse_boundary_condition::fix_u_x_down;
					domain_boundary_condition_v_x_down = inse_boundary_condition::fix_v_x_down;
					domain_boundary_condition_w_x_down = inse_boundary_condition::fix_w_x_down;
				}
				if (InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.BC_Down_X.pressure", pressure_boundary_condition[Boundary::DOWN_X], infile_debug)) {
					domain_boundary_condition_pressure_x_down = inse_boundary_condition::fix_pressure_x_down;
				}
			}
			domain_boundary_condition_pressure_y_down = inse_boundary_condition::default_boundary_condition_node;
			domain_boundary_condition_u_y_down = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_v_y_down = inse_boundary_condition::free_v_y_down;
			domain_boundary_condition_w_y_down = inse_boundary_condition::default_boundary_condition_stagered;
			if (phaseMesh._bc_y_down != BoundaryCondition::PERIODIC) {
				string domain_boundary_key = "Postprocess.FluidDynamics.PressureCorrection.BC_Down_Y.fix_velocity", domain_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(domain_boundary_key, domain_boundary_input, infile_debug)) {
					vector<input_value> domain_boundary_val = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, domain_boundary_key, domain_boundary_input, infile_debug);
					flow_boundary_condition[Boundary::DOWN_Y] = Vector3(domain_boundary_val[0].double_value, domain_boundary_val[1].double_value, domain_boundary_val[2].double_value);
					domain_boundary_condition_u_y_down = inse_boundary_condition::fix_u_y_down;
					domain_boundary_condition_v_y_down = inse_boundary_condition::fix_v_y_down;
					domain_boundary_condition_w_y_down = inse_boundary_condition::fix_w_y_down;
				}
				if (InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.BC_Down_Y.pressure", pressure_boundary_condition[Boundary::DOWN_Y], infile_debug)) {
					domain_boundary_condition_pressure_y_down = inse_boundary_condition::fix_pressure_y_down;
				}
			}
			domain_boundary_condition_pressure_z_down = inse_boundary_condition::default_boundary_condition_node;
			domain_boundary_condition_u_z_down = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_v_z_down = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_w_z_down = inse_boundary_condition::free_w_z_down;
			if (phaseMesh._bc_z_down != BoundaryCondition::PERIODIC) {
				string domain_boundary_key = "Postprocess.FluidDynamics.PressureCorrection.BC_Down_Z.fix_velocity", domain_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(domain_boundary_key, domain_boundary_input, infile_debug)) {
					vector<input_value> domain_boundary_val = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, domain_boundary_key, domain_boundary_input, infile_debug);
					flow_boundary_condition[Boundary::DOWN_Z] = Vector3(domain_boundary_val[0].double_value, domain_boundary_val[1].double_value, domain_boundary_val[2].double_value);
					domain_boundary_condition_u_z_down = inse_boundary_condition::fix_u_z_down;
					domain_boundary_condition_v_z_down = inse_boundary_condition::fix_v_z_down;
					domain_boundary_condition_w_z_down = inse_boundary_condition::fix_w_z_down;
				}
				if (InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.BC_Down_Z.pressure", pressure_boundary_condition[Boundary::DOWN_Z], infile_debug)) {
					domain_boundary_condition_pressure_z_down = inse_boundary_condition::fix_pressure_z_down;
				}
			}
			domain_boundary_condition_pressure_x_up = inse_boundary_condition::default_boundary_condition_node;
			domain_boundary_condition_u_x_up = inse_boundary_condition::free_u_x_up;
			domain_boundary_condition_v_x_up = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_w_x_up = inse_boundary_condition::default_boundary_condition_stagered;
			if (phaseMesh._bc_x_up != BoundaryCondition::PERIODIC) {
				string domain_boundary_key = "Postprocess.FluidDynamics.PressureCorrection.BC_Up_X.fix_velocity", domain_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(domain_boundary_key, domain_boundary_input, infile_debug)) {
					vector<input_value> domain_boundary_val = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, domain_boundary_key, domain_boundary_input, infile_debug);
					flow_boundary_condition[Boundary::UP_X] = Vector3(domain_boundary_val[0].double_value, domain_boundary_val[1].double_value, domain_boundary_val[2].double_value);
					domain_boundary_condition_u_x_up = inse_boundary_condition::fix_u_x_up;
					domain_boundary_condition_v_x_up = inse_boundary_condition::fix_v_x_up;
					domain_boundary_condition_w_x_up = inse_boundary_condition::fix_w_x_up;
				}
				if (InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.BC_Up_X.pressure", pressure_boundary_condition[Boundary::UP_X], infile_debug)) {
					domain_boundary_condition_pressure_x_up = inse_boundary_condition::fix_pressure_x_up;
				}
			}
			domain_boundary_condition_pressure_y_up = inse_boundary_condition::default_boundary_condition_node;
			domain_boundary_condition_u_y_up = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_v_y_up = inse_boundary_condition::free_v_y_up;
			domain_boundary_condition_w_y_up = inse_boundary_condition::default_boundary_condition_stagered;
			if (phaseMesh._bc_y_up != BoundaryCondition::PERIODIC) {
				string domain_boundary_key = "Postprocess.FluidDynamics.PressureCorrection.BC_Up_Y.fix_velocity", domain_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(domain_boundary_key, domain_boundary_input, infile_debug)) {
					vector<input_value> domain_boundary_val = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, domain_boundary_key, domain_boundary_input, infile_debug);
					flow_boundary_condition[Boundary::UP_Y] = Vector3(domain_boundary_val[0].double_value, domain_boundary_val[1].double_value, domain_boundary_val[2].double_value);
					domain_boundary_condition_u_y_up = inse_boundary_condition::fix_u_y_up;
					domain_boundary_condition_v_y_up = inse_boundary_condition::fix_v_y_up;
					domain_boundary_condition_w_y_up = inse_boundary_condition::fix_w_y_up;
				}
				if (InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.BC_Up_Y.pressure", pressure_boundary_condition[Boundary::UP_Y], infile_debug)) {
					domain_boundary_condition_pressure_y_up = inse_boundary_condition::fix_pressure_y_up;
				}
			}
			domain_boundary_condition_pressure_z_up = inse_boundary_condition::default_boundary_condition_node;
			domain_boundary_condition_u_z_up = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_v_z_up = inse_boundary_condition::default_boundary_condition_stagered;
			domain_boundary_condition_w_z_up = inse_boundary_condition::free_w_z_up;
			if (phaseMesh._bc_z_up != BoundaryCondition::PERIODIC) {
				string domain_boundary_key = "Postprocess.FluidDynamics.PressureCorrection.BC_Up_Z.fix_velocity", domain_boundary_input = "(0,0,0)";
				if (InputFileReader::get_instance()->read_string_value(domain_boundary_key, domain_boundary_input, infile_debug)) {
					vector<input_value> domain_boundary_val = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, domain_boundary_key, domain_boundary_input, infile_debug);
					flow_boundary_condition[Boundary::UP_Z] = Vector3(domain_boundary_val[0].double_value, domain_boundary_val[1].double_value, domain_boundary_val[2].double_value);
					domain_boundary_condition_u_z_up = inse_boundary_condition::fix_u_z_up;
					domain_boundary_condition_v_z_up = inse_boundary_condition::fix_v_z_up;
					domain_boundary_condition_w_z_up = inse_boundary_condition::fix_w_z_up;
				}
				if (InputFileReader::get_instance()->read_double_value("Postprocess.FluidDynamics.PressureCorrection.BC_Up_Z.pressure", pressure_boundary_condition[Boundary::UP_Z], infile_debug)) {
					domain_boundary_condition_pressure_z_up = inse_boundary_condition::fix_pressure_z_up;
				}
			}
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream output;
			fluid_solver.assign_velocity_to_fluid_domain();
			if (debug_solver) {
				output << "> Fluid field solver debug:" << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			}
			vector<double> MAX_ABS_PRESSURE;
			MAX_ABS_PRESSURE.resize(2);
			Vector3 MAX_ABS_dVelocity;
			fluid_solver.cal_parameters_before_calculation();
			for (int istep = 1; istep <= max_solver_iterate_steps; istep++) {
				MAX_ABS_PRESSURE.clear();

				MAX_ABS_dVelocity = fluid_solver.evolve_momentum_equation(fluid_dt);

				MAX_ABS_PRESSURE = fluid_solver.do_pressure_correction(fluid_solid_interface_threshold, pressure_accuracy, max_pressure_correction_steps, false, debug_output_step);

				fluid_solver.correcting_velocity_field(fluid_dt);

				if (debug_solver && (istep % debug_output_step == 0 || istep == 1)) {
					output.str("");
					output << "> Fluid field iterate step:" << istep << endl
						<< "		ABS_dPressure      = " << MAX_ABS_PRESSURE[0] << endl
						<< "		Pressure_Corrects  = " << int(MAX_ABS_PRESSURE[1] + 0.5) << endl
						<< "		MAX_ABS_dU		   = " << MAX_ABS_dVelocity[0] << endl
						<< "		MAX_ABS_dV		   = " << MAX_ABS_dVelocity[1] << endl
						<< "		MAX_ABS_dW		   = " << MAX_ABS_dVelocity[2] << endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
				}
				if (MAX_ABS_PRESSURE[0] < pressure_accuracy && int(MAX_ABS_PRESSURE[1] + 0.5) <= correct_steps)
					break;
			}
			output.str("");
			output << "> Fluid field solver:" << endl
				<< "		ABS_dPressure      = " << MAX_ABS_PRESSURE[0] << endl
				<< "		Pressure_Corrects  = " << int(MAX_ABS_PRESSURE[1] + 0.5) << endl
				<< "		MAX_ABS_dU		   = " << MAX_ABS_dVelocity[0] << endl
				<< "		MAX_ABS_dV		   = " << MAX_ABS_dVelocity[1] << endl
				<< "		MAX_ABS_dW		   = " << MAX_ABS_dVelocity[2] << endl;
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
			fluid_solver.assign_velocity_to_main_domain();
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream report;
			vector<double> MAX_ABS_PRESSURE;
			MAX_ABS_PRESSURE.resize(2);
			Vector3 MAX_ABS_dVelocity;
			fluid_solver.cal_parameters_before_calculation();
			for (int istep = 1; istep <= max_solver_iterate_steps; istep++) {
				MAX_ABS_PRESSURE.clear();

				MAX_ABS_dVelocity = fluid_solver.evolve_momentum_equation(fluid_dt);

				MAX_ABS_PRESSURE = fluid_solver.do_pressure_correction(fluid_solid_interface_threshold, pressure_accuracy, max_pressure_correction_steps, false, debug_output_step);

				fluid_solver.correcting_velocity_field(fluid_dt);

				if (MAX_ABS_PRESSURE[0] < pressure_accuracy && int(MAX_ABS_PRESSURE[1] + 0.5) <= correct_steps)
					break;
			}
			report << "> Fluid field solver:" << endl
				<< "		ABS_dPressure      = " << MAX_ABS_PRESSURE[0] << endl
				<< "		Pressure_Corrects  = " << int(MAX_ABS_PRESSURE[1] + 0.5) << endl
				<< "		MAX_ABS_dU		   = " << MAX_ABS_dVelocity[0] << endl
				<< "		MAX_ABS_dV		   = " << MAX_ABS_dVelocity[1] << endl
				<< "		MAX_ABS_dW		   = " << MAX_ABS_dVelocity[2] << endl;
			fluid_solver.assign_velocity_to_main_domain();
			return report.str();
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			fluid_solver.free();
			solid_phase_properties.clear();
			flow_boundary_condition.clear();
			pressure_boundary_condition.clear();
			domain_boundary_condition_pressure_x_down = nullptr;
			domain_boundary_condition_pressure_x_up = nullptr;
			domain_boundary_condition_pressure_y_down = nullptr;
			domain_boundary_condition_pressure_y_up = nullptr;
			domain_boundary_condition_pressure_z_down = nullptr;
			domain_boundary_condition_pressure_z_up = nullptr;
			domain_boundary_condition_u_x_down = nullptr;
			domain_boundary_condition_u_x_up = nullptr;
			domain_boundary_condition_u_y_down = nullptr;
			domain_boundary_condition_u_y_up = nullptr;
			domain_boundary_condition_u_z_down = nullptr;
			domain_boundary_condition_u_z_up = nullptr;
			domain_boundary_condition_v_x_down = nullptr;
			domain_boundary_condition_v_x_up = nullptr;
			domain_boundary_condition_v_y_down = nullptr;
			domain_boundary_condition_v_y_up = nullptr;
			domain_boundary_condition_v_z_down = nullptr;
			domain_boundary_condition_v_z_up = nullptr;
			domain_boundary_condition_w_x_down = nullptr;
			domain_boundary_condition_w_x_up = nullptr;
			domain_boundary_condition_w_y_down = nullptr;
			domain_boundary_condition_w_y_up = nullptr;
			domain_boundary_condition_w_z_down = nullptr;
			domain_boundary_condition_w_z_up = nullptr;
			smooth_boundary_condition_pressure = nullptr;
			smooth_boundary_condition_u = nullptr;
			smooth_boundary_condition_v = nullptr;
			smooth_boundary_condition_w = nullptr;
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_pressure" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customValues[ExternalFields::FLUID_pressure] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_density" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customValues[ExternalFields::FLUID_mass] << endl;
					}
			fout << "</DataArray>" << endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_speed" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity].abs() << endl;
					}
			fout << "</DataArray>" << endl;
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			string name;
			name = "\"fluid_velocity\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity][0] << " "
							<< phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity][1] << " "
							<< phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_velocity][2] << endl;
					}
			fout << "</DataArray>" << endl;
			name = "\"fluid_volume_force\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < phaseMesh.limit_z; k++)
				for (int j = 0; j < phaseMesh.limit_y; j++)
					for (int i = 0; i < phaseMesh.limit_x; i++) {
						fout << phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_volume_force][0] << " "
							<< phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_volume_force][1] << " "
							<< phaseMesh(i, j, k).customVec3s[ExternalFields::FLUID_volume_force][2] << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
}