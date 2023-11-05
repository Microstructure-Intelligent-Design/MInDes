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
#include"base.h"

namespace pf {
	namespace fluid_boundary_condition_funcs {
		static void boundary_condition_for_domain_U(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			return;
		}
		static void boundary_condition_for_domain_V(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			return;
		}
		static void boundary_condition_for_domain_W(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			return;
		}
		static void boundary_condition_for_pressure(PhaseNode& node, int Nx, int Ny, int Nz) {
			return;
		}
		static void cal_parameters_for_main_domain(PhaseNode& node, int Nx, int Ny, int Nz) {
			node.customValues[ExternalFields::FLUID_viscosity] = 1.0;
			node.customValues[ExternalFields::FLUID_mass] = 1.0;
			node.customVec3s[ExternalFields::FLUID_volume_force].set_to_zero();
			return;
		}
	};
	class FluidField
	{
	public:
		FluidField(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition x_up_bc,
			BoundaryCondition x_down_bc, BoundaryCondition y_up_bc, BoundaryCondition y_down_bc, BoundaryCondition z_up_bc, BoundaryCondition z_down_bc) {
			init(_phaseMesh, x_up_bc, x_down_bc, y_up_bc, y_down_bc, z_up_bc, z_down_bc);
		};
		FluidField() {};
		~FluidField() { 
			free();
		};
		void init(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition x_up_bc,
			BoundaryCondition x_down_bc, BoundaryCondition y_up_bc, BoundaryCondition y_down_bc, BoundaryCondition z_up_bc, BoundaryCondition z_down_bc) {
			phaseMesh = &_phaseMesh;
			if (x_up_bc == PERIODIC && x_down_bc == PERIODIC) {
				U.init(phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z,
					phaseMesh->dr, x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc);
			}
			else {
				U.init(phaseMesh->limit_x + 1, phaseMesh->limit_y, phaseMesh->limit_z,
					phaseMesh->dr, x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc);
			}
			if (y_up_bc == PERIODIC && y_down_bc == PERIODIC) {
				V.init(phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z,
					phaseMesh->dr, x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc);
			}
			else {
				V.init(phaseMesh->limit_x, phaseMesh->limit_y + 1, phaseMesh->limit_z,
					phaseMesh->dr, x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc);
			}
			if (z_up_bc == PERIODIC && z_down_bc == PERIODIC) {
				W.init(phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z,
					phaseMesh->dr, x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc);
			}
			else {
				W.init(phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z + 1,
					phaseMesh->dr, x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc);
			}
			for (auto vecNode = U._mesh.begin(); vecNode < U._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value: 0
				vecNode->vals.push_back(0.0); // increment: 1
			}
			for (auto vecNode = V._mesh.begin(); vecNode < V._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value
				vecNode->vals.push_back(0.0); // increment
			}
			for (auto vecNode = W._mesh.begin(); vecNode < W._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value
				vecNode->vals.push_back(0.0); // increment
			}
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				node->customValues.add_double(ExternalFields::FLUID_Fluid_Domain, 0.0);
				node->customValues.add_double(ExternalFields::FLUID_mass, 0.0);
				node->customValues.add_double(ExternalFields::FLUID_pressure, 0.0);
				node->customValues.add_double(ExternalFields::FLUID_viscosity, 0.0);
				node->customVec3s.add_vec(ExternalFields::FLUID_volume_force, Vector3(0.0, 0.0, 0.0));
				node->customVec3s.add_vec(ExternalFields::FLUID_velocity, Vector3(0.0, 0.0, 0.0));
			}
			boundary_condition_for_domain_U = fluid_boundary_condition_funcs::boundary_condition_for_domain_U;
			boundary_condition_for_domain_V = fluid_boundary_condition_funcs::boundary_condition_for_domain_V;
			boundary_condition_for_domain_W = fluid_boundary_condition_funcs::boundary_condition_for_domain_W;
			boundary_condition_for_pressure = fluid_boundary_condition_funcs::boundary_condition_for_pressure;
			cal_parameters_for_main_domain = fluid_boundary_condition_funcs::cal_parameters_for_main_domain;
		}

		void set_u_in_velocity_field(int x, int y, int z, double u);
		void set_v_in_velocity_field(int x, int y, int z, double v);
		void set_w_in_velocity_field(int x, int y, int z, double w);

		void define_funcs_for_fluid(void(*boundary_U)(pf::VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int) = fluid_boundary_condition_funcs::boundary_condition_for_domain_U,
			void(*boundary_V)(pf::VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int) = fluid_boundary_condition_funcs::boundary_condition_for_domain_V,
			void(*boundary_W)(pf::VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int) = fluid_boundary_condition_funcs::boundary_condition_for_domain_W,
			void(*boundary_main_node)(pf::PhaseNode&, int, int, int) = fluid_boundary_condition_funcs::boundary_condition_for_pressure,
			void(*cal_parameters)(pf::PhaseNode&, int, int, int) = fluid_boundary_condition_funcs::cal_parameters_for_main_domain) {
			boundary_condition_for_domain_U = boundary_U;
			boundary_condition_for_domain_V = boundary_V;
			boundary_condition_for_domain_W = boundary_W;
			boundary_condition_for_pressure = boundary_main_node;
			cal_parameters_for_main_domain = cal_parameters;
		}

		void cal_parameters_before_calculation();

		//return [MAX_ABS_dU, MAX_ABS_dV, MAX_ABS_dW]
		Vector3 evolve_momentum_equation(double fluid_dt);

		//return [MAX_abs_dPressure, pressure_iterate_times]
		vector<double> do_pressure_correction(double fluid_threshold = 0.5, double accuracy = 1e-6, int max_iterate_steps = 1000
			, bool debug_solver = false, int output_step = 1000);

		void correcting_velocity_field(double fluid_dt);

		void assign_velocity_to_main_domain();
		void assign_velocity_to_fluid_domain();

		void free() {
			phaseMesh = nullptr;
			boundary_condition_for_domain_U = nullptr;
			boundary_condition_for_domain_V = nullptr;
			boundary_condition_for_domain_W = nullptr;
			boundary_condition_for_pressure = nullptr;
			cal_parameters_for_main_domain = nullptr;
			U.free();
			V.free();
			W.free();
		}

		void(*boundary_condition_for_domain_U)(VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int);
		void(*boundary_condition_for_domain_V)(VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int);
		void(*boundary_condition_for_domain_W)(VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int);
		void(*boundary_condition_for_pressure)(PhaseNode&, int, int, int);
		void(*cal_parameters_for_main_domain)(PhaseNode&, int, int, int);
	private:
		FieldStorage_forPhaseNode* phaseMesh;
		FieldStorage_forVector U;
		FieldStorage_forVector V;
		FieldStorage_forVector W;
	};

	enum LBM_LATTICE_MODEL { LBM_D2Q9, LBM_D3Q19 };
	namespace lbm_funcs {
		static void init_distribution_functions(pf::PhaseNode& node, int LBM_F_INDEX) {
			return;
		}
		static void collision(pf::PhaseNode & node, int LBM_F_INDEX) {
			return;
		}
		static void cal_macro_variables(pf::PhaseNode& node, int LBM_F_INDEX) {
			return;
		}
		static void boundary_condition(pf::PhaseNode& node, int LBM_F_INDEX) {
			return;
		}
	}
	class LBM
	{
	public:
		LBM() {};
		LBM(FieldStorage_forPhaseNode& _phaseMesh, int _LBM_F_INDEX = ExternalFields::LBM_Symbols_INDEX_0) {
			init(_phaseMesh, _LBM_F_INDEX);
		}
		~LBM() {
			free();
		}
		void init(FieldStorage_forPhaseNode& _phaseMesh, int _LBM_F_INDEX = ExternalFields::LBM_Symbols_INDEX_0) {
			phaseMesh = &_phaseMesh;
			if (phaseMesh->_dimention == Dimension::Two_Dimension) {
				lbm_lattice_model = LBM_LATTICE_MODEL::LBM_D2Q9;

			}
			else if (phaseMesh->_dimention == Dimension::Three_Dimension) {
				lbm_lattice_model = LBM_LATTICE_MODEL::LBM_D3Q19;

			}
			else {
				cout << "> error, LBM don't work in one-dimension in MID software !" << endl;
				SYS_PROGRAM_STOP;
			}
			LBM_F_INDEX = _LBM_F_INDEX;
			solver_name = "DistributionFunction";
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				node->customValues.add_double(ExternalFields::FLUID_pressure, 0.0);
				node->customValues.add_double(ExternalFields::FLUID_mass, 0.0);
				node->customVec3s.add_vec(ExternalFields::FLUID_velocity, Vector3(0.0, 0.0, 0.0));
				// distribution function
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_0, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_1, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_2, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_3, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_4, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_5, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_6, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_7, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_8, 0.0); 
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_0, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_1, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_2, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_3, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_4, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_5, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_6, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_7, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_8, 0.0);
				node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_macro, 0.0);
				node->customVec3s.add_vec(LBM_F_INDEX + LBM_Symbols::LBM_fv_macro, Vector3(0.0, 0.0, 0.0));
				if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_9 , 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_10, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_11, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_12, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_13, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_14, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_15, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_16, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_17, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_f_18, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_9, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_10, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_11, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_12, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_13, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_14, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_15, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_16, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_17, 0.0);
					node->customValues.add_double(LBM_F_INDEX + LBM_Symbols::LBM_m_18, 0.0);
				}
			}
			_init_distribution_functions = lbm_funcs::init_distribution_functions;
			_collision = lbm_funcs::collision;
			_boundary_condition = lbm_funcs::boundary_condition;
			_cal_macro_variables = lbm_funcs::cal_macro_variables;
		}

		void free() {
			phaseMesh = nullptr;
			_init_distribution_functions = nullptr;
			_collision = nullptr;
			_boundary_condition = nullptr;
			_cal_macro_variables = nullptr;
		}

		void init_distribution_functions();


		/*
		// LBM_D2Q9 
			f0  c( 0,  0)
			f1  c( 1,  0)
			f2  c( 0,  1)
			f3  c(-1,  0)
			f4  c( 0, -1)
			f5  c( 1,  1)
			f6  c(-1,  1)
			f7  c(-1, -1)
			f8  c( 1, -1)
		//-----------------------------------------------------------------------------------------------
		// LBM_D3Q19 
			f0   c( 0,  0,  0)
			f1   c( 1,  0,  0)
			f2   c(-1,  0,  0)
			f3   c( 1,  0,  0)
			f4   c(-1,  0,  0)
			f5   c( 1,  0,  0)
			f6   c(-1,  0,  0)
			f7   c( 1,  1,  0)
			f8   c(-1,  1,  0)
			f9   c(-1, -1,  0)
			f10  c( 1, -1,  0)
			f11  c( 1,  0,  1)
			f12  c(-1,  0,  1)
			f13  c(-1,  0, -1)
			f14  c( 1,  0, -1)
			f15  c( 0,  1,  1)
			f16  c( 0, -1,  1)
			f17  c( 0, -1, -1)
			f18  c( 0,  1, -1)
		*/
		void collision();

		// return MAX_VARIATION
		void streaming();

		void boundary_condition();

		Vector3 cal_macro_variables();

		void cal_macro_variables(double& F_MACRO_MAX_VARIATION);

		void (*_init_distribution_functions)(pf::PhaseNode&, int);

		void (*_collision)(pf::PhaseNode&, int);

		void (*_boundary_condition)(pf::PhaseNode&, int);

		void (*_cal_macro_variables)(pf::PhaseNode&, int);

		LBM_LATTICE_MODEL lbm_lattice_model;
		int LBM_F_INDEX;
		string solver_name;
	private:
		FieldStorage_forPhaseNode* phaseMesh;
	};

}