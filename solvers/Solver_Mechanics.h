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
	enum VelocityDomainIndex { VDIndex_OLD, VDIndex_NOW, VDIndex_FUTURE };
	namespace mechanical_boundary_condition_funcs {
		static void cal_parameters_for_main_domain(PhaseNode& node, int Nx, int Ny, int Nz) {
			node.customVec6s[ExternalFields::MECH_eigen_strain].set_to_zero();
			node.customMatrix6x6s[ExternalFields::MECH_stiffness].set_to_zero();
			return;
		}
	};
	class MechanicalField_Explicit
	{
	public:
		MechanicalField_Explicit(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition x_up_bc,
			BoundaryCondition x_down_bc, BoundaryCondition y_up_bc, BoundaryCondition y_down_bc, BoundaryCondition z_up_bc, BoundaryCondition z_down_bc) {
			init(_phaseMesh, x_up_bc, x_down_bc, y_up_bc, y_down_bc, z_up_bc, z_down_bc);
		};
		MechanicalField_Explicit() {};
		~MechanicalField_Explicit() {
			free();
		};
		void init(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition x_up_bc,
			BoundaryCondition x_down_bc, BoundaryCondition y_up_bc, BoundaryCondition y_down_bc, BoundaryCondition z_up_bc, BoundaryCondition z_down_bc) {
			phaseMesh = &_phaseMesh;
			AvgStrainMask.resize(3);
			LoadStressMask.resize(3);
			AppStrainMask.resize(3);
			for (int index = 0; index < 3; index++) {
				AvgStrainMask[index] = 1;
				LoadStressMask[index] = 0;
				AppStrainMask[index] = 0;
			}
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
				vecNode->vals.push_back(0.0); // value1: 0
				vecNode->vals.push_back(0.0); // value2: 1
				vecNode->vals.push_back(0.0); // increment: 2
			}
			for (auto vecNode = V._mesh.begin(); vecNode < V._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value1
				vecNode->vals.push_back(0.0); // value2
				vecNode->vals.push_back(0.0); // increment
			}
			for (auto vecNode = W._mesh.begin(); vecNode < W._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value1
				vecNode->vals.push_back(0.0); // value2
				vecNode->vals.push_back(0.0); // increment
			}
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				Vector6 vec6; Matrix6x6 matrix;
				node->customVec6s.add_vec(ExternalFields::MECH_stress, vec6);
				node->customVec6s.add_vec(ExternalFields::MECH_strain, vec6);
				node->customVec6s.add_vec(ExternalFields::MECH_eigen_strain, vec6);
				node->customMatrix6x6s.add_matrix(ExternalFields::MECH_stiffness, matrix);
			}
			cal_parameters_for_main_domain = mechanical_boundary_condition_funcs::cal_parameters_for_main_domain;
			velocity_doamin_index.push_back(VelocityDomainIndex::VDIndex_OLD);
			velocity_doamin_index.push_back(VelocityDomainIndex::VDIndex_NOW);
			velocity_doamin_index.push_back(VelocityDomainIndex::VDIndex_FUTURE);
		}

		void define_funcs_for_mechanics(vector<bool> _avgStrainMask, vector<bool> _loadStressMask, vector<bool> _appStrainMask,
			vStress _applied_stress, vStrain _applied_strain, 
			void(*cal_parameters)(pf::PhaseNode&, int, int, int) = mechanical_boundary_condition_funcs::cal_parameters_for_main_domain) {
			cal_parameters_for_main_domain = cal_parameters;
		}

		void init_velocity_field();

		void cal_stress(vStress& AverageStress, bool is_plastic);

		void cal_parameters_before_calculation();

		double boundary_condition(vStress& AverageStress, vStrain& AverageStrain, double incre_rate);

		// return MAX_abs_dstrain in information
		double evolve_momentum_equation(double mass_density, double mechanic_dt, vStrain& AverageStrain);

		void free() {
			phaseMesh = nullptr;
			cal_parameters_for_main_domain = nullptr;
			velocity_doamin_index.clear();
			U.free();
			V.free();
			W.free();
		}

		double get_u_main_node(int _x, int _y, int _z);
		double get_v_main_node(int _x, int _y, int _z);
		double get_w_main_node(int _x, int _y, int _z);

		//----------------------------------------------------------------- settings
		vStress applied_stress;
		vStrain applied_strain;

		vector<bool>     AvgStrainMask;
		vector<bool>     LoadStressMask;
		vector<bool>     AppStrainMask;
		void(*cal_parameters_for_main_domain)(PhaseNode&, int, int, int);
		//----------------------------------------------------------------- settings
	private:
		Matrix6x6   average_stiffness;
		Matrix6x6   average_compliences;
		FieldStorage_forPhaseNode* phaseMesh;
		FieldStorage_forVector U;
		FieldStorage_forVector V;
		FieldStorage_forVector W;
		vector<int> velocity_doamin_index;
	};

	class MechanicalField_Implicit
	{
	public:
		MechanicalField_Implicit() {};
		MechanicalField_Implicit(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition _x_bc, BoundaryCondition _y_bc,  BoundaryCondition _z_bc) {
			init(_phaseMesh, _x_bc, _y_bc, _z_bc);
		}
		~MechanicalField_Implicit() {
			mechanicalField.free();
			phaseMesh = nullptr;
			for (int n = 0; n < 6; n++)
			{
				fftw_destroy_plan(ForwardPlanRHS[n]);

				delete[] rlRHSide[n];
				delete[] rcRHSide[n];
			}
			for (int n = 0; n < 3; n++)
			{
				fftw_destroy_plan(BackwardPlanU[n]);

				delete[] rlU[n];
				delete[] rcU[n];
				delete[] Q[n];
			}
			for (int n = 0; n < 9; n++)
			{
				fftw_destroy_plan(BackwardPlanDefGrad[n]);
				delete[] rcDefGrad[n];
				delete[] rlDefGrad[n];
			}
		}
		void init(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition _x_bc, BoundaryCondition _y_bc, BoundaryCondition _z_bc) {
			phaseMesh = &_phaseMesh;
			AvgStrainMask.resize(3);
			LoadStressMask.resize(3);
			AppStrainMask.resize(3);
			for (int index = 0; index < 3; index++) {
				AvgStrainMask[index] = 1;
				LoadStressMask[index] = 0;
				AppStrainMask[index] = 0;
			}

			mechanicalField.init(_phaseMesh.limit_x, _phaseMesh.limit_y, _phaseMesh.limit_z, _phaseMesh.dr, _x_bc, _y_bc, _z_bc);

			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				Vector6 vec6; Matrix6x6 matrix;
				node->customVec6s.add_vec(ExternalFields::MECH_stress, vec6);
				node->customVec6s.add_vec(ExternalFields::MECH_strain, vec6);
				node->customVec6s.add_vec(ExternalFields::MECH_eigen_strain, vec6);
				node->customMatrix6x6s.add_matrix(ExternalFields::MECH_stiffness, matrix);
			}

			Nz2 = (mechanicalField.Nz) / 2 + 1;
			rlSIZE = mechanicalField.Nx * mechanicalField.Ny * mechanicalField.Nz;
			rcSIZE = mechanicalField.Nx * mechanicalField.Ny * Nz2;

			DPi_Nx = 2.0 * PI / double(mechanicalField.Nx);
			DPi_Ny = 2.0 * PI / double(mechanicalField.Ny);
			DPi_Nz = 2.0 * PI / double(mechanicalField.Nz);

			MAX_ElasticConstants.set_to_zero();

			Norm = 1.0 / double(rlSIZE);
			// Arrays allocation:
			for (int n = 0; n < 6; n++)
			{
				rlRHSide[n] = new double[rlSIZE]();
				rcRHSide[n] = new complex<double>[rcSIZE]();
			}
			for (int n = 0; n < 3; n++)
			{
				rlU[n] = new double[rlSIZE]();
				rcU[n] = new complex<double>[rcSIZE]();
				Q[n] = new double[rcSIZE]();
			}
			for (int n = 0; n < 9; n++)
			{
				rlDefGrad[n] = new double[rlSIZE]();
				rcDefGrad[n] = new complex<double>[rcSIZE]();
			}

			// set Q
			for (int i = 0; i < mechanicalField.Nx; i++)
				for (int j = 0; j < mechanicalField.Ny; j++)
					for (int k = 0; k < Nz2; k++)
					{
						int XYZ = k + Nz2 * (j + mechanicalField.Ny * i);

						Q[0][XYZ] = DPi_Nx * (i * (i <= mechanicalField.Nx / 2) - (mechanicalField.Nx - i) * (i > mechanicalField.Nx / 2)) / mechanicalField.dx;
						Q[1][XYZ] = DPi_Ny * (j * (j <= mechanicalField.Ny / 2) - (mechanicalField.Ny - j) * (j > mechanicalField.Ny / 2)) / mechanicalField.dx;
						Q[2][XYZ] = DPi_Nz * (k * (k <= mechanicalField.Nz / 2) - (mechanicalField.Nz - k) * (k > mechanicalField.Nz / 2)) / mechanicalField.dx;
					}

			for (int n = 0; n < 6; n++)
			{
				ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
				(mechanicalField.Nx, mechanicalField.Ny, mechanicalField.Nz, rlRHSide[n],
					reinterpret_cast<fftw_complex*> (rcRHSide[n]),
					FFTW_PATIENT);
			}

			for (int n = 0; n < 9; n++)
			{
				BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
				(mechanicalField.Nx, mechanicalField.Ny, mechanicalField.Nz,
					reinterpret_cast<fftw_complex*> (rcDefGrad[n]),
					rlDefGrad[n],
					FFTW_PATIENT);
			}

			for (int n = 0; n < 3; n++)
			{
				BackwardPlanU[n] = fftw_plan_dft_c2r_3d
				(mechanicalField.Nx, mechanicalField.Ny, mechanicalField.Nz,
					reinterpret_cast<fftw_complex*> (rcU[n]),
					rlU[n],
					FFTW_PATIENT);
			}
		}
		void define_funcs_for_mechanics(vector<bool> _avgStrainMask, vector<bool> _loadStressMask, vector<bool> _appStrainMask,
			vStress _applied_stress, vStrain _applied_strain, 
			void(*cal_parameters)(pf::PhaseNode&, int, int, int) = mechanical_boundary_condition_funcs::cal_parameters_for_main_domain) {
			AvgStrainMask = _avgStrainMask;
			LoadStressMask = _loadStressMask;
			AppStrainMask = _appStrainMask;
			applied_stress = _applied_stress;
			applied_strain = _applied_strain;
			cal_parameters_for_main_domain = cal_parameters;
		}
		void free() {
			mechanicalField.free();
			phaseMesh = nullptr;
			for (int n = 0; n < 6; n++)
			{
				fftw_destroy_plan(ForwardPlanRHS[n]);

				delete[] rlRHSide[n];
				delete[] rcRHSide[n];
			}

			for (int n = 0; n < 3; n++)
			{
				fftw_destroy_plan(BackwardPlanU[n]);

				delete[] rlU[n];
				delete[] rcU[n];
				delete[] Q[n];
			}
			for (int n = 0; n < 9; n++)
			{
				fftw_destroy_plan(BackwardPlanDefGrad[n]);
				delete[] rcDefGrad[n];
				delete[] rlDefGrad[n];
			}
		}

		void cal_parameters_before_calculation(bool is_plastic);
		void recal_eigenstrain_with_plasticity();
		void	SetMAXElasticConstants(vector<Matrix6x6> Cijs);
		void	SetMAXElasticConstants(tensor1_matrix6 Cijs);

		double get_u_main_node(int _x, int _y, int _z);
		double get_v_main_node(int _x, int _y, int _z);
		double get_w_main_node(int _x, int _y, int _z);

		// Ingo Steinbach Method
		void initStrainIncrements();
		string Solve(double StrainAccuracy, int MAXIterations, WriteToFile& writer, double incre_rate = 1.0, bool is_dvStraindt_output = false, bool getU = false);

		// Armen G. Khachaturyan Method
		void initVirtualEigenstrain();
		string Solve2(double StrainAccuracy, int MAXIterations, double iterate_rate, WriteToFile& writer, bool is_dvStraindt_output = false, bool getU = false);

//----------------------------------------------------------------- settings
		vStress applied_stress;
		vStrain applied_strain;

		vector<bool>     AvgStrainMask;
		vector<bool>     LoadStressMask;
		vector<bool>     AppStrainMask;

		void(*cal_parameters_for_main_domain)(PhaseNode&, int, int, int);
//----------------------------------------------------------------- settings

	private:
		PhaseNode& get_phaseNode_in_elasticField(int elas_x, int elas_y, int elas_z);
		// Ingo Steinbach
		void CalculateRHS(Matrix6x6 Cij);
		void ExecuteForwardFFT();
		void CalculateFourierSolution(Matrix6x6 Cij, bool getU);
		void ExecuteBackwardFFT(bool getU);
		void SetElasticProperties1(double& MAXStrainDifference, vStress& AverageStress);
		void SetElasticProperties2(double& MAXStrainDifference);
		void SetElasticBoundaryConditions(vStress TargetStress);
		void assignment_to_phaseField();
		// Armen G. Khachaturyan
		void CalculateRHS2(Matrix6x6& Cij);
		void evaluate_virtualEigenstrain(Matrix6x6 Cij, Matrix6x6 Sij, double& MAXvStrainDifference, double iterate_rate);
		void SetElasticBoundaryConditions2(Matrix6x6& Cij, Matrix6x6& Sij);
		FieldStorage_forMechanicNode mechanicalField;

		vStrain average_strain;
		vStrain average_virtual_strain;

		//Elasticity Tensors:
		Matrix6x6   average_stiffness;
		Matrix6x6   average_compliences;
		//Matrix6x6   MINCompliences;

		FieldStorage_forPhaseNode* phaseMesh;

		int Nz2;                                              

		int rlSIZE;                                           
		int rcSIZE;                                           

		double Norm;                                          
		double DPi_Nx;                                        
		double DPi_Ny;                                        
		double DPi_Nz;                                        

		//Matrix6x6 C0;
		//Matrix6x6 C0inverse;

		Matrix6x6 MAX_ElasticConstants;

		double* Q[3];

		double* rlRHSide[6];                                  
		std::complex<double>* rcRHSide[6];                    

		double* rlU[3];
		std::complex<double>* rcU[3];							

		double* rlDefGrad[9];                                 
		std::complex<double>* rcDefGrad[9];                   

		fftw_plan               ForwardPlanRHS[6];            
		fftw_plan               BackwardPlanDefGrad[9];       
		fftw_plan               BackwardPlanU[3];             
	};
}
