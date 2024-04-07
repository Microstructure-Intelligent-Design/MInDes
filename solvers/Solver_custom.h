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
#include "base.h"
using namespace std;
namespace pf {
	enum CUSTOM_SOLVER { SOLVER_ALLEN_CAHN = 10000, SOLVER_CAHN_HILLIARD = 11000, BUFF_FOR_SOLVER_CAHN_HILLIARD = 12000 };
	namespace allenCahnEquationSolver {
		// grain_index from 0 to grain_number - 1
		static double dF_dphi_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
			return 0.0;
		}
		static double L_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
			return 0.0;
		}
		static double Source_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
			return 0.0;
		}
		static void BoundaryCondition(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return;
		}
	}

	// partial(phi_i) / partial(t) = - Mobility_i * variation(F) / variation(phi_i) + Source
	class CustomAllenCahnSolver {
	public:
		CustomAllenCahnSolver() {};
		CustomAllenCahnSolver(FieldStorage_forPhaseNode& _phaseMesh, int _grains_number = 0, int _grains_start_index = SOLVER_ALLEN_CAHN, string _solver_name = "AllenCahnSolver") {
			init(_phaseMesh, _grains_number, _grains_start_index, _solver_name);
		}
		~CustomAllenCahnSolver() {
			clear();
		};
		void init(FieldStorage_forPhaseNode& _phaseMesh, int _grains_number = 0, int _grains_start_index = SOLVER_ALLEN_CAHN, string _solver_name = "AllenCahnSolver") {
			solver_name = _solver_name;
			phaseMesh = &_phaseMesh;
			grains_number = _grains_number;
			grains_start_index = _grains_start_index;
			dF_dphi = allenCahnEquationSolver::dF_dphi_cal;
			L = allenCahnEquationSolver::L_cal;
			Source = allenCahnEquationSolver::Source_cal;
			BoundaryCondition = allenCahnEquationSolver::BoundaryCondition;
		}
		void define_funcs_for_AC_solver(double(*dF_dphi_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::dF_dphi_cal,
			double(*L_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::L_cal,
			double(*Source_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::Source_cal,
			void(*BoundaryCondition_cal)(pf::PhaseNode&, int, int) = allenCahnEquationSolver::BoundaryCondition) {
			dF_dphi = dF_dphi_cal;
			L = L_cal;
			Source = Source_cal;
			BoundaryCondition = BoundaryCondition_cal;
		}
		void init_mesh() {
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				for (int grain = 0; grain < grains_number; grain++) {
					node->customValues.add_double((grain + grains_start_index), 0.0);
					node->customValues.add_double(-(grain + grains_start_index), 0.0);
				}
			}
		}
		string print_model() {
			stringstream rep;
			rep << "External Allen Cahn Model : dphi_dt = - L * df_dphi + S " << endl;
			return rep.str();
		}
		// retuan MAX_PHI_VARIATION;
		double solve_one_step(double dt, bool adjust_phi_0_1 = false);
		string solver_name;
		int grains_number;
		int grains_start_index;
		double(*dF_dphi)(pf::PhaseNode&, int, int);
		double(*L)(pf::PhaseNode&, int, int);
		double(*Source)(pf::PhaseNode&, int, int);
		void(*BoundaryCondition)(pf::PhaseNode&, int, int);
	private:
		FieldStorage_forPhaseNode* phaseMesh;
		void clear() {
			phaseMesh = nullptr;
			dF_dphi = nullptr;
			L = nullptr;
			Source = nullptr;
			BoundaryCondition = nullptr;
		}
	};

	namespace cahnHilliardEquationSolver {
		// grain_index from 0 to grain_number - 1
		static double dF_dc_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return 0.0;
		}
		static double Mobility_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return 0.0;
		}
		static double Source_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return 0.0;
		}
		static void BoundaryCondition(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return;
		}
	}

	// partial(c_i) / partial(t) = delt(M_i * delt(variation(F) / variation(c_i))) + Source
	class CustomCahnHilliardSolver {
	public:
		CustomCahnHilliardSolver() { };
		CustomCahnHilliardSolver(FieldStorage_forPhaseNode& _phaseMesh, int _components_number = 0, int _comps_start_index = SOLVER_CAHN_HILLIARD, string _solver_name = "CahnHilliardSolver") {
			init(_phaseMesh, _components_number, _comps_start_index, _solver_name);
		}
		~CustomCahnHilliardSolver() {
			clear();
		};
		void init(FieldStorage_forPhaseNode& _phaseMesh, int _components_number = 0, int _comps_start_index = SOLVER_CAHN_HILLIARD, string _solver_name = "CahnHilliardSolver") {
			solver_name = _solver_name;
			phaseMesh = &_phaseMesh;
			components_number = _components_number;
			comps_start_index = _comps_start_index;
			dF_dc = cahnHilliardEquationSolver::dF_dc_cal;
			Mobility = cahnHilliardEquationSolver::Mobility_cal;
			Source = cahnHilliardEquationSolver::Source_cal;
			BoundaryCondition = cahnHilliardEquationSolver::BoundaryCondition;
		}
		void init_mesh() {
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				for (int comp = 0; comp < components_number; comp++) {
					node->customValues.add_double(comp + comps_start_index, 0.0);
					node->customValues.add_double(-(comp + comps_start_index), 0.0);
					node->customValues.add_double(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD, 0.0);
					node->customValues.add_double(-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD), 0.0);
				}
			}
		}
		void define_funcs_for_AC_solver(double(*dF_dc_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::dF_dc_cal,
			double(*Mobility_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::Mobility_cal,
			double(*Source_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::Source_cal,
			void(*BoundaryCondition_cal)(pf::PhaseNode&, int, int) = cahnHilliardEquationSolver::BoundaryCondition) {
			dF_dc = dF_dc_cal;
			Mobility = Mobility_cal;
			Source = Source_cal;
			BoundaryCondition = BoundaryCondition_cal;
		}
		string print_model() {
			stringstream rep;
			rep << "External Cahn Hilliard Model : dc_dt = \\labla{ M * \\labla( df_dc ) } + S " << endl;
			return rep.str();
		}
		// return MAX_COMP_VARIATION
		double solve_one_step(double dt, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT, bool adjust_phi_0_1 = false);
		string solver_name;
		int components_number;
		int comps_start_index;
	private:
		FieldStorage_forPhaseNode* phaseMesh;
		double(*dF_dc)(pf::PhaseNode&, int, int);
		double(*Mobility)(pf::PhaseNode&, int, int);
		double(*Source)(pf::PhaseNode&, int, int);
		void(*BoundaryCondition)(pf::PhaseNode&, int, int);
		void clear() {
			phaseMesh = nullptr;
			dF_dc = nullptr;
			Mobility = nullptr;
			Source = nullptr;
			BoundaryCondition = nullptr;
		}
	};

	namespace fourier_transform_solver_funcs {
		static void fill_node_real_space(vector<std::complex<double>>& real_space, PhaseNode& node, int real_x, int real_y, int real_z) {
			for (int index = 0; index < real_space.size(); index++)
				real_space[index]._Val[FFTW_REAL] = 0.0;
		}
		static std::complex<double> dynamic_equation_fourier_space(vector<std::complex<double>> fourier_space, PhaseNode& node, double Q2, double Q4) {
			return { fourier_space[0]._Val[FFTW_REAL], fourier_space[0]._Val[FFTW_IMAG] };
		}
		static void boundary_condition_real_space(std::complex<double>& basic_real_space, PhaseNode& node, int real_x, int real_y, int real_z) {
			return;
		}
	}

	class FourierTransformSolver
	{
	public:
		FourierTransformSolver() {};
		FourierTransformSolver(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition _x_bc, BoundaryCondition _y_bc, BoundaryCondition _z_bc) {
			init(_phaseMesh, _x_bc, _y_bc, _z_bc);
		};
		~FourierTransformSolver() {
			mirrorMesh.free();
			phaseMesh = nullptr;
			fill_node_real_space = nullptr;
			dynamic_equation_fourier_space = nullptr;
			boundary_condition_real_space = nullptr;
			for (int n = 0; n < ForwardPlan.size(); n++)
			{
				fftw_destroy_plan(ForwardPlan[n]);
			}
			fftw_destroy_plan(BackwardPlan);

			for (int n = 0; n < 3; n++)
			{
				delete[] Q[n];
			}
			for (int n = 0; n < rlDomain.size(); n++)
			{
				delete[] rlDomain[n];
			}
			for (int n = 0; n < rcDomain.size(); n++)
			{
				delete[] rcDomain[n];
			}
		}

		void init(FieldStorage_forPhaseNode& _phaseMesh, BoundaryCondition _x_bc, BoundaryCondition _y_bc, BoundaryCondition _z_bc) {
			phaseMesh = &_phaseMesh;

			mirrorMesh.init(_phaseMesh.limit_x, _phaseMesh.limit_y, _phaseMesh.limit_z, _phaseMesh.dr, _x_bc, _y_bc, _z_bc);

			//Nz2 = (mirrorMesh.Nz) / 2 + 1;
			rlcSIZE = mirrorMesh.Nx * mirrorMesh.Ny * mirrorMesh.Nz;
			//rcSIZE = mirrorMesh.Nx * mirrorMesh.Ny * Nz2;

			DPi_Nx = 2.0 * PI / double(mirrorMesh.Nx);
			DPi_Ny = 2.0 * PI / double(mirrorMesh.Ny);
			DPi_Nz = 2.0 * PI / double(mirrorMesh.Nz);

			Norm = 1.0 / double(rlcSIZE);
			// Arrays allocation:
			for (int n = 0; n < 3; n++)
			{
				Q[n] = new double[rlcSIZE]();
			}

			std::complex<double>* rlBaseDomain = new std::complex<double>[rlcSIZE]();
			std::complex<double>* rcBaseDomain = new std::complex<double>[rlcSIZE]();
			rlDomain.push_back(rlBaseDomain);
			rcDomain.push_back(rcBaseDomain);
			rlc_space_size = 1;

			// set Q
			for (int i = 0; i < mirrorMesh.Nx; i++)
				for (int j = 0; j < mirrorMesh.Ny; j++)
					for (int k = 0; k < mirrorMesh.Nz; k++)
					{
						int XYZ = k + mirrorMesh.Nz * (j + mirrorMesh.Ny * i);
						Q[0][XYZ] = DPi_Nx * (i * (i <= mirrorMesh.Nx / 2) - (mirrorMesh.Nx - i) * (i > mirrorMesh.Nx / 2)) / mirrorMesh.dx;
						Q[1][XYZ] = DPi_Ny * (j * (j <= mirrorMesh.Ny / 2) - (mirrorMesh.Ny - j) * (j > mirrorMesh.Ny / 2)) / mirrorMesh.dx;
						Q[2][XYZ] = DPi_Nz * (k * (k <= mirrorMesh.Nz / 2) - (mirrorMesh.Nz - k) * (k > mirrorMesh.Nz / 2)) / mirrorMesh.dx;
					}

			// fftw_complex* rlDomain;
			// fftw_plan_dft_3d(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, rlDomain[index], rcDomain[index], FFTW_FORWARD, FFTW_MEASURE);
			// fftw_plan_dft_3d(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, rlDomain[index], rcDomain[index], FFTW_BACKWARD, FFTW_MEASURE);

			ForwardPlan.push_back(fftw_plan_dft_3d(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, reinterpret_cast<fftw_complex*>(rlDomain[0]), reinterpret_cast<fftw_complex*>(rcDomain[0]), FFTW_FORWARD, FFTW_MEASURE));
			/*ForwardPlan[0] = fftw_plan_dft_r2c_3d
			(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, rlDomain[0],
				reinterpret_cast<fftw_complex*> (rcDomain[0]),
				FFTW_PATIENT);*/

			BackwardPlan = fftw_plan_dft_3d(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, reinterpret_cast<fftw_complex*>(rcDomain[0]), reinterpret_cast<fftw_complex*>(rlDomain[0]), FFTW_FORWARD, FFTW_MEASURE);
			/*BackwardPlan = fftw_plan_dft_c2r_3d
			(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz,
				reinterpret_cast<fftw_complex*> (rcDomain[0]),
				rlDomain[0],
				FFTW_PATIENT);*/

			// rcDomain[0][1][FFTW_REAL] = rcDomain[0][1][FFTW_REAL];
			// rcDomain[0][1][FFTW_IMAG] = rcDomain[0][1][FFTW_IMAG];

			init_real_space = fourier_transform_solver_funcs::boundary_condition_real_space;

			fill_node_real_space = fourier_transform_solver_funcs::fill_node_real_space;

			dynamic_equation_fourier_space = fourier_transform_solver_funcs::dynamic_equation_fourier_space;

			boundary_condition_real_space = fourier_transform_solver_funcs::boundary_condition_real_space;

		}

		int new_one_real_and_fourier_space() {
			int index = int(rlDomain.size());

			std::complex<double>* rlBaseDomain = new std::complex<double>[rlcSIZE]();
			std::complex<double>* rcBaseDomain = new std::complex<double>[rlcSIZE]();
			rlDomain.push_back(rlBaseDomain);
			rcDomain.push_back(rcBaseDomain);

			ForwardPlan.push_back(fftw_plan_dft_3d(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, reinterpret_cast<fftw_complex*>(rlDomain[index]), reinterpret_cast<fftw_complex*>(rcDomain[index]), FFTW_FORWARD, FFTW_MEASURE));
			/*ForwardPlan[0] = fftw_plan_dft_r2c_3d
			(mirrorMesh.Nx, mirrorMesh.Ny, mirrorMesh.Nz, rlDomain[index],
				reinterpret_cast<fftw_complex*> (rcDomain[index]),
				FFTW_PATIENT);*/

			rlc_space_size = int(rlDomain.size());

			return index;
		}

		double inline Q2(int x, int y, int z) {
			int XYZ = z + mirrorMesh.Nz * (y + mirrorMesh.Ny * x);
			return Q[0][XYZ] * Q[0][XYZ] + Q[1][XYZ] * Q[1][XYZ] + Q[2][XYZ] * Q[2][XYZ];
		}

		double inline Q2(int XYZ) {
			return Q[0][XYZ] * Q[0][XYZ] + Q[1][XYZ] * Q[1][XYZ] + Q[2][XYZ] * Q[2][XYZ];
		}

		double inline Q4(int x, int y, int z) {
			double q2 = Q2(x, y, z);
			return q2 * q2;
		}

		double inline Q4(int XYZ) {
			double q2 = Q2(XYZ);
			return q2 * q2;
		}

		void init_basic_real_space();

		void solve_one_step();

		// this function will treatment the basic real space, for rlNode[0] & main node
		void (*init_real_space)(std::complex<double>& basic_real_space, PhaseNode& node, int real_x, int real_y, int real_z);

		// this function will not fill the basic real space, for rlNode[0]
		void (*fill_node_real_space)(vector<std::complex<double>>& real_space, PhaseNode& node, int real_x, int real_y, int real_z);

		// this function will calculate basic fourier space, for rcNode[0]
		std::complex<double>(*dynamic_equation_fourier_space)(vector<std::complex<double>> fourier_space, PhaseNode& node, double Q2, double Q4);

		// this function will treatment the basic real space, for rlNode[0] & main node
		void (*boundary_condition_real_space)(std::complex<double>& basic_real_space, PhaseNode& node, int real_x, int real_y, int real_z);

	private:

		FieldStorage_forCustomNode mirrorMesh;

		FieldStorage_forPhaseNode* phaseMesh;

		PhaseNode& get_phaseNode_by_mirrorMesh_XYZ(int mirror_x, int mirror_y, int mirror_z) {
			if (mirror_x >= phaseMesh->limit_x) {
				mirror_x = mirrorMesh.Nx - 1 - mirror_x;
			}
			if (mirror_y >= phaseMesh->limit_y) {
				mirror_y = mirrorMesh.Ny - 1 - mirror_y;
			}
			if (mirror_z >= phaseMesh->limit_z) {
				mirror_z = mirrorMesh.Nz - 1 - mirror_z;
			}
			return (*phaseMesh)(mirror_x, mirror_y, mirror_z);
		}

		//int Nz2;

		int rlcSIZE;
		//int rcSIZE;

		double Norm;
		double DPi_Nx;
		double DPi_Ny;
		double DPi_Nz;

		double* Q[3];

		int rlc_space_size;

		vector<std::complex<double>*> rlDomain;
		vector<std::complex<double>*> rcDomain;

		vector<fftw_plan>             ForwardPlan;
		fftw_plan                     BackwardPlan;
	};

}