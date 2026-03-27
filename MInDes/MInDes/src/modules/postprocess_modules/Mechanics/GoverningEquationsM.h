#pragma once
#include "../../base/Mesh_0.h"
#include "../../base/VectorMatrix.h"
#include "../../input_modules/ioFiles_Params.h"
//x64
#include "../../../../lib/x64/fftw3.h"
#include <complex>
#include "MechanicalPoint.h"
const std::complex<double> I(0.0, 1.0);
namespace pf {
	enum VelocityDomainIndex { VDIndex_OLD, VDIndex_NOW, VDIndex_FUTURE };
	namespace mechanical_boundary_condition_funcs {
		inline Matrix6x6 cal_stiffness(int Nx, int Ny, int Nz) {
			return Matrix6x6();
		}
		inline vStrain cal_eigenstrain(int Nx, int Ny, int Nz) {
			return vStrain();
		}
	};
	class MechanicalField_Implicit
	{
	public:
		MechanicalField_Implicit() {};
		~MechanicalField_Implicit() {
			free();
		}
		void init(int _Nx, int _Ny, int _Nz, BoundaryCondition _x_bc, BoundaryCondition _y_bc, BoundaryCondition _z_bc, Mesh<ElasticPoint>& _elastic_field) {
			AvgStrainMask.resize(3);
			LoadStressMask.resize(3);
			AppStrainMask.resize(3);
			for (int index = 0; index < 3; index++) {
				AvgStrainMask[index] = 1;
				LoadStressMask[index] = 0;
				AppStrainMask[index] = 0;
			}
			Nx = _Nx; mech_Nx = _Nx;
			Ny = _Ny; mech_Ny = _Ny;
			Nz = _Nz; mech_Nz = _Nz;
			if (_x_bc != BoundaryCondition::PERIODIC)
				mech_Nx = _Nx * 2;
			if (_y_bc != BoundaryCondition::PERIODIC)
				mech_Ny = _Ny * 2;
			if (_z_bc != BoundaryCondition::PERIODIC)
				mech_Nz = _Nz * 2;
			elastic_field = &_elastic_field;
			elastic_field->init(mech_Nx, mech_Ny, mech_Nz, 1);

			Nz2 = (mech_Nz) / 2 + 1;
			rlSIZE = mech_Nx * mech_Ny * mech_Nz;
			rcSIZE = mech_Nx * mech_Ny * Nz2;

			DPi_Nx = 2.0 * PI / double(mech_Nx);
			DPi_Ny = 2.0 * PI / double(mech_Ny);
			DPi_Nz = 2.0 * PI / double(mech_Nz);

			MAX_ElasticConstants.set_to_zero();

			Norm = 1.0 / double(rlSIZE);
			// Arrays allocation:
			for (int n = 0; n < 6; n++)
			{
				rlRHSide[n] = new double[rlSIZE]();
				rcRHSide[n] = new std::complex<double>[rcSIZE]();
			}
			for (int n = 0; n < 3; n++)
			{
				rlU[n] = new double[rlSIZE]();
				rcU[n] = new std::complex<double>[rcSIZE]();
				Q[n] = new double[rcSIZE]();
			}
			for (int n = 0; n < 9; n++)
			{
				rlDefGrad[n] = new double[rlSIZE]();
				rcDefGrad[n] = new std::complex<double>[rcSIZE]();
			}

			// set Q
			for (int i = 0; i < mech_Nx; i++)
				for (int j = 0; j < mech_Ny; j++)
					for (int k = 0; k < Nz2; k++)
					{
						int XYZ = k + Nz2 * (j + mech_Ny * i);

						Q[0][XYZ] = DPi_Nx * (i * (i <= mech_Nx / 2) - (mech_Nx - i) * (i > mech_Nx / 2)) / 1.0;
						Q[1][XYZ] = DPi_Ny * (j * (j <= mech_Ny / 2) - (mech_Ny - j) * (j > mech_Ny / 2)) / 1.0;
						Q[2][XYZ] = DPi_Nz * (k * (k <= mech_Nz / 2) - (mech_Nz - k) * (k > mech_Nz / 2)) / 1.0;
					}

			for (int n = 0; n < 6; n++)
			{
				ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
				(mech_Nx, mech_Ny, mech_Nz, rlRHSide[n],
					reinterpret_cast<fftw_complex*>(rcRHSide[n]),
					FFTW_PATIENT);
			}

			for (int n = 0; n < 9; n++)
			{
				BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
				(mech_Nx, mech_Ny, mech_Nz,
					reinterpret_cast<fftw_complex*>(rcDefGrad[n]),
					rlDefGrad[n],
					FFTW_PATIENT);
			}

			for (int n = 0; n < 3; n++)
			{
				BackwardPlanU[n] = fftw_plan_dft_c2r_3d
				(mech_Nx, mech_Ny, mech_Nz,
					reinterpret_cast<fftw_complex*>(rcU[n]),
					rlU[n],
					FFTW_PATIENT);
			}
		}
		void define_funcs_for_mechanics(std::vector<bool> _avgStrainMask, std::vector<bool> _loadStressMask, std::vector<bool> _appStrainMask,
			vStress _applied_stress, vStrain _applied_strain,
			Matrix6x6(*_cal_stiffness)(int, int, int) = mechanical_boundary_condition_funcs::cal_stiffness,
			vStrain(*_cal_eigenstrain)(int, int, int) = mechanical_boundary_condition_funcs::cal_eigenstrain) {
			AvgStrainMask = _avgStrainMask;
			LoadStressMask = _loadStressMask;
			AppStrainMask = _appStrainMask;
			applied_stress = _applied_stress;
			applied_strain = _applied_strain;
			cal_stiffness = _cal_stiffness;
			cal_eigenstrain = _cal_eigenstrain;
		}
		void free() {
			if(elastic_field)
				elastic_field->clear();
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

		void cal_eigenstrain_stiffness();
		void recal_eigenstrain();
		void SetMAXElasticConstants(std::vector<Matrix6x6> Cijs);

		double get_u_main_node(size_t _x, size_t _y, size_t _z);
		double get_v_main_node(size_t _x, size_t _y, size_t _z);
		double get_w_main_node(size_t _x, size_t _y, size_t _z);

		// Ingo Steinbach Method
		void initStrainIncrements();
		std::string Solve(double StrainAccuracy, int MAXIterations, double incre_rate = 1.0, bool is_dvStraindt_output = false, bool getU = false);

		// Armen G. Khachaturyan Method
		void initVirtualEigenstrain();
		std::string Solve2(double StrainAccuracy, int MAXIterations, double iterate_rate, bool is_dvStraindt_output = false, bool getU = false);

		//----------------------------------------------------------------- settings
		vStress applied_stress;
		vStrain applied_strain;

		std::vector<bool>     AvgStrainMask;
		std::vector<bool>     LoadStressMask;
		std::vector<bool>     AppStrainMask;

		Matrix6x6(*cal_stiffness)(int x, int y, int z);
		vStrain(*cal_eigenstrain)(int x, int y, int z);
		//----------------------------------------------------------------- settings

	private:
		// Ingo Steinbach
		void CalculateRHS(Matrix6x6 Cij);
		void ExecuteForwardFFT();
		void CalculateFourierSolution(Matrix6x6 Cij, bool getU);
		void ExecuteBackwardFFT(bool getU);
		void SetElasticProperties1(double& MAXStrainDifference, vStress& AverageStress);
		void SetElasticProperties2(double& MAXStrainDifference);
		void SetElasticBoundaryConditions(vStress TargetStress);
		// Armen G. Khachaturyan
		void CalculateRHS2(Matrix6x6& Cij);
		void evaluate_virtualEigenstrain(Matrix6x6 Cij, Matrix6x6 Sij, double& MAXvStrainDifference, double iterate_rate);
		void SetElasticBoundaryConditions2(Matrix6x6& Cij, Matrix6x6& Sij);
		int Nx;
		int Ny;
		int Nz;
		int mech_Nx;
		int mech_Ny;
		int mech_Nz;
		vStrain average_strain;
		vStrain average_virtual_strain;
		Mesh<ElasticPoint>* elastic_field;
		//Elasticity Tensors:
		Matrix6x6   average_stiffness;
		Matrix6x6   average_compliences;
		//Matrix6x6   MINCompliences;

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
