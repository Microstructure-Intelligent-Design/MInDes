#pragma once
#include "../../base/MACRO_DEF.h"
#include "../../Modules_Params.h"
namespace pf {
	struct MACRO_MAX_VARIATION {
		Vector3 FV_MACRO_MAX_VARIATION;
		double F_MACRO_MAX_VARIATION;
		MACRO_MAX_VARIATION() {
			F_MACRO_MAX_VARIATION = 0.0;
		}
		void operator=(const MACRO_MAX_VARIATION& n) {
			FV_MACRO_MAX_VARIATION = n.FV_MACRO_MAX_VARIATION;
			F_MACRO_MAX_VARIATION = n.F_MACRO_MAX_VARIATION;
		}
	};
	namespace lbm_funcs {
		static void collision(long long x, long long y, long long z) {
			return;
		}
		static void cal_macro_variables(long long x, long long y, long long z) {
			return;
		}
		static void boundary_condition(long long x, long long y, long long z) {
			return;
		}
	}
	class LBM
	{
	public:
		LBM() {};
		~LBM() {
			free();
		}
		void init(Mesh_Boundry<LBMPoint>& _lbm_field, size_t Nx, size_t Ny, size_t Nz, REAL dr, BoundaryCondition x_down, BoundaryCondition x_up
			, BoundaryCondition y_down, BoundaryCondition y_up, BoundaryCondition z_down, BoundaryCondition z_up) {
			lbm_field = &_lbm_field; 
			lbm_field->init(Nx, Ny, Nz, dr, x_down, x_up, y_down, y_up, z_down, z_up);
			lbm_lattice_model = LBM_LATTICE_MODEL::LBM_D3Q19;
			if (mesh_parameters::MESH_NX <= 1 || mesh_parameters::MESH_NY <= 1 || mesh_parameters::MESH_NZ <= 1)
				lbm_lattice_model = LBM_LATTICE_MODEL::LBM_D2Q9;
			for (long long z = 0; z < lbm_field->Nz(); z++)
				for (long long y = 0; y < lbm_field->Ny(); y++)
					for (long long x = 0; x < lbm_field->Nx(); x++)
						lbm_field->at(x, y, z).init(lbm_lattice_model);
			_collision = lbm_funcs::collision;
			_boundary_condition = lbm_funcs::boundary_condition;
			_cal_macro_variables = lbm_funcs::cal_macro_variables;
		}

		void free() {
			lbm_field = nullptr;
			_collision = nullptr;
			_boundary_condition = nullptr;
			_cal_macro_variables = nullptr;
		}

		void do_collision();

		void do_streaming();

		void do_boundary_condition();

		MACRO_MAX_VARIATION cal_macro_variables();

		void (*_collision)(long long x, long long y, long long z);

		void (*_boundary_condition)(long long x, long long y, long long z);

		void (*_cal_macro_variables)(long long x, long long y, long long z);

		Mesh_Boundry<LBMPoint>* lbm_field;

		LBM_LATTICE_MODEL lbm_lattice_model;
	};

}