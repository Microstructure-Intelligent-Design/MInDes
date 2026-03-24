#pragma once
#include "../../base/VectorMatrix.h"
namespace pf {
	enum LBM_LATTICE_MODEL { LBM_D2Q9, LBM_D3Q19 };
	enum LBM_LATTICE_ELEMENT {
		LBM_0, LBM_1, LBM_2, LBM_3, LBM_4, LBM_5, LBM_6, LBM_7, LBM_8, LBM_9, 
		LBM_10, LBM_11, LBM_12, LBM_13, LBM_14, LBM_15, LBM_16, LBM_17, LBM_18
	};
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
	*/
	const std::vector<Vector3> LBM_d2q9_w_vec = {
		Vector3(0.0, 0.0, 0.0) ,	 // f0  c( 0,  0)
		Vector3(1.0, 0.0, 0.0),	     // f1  c( 1,  0)
		Vector3(0.0, 1.0, 0.0),	     // f2  c( 0,  1)
		Vector3(-1.0, 0.0, 0.0),	 // f3  c(-1,  0)
		Vector3(0.0, -1.0, 0.0),	 // f4  c( 0, -1)
		Vector3(1.0, 1.0, 0.0),	     // f5  c( 1,  1)
		Vector3(-1.0, 1.0, 0.0),	 // f6  c(-1,  1)
		Vector3(-1.0, -1.0, 0.0),	 // f7  c(-1, -1)
		Vector3(1.0, -1.0, 0.0)	     // f8  c( 1, -1)
	};
	const std::vector<REAL> LBM_d2q9_w = {
		4.0 / 9.0,
		1.0 / 9.0,
		1.0 / 9.0,
		1.0 / 9.0,
		1.0 / 9.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0
	};
	/*
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
	const std::vector<Vector3> LBM_d3q19_w_vec = {
		Vector3(0.0, 0.0, 0.0),	     // f0   c( 0,  0,  0)
		Vector3(1.0, 0.0, 0.0),	     // f1   c( 1,  0,  0)
		Vector3(-1.0, 0.0, 0.0),	 // f2   c(-1,  0,  0)
		Vector3(0.0, 1.0, 0.0),	     // f3   c( 0,  1,  0)
		Vector3(0.0, -1.0, 0.0),	 // f4   c( 0, -1,  0)
		Vector3(0.0, 0.0, 1.0),	     // f5   c( 0,  0,  1)
		Vector3(0.0, 0.0, -1.0),	 // f6   c( 0,  0, -1)
		Vector3(1.0, 1.0, 0.0),	     // f7   c( 1,  1,  0)
		Vector3(-1.0, 1.0, 0.0),	 // f8   c(-1,  1,  0)
		Vector3(-1.0, -1.0, 0.0),    // f9   c(-1, -1,  0)
		Vector3(1.0, -1.0, 0.0),	 // f10  c( 1, -1,  0)
		Vector3(1.0, 0.0, 1.0),	     // f11  c( 1,  0,  1)
		Vector3(-1.0, 0.0, 1.0),	 // f12  c(-1,  0,  1)
		Vector3(-1.0, 0.0, -1.0),    // f13  c(-1,  0, -1)
		Vector3(1.0, 0.0, -1.0),	 // f14  c( 1,  0, -1)
		Vector3(0.0, 1.0, 1.0),	     // f15  c( 0,  1,  1)
		Vector3(0.0, -1.0, 1.0),	 // f16  c( 0, -1,  1)
		Vector3(0.0, -1.0, -1.0),    // f17  c( 0, -1, -1)
		Vector3(0.0, 1.0, -1.0) 	 // f18  c( 0,  1, -1)
	};
	const std::vector<REAL> LBM_d3q19_w = {
		12.0 / 36.0,
		2.0 / 36.0,
		2.0 / 36.0,
		2.0 / 36.0,
		2.0 / 36.0,
		2.0 / 36.0,
		2.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0,
		1.0 / 36.0
	};
	const REAL LBM_Cs2 = 1.0 / 3.0;
	const REAL LBM_Cs4 = 1.0 / 9.0;

	class LBMPoint {
	public:
		std::vector<REAL>       F;                                        ///< Storage for distribution functions F
		std::vector<REAL>       M;                                        ///< Storage for collision M
		REAL F_MACRO;
		Vector3 FV_MACRO;
		REAL pressure;
		REAL mass;
		Vector3 velocity;
		REAL fluid_region;             // fluid region: 1 , solid region : 0
		void operator=(const LBMPoint& n) {
			F = n.F;
			M = n.M;
			F_MACRO = n.F_MACRO;
			FV_MACRO = n.FV_MACRO;
			pressure = n.pressure;
			mass = n.mass;
			velocity = n.velocity;
			fluid_region = n.fluid_region;
		}
		LBMPoint() {
			F_MACRO = 0;
			pressure = 0;
			mass = 0;
			fluid_region = 0;
		};
		void init(LBM_LATTICE_MODEL model) {
			if (model == LBM_LATTICE_MODEL::LBM_D2Q9) {
				F.resize(9, 0);
				M.resize(9, 0);
			}
			else if (model == LBM_LATTICE_MODEL::LBM_D3Q19) {
				F.resize(19, 0);
				M.resize(19, 0);
			}
			F_MACRO = 0;
			FV_MACRO = Vector3(0, 0, 0);
			pressure = 0;
			mass = 0;
			velocity = Vector3(0, 0, 0);
			fluid_region = 0;
		}
	};
}