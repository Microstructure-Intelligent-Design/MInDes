#include "GoverningEquationsF.h"
using namespace std;

namespace pf {
	void LBM::do_collision() {
#pragma omp parallel for
		for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
			for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
				for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++)
					_collision(x, y, z);
		if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
#pragma omp parallel for
			for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
				for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
					for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
						LBMPoint& point = main_field::lbm_field(x, y, z);
						for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_8; index++)
							point.F[index] += point.M[index];
					}
		}
		else if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D3Q19) {
#pragma omp parallel for
			for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
				for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
					for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
						LBMPoint& point = main_field::lbm_field(x, y, z);
						for (size_t index = LBM_LATTICE_ELEMENT::LBM_0; index <= LBM_LATTICE_ELEMENT::LBM_18; index++)
							point.F[index] += point.M[index];
					}
		}
	}

	void LBM::do_streaming() {
		if (lbm_lattice_model == LBM_LATTICE_MODEL::LBM_D2Q9) {
			long long z = 0;
#pragma omp parallel sections// OMP BEGIN
			{
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_1] = 
								main_field::lbm_field(x - 1, y, z).F[LBM_LATTICE_ELEMENT::LBM_1];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_2] = 
								main_field::lbm_field(x, y - 1, z).F[LBM_LATTICE_ELEMENT::LBM_2];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_3] = 
								main_field::lbm_field(x + 1, y, z).F[LBM_LATTICE_ELEMENT::LBM_3];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_4] = 
								main_field::lbm_field(x, y + 1, z).F[LBM_LATTICE_ELEMENT::LBM_4];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_5] = 
								main_field::lbm_field(x - 1, y - 1, z).F[LBM_LATTICE_ELEMENT::LBM_5];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_6] = 
								main_field::lbm_field(x + 1, y - 1, z).F[LBM_LATTICE_ELEMENT::LBM_6];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_7] = 
								main_field::lbm_field(x + 1, y + 1, z).F[LBM_LATTICE_ELEMENT::LBM_7];
						}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++) {
							main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_8] = 
								main_field::lbm_field(x - 1, y + 1, z).F[LBM_LATTICE_ELEMENT::LBM_8];
						}
				}
			}//OMP END
		}
		else {
#pragma omp parallel sections// OMP BEGIN
			{
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_1] = 
									main_field::lbm_field(x - 1, y, z).F[LBM_LATTICE_ELEMENT::LBM_1];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_2] = 
									main_field::lbm_field(x + 1, y, z).F[LBM_LATTICE_ELEMENT::LBM_2];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_3] = 
									main_field::lbm_field(x, y - 1, z).F[LBM_LATTICE_ELEMENT::LBM_3];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_4] = 
									main_field::lbm_field(x, y + 1, z).F[LBM_LATTICE_ELEMENT::LBM_4];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_END(); z >= main_field::lbm_field.COMP_Z_BGN(); z--) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_5] = 
									main_field::lbm_field(x, y, z - 1).F[LBM_LATTICE_ELEMENT::LBM_5];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_6] = 
									main_field::lbm_field(x, y, z + 1).F[LBM_LATTICE_ELEMENT::LBM_6];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_7] = 
									main_field::lbm_field(x - 1, y - 1, z).F[LBM_LATTICE_ELEMENT::LBM_7];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_8] = 
									main_field::lbm_field(x + 1, y - 1, z).F[LBM_LATTICE_ELEMENT::LBM_8];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_9] = 
									main_field::lbm_field(x + 1, y + 1, z).F[LBM_LATTICE_ELEMENT::LBM_9];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_10] = 
									main_field::lbm_field(x - 1, y + 1, z).F[LBM_LATTICE_ELEMENT::LBM_10];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_END(); z >= main_field::lbm_field.COMP_Z_BGN(); z--) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_11] = 
									main_field::lbm_field(x - 1, y, z - 1).F[LBM_LATTICE_ELEMENT::LBM_11];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_END(); z >= main_field::lbm_field.COMP_Z_BGN(); z--) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_12] = 
									main_field::lbm_field(x + 1, y, z - 1).F[LBM_LATTICE_ELEMENT::LBM_12];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_13] = 
									main_field::lbm_field(x + 1, y, z + 1).F[LBM_LATTICE_ELEMENT::LBM_13];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_END(); x >= main_field::lbm_field.COMP_X_BGN(); x--)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_14] = 
									main_field::lbm_field(x - 1, y, z + 1).F[LBM_LATTICE_ELEMENT::LBM_14];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--)
							for (long long z = main_field::lbm_field.COMP_Z_END(); z >= main_field::lbm_field.COMP_Z_BGN(); z--) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_15] = 
									main_field::lbm_field(x, y - 1, z - 1).F[LBM_LATTICE_ELEMENT::LBM_15];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_END(); z >= main_field::lbm_field.COMP_Z_BGN(); z--) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_16] = 
									main_field::lbm_field(x, y + 1, z - 1).F[LBM_LATTICE_ELEMENT::LBM_16];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_17] = 
									main_field::lbm_field(x, y + 1, z + 1).F[LBM_LATTICE_ELEMENT::LBM_17];
							}
				}
#pragma omp section
				{
					for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
						for (long long y = main_field::lbm_field.COMP_Y_END(); y >= main_field::lbm_field.COMP_Y_BGN(); y--)
							for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
								main_field::lbm_field(x, y, z).F[LBM_LATTICE_ELEMENT::LBM_18] = 
									main_field::lbm_field(x, y - 1, z + 1).F[LBM_LATTICE_ELEMENT::LBM_18];
							}
				}
			}//OMP END
		}
	}

	void LBM::do_boundary_condition() {
#pragma omp parallel for
		for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
			for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
				for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++)
					_boundary_condition(x, y, z);
	}

	MACRO_MAX_VARIATION LBM::cal_macro_variables() {
		MACRO_MAX_VARIATION MAX_VARIATION;
#pragma omp parallel for
		for (long long x = main_field::lbm_field.COMP_X_BGN(); x <= main_field::lbm_field.COMP_X_END(); x++)
			for (long long y = main_field::lbm_field.COMP_Y_BGN(); y <= main_field::lbm_field.COMP_Y_END(); y++)
				for (long long z = main_field::lbm_field.COMP_Z_BGN(); z <= main_field::lbm_field.COMP_Z_END(); z++) {
					LBMPoint& point = main_field::lbm_field(x, y, z);
					Vector3 old_momuntum = point.FV_MACRO;
					double old_f_macro = point.F_MACRO;
					_cal_macro_variables(x, y, z);
					Vector3 MOMENTUM_VARIATION = old_momuntum - point.FV_MACRO;
					double VARIATION = old_f_macro - point.F_MACRO;
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (MAX_VARIATION.FV_MACRO_MAX_VARIATION[0] < abs(MOMENTUM_VARIATION[0]))
							MAX_VARIATION.FV_MACRO_MAX_VARIATION[0] = abs(MOMENTUM_VARIATION[0]);
						if (MAX_VARIATION.FV_MACRO_MAX_VARIATION[1] < abs(MOMENTUM_VARIATION[1]))
							MAX_VARIATION.FV_MACRO_MAX_VARIATION[1] = abs(MOMENTUM_VARIATION[1]);
						if (MAX_VARIATION.FV_MACRO_MAX_VARIATION[2] < abs(MOMENTUM_VARIATION[2]))
							MAX_VARIATION.FV_MACRO_MAX_VARIATION[2] = abs(MOMENTUM_VARIATION[2]);
						if (MAX_VARIATION.F_MACRO_MAX_VARIATION < abs(VARIATION))
							MAX_VARIATION.F_MACRO_MAX_VARIATION = abs(VARIATION);
					}
				}
		return MAX_VARIATION;
	}

}
