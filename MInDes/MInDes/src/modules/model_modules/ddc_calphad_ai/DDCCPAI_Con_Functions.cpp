#include "DDCCPAI_Con_Functions.h"
#include <cmath>
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace concentration_field_functions {
			void init_con_in_moving_region_default(size_t x, size_t y, size_t z, size_t region_index) {
				std::vector<FIELD_Con*> point_box;
				REAL sum_rg_phi = 0;
				for(int relate_x = -1; relate_x <= 1; relate_x++)
					for (int relate_y = -1; relate_y <= 1; relate_y++)
						for (int relate_z = -1; relate_z <= 1; relate_z++) {
							if (relate_x == 0 && relate_y == 0 && relate_z == 0)
								continue;
							FIELD_Con* field_var_ave = &parameters::Con_field(x + relate_x, y + relate_y, z + relate_z);
							if (field_var_ave->new_region[region_index] > parameters::PhiCon_Cut_Off
								&& field_var_ave->old_region[region_index] > parameters::PhiCon_Cut_Off) {
								point_box.push_back(field_var_ave);
								sum_rg_phi += field_var_ave->new_region[region_index];
							}
						}
				size_t rg_con_number = ConRegions::instance().region_con_number(region_index);
				if (sum_rg_phi < parameters::PhiCon_Cut_Off) {
					Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
					for (size_t index = 0; index < rg_con_number; index++) {
						size_t con_index = ConRegions::instance().region_con(region_index, index);
						con[con_index] = parameters::PhiCon_Cut_Off_R / rg_con_number;
					}
				}
				else {
					Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
					for (size_t index = 0; index < rg_con_number; index++) {
						size_t con_index = ConRegions::instance().region_con(region_index, index);
						con[con_index] = 0;
						for (auto point = point_box.begin(); point < point_box.end(); point++)
							con[con_index] += (*point)->new_region[region_index] / sum_rg_phi * (*point)->new_con[con_index];
					}
				}
			}

			void deinit_con_in_moving_region_default(size_t x, size_t y, size_t z, size_t region_index) {
				Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
				size_t rg_con_number = ConRegions::instance().region_con_number(region_index);
				for (size_t index = 0; index < rg_con_number; index++) {
					size_t con_index = ConRegions::instance().region_con(region_index, index);
					con[con_index] = 0;
				}
			}

			void cal_mob_miu_grad_lap_7P(size_t x, size_t y, size_t z, size_t region_index,
				Vector3& grad_region, std::vector<Vector3>& grad_miu, std::vector<Vector3>& grad_mob, std::vector<REAL>& lap_miu) {
				FIELD_Con* point = &parameters::Con_field(x, y, z);
				FIELD_Con* point_upx = &parameters::Con_field(x + 1, y, z);
				FIELD_Con* point_downx = &parameters::Con_field(x - 1, y, z);
				FIELD_Con* point_upy = &parameters::Con_field(x, y + 1, z);
				FIELD_Con* point_downy = &parameters::Con_field(x, y - 1, z);
				FIELD_Con* point_upz = &parameters::Con_field(x, y, z + 1);
				FIELD_Con* point_downz = &parameters::Con_field(x, y, z - 1);
				if (point_upx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upx = point;
				if (point_downx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downx = point;
				if (point_upy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upy = point;
				if (point_downy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downy = point;
				if (point_upz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upz = point;
				if (point_downz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downz = point;
				grad_region[0] = (point_upx->new_region[region_index] - point_downx->new_region[region_index]) / 2 / mesh_parameters::delt_r;
				grad_region[1] = (point_upy->new_region[region_index] - point_downy->new_region[region_index]) / 2 / mesh_parameters::delt_r;
				grad_region[2] = (point_upz->new_region[region_index] - point_downz->new_region[region_index]) / 2 / mesh_parameters::delt_r;
				size_t rg_con_number = ConRegions::instance().region_con_number(region_index);
				for (size_t index = 0; index < rg_con_number; index++) {
					size_t con_index = ConRegions::instance().region_con(region_index, index);
					// - 
					grad_miu[con_index][0] = (point_upx->miu_con[con_index] - point_downx->miu_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_miu[con_index][1] = (point_upy->miu_con[con_index] - point_downy->miu_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_miu[con_index][2] = (point_upz->miu_con[con_index] - point_downz->miu_con[con_index]) / 2 / mesh_parameters::delt_r;
					lap_miu[con_index] =
						(point_upx->miu_con[con_index] + point_downx->miu_con[con_index]
							+ point_upy->miu_con[con_index] + point_downy->miu_con[con_index]
							+ point_upz->miu_con[con_index] + point_downz->miu_con[con_index] - 6 * point->miu_con[con_index])
						/ mesh_parameters::delt_r / mesh_parameters::delt_r;
					// - 
					grad_mob[con_index][0] = (point_upx->mob_con[con_index] - point_downx->mob_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_mob[con_index][1] = (point_upy->mob_con[con_index] - point_downy->mob_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_mob[con_index][2] = (point_upz->mob_con[con_index] - point_downz->mob_con[con_index]) / 2 / mesh_parameters::delt_r;
				}
			}

			void cal_mob_miu_grad_lap_19P(size_t x, size_t y, size_t z, size_t region_index,
				Vector3& grad_region, std::vector<Vector3>& grad_miu, std::vector<Vector3>& grad_mob, std::vector<REAL>& lap_miu) {
				FIELD_Con* point = &parameters::Con_field(x, y, z);
				FIELD_Con* point_upx = &parameters::Con_field(x + 1, y, z);
				FIELD_Con* point_downx = &parameters::Con_field(x - 1, y, z);
				FIELD_Con* point_upy = &parameters::Con_field(x, y + 1, z);
				FIELD_Con* point_downy = &parameters::Con_field(x, y - 1, z);
				FIELD_Con* point_upz = &parameters::Con_field(x, y, z + 1);
				FIELD_Con* point_downz = &parameters::Con_field(x, y, z - 1);
				if (point_upx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upx = point;
				if (point_downx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downx = point;
				if (point_upy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upy = point;
				if (point_downy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downy = point;
				if (point_upz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upz = point;
				if (point_downz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downz = point;
				FIELD_Con* point_upxupy = &parameters::Con_field(x + 1, y + 1, z);
				FIELD_Con* point_downxdowny = &parameters::Con_field(x - 1, y - 1, z);
				FIELD_Con* point_upydownx = &parameters::Con_field(x - 1, y + 1, z);
				FIELD_Con* point_downyupx = &parameters::Con_field(x + 1, y - 1, z);
				FIELD_Con* point_upxupz = &parameters::Con_field(x + 1, y, z + 1);
				FIELD_Con* point_downxdownz = &parameters::Con_field(x - 1, y, z - 1);
				FIELD_Con* point_upzdownx = &parameters::Con_field(x - 1, y, z + 1);
				FIELD_Con* point_downzupx = &parameters::Con_field(x + 1, y, z - 1);
				FIELD_Con* point_upzupy = &parameters::Con_field(x, y + 1, z + 1);
				FIELD_Con* point_downzdowny = &parameters::Con_field(x, y - 1, z - 1);
				FIELD_Con* point_upydownz = &parameters::Con_field(x, y + 1, z - 1);
				FIELD_Con* point_downyupz = &parameters::Con_field(x, y - 1, z + 1);
				if (point_upxupy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upxupy = point;
				if (point_downxdowny->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downxdowny = point;
				if (point_upydownx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upydownx = point;
				if (point_downyupx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downyupx = point;
				if (point_upxupz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upxupz = point;
				if (point_downxdownz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downxdownz = point;
				if (point_upzdownx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upzdownx = point;
				if (point_downzupx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downzupx = point;
				if (point_upzupy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upzupy = point;
				if (point_downzdowny->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downzdowny = point;
				if (point_upydownz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upydownz = point;
				if (point_downyupz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downyupz = point;
				grad_region[0] = (point_upx->new_region[region_index] - point_downx->new_region[region_index]) / 2 / mesh_parameters::delt_r;
				grad_region[1] = (point_upy->new_region[region_index] - point_downy->new_region[region_index]) / 2 / mesh_parameters::delt_r;
				grad_region[2] = (point_upz->new_region[region_index] - point_downz->new_region[region_index]) / 2 / mesh_parameters::delt_r;
				size_t rg_con_number = ConRegions::instance().region_con_number(region_index);
				for (size_t index = 0; index < rg_con_number; index++) {
					size_t con_index = ConRegions::instance().region_con(region_index, index);
					// - 
					grad_miu[con_index][0] = (point_upx->miu_con[con_index] - point_downx->miu_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_miu[con_index][1] = (point_upy->miu_con[con_index] - point_downy->miu_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_miu[con_index][2] = (point_upz->miu_con[con_index] - point_downz->miu_con[con_index]) / 2 / mesh_parameters::delt_r;
					lap_miu[con_index] =
						(4 * (point_upx->miu_con[con_index] + point_downx->miu_con[con_index]
							+ point_upy->miu_con[con_index] + point_downy->miu_con[con_index]
							+ point_upz->miu_con[con_index] + point_downz->miu_con[con_index])
							+ point_upxupy->miu_con[con_index] + point_downxdowny->miu_con[con_index]
							+ point_upydownx->miu_con[con_index] + point_downyupx->miu_con[con_index]
							+ point_upxupz->miu_con[con_index] + point_downxdownz->miu_con[con_index]
							+ point_upzdownx->miu_con[con_index] + point_downzupx->miu_con[con_index]
							+ point_upzupy->miu_con[con_index] + point_downzdowny->miu_con[con_index]
							+ point_upydownz->miu_con[con_index] + point_downyupz->miu_con[con_index]
							- 36 * point->miu_con[con_index]) / 6 / mesh_parameters::delt_r / mesh_parameters::delt_r;
					// - 
					grad_mob[con_index][0] = (point_upx->mob_con[con_index] - point_downx->mob_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_mob[con_index][1] = (point_upy->mob_con[con_index] - point_downy->mob_con[con_index]) / 2 / mesh_parameters::delt_r;
					grad_mob[con_index][2] = (point_upz->mob_con[con_index] - point_downz->mob_con[con_index]) / 2 / mesh_parameters::delt_r;
				}
			}

			REAL dfdcon_con(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				return parameters::Con_field(x, y, z).old_con[con_index];
			}

			REAL dfdcon_lap_con_7P(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				FIELD_Con* point = &parameters::Con_field(x, y, z);
				FIELD_Con* point_upx = &parameters::Con_field(x + 1, y, z);
				FIELD_Con* point_downx = &parameters::Con_field(x - 1, y, z);
				FIELD_Con* point_upy = &parameters::Con_field(x, y + 1, z);
				FIELD_Con* point_downy = &parameters::Con_field(x, y - 1, z);
				FIELD_Con* point_upz = &parameters::Con_field(x, y, z + 1);
				FIELD_Con* point_downz = &parameters::Con_field(x, y, z - 1);
				if (point_upx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upx = point;
				if (point_downx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downx = point;
				if (point_upy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upy = point;
				if (point_downy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downy = point;
				if (point_upz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upz = point;
				if (point_downz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downz = point;
				REAL lap_con = (point_upx->old_con[con_index] + point_downx->old_con[con_index]
					+ point_upy->old_con[con_index] + point_downy->old_con[con_index]
					+ point_upz->old_con[con_index] + point_downz->old_con[con_index] - 6 * point->old_con[con_index])
					/ mesh_parameters::delt_r / mesh_parameters::delt_r;
				return -parameters::Ki[con_index] * lap_con;
			}

			REAL dfdcon_lap_con_19P(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				FIELD_Con* point = &parameters::Con_field(x, y, z);
				FIELD_Con* point_upx = &parameters::Con_field(x + 1, y, z);
				FIELD_Con* point_downx = &parameters::Con_field(x - 1, y, z);
				FIELD_Con* point_upy = &parameters::Con_field(x, y + 1, z);
				FIELD_Con* point_downy = &parameters::Con_field(x, y - 1, z);
				FIELD_Con* point_upz = &parameters::Con_field(x, y, z + 1);
				FIELD_Con* point_downz = &parameters::Con_field(x, y, z - 1);
				if (point_upx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upx = point;
				if (point_downx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downx = point;
				if (point_upy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upy = point;
				if (point_downy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downy = point;
				if (point_upz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upz = point;
				if (point_downz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downz = point;
				FIELD_Con* point_upxupy = &parameters::Con_field(x + 1, y + 1, z);
				FIELD_Con* point_downxdowny = &parameters::Con_field(x - 1, y - 1, z);
				FIELD_Con* point_upydownx = &parameters::Con_field(x - 1, y + 1, z);
				FIELD_Con* point_downyupx = &parameters::Con_field(x + 1, y - 1, z);
				FIELD_Con* point_upxupz = &parameters::Con_field(x + 1, y, z + 1);
				FIELD_Con* point_downxdownz = &parameters::Con_field(x - 1, y, z - 1);
				FIELD_Con* point_upzdownx = &parameters::Con_field(x - 1, y, z + 1);
				FIELD_Con* point_downzupx = &parameters::Con_field(x + 1, y, z - 1);
				FIELD_Con* point_upzupy = &parameters::Con_field(x, y + 1, z + 1);
				FIELD_Con* point_downzdowny = &parameters::Con_field(x, y - 1, z - 1);
				FIELD_Con* point_upydownz = &parameters::Con_field(x, y + 1, z - 1);
				FIELD_Con* point_downyupz = &parameters::Con_field(x, y - 1, z + 1);
				if (point_upxupy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upxupy = point;
				if (point_downxdowny->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downxdowny = point;
				if (point_upydownx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upydownx = point;
				if (point_downyupx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downyupx = point;
				if (point_upxupz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upxupz = point;
				if (point_downxdownz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downxdownz = point;
				if (point_upzdownx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upzdownx = point;
				if (point_downzupx->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downzupx = point;
				if (point_upzupy->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upzupy = point;
				if (point_downzdowny->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downzdowny = point;
				if (point_upydownz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_upydownz = point;
				if (point_downyupz->new_region[region_index] <= parameters::PhiCon_Cut_Off)
					point_downyupz = point;
				REAL lap_con = (4 * (point_upx->old_con[con_index] + point_downx->old_con[con_index]
					+ point_upy->old_con[con_index] + point_downy->old_con[con_index]
					+ point_upz->old_con[con_index] + point_downz->old_con[con_index])
					+ point_upxupy->old_con[con_index] + point_downxdowny->old_con[con_index]
					+ point_upydownx->old_con[con_index] + point_downyupx->old_con[con_index]
					+ point_upxupz->old_con[con_index] + point_downxdownz->old_con[con_index]
					+ point_upzdownx->old_con[con_index] + point_downzupx->old_con[con_index]
					+ point_upzupy->old_con[con_index] + point_downzdowny->old_con[con_index]
					+ point_upydownz->old_con[con_index] + point_downyupz->old_con[con_index]
					- 36 * point->old_con[con_index]) / 6 / mesh_parameters::delt_r / mesh_parameters::delt_r;
				return -parameters::Ki[con_index] * lap_con;
			}

			REAL surface_reaction_flux(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				REAL surface_reaction_flux = 0;
				for (auto reaction : parameters::surface_reaction)
					surface_reaction_flux += reaction(x, y, z, region_index, con_index);
				return surface_reaction_flux;
			}

			REAL bulk_reaction_flux(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				REAL bulk_reaction_flux = 0;
				for (auto reaction : parameters::bulk_reaction)
					bulk_reaction_flux += reaction(x, y, z, region_index, con_index);
				return bulk_reaction_flux;
			}

			REAL bulk_diffusion_mobility(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
				REAL Mii = 0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_var.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_index = field_var.active_index[index];
					if (ConRegions::instance().phi_region(phi_index) == region_index)
						Mii += field_var.new_phi[phi_index] * parameters::Mii(PhiProperties::instance()[phi_index], con_index);
				}
				return Mii;
			}

			REAL bulk_diffusion_mobility_with_temperature(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
				REAL Mii = 0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_var.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_index = field_var.active_index[index];
					if (ConRegions::instance().phi_region(phi_index) == region_index)
						Mii += field_var.new_phi[phi_index] 
						* (parameters::Mii(PhiProperties::instance()[phi_index], con_index) 
							* std::exp(-parameters::Qii(PhiProperties::instance()[phi_index], con_index) / parameters::R / field_var.new_temp));
				}
				return Mii;
			}

			REAL surface_diffusion_mobility(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				REAL region_phi = parameters::Con_field(x, y, z).new_region[region_index];
				return region_phi * (1 - region_phi) * parameters::Mii_surf(region_index, con_index);
			}

			REAL interphase_diffusion_mobility(size_t x, size_t y, size_t z, size_t region_index, size_t con_index) {
				FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
				REAL Mii = 0;
				for (size_t index1 = 0; index1 < parameters::PHI_ACC_NUMBER && field_var.active_index[index1] != parameters::PAIRWISE_ACC_STOP; index1++) {
					size_t alpha_index = field_var.active_index[index1];
					if (ConRegions::instance().phi_region(alpha_index) == region_index)
						for (size_t index2 = index1 + 1; index2 < parameters::PHI_ACC_NUMBER && field_var.active_index[index2] != parameters::PAIRWISE_ACC_STOP; index2++) {
							size_t beta_index = field_var.active_index[index2];
							if (ConRegions::instance().phi_region(beta_index) == region_index)
								Mii += field_var.new_phi[alpha_index] * field_var.new_phi[beta_index] 
								* parameters::Mii_grain(PhiProperties::instance()[alpha_index], PhiProperties::instance()[beta_index], con_index);
						}
				}
				return Mii;
			}

			//========================================================================================================================
			// - evolution equations
			void init_total_concentration() {
				for (long long x = 0; x < main_field::concentration_field.Nx(); x++)
					for (long long y = 0; y < main_field::concentration_field.Ny(); y++)
						for (long long z = 0; z < main_field::concentration_field.Nz(); z++) {
							Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
							FIELD_PhiTemp& field_phi = parameters::PhiTemp_field(x, y, z);
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							field_con.init_con(main_field::con_number, ConRegions::instance().region_number());
							for (size_t index = 0; index < main_field::phi_number; index++) {
								// - calculate after init phase field
								if ((field_phi.intflag[index] || field_phi.new_phi[index] > parameters::PhiCon_Cut_Off))
									field_con.new_region[ConRegions::instance().phi_region(index)] += field_phi.new_phi[index];
							}
							for (size_t index = 0; index < main_field::con_number; index++) {
								if (field_con.new_region[ConRegions::instance().con_region(index)] < parameters::PhiCon_Cut_Off) {
									con[index] = 0;
								}
								else {
									field_con.old_con[index] = con[index];
									field_con.new_con[index] = con[index];
								}
							}
						}
				chemical_energy_functions::init_phase_concentration();
				chemical_energy_functions::calculation_energy_minimazation_pre();
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							// - 
							for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
								if (field_con.new_region[index] > parameters::PhiCon_Cut_Off) {
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
										size_t con_index = ConRegions::instance().region_con(index, index2);
										// - 
										// - diffusion driving force calculation
										for (auto dfdcon : parameters::delt_Fbulk_delt_con)
											field_con.miu_con[con_index] += dfdcon(x, y, z, index, con_index);
										// - mobility calculation
										for (auto mobi : parameters::mobility)
											field_con.mob_con[con_index] += mobi(x, y, z, index, con_index);
									}
								}
							}
						}
				parameters::PhaseCon_field.init_boundary_condition();
				parameters::Con_field.init_boundary_condition();
			}

			void pre_calculation_total_concentration() {
				// - region
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_phi = parameters::PhiTemp_field(x, y, z);
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
								field_con.old_region[index] = field_con.new_region[index];
								field_con.new_region[index] = 0;
							}
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_phi.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
								size_t phi_index = field_phi.active_index[index];
								field_con.new_region[ConRegions::instance().phi_region(phi_index)] += field_phi.new_phi[phi_index];
							}
						}
				// - moving region method
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
								if (field_con.new_region[index] > parameters::PhiCon_Cut_Off
									&& field_con.old_region[index] < parameters::PhiCon_Cut_Off)
									parameters::init_con_in_moving_region(x, y, z, index);
								else if(field_con.new_region[index] < parameters::PhiCon_Cut_Off
									&& field_con.old_region[index] > parameters::PhiCon_Cut_Off)
									parameters::deinit_con_in_moving_region(x, y, z, index);
							}
						}
				// - pre_calculation
				parameters::MAX_MINIMIZATION_RESULTS = { 0, 0 };
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							// - local concentration redistribution
							std::pair<REAL, size_t> con_results = { 0, 0 };
							if (parameters::local_concentration_redistribution)
								con_results = parameters::local_concentration_redistribution(x, y, z);
							// - 
							for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
								if (field_con.new_region[index] > parameters::PhiCon_Cut_Off) {
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
										size_t con_index = ConRegions::instance().region_con(index, index2);
										// - 
										field_con.old_con[con_index] = con[con_index];
										field_con.new_con[con_index] = 0;
										field_con.mob_con[con_index] = 0;
										field_con.miu_con[con_index] = 0;
										// - diffusion driving force calculation
										for (auto dfdcon : parameters::delt_Fbulk_delt_con)
											field_con.miu_con[con_index] += dfdcon(x, y, z, index, con_index);
										// - mobility calculation
										for (auto mobi : parameters::mobility)
											field_con.mob_con[con_index] += mobi(x, y, z, index, con_index);
									}
								}
							}
#ifdef _OPENMP
#pragma omp critical
#endif
							{
								if (con_results.first > parameters::MAX_MINIMIZATION_RESULTS.first)
									parameters::MAX_MINIMIZATION_RESULTS.first = con_results.first;
								if (con_results.second > parameters::MAX_MINIMIZATION_RESULTS.second)
									parameters::MAX_MINIMIZATION_RESULTS.second = con_results.second;
							}
						}
				parameters::PhaseCon_field.do_boundary_condition();
				parameters::Con_field.do_boundary_condition();
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
								if (field_con.new_region[index] > parameters::PhiCon_Cut_Off) {
									Vector3 grad_region(0, 0, 0);
									std::vector<Vector3> grad_miu(main_field::con_number, Vector3(0, 0, 0));
									std::vector<Vector3> grad_mob(main_field::con_number, Vector3(0, 0, 0));
									std::vector<REAL> lap_miu(main_field::con_number, 0);
									parameters::cal_mob_miu_grad_lap(x, y, z, index, grad_region, grad_miu, grad_mob, lap_miu);
									REAL abs_grad_region = grad_region.abs();
									// diffusion flux
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
										size_t con_index = ConRegions::instance().region_con(index, index2);
										REAL diffusion_flux = grad_mob[con_index] * grad_miu[con_index] + field_con.mob_con[con_index] * lap_miu[con_index];
#ifdef _DEBUG
										if (std::isnan(diffusion_flux)) {
											std::cout << "DEBUG: Diffusion flux on position x = " << x << " ,y = " << y << " , z = " << z << 
												" , con index = " << con_index << " is NaN !" << std::endl;
											SYS_PROGRAM_STOP;
										}
#endif
										// surface reaction flux
										REAL reaction_flux1 = abs_grad_region * surface_reaction_flux(x, y, z, index, con_index);
#ifdef _DEBUG
										if (std::isnan(reaction_flux1)) {
											std::cout << "DEBUG: Surface reaction flux on position x = " << x << " ,y = " << y << " , z = " << z <<
												" , con index = " << con_index << " is NaN !" << std::endl;
											SYS_PROGRAM_STOP;
										}
#endif
										// bulk reaction flux
										REAL reaction_flux2 = field_con.new_region[index] * bulk_reaction_flux(x, y, z, index, con_index);
#ifdef _DEBUG
										if (std::isnan(reaction_flux2)) {
											std::cout << "DEBUG: Bulk reaction flux on position x = " << x << " ,y = " << y << " , z = " << z <<
												" , con index = " << con_index << " is NaN !" << std::endl;
											SYS_PROGRAM_STOP;
										}
#endif
										// moving surface term
										REAL moving_surface_flux = -field_con.old_con[con_index]
											* (field_con.new_region[index] - field_con.old_region[index]) / time_parameters::delt_t;
#ifdef _DEBUG
										if (std::isnan(moving_surface_flux)) {
											std::cout << "DEBUG: Moving surface flux on position x = " << x << " ,y = " << y << " , z = " << z <<
												" , con index = " << con_index << " is NaN !" << std::endl;
											SYS_PROGRAM_STOP;
										}
#endif
										field_con.new_con[con_index] += (diffusion_flux + reaction_flux1 + reaction_flux2 + moving_surface_flux)
											/ field_con.new_region[index];
									}
								}
							}
						}
			}

			REAL solve_total_concentration() {
				REAL MAX_VARIATION = 0;
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_Con& field_con = parameters::Con_field(x, y, z);
							Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
							// assignment
							for (size_t index = 0; index < ConRegions::instance().region_number(); index++) {
								if (field_con.new_region[index] > parameters::PhiCon_Cut_Off) {
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
										size_t con_index = ConRegions::instance().region_con(index, index2);
										field_con.new_con[con_index] = field_con.old_con[con_index] +
											time_parameters::delt_t * field_con.new_con[con_index];
										if (parameters::is_con_normalized) {
											if (field_con.new_con[con_index] < SYS_EPSILON)
												field_con.new_con[con_index] = SYS_EPSILON;
											else if (field_con.new_con[con_index] > SYS_EPSILON_R)
												field_con.new_con[con_index] = SYS_EPSILON_R;
										}
									}
									REAL sum_con = 0;
									if (parameters::is_con_normalized) {
										for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
											size_t con_index = ConRegions::instance().region_con(index, index2);
											sum_con += field_con.new_con[con_index];
										}
										if(sum_con > SYS_EPSILON_R)
											for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
												size_t con_index = ConRegions::instance().region_con(index, index2);
												field_con.new_con[con_index] *= SYS_EPSILON_R / sum_con;
											}
									}
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(index); index2++) {
										size_t con_index = ConRegions::instance().region_con(index, index2);
										con[con_index] = field_con.new_con[con_index];
#ifdef _OPENMP
#pragma omp critical
#endif
										{
											if (abs(field_con.new_con[con_index] - field_con.old_con[con_index]) > MAX_VARIATION) // MAX_X_INCREMENT
												MAX_VARIATION = abs(field_con.new_con[con_index] - field_con.old_con[con_index]);
										}
									}
								}
							}
						}
				main_field::concentration_field.do_boundary_condition();
				return MAX_VARIATION;
			}

			void write_scalar_con_all(std::ofstream& fout) {
				for (size_t cindex = 0; cindex < main_field::con_number; cindex++) {
					size_t region_index = ConRegions::instance().con_region(cindex);
					std::string con_name = "ddc_con_" + std::to_string(cindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << con_name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								FIELD_Con& field_con = parameters::Con_field(i, j, k);
								Matrix1D<REAL>& con = main_field::concentration_field(i, j, k);
								fout << con[cindex] * field_con.new_region[region_index] << std::endl;
							}
					fout << "</DataArray>" << std::endl;
				}
			}

			//========================================================================================================================
			// - output

		}
	}
}
