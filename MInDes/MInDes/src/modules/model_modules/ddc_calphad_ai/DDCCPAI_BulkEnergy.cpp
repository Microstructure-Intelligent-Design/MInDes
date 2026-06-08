#include "DDCCPAI_BulkEnergy.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace chemical_energy_functions {
			REAL fchem_polynomial(std::vector<REAL> con, size_t phi_property) {
				size_t region_index = ConRegions::instance().phi_property_region(phi_property), con_index = 0;
				REAL fchem = 0, term = 0;
				for (size_t term_index = 0; term_index < parameters::terms_number[phi_property]; term_index++) {
					term = parameters::params(phi_property, term_index);
					for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
						con_index = ConRegions::instance().region_con(region_index, index);
						term *= std::pow(con[con_index], parameters::orders(phi_property, term_index, con_index));
					}
					fchem += term;
				}
				return fchem;
			}

			REAL miu_polynomial(std::vector<REAL> con, size_t phi_property, size_t _con_index) {
				size_t region_index = ConRegions::instance().phi_property_region(phi_property), con_index = 0;
				REAL fchem = 0, term = 0;
				for (size_t term_index = 0; term_index < parameters::terms_number[phi_property]; term_index++) {
					if (parameters::orders(phi_property, term_index, _con_index) > 0)
						term = 0;
					else {
						term = parameters::params(phi_property, term_index);
						for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
							con_index = ConRegions::instance().region_con(region_index, index);
							if (con_index == _con_index)
								term *= std::pow(con[con_index], parameters::orders(phi_property, term_index, con_index) - 1) * parameters::orders(phi_property, term_index, con_index);
							else
								term *= std::pow(con[con_index], parameters::orders(phi_property, term_index, con_index));
						}
					}
					fchem += term;
				}
				return fchem;
			}

			REAL delt_Fchem_delt_phi_polynomial(size_t x, size_t y, size_t z, size_t phi_index) {
				size_t phi_property = PhiProperties::instance()[phi_index], con_index = 0,
					region_index = ConRegions::instance().phi_property_region(phi_property);
				FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
				REAL dFdphi = fchem_polynomial(phase_con.phase_con[phi_property], phi_property);
				for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
					con_index = ConRegions::instance().region_con(region_index, index);
					dFdphi -= phase_con.phase_con[phi_property][con_index] * phase_con.phase_miu[phi_property][con_index];
				}
				return dFdphi;
			}

			REAL delt_Fchem_delt_con_polynomial(size_t x, size_t y, size_t z, size_t region_index, size_t _con_index) {
				REAL miu = 0, sum_phi = 0;
				FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
				for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(region_index); index++) {
					size_t phi_property = ConRegions::instance().region_phi_property(region_index, index);
					if (phase_con.phase_phi[phi_property] > parameters::PhiCon_Cut_Off) {
						sum_phi += phase_con.phase_phi[phi_property];
						miu += phase_con.phase_phi[phi_property] * phase_con.phase_miu[phi_property][_con_index];
					}
				}
				if (sum_phi > parameters::PhiCon_Cut_Off)
					return miu / sum_phi;
				else
					return 0;
			}

			void init_phase_concentration() {
				parameters::PhaseCon_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, mesh_parameters::delt_r, mesh_parameters::x_down,
					mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
				for (long long x = parameters::PhaseCon_field.COMP_X_BGN(); x <= parameters::PhaseCon_field.COMP_X_END(); x++)
					for (long long y = parameters::PhaseCon_field.COMP_Y_BGN(); y <= parameters::PhaseCon_field.COMP_Y_END(); y++)
						for (long long z = parameters::PhaseCon_field.COMP_Z_BGN(); z <= parameters::PhaseCon_field.COMP_Z_END(); z++) {
							FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
							phase_con.init_con(PhiProperties::instance().phi_property_number(), main_field::con_number);
						}
			}

			std::pair<REAL, size_t> energy_minimazation(FIELD_PhaseCon& phase_con, size_t region_index, std::vector<size_t> active_phase,
				size_t active_phase_number, REAL active_phi) {
				REAL variation = 0, dxdt = 0;
				std::vector<REAL> miu_ave(main_field::con_number, 0);
				size_t istep = 0, index = 0, index2 = 0, phi_property = 0, con_index = 0, con_number = ConRegions::instance().region_con_number(region_index);
				for (istep = 1; istep <= parameters::max_phasecon_variation_step; istep++) {
					variation = 0;
					// - cal diffusion potential
					for (index2 = 0; index2 < con_number; index2++) {
						con_index = ConRegions::instance().region_con(region_index, index2);
						miu_ave[con_index] = 0;
					}
					for (index = 0; index < active_phase_number; index++) {
						phi_property = active_phase[index];
						for (index2 = 0; index2 < con_number; index2++) {
							con_index = ConRegions::instance().region_con(region_index, index2);
							phase_con.phase_miu[phi_property][con_index] = parameters::miu(phase_con.phase_con[phi_property], phi_property, con_index);
							miu_ave[con_index] += phase_con.phase_miu[phi_property][con_index] / active_phase_number;
						}
					}
					// - do variation
					for (index = 0; index < active_phase_number; index++) {
						phi_property = active_phase[index];
						for (index2 = 0; index2 < con_number; index2++) {
							con_index = ConRegions::instance().region_con(region_index, index2);
							dxdt = -parameters::L_phasecon / phase_con.phase_phi[phi_property] * active_phi 
								* (phase_con.phase_miu[phi_property][con_index] - miu_ave[con_index]);
							phase_con.phase_con[phi_property][con_index] += dxdt;
							if (abs(dxdt) > variation)
								variation = abs(dxdt);
						}
					}
					if(variation < parameters::max_phasecon_variation)
						return { variation , istep };
				}
				return { variation , istep };
			}

			void calculation_energy_minimazation_pre() {
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
							FIELD_PhiTemp& phi_temp = parameters::PhiTemp_field(x, y, z);
							FIELD_Con& conf = parameters::Con_field(x, y, z);
							for (size_t index = 0; index < PhiProperties::instance().phi_property_number(); index++)
								phase_con.phase_phi[index] = 0;
							for (size_t index = 0; index < main_field::phi_number; index++)
								phase_con.phase_phi[PhiProperties::instance()[index]] += phi_temp.new_phi[index];
							for (size_t phi_property = 0; phi_property < PhiProperties::instance().phi_property_number(); phi_property++) {
								phase_con.phase_old_phi[phi_property] = phase_con.phase_phi[phi_property];
								size_t region_index = ConRegions::instance().phi_property_region(phi_property);
								if (parameters::is_energy_minimization[region_index]) {
									for (size_t con_index = 0; con_index < main_field::con_number; con_index++) {
										phase_con.phase_con[phi_property][con_index] = 0;
										phase_con.phase_miu[phi_property][con_index] = 0;
									}
									if (phase_con.phase_phi[phi_property] < parameters::PhiCon_Cut_Off)
										continue;
									for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
										size_t con_index = ConRegions::instance().region_con(region_index, index);
										phase_con.phase_con[phi_property][con_index] = conf.new_con[con_index];
									}
								}
							}
							for (size_t region_index = 0; region_index < ConRegions::instance().region_number(); region_index++) {
								if (parameters::is_energy_minimization[region_index]) {

								}
							}
						}
			}

		}
	}
}