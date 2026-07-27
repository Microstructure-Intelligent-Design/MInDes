#include "DDCCPAI_BulkEnergy.h"
#include <cmath>
#include <functional>
#include <iomanip>
#include <locale>
namespace pf {
	namespace ddc_calphad_ai_model {
		// - default functions
		namespace chemical_energy_functions {
			REAL fchem_polynomial(std::vector<REAL> con, REAL temperature, size_t phi_property) {
				size_t region_index = ConRegions::instance().phi_property_region(phi_property), con_index = 0;
				REAL fchem = 0, term = 0;
				for (size_t term_index = 0; term_index < parameters::terms_number[phi_property]; term_index++) {
					term = 0;
					for(size_t index = 0; index < parameters::temp_terms_number[phi_property][term_index]; index++)
						term += parameters::params[phi_property](term_index, index) * std::pow(temperature, parameters::temp_orders[phi_property](term_index, index));
					for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
						con_index = ConRegions::instance().region_con(region_index, index);
						term *= std::pow(con[con_index], parameters::con_orders[phi_property](term_index, con_index));
					}
					fchem += term;
				}
				return fchem;
			}
			
			REAL miu_polynomial(std::vector<REAL> con, REAL temperature, size_t phi_property, size_t _con_index) {
				size_t region_index = ConRegions::instance().phi_property_region(phi_property), con_index = 0;
				REAL fchem = 0, term = 0;
				for (size_t term_index = 0; term_index < parameters::terms_number[phi_property]; term_index++) {
					term = 0;
					if (parameters::con_orders[phi_property](term_index, _con_index) != 0) {
						for (size_t index = 0; index < parameters::temp_terms_number[phi_property][term_index]; index++)
							term += parameters::params[phi_property](term_index, index) * std::pow(temperature, parameters::temp_orders[phi_property](term_index, index));
						for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
							con_index = ConRegions::instance().region_con(region_index, index);
							if (con_index == _con_index)
								term *= std::pow(con[con_index], parameters::con_orders[phi_property](term_index, con_index) - 1) * parameters::con_orders[phi_property](term_index, con_index);
							else
								term *= std::pow(con[con_index], parameters::con_orders[phi_property](term_index, con_index));
						}
					}
					fchem += term;
				}
				return fchem;
			}

			REAL delt_Fchem_delt_phi(size_t x, size_t y, size_t z, size_t phi_index) {
				size_t phi_property = PhiProperties::instance()[phi_index], con_index = 0,
					region_index = ConRegions::instance().phi_property_region(phi_property);
				FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
				REAL& temp = parameters::PhiTemp_field(x, y, z).new_temp;
				REAL dFdphi = parameters::fchem(phase_con.phase_con[phi_property], temp, phi_property);
				for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
					con_index = ConRegions::instance().region_con(region_index, index);
					dFdphi -= phase_con.phase_con[phi_property][con_index] * phase_con.phase_miu[phi_property][con_index];
				}
				return dFdphi;
			}

			REAL delt_Fchem_delt_con(size_t x, size_t y, size_t z, size_t region_index, size_t _con_index) {
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
				for (long long x = 0; x < parameters::PhaseCon_field.Nx(); x++)
					for (long long y = 0; y < parameters::PhaseCon_field.Ny(); y++)
						for (long long z = 0; z < parameters::PhaseCon_field.Nz(); z++) {
							FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
							phase_con.init_con(PhiProperties::instance().phi_property_number(), main_field::con_number);
						}
			}

			std::pair<REAL, size_t> energy_minimazation(std::vector<REAL>& phase_phi, std::vector<std::vector<REAL>>& phase_con, std::vector<std::vector<REAL>>& phase_miu, std::vector<size_t> active_phase,
				size_t active_phase_number, REAL active_phi, REAL temperature, size_t region_index) {
				REAL max_variation = 0, dxdt = 0;
				std::vector<REAL> miu_ave(main_field::con_number, 0);
				size_t istep = 0, index = 0, index2 = 0, phi_property = 0, con_index = 0, con_number = ConRegions::instance().region_con_number(region_index);
				for (istep = 1; istep <= parameters::max_phasecon_variation_step; istep++) {
					REAL variation = 0;
					// - cal diffusion potential
					for (index2 = 0; index2 < con_number; index2++) {
						con_index = ConRegions::instance().region_con(region_index, index2);
						miu_ave[con_index] = 0;
					}
					for (index = 0; index < active_phase_number; index++) {
						phi_property = active_phase[index];
						for (index2 = 0; index2 < con_number; index2++) {
							con_index = ConRegions::instance().region_con(region_index, index2);
							phase_miu[phi_property][con_index] = parameters::miu(phase_con[phi_property], temperature, phi_property, con_index);
							miu_ave[con_index] += phase_miu[phi_property][con_index] / active_phase_number;
						}
					}
					// - do variation
					for (index = 0; index < active_phase_number; index++) {
						phi_property = active_phase[index];
						for (index2 = 0; index2 < con_number; index2++) {
							con_index = ConRegions::instance().region_con(region_index, index2);
							dxdt = -parameters::L_phasecon / phase_phi[phi_property] * active_phi
								* (phase_miu[phi_property][con_index] - miu_ave[con_index]);
							phase_con[phi_property][con_index] += dxdt;
							if (abs(dxdt) > variation)
								variation = abs(dxdt);
						}
					}
					if (variation > max_variation)
						max_variation = variation;
					if (variation < parameters::phasecon_epsilon)
						return { max_variation , istep };
				}
				return { max_variation , istep };
			}

			std::pair<REAL, size_t> local_concentration_redistribution(size_t x, size_t y, size_t z) {
				// { MAX_VARIATION , MAX_ITERATION_STEP }
				std::pair<REAL, size_t> max_results = { 0, 0 };
				FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
				FIELD_PhiTemp& phi_temp = parameters::PhiTemp_field(x, y, z);
				FIELD_Con& conf = parameters::Con_field(x, y, z);
				for (size_t phi_property = 0; phi_property < PhiProperties::instance().phi_property_number(); phi_property++) {
					phase_con.phase_old_phi[phi_property] = phase_con.phase_phi[phi_property];
					phase_con.phase_phi[phi_property] = 0;
				}
				for (size_t index = 0; index < main_field::phi_number; index++)
					phase_con.phase_phi[PhiProperties::instance()[index]] += phi_temp.new_phi[index];
				for (size_t region_index = 0; region_index < ConRegions::instance().region_number(); region_index++) {
					if (parameters::is_energy_minimization[region_index]) {
						std::vector<size_t> active_phase(PhiProperties::instance().phi_property_number(), 0);
						std::vector<REAL> sum_con(main_field::con_number, 0);
						size_t active_phase_number = 0;
						REAL active_phi = 0;
						for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(region_index); index++) {
							size_t phi_property = ConRegions::instance().region_phi_property(region_index, index);
							if (phase_con.phase_phi[phi_property] > parameters::PhiCon_Cut_Off) {
								active_phase[active_phase_number] = phi_property;
								active_phase_number++;
								active_phi += phase_con.phase_phi[phi_property];
								if (phase_con.phase_old_phi[phi_property] <= parameters::PhiCon_Cut_Off) {
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
										size_t con_index = ConRegions::instance().region_con(region_index, index2);
										phase_con.phase_con[phi_property][con_index] = conf.new_con[con_index];
									}
								}
								for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
									size_t con_index = ConRegions::instance().region_con(region_index, index2);
									sum_con[con_index] += phase_con.phase_con[phi_property][con_index] * phase_con.phase_phi[phi_property];
								}
							}
							else {
								if (phase_con.phase_old_phi[phi_property] > parameters::PhiCon_Cut_Off) {
									for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
										size_t con_index = ConRegions::instance().region_con(region_index, index2);
										phase_con.phase_con[phi_property][con_index] = 0;
										phase_con.phase_miu[phi_property][con_index] = 0;
									}
								}
							}
						}
						if (active_phase_number == 1) {
							size_t phi_property = active_phase[0];
							for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
								size_t con_index = ConRegions::instance().region_con(region_index, index);
								phase_con.phase_con[phi_property][con_index] = conf.new_con[con_index];
							}
							for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
								size_t con_index = ConRegions::instance().region_con(region_index, index);
								phase_con.phase_miu[phi_property][con_index] =
									parameters::miu(phase_con.phase_con[phi_property], phi_temp.new_temp, phi_property, con_index);
							}
						}
						else {
							for (size_t index = 0; index < active_phase_number; index++) {
								size_t phi_property = active_phase[index];
								for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
									size_t con_index = ConRegions::instance().region_con(region_index, index2);
									phase_con.phase_con[phi_property][con_index] *= conf.new_con[con_index] * conf.new_region[region_index] / sum_con[con_index];
								}
							}
							// { VARIATION , ITERATION_STEP } for each region
							std::pair<REAL, size_t> results =
								energy_minimazation(phase_con.phase_phi, phase_con.phase_con, phase_con.phase_miu, active_phase, active_phase_number, active_phi, phi_temp.new_temp, region_index);
							if (results.first > max_results.first)
								max_results.first = results.first;
							if (results.second > max_results.second)
								max_results.second = results.second;
						}
					}
				}
				return max_results;
			}

			void calculation_energy_minimazation_pre() {
				WriteLog("> Do chemical energy minimazation : ");
#pragma omp parallel for
				for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
					for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
						for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
							FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(x, y, z);
							FIELD_PhiTemp& phi_temp = parameters::PhiTemp_field(x, y, z);
							Matrix1D<REAL>& conf = main_field::concentration_field(x, y, z);
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
										phase_con.phase_con[phi_property][con_index] = conf[con_index];
									}
								}
							}
							for (size_t region_index = 0; region_index < ConRegions::instance().region_number(); region_index++) {
								if (parameters::is_energy_minimization[region_index]) {
									std::pair<REAL, size_t> max_results = { 0, 0 };
									std::vector<size_t> active_phase(PhiProperties::instance().phi_property_number(), false);
									size_t active_phase_number = 0;
									REAL active_phi = 0;
									for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(region_index); index++) {
										size_t phi_property = ConRegions::instance().region_phi_property(region_index, index);
										if (phase_con.phase_phi[phi_property] > parameters::PhiCon_Cut_Off) {
											active_phase[active_phase_number] = phi_property;
											active_phase_number++;
											active_phi += phase_con.phase_phi[phi_property];
											for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
												size_t con_index = ConRegions::instance().region_con(region_index, index2);
												phase_con.phase_con[phi_property][con_index] = conf[con_index];
											}
										}
										else {
											for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
												size_t con_index = ConRegions::instance().region_con(region_index, index2);
												phase_con.phase_con[phi_property][con_index] = 0;
												phase_con.phase_miu[phi_property][con_index] = 0;
											}
										}
									}
									if (active_phase_number == 1) {
										size_t phi_property = active_phase[0];
										for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
											size_t con_index = ConRegions::instance().region_con(region_index, index);
											phase_con.phase_miu[phi_property][con_index] =
												parameters::miu(phase_con.phase_con[phi_property], phi_temp.new_temp, phi_property, con_index);
										}
									}
									else if (active_phase_number >= 2) {
										std::pair<REAL, size_t> results = 
											energy_minimazation(phase_con.phase_phi, phase_con.phase_con, phase_con.phase_miu, active_phase, active_phase_number, active_phi, phi_temp.new_temp, region_index);
										if (results.first > max_results.first)
											max_results.first = results.first;
										if (results.second > max_results.second)
											max_results.second = results.second;
										stringstream log;
										log << std::endl;
										log << "     position ( " << x << " , " << y << " , " << z << " ) , active region : " << ConRegions::instance().region_name(region_index) << " , temperature : " << phi_temp.new_temp << std::endl;
										log << "     phase volume : ";
										for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(region_index); index++) {
											size_t phi_property = ConRegions::instance().region_phi_property(region_index, index);
											log << "phi_" << PhiProperties::instance().phi_property_name(phi_property) << " = " << phase_con.phase_phi[phi_property] << " , ";
										}
										log << std::endl;
										log << "     total concentration : ";
										for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
											size_t con_index = ConRegions::instance().region_con(region_index, index);
											log << "con_" << ConRegions::instance().con_name(con_index) << " = " << conf[con_index] << " , ";
										}
										log << std::endl;
										log << "     max variation = " << max_results.first << ", max iteration step = " << max_results.second << endl;
										log << "     results : " << endl;
										for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(region_index); index++) {
											size_t phi_property = ConRegions::instance().region_phi_property(region_index, index);
											if (phase_con.phase_phi[phi_property] > parameters::PhiCon_Cut_Off) {
												log << "       phi_" << PhiProperties::instance().phi_property_name(phi_property) << " = " << phase_con.phase_phi[phi_property] << " , ";
												for (size_t index2 = 0; index2 < ConRegions::instance().region_con_number(region_index); index2++) {
													size_t con_index = ConRegions::instance().region_con(region_index, index2);
													log << "phase_con_" << ConRegions::instance().con_name(con_index) << " = " << phase_con.phase_con[phi_property][con_index] << " , ";
													log << "phase_miu_" << ConRegions::instance().con_name(con_index) << " = " << phase_con.phase_miu[phi_property][con_index] << " , ";
												}
												log << std::endl;
											}
										}
										log << std::endl;

#ifdef _OPENMP
#pragma omp critical
#endif
										{
											WriteLog(log.str());
										}
									}
								}
							}
							parameters::PhaseCon_field.init_boundary_condition();
						}
			}
			// - 
			void write_scalar_phasecon(std::ofstream& fout) {
				for (size_t index = 0; index < parameters::is_write_phase_con.size(); index++) {
					size_t phi_property = parameters::is_write_phase_con[index].first, con_index = parameters::is_write_phase_con[index].second;
					std::string name = "con_" + PhiProperties::instance().phi_property_name(phi_property) + "_" + ConRegions::instance().con_name(con_index);
					fout << "<DataArray type = \"Float64\" Name = \"" << name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								FIELD_PhaseCon& field = parameters::PhaseCon_field(i, j, k);
								fout << field.phase_con[phi_property][con_index] << std::endl;
							}
					fout << "</DataArray>" << std::endl;
				}
			}
			void write_scalar_con(std::ofstream& fout) {
				for (size_t index = 0; index < parameters::is_write_con.size(); index++) {
					size_t con_index = parameters::is_write_con[index];
					std::string name = "con_" + ConRegions::instance().con_name(con_index);
					fout << "<DataArray type = \"Float64\" Name = \"" << name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								Matrix1D<REAL>& con = main_field::concentration_field(i, j, k);
								fout << con[con_index] << std::endl;
							}
					fout << "</DataArray>" << std::endl;
				}
			}
			void write_scalar_phasemiu(std::ofstream& fout) {
				for (size_t index = 0; index < parameters::is_write_phase_miu.size(); index++) {
					size_t phi_property = parameters::is_write_phase_miu[index].first, con_index = parameters::is_write_phase_miu[index].second;
					std::string name = "miu_" + PhiProperties::instance().phi_property_name(phi_property) + "_" + ConRegions::instance().con_name(con_index);
					fout << "<DataArray type = \"Float64\" Name = \"" << name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								FIELD_PhaseCon& field = parameters::PhaseCon_field(i, j, k);
								fout << field.phase_miu[phi_property][con_index] << std::endl;
							}
					fout << "</DataArray>" << std::endl;
				}
			}
			void write_scalar_miu(std::ofstream& fout) {
				for (size_t index = 0; index < parameters::is_write_miu.size(); index++) {
					size_t con_index = parameters::is_write_miu[index];
					std::string name = "miu_" + ConRegions::instance().con_name(con_index);
					fout << "<DataArray type = \"Float64\" Name = \"" << name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								fout << delt_Fchem_delt_con(i, j, k, ConRegions::instance().con_region(con_index), con_index) << std::endl;
							}
					fout << "</DataArray>" << std::endl;
				}
			}
			void write_scalar_fchem(std::ofstream& fout) {
				for (size_t index = 0; index < parameters::is_write_fchem.size(); index++) {
					size_t phi_property = parameters::is_write_fchem[index];
					std::string name = "fchem_" + PhiProperties::instance().phi_property_name(phi_property);
					fout << "<DataArray type = \"Float64\" Name = \"" << name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								FIELD_PhaseCon& phase_con = parameters::PhaseCon_field(i, j, k);
								REAL& temp = parameters::PhiTemp_field(i, j, k).new_temp;
								if (phase_con.phase_phi[phi_property] > parameters::PhiCon_Cut_Off)
									fout << fchem_polynomial(phase_con.phase_con[phi_property], temp, phi_property) << std::endl;
								else
									fout << 0 << std::endl;
							}
					fout << "</DataArray>" << std::endl;
				}
			}

			void write_csv_energy_results() {
				if (parameters::thermo_calc_energy_phases.empty())
					return;
				if (!main_field::is_con_field_on) {
					WriteLog("> DDCCPAI chemical-energy CSV output skipped: concentration field is disabled.\n");
					return;
				}

				auto scan_values = [](const ThermoCalcScanRange& range) {
					std::vector<REAL> values;
					if (range.begin == range.end) {
						values.push_back(range.begin);
						return values;
					}
					const long double begin = static_cast<long double>(range.begin);
					const long double end = static_cast<long double>(range.end);
					const long double step = static_cast<long double>(range.step);
					const long double scale = std::abs(end) > 1.0L ? std::abs(end) : 1.0L;
					const long double tolerance = std::numeric_limits<REAL>::epsilon() * scale * 16.0L;
					for (size_t index = 0; ; index++) {
						long double raw_value = begin + static_cast<long double>(index) * step;
						if (raw_value > end + tolerance)
							break;
						if (std::abs(raw_value - end) <= tolerance)
							raw_value = end;
						values.push_back(static_cast<REAL>(raw_value));
					}
					return values;
				};

				std::vector<REAL> temperature_axis = scan_values(parameters::thermo_calc_energy_temp_range);
				for (size_t phi_property : parameters::thermo_calc_energy_phases) {
					const size_t region_index = ConRegions::instance().phi_property_region(phi_property);
					if (region_index >= parameters::is_energy_minimization.size() ||
						!parameters::is_energy_minimization[region_index]) {
						WriteLog("> DDCCPAI chemical-energy CSV output skipped for an invalid phase selection.\n");
						continue;
					}

					std::vector<size_t> con_indices;
					std::vector<std::vector<REAL>> con_axes;
					for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
						size_t con_index = ConRegions::instance().region_con(region_index, index);
						con_indices.push_back(con_index);
						con_axes.push_back(scan_values(parameters::thermo_calc_energy_con_ranges[con_index]));
					}

					const std::string phase_name = PhiProperties::instance().phi_property_name(phi_property);
					std::string csv_name = parameters::thermo_calc_energy_csv_name;
					size_t extension_position = csv_name.find_last_of('.');
					if (extension_position == std::string::npos)
						csv_name += "_" + phase_name;
					else
						csv_name.insert(extension_position, "_" + phase_name);
					const std::string csv_path = input_output_files_parameters::WorkingFolder_Path + dirSeparator + csv_name;
					std::ofstream fout(csv_path, std::ios::out | std::ios::trunc);
					if (!fout.is_open()) {
						WriteLog("> ERROR: Cannot create DDCCPAI chemical-energy CSV file: " + csv_path + "\n");
						SYS_PROGRAM_STOP;
					}
					fout.imbue(std::locale::classic());
					fout << std::setprecision(std::numeric_limits<REAL>::max_digits10);
					fout << "sample_id,temperature";
					for (size_t con_index : con_indices)
						fout << ",con_" << ConRegions::instance().con_name(con_index);
					fout << ",fchem_" << phase_name << '\n';

					std::vector<REAL> scanned_con(con_indices.size(), 0);
					std::vector<REAL> phase_con(main_field::con_number, 0);
					size_t sample_id = 0;
					size_t skipped_composition = 0;
					auto calculate_and_write = [&](REAL temperature) {
						REAL con_sum = 0;
						for (REAL concentration : scanned_con)
							con_sum += concentration;
						if (con_sum > REAL(1)) {
							skipped_composition++;
							return;
						}
						for (size_t index = 0; index < con_indices.size(); index++)
							phase_con[con_indices[index]] = scanned_con[index];
						sample_id++;
						fout << sample_id << ',' << temperature;
						for (REAL concentration : scanned_con)
							fout << ',' << concentration;
						fout << ',' << parameters::fchem(phase_con, temperature, phi_property) << '\n';
					};

					std::function<void(size_t)> scan_con;
					scan_con = [&](size_t axis) {
						if (axis == con_axes.size()) {
							for (REAL temperature : temperature_axis)
								calculate_and_write(temperature);
							return;
						}
						for (REAL value : con_axes[axis]) {
							scanned_con[axis] = value;
							scan_con(axis + 1);
						}
					};

					WriteLog("> Write DDCCPAI chemical-energy CSV data for phase " + phase_name + " to " + csv_path + "\n");
					scan_con(0);
					fout.close();
					WriteLog("> DDCCPAI chemical-energy CSV output for phase " + phase_name + " finished: " +
						std::to_string(sample_id) + " samples, " + std::to_string(skipped_composition) +
						" composition combinations skipped.\n");
				}
			}

			void write_csv_energy_minimization_results() {
				if (!parameters::is_write_thermo_calc_csv)
					return;
				if (!main_field::is_con_field_on) {
					WriteLog("> DDCCPAI energy-minimization CSV output skipped: concentration field is disabled.\n");
					return;
				}
				if (parameters::thermo_calc_region >= parameters::is_energy_minimization.size() ||
					!parameters::is_energy_minimization[parameters::thermo_calc_region]) {
					WriteLog("> DDCCPAI energy-minimization CSV output skipped: selected region is not enabled for ThermoCalc.\n");
					return;
				}
				if (parameters::thermo_calc_phases.empty()) {
					WriteLog("> DDCCPAI energy-minimization CSV output skipped: selected region has no phase.\n");
					return;
				}

				const size_t region_index = parameters::thermo_calc_region;
				std::vector<size_t> con_indices;
				for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
					con_indices.push_back(ConRegions::instance().region_con(region_index, index));

				auto scan_values = [](const ThermoCalcScanRange& range) {
					std::vector<REAL> values;
					if (range.begin == range.end) {
						values.push_back(range.begin);
						return values;
					}
					const long double begin = static_cast<long double>(range.begin);
					const long double end = static_cast<long double>(range.end);
					const long double step = static_cast<long double>(range.step);
					const long double scale = std::abs(end) > 1.0L ? std::abs(end) : 1.0L;
					const long double tolerance = std::numeric_limits<REAL>::epsilon() * scale * 16.0L;
					for (size_t index = 0; ; index++) {
						long double raw_value = begin + static_cast<long double>(index) * step;
						if (raw_value > end + tolerance)
							break;
						if (std::abs(raw_value - end) <= tolerance)
							raw_value = end;
						values.push_back(static_cast<REAL>(raw_value));
					}
					return values;
				};

				std::vector<std::vector<REAL>> con_axes;
				for (size_t con_index : con_indices)
					con_axes.push_back(scan_values(parameters::thermo_calc_con_ranges[con_index]));
				std::vector<std::vector<REAL>> phi_axes;
				for (size_t index = 0; index < parameters::thermo_calc_phases.size(); index++)
					phi_axes.push_back(scan_values(parameters::thermo_calc_phi_ranges[parameters::thermo_calc_phases[index]]));
				std::vector<REAL> temperature_axis = scan_values(parameters::thermo_calc_temp_range);

				const std::string csv_path = input_output_files_parameters::WorkingFolder_Path + dirSeparator + parameters::thermo_calc_equi_csv_name;
				std::ofstream fout(csv_path, std::ios::out | std::ios::trunc);
				if (!fout.is_open()) {
					WriteLog("> ERROR: Cannot create DDCCPAI energy-minimization CSV file: " + csv_path + "\n");
					SYS_PROGRAM_STOP;
				}
				fout.imbue(std::locale::classic());
				fout << std::setprecision(std::numeric_limits<REAL>::max_digits10);
				fout << "sample_id,temperature";
				for (size_t con_index : con_indices)
					fout << ",con_" << ConRegions::instance().con_name(con_index);
				for (size_t phi_property : parameters::thermo_calc_phases)
					fout << ",phi_" << PhiProperties::instance().phi_property_name(phi_property);
				for (size_t phi_property : parameters::thermo_calc_phases)
					for (size_t con_index : con_indices)
						fout << ",phase_con_" << PhiProperties::instance().phi_property_name(phi_property)
							<< "_" << ConRegions::instance().con_name(con_index);
				for (size_t phi_property : parameters::thermo_calc_phases)
					for (size_t con_index : con_indices)
						fout << ",phase_miu_" << PhiProperties::instance().phi_property_name(phi_property)
							<< "_" << ConRegions::instance().con_name(con_index);
				for (size_t phi_property : parameters::thermo_calc_phases)
					fout << ",fchem_" << PhiProperties::instance().phi_property_name(phi_property);
				fout << ",variation,iteration_step\n";

				std::vector<REAL> total_con(con_indices.size(), 0);
				std::vector<REAL> selected_phi(parameters::thermo_calc_phases.size(), 0);
				size_t sample_id = 0;
				size_t skipped_composition = 0;
				size_t skipped_phase_fraction = 0;

				auto calculate_and_write = [&](REAL temperature) {
					REAL con_sum = 0;
					for (REAL concentration : total_con)
						con_sum += concentration;
					if (!(con_sum <= REAL(1))) {
						skipped_composition++;
						return;
					}

					REAL active_phi = 0;
					for (REAL phi : selected_phi) {
						if (!(phi > parameters::PhiCon_Cut_Off)) {
							skipped_phase_fraction++;
							return;
						}
						active_phi += phi;
					}
					if (!(active_phi <= REAL(1))) {
						skipped_phase_fraction++;
						return;
					}

					std::vector<REAL> phase_phi(PhiProperties::instance().phi_property_number(), 0);
					std::vector<std::vector<REAL>> phase_con(PhiProperties::instance().phi_property_number(),
						std::vector<REAL>(main_field::con_number, 0));
					std::vector<std::vector<REAL>> phase_miu(PhiProperties::instance().phi_property_number(),
						std::vector<REAL>(main_field::con_number, 0));
					std::vector<size_t> active_phase;
					for (size_t index = 0; index < parameters::thermo_calc_phases.size(); index++) {
						size_t phi_property = parameters::thermo_calc_phases[index];
						phase_phi[phi_property] = selected_phi[index];
						active_phase.push_back(phi_property);
						for (size_t con_position = 0; con_position < con_indices.size(); con_position++)
							phase_con[phi_property][con_indices[con_position]] = total_con[con_position];
					}

					std::pair<REAL, size_t> result = { 0, 0 };
					if (active_phase.size() == 1) {
						size_t phi_property = active_phase[0];
						for (size_t con_index : con_indices)
							phase_miu[phi_property][con_index] =
								parameters::miu(phase_con[phi_property], temperature, phi_property, con_index);
					}
					else {
						result = energy_minimazation(phase_phi, phase_con, phase_miu, active_phase,
							active_phase.size(), active_phi, temperature, region_index);
					}

					sample_id++;
					fout << sample_id << ',' << temperature;
					for (REAL concentration : total_con)
						fout << ',' << concentration;
					for (REAL phi : selected_phi)
						fout << ',' << phi;
					for (size_t phi_property : parameters::thermo_calc_phases)
						for (size_t con_index : con_indices)
							fout << ',' << phase_con[phi_property][con_index];
					for (size_t phi_property : parameters::thermo_calc_phases)
						for (size_t con_index : con_indices)
							fout << ',' << phase_miu[phi_property][con_index];
					for (size_t phi_property : parameters::thermo_calc_phases) {
						if (phase_phi[phi_property] > 0)
							fout << ',' << parameters::fchem(phase_con[phi_property], temperature, phi_property);
						else
							fout << ",0";
					}
					fout << ',' << result.first << ',' << result.second << '\n';

					std::stringstream log;
					log << std::setprecision(std::numeric_limits<REAL>::max_digits10);
					log << std::endl << "     energy-minimization CSV sample " << sample_id
						<< " , active region : " << ConRegions::instance().region_name(region_index) << std::endl;
					log << "     phase volume : ";
					for (size_t index = 0; index < parameters::thermo_calc_phases.size(); index++)
						log << "phi_" << PhiProperties::instance().phi_property_name(parameters::thermo_calc_phases[index])
							<< " = " << selected_phi[index] << " , ";
					log << "phi_S = " << active_phi << " , ";
					log << std::endl << "     total concentration : ";
					for (size_t index = 0; index < con_indices.size(); index++)
						log << "con_" << ConRegions::instance().con_name(con_indices[index]) << " = " << total_con[index] << " , ";
					log << std::endl << "     temperature = " << temperature;
					log << std::endl << "     max variation = " << result.first << ", max iteration step = " << result.second << std::endl;
					log << "     results : " << std::endl;
					for (size_t phi_property : parameters::thermo_calc_phases) {
						if (phase_phi[phi_property] > 0) {
							log << "       phi_" << PhiProperties::instance().phi_property_name(phi_property)
								<< " = " << phase_phi[phi_property] << " , ";
							for (size_t con_index : con_indices) {
								log << "phase_con_" << ConRegions::instance().con_name(con_index) << " = "
									<< phase_con[phi_property][con_index] << " , ";
								log << "phase_miu_" << ConRegions::instance().con_name(con_index) << " = "
									<< phase_miu[phi_property][con_index] << " , ";
							}
							log << std::endl;
						}
					}
					WriteLog(log.str());
				};

				std::function<void(size_t)> scan_phi;
				std::function<void(size_t)> scan_con;
				scan_phi = [&](size_t axis) {
					if (axis == phi_axes.size()) {
						for (REAL temperature : temperature_axis)
							calculate_and_write(temperature);
						return;
					}
					for (REAL value : phi_axes[axis]) {
						selected_phi[axis] = value;
						scan_phi(axis + 1);
					}
				};
				scan_con = [&](size_t axis) {
					if (axis == con_axes.size()) {
						scan_phi(0);
						return;
					}
					for (REAL value : con_axes[axis]) {
						total_con[axis] = value;
						scan_con(axis + 1);
					}
				};

				WriteLog("> Write DDCCPAI energy-minimization CSV data to " + csv_path + "\n");
				scan_con(0);
				fout.close();
				WriteLog("> DDCCPAI energy-minimization CSV output finished: " + std::to_string(sample_id) +
					" samples, " + std::to_string(skipped_composition) + " composition combinations skipped, " +
					std::to_string(skipped_phase_fraction) + " phase-fraction combinations skipped.\n");
			}
		}
	}
}
