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
#include "../Base.h"

using namespace std;
namespace pf {
	namespace micro_structure_init {
		NucleationBox nucleation_box;
		static double voronoi_point_distance = -1.0;
		bool is_voronoi_rand = true;
		static int rand_seed = 0;
		static bool is_datafile_init = false;
		static pf::Data_report datafile_report;
		//-
		static bool is_init_by_datefile() {
			return is_datafile_init;
		}
		static pf::Data_report get_datafile_info() {
			return datafile_report;
		}
		//-

		static bool add_new_phi_by_index(int phi_index, int phi_property, double phi_value = 0.0) {
			//int componentFluxDrivingForce = information->settings.details_settings.flux_model;
			Info_Phase phases = Solvers::get_instance()->parameters.Phases[phi_property];
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
			for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++)
				if (phase->index == phi_index)
					return false;
			phaseMesh.info_node.add_phase(phi_index, phi_property);
			for (auto x = phases.x.begin(); x < phases.x.end(); x++) {
				phaseMesh.info_node[phi_index].x.add_con(x->index, 0.0);
				phaseMesh.info_node[phi_index].potential.add_con(x->index, 0.0);
			}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.add_phase(phi_index, phi_property, pf_BULK, phi_value);
						PhaseEntry& new_phase = node[phi_index];
						new_phase.old_phi = new_phase.phi;
						for (auto x = phases.x.begin(); x < phases.x.end(); x++) {
							new_phase.x.add_con(x->index, 0.0);
							new_phase.potential.add_con(x->index, 0.0);
						}
						for (auto comp1 = new_phase.x.begin(); comp1 < new_phase.x.end(); comp1++)
							for (auto comp2 = new_phase.x.begin(); comp2 < new_phase.x.end(); comp2++)
								new_phase.kinetics_coeff.set(comp1->index, comp2->index, 0.0);
					}
			return true;
		}
		static bool add_new_phi_by_index_with_index(int phi_index, int phi_property, int with_phi_index, double phi_value = 0.0) {
			//int componentFluxDrivingForce = information->settings.details_settings.flux_model;
			Info_Phase phases = Solvers::get_instance()->parameters.Phases[phi_property];
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
			for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++)
				if (phase->index == phi_index)
					return false;
			phaseMesh.info_node.add_phase(phi_index, phi_property);
			for (auto x = phases.x.begin(); x < phases.x.end(); x++) {
				phaseMesh.info_node[phi_index].x.add_con(x->index, 0.0);
				phaseMesh.info_node[phi_index].potential.add_con(x->index, 0.0);
			}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						bool is_phi = false;
						for (auto phase = node.begin(); phase < node.end(); phase++)
							if (phase->index == with_phi_index)
								is_phi = true;
						if (is_phi) {
							node.add_phase(phi_index, phi_property, pf_BULK, phi_value);
							PhaseEntry& new_phase = node[phi_index];
							new_phase.old_phi = new_phase.phi;
							for (auto x = phases.x.begin(); x < phases.x.end(); x++) {
								new_phase.x.add_con(x->index, 0.0);
								new_phase.potential.add_con(x->index, 0.0);
							}
							for (auto comp1 = new_phase.x.begin(); comp1 < new_phase.x.end(); comp1++)
								for (auto comp2 = new_phase.x.begin(); comp2 < new_phase.x.end(); comp2++)
									new_phase.kinetics_coeff.set(comp1->index, comp2->index, 0.0);
						}
					}
			return true;
		}
		static void matrix_phase_component(int phi_index, int phi_property, vector<input_value> comp_value) {
			Info_Node phase_x = Solvers::get_instance()->parameters.Phases[phi_property].x;
			if (phase_x.size() != comp_value.size()) {
				string error_report = "> define phase_x and their values error, size of phase component in matrix mismatch !\n";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				std::exit(0);
			}
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseEntry& phase = phaseMesh(x, y, z)[phi_index];
						int index = 0;
						for (auto x = phase_x.begin(); x < phase_x.end(); x++) {
							phase.x.add_con(x->index, comp_value[index].double_value);
							index++;
						}
					}
		}
		static void total_component(vector<input_value> comp_value) {
			Info_Node total_x = Solvers::get_instance()->parameters.Components;
			if (total_x.size() != comp_value.size()) {
				string error_report = "> define total_x and their values error, size of total component in matrix mismatch !\n";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				std::exit(0);
			}
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
			int index = 0;
			for (auto x = total_x.begin(); x < total_x.end(); x++) {
				phaseMesh.info_node.x.add_con(x->index, comp_value[index].double_value);
				phaseMesh.info_node.potential.add_con(x->index, 0.0);
				index++;
			}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						index = 0;
						for (auto x = total_x.begin(); x < total_x.end(); x++) {
							node.x.add_con(x->index, comp_value[index].double_value);
							node.potential.add_con(x->index, 0.0);
							index++;
						}
						for (auto comp1 = total_x.begin(); comp1 < total_x.end(); comp1++)
							for (auto comp2 = total_x.begin(); comp2 < total_x.end(); comp2++)
								node.kinetics_coeff.set(comp1->index, comp2->index, 0.0);
					}
		}
		static void temperature(double temperature) {
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.temperature.T = temperature;
					}
		}
		static void free_nonexistence_phi(FieldStorage_forPhaseNode& mesh) {
			PhaseNode& inf_node = Solvers::get_instance()->statistics_information_in_phaseMesh();
			for (auto p = inf_node.begin(); p < inf_node.end(); p++)
				if (p->phi < SYS_EPSILON) {
					for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++)
						node->erase((*node)[p->index]);
				}
		}
		static void definiteNucleation() {
			bool _creatNewStorage = true;
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
			for (auto nucleus = nucleation_box.nucleus_box.begin(); nucleus < nucleation_box.nucleus_box.end();) {
				if (nucleus->generate_step == Solvers::get_instance()->current_istep && nucleus->phasefraction < (1.0 + SYS_EPSILON) && nucleus->phasefraction > Phi_Num_Cut_Off) {
					bool checkk = false;

					_creatNewStorage = true;
					for (auto phase = phaseMesh(0, 0, 0).begin(); phase < phaseMesh(0, 0, 0).end(); phase++)
						if (phase->index == nucleus->phaseIndex) {
							if (phase->property != nucleus->phaseProperty) {
								string error_report = "> init microstructure error, the phase name of new generate structure should be the same as the old one as they are in one grain (phi index) !\n";
								InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
								std::exit(0);
							}
							_creatNewStorage = false;
						}
					if (_creatNewStorage)
						add_new_phi_by_index(nucleus->phaseIndex, nucleus->phaseProperty, 0.0);

					int x = double_to_int(nucleus->core.x), y = double_to_int(nucleus->core.y), z = double_to_int(nucleus->core.z);
					PhaseNode& node = phaseMesh(x, y, z);

					for (auto cc = node[nucleus->phaseIndex].x.begin(); cc < node[nucleus->phaseIndex].x.end(); cc++)
						cc->value = nucleus->x[cc->index].value;
					for (auto cc = node.x.begin(); cc < node.x.end(); cc++)
						for (auto xx = nucleus->x.begin(); xx < nucleus->x.end(); xx++)
							if (cc->index == xx->index)
								cc->value = xx->value;
					node[nucleus->phaseIndex].phi = nucleus->phasefraction;

					stringstream report;
					report << "> A new nucleus for grain : " << to_string(nucleus->phaseIndex) << " phase: " << Solvers::get_instance()->parameters.Phases[nucleus->phaseProperty].phi_name << " has been initialized on position_x:" << nucleus->core.x
						<< " position_y:" << nucleus->core.y << " position_z:" << nucleus->core.z << " ," << std::endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
					nucleus = nucleation_box.nucleus_box.erase(nucleus);
				}
				else if (nucleus->phasefraction > 1.0 || nucleus->phasefraction < Phi_Num_Cut_Off) {
					nucleus = nucleation_box.nucleus_box.erase(nucleus);
				}
				else {
					nucleus++;
				}
			}
			for (auto geo = nucleation_box.geometry_box.begin(); geo < nucleation_box.geometry_box.end();) {
				if (geo->generate_step == Solvers::get_instance()->current_istep) {
					_creatNewStorage = true;
					for (auto phase = phaseMesh(0, 0, 0).begin(); phase < phaseMesh(0, 0, 0).end(); phase++)
						if (phase->index == geo->phaseIndex) {
							if (phase->property != geo->phaseProperty) {
								string error_report = "> init microstructure error, the phase name of new generate structure should be the same as the old one as they are in one grain (phi index) !\n";
								InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
								std::exit(0);
							}
							_creatNewStorage = false;
						}
					if (_creatNewStorage)
						add_new_phi_by_index(geo->phaseIndex, geo->phaseProperty, 0.0);

					if (geo->isNormalized) {
						if (geo->phi > 1.0)
							geo->phi = 1.0;
						else if(geo->phi < 1.0)
							geo->phi = 0.0;
					}

					if (geo->geometryProperty == Geometry::Geo_Ellipsoid) {
						for (int z = 0; z < phaseMesh.limit_z; z++)
							for (int y = 0; y < phaseMesh.limit_y; y++)
								for (int x = 0; x < phaseMesh.limit_x; x++) {
									PhaseNode& node = phaseMesh(x, y, z);
									// Vector3 p(x, y, z);
									Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
									p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->ellipSolid.radian_x, geo->ellipSolid.radian_y, geo->ellipSolid.radian_z),
										geo->ellipSolid.rotationGauge));
									p[0] += geo->ellipSolid.core.x;
									p[1] += geo->ellipSolid.core.y;
									p[2] += geo->ellipSolid.core.z;
									Point po(p[0], p[1], p[2]);
									bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
									if (!geo->isReverseRegion && check0) {
										if (geo->isNormalized) {
											double sum_phis = 0.0;
											for (auto p = node.begin(); p < node.end(); p++)
												if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
													sum_phis += p->phi;
											if (sum_phis > SYS_EPSILON) {
												for (auto p = node.begin(); p < node.end(); p++) {
													if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
														p->phi *= (1.0 - geo->phi) / sum_phis;
													if (p->phi < Phi_Num_Cut_Off) {
														for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
															comp->value = 0.0;
													}
												}
											}
											else {
												for (auto p = node.begin(); p < node.end(); p++) {
													for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
														comp->value = 0.0;
												}
											}
										}
										node[geo->phaseIndex].phi = geo->phi;
										if (node[geo->phaseIndex].phi > Phi_Num_Cut_Off) {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = geo->x[cc->index].value;
										}
										else {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = 0.0;
										}
										for (auto cc = node.x.begin(); cc < node.x.end(); cc++)
											for (auto xx = geo->x.begin(); xx < geo->x.end(); xx++)
												if (cc->index == xx->index)
													cc->value = xx->value;
										for (auto value = geo->customValues.begin(); value < geo->customValues.end(); value++)
											node.customValues.add_double(value->index, value->value);
										for (auto flag = geo->customFlags.begin(); flag < geo->customFlags.end(); flag++)
											node.customFlags.add_int(flag->index, flag->value);
										for (auto vec = geo->customVec3s.begin(); vec < geo->customVec3s.end(); vec++)
											node.customVec3s.add_vec(vec->index, vec->vec);
										for (auto vec = geo->customVec6s.begin(); vec < geo->customVec6s.end(); vec++)
											node.customVec6s.add_vec(vec->index, vec->vec);
										node.temperature.T = geo->temperature;
									}
									else if (geo->isReverseRegion && !check0) {
										if (geo->isNormalized) {
											double sum_phis = 0.0;
											for (auto p = node.begin(); p < node.end(); p++)
												if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
													sum_phis += p->phi;
											if (sum_phis > SYS_EPSILON) {
												for (auto p = node.begin(); p < node.end(); p++) {
													if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
														p->phi *= (1.0 - geo->phi) / sum_phis;
													if (p->phi < Phi_Num_Cut_Off) {
														for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
															comp->value = 0.0;
													}
												}
											}
											else {
												for (auto p = node.begin(); p < node.end(); p++) {
													for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
														comp->value = 0.0;
												}
											}
										}
										node[geo->phaseIndex].phi = geo->phi;
										if (node[geo->phaseIndex].phi > Phi_Num_Cut_Off) {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = geo->x[cc->index].value;
										}
										else {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = 0.0;
										}
										for (auto cc = node.x.begin(); cc < node.x.end(); cc++)
											for (auto xx = geo->x.begin(); xx < geo->x.end(); xx++)
												if (cc->index == xx->index)
													cc->value = xx->value;
										for (auto value = geo->customValues.begin(); value < geo->customValues.end(); value++)
											node.customValues.add_double(value->index, value->value);
										for (auto flag = geo->customFlags.begin(); flag < geo->customFlags.end(); flag++)
											node.customFlags.add_int(flag->index, flag->value);
										for (auto vec = geo->customVec3s.begin(); vec < geo->customVec3s.end(); vec++)
											node.customVec3s.add_vec(vec->index, vec->vec);
										for (auto vec = geo->customVec6s.begin(); vec < geo->customVec6s.end(); vec++)
											node.customVec6s.add_vec(vec->index, vec->vec);
										node.temperature.T = geo->temperature;
									}
								}
								stringstream report;
								report << "> A new Ellipsoid for grain : " << to_string(geo->phaseIndex) << " phase: " << Solvers::get_instance()->parameters.Phases[geo->phaseProperty].phi_name << " has been initialized at position ( " 
									<< geo->ellipSolid.core.x << ", " << geo->ellipSolid.core.y << ", " << geo->ellipSolid.core.z << " )." << std::endl;
								Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
		
					}
					else if (geo->geometryProperty == Geometry::Geo_Polyhedron) {
						for (int z = 0; z < phaseMesh.limit_z; z++)
							for (int y = 0; y < phaseMesh.limit_y; y++)
								for (int x = 0; x < phaseMesh.limit_x; x++) {
									PhaseNode& node = phaseMesh(x, y, z);
									Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
									pv.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
										geo->polyhedron.rotationGauge));
									pv[0] += geo->polyhedron.point_inside_polyhedron.x;
									pv[1] += geo->polyhedron.point_inside_polyhedron.y;
									pv[2] += geo->polyhedron.point_inside_polyhedron.z;
									Point p(pv[0], pv[1], pv[2]);
									//Point p(x, y, z);
									p.do_boundary(phaseMesh._bc_x_up, phaseMesh._bc_y_up, phaseMesh._bc_z_up, phaseMesh._bc_x_down, phaseMesh._bc_y_down, phaseMesh._bc_z_down, 
										phaseMesh.limit_x, phaseMesh.limit_y, phaseMesh.limit_z);
									bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
									if (!geo->isReverseRegion && check0) {
										if (geo->isNormalized) {
											double sum_phis = 0.0;
											for (auto p = node.begin(); p < node.end(); p++)
												if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
													sum_phis += p->phi;
											if (sum_phis > SYS_EPSILON) {
												for (auto p = node.begin(); p < node.end(); p++) {
													if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
														p->phi *= (1.0 - geo->phi) / sum_phis;
													if (p->phi < Phi_Num_Cut_Off) {
														for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
															comp->value = 0.0;
													}
												}
											}
											else {
												for (auto p = node.begin(); p < node.end(); p++) {
													for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
														comp->value = 0.0;
												}
											}
										}
										node[geo->phaseIndex].phi = geo->phi;
										if (node[geo->phaseIndex].phi > Phi_Num_Cut_Off) {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = geo->x[cc->index].value;
										}
										else {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = 0.0;
										}
										for (auto cc = node.x.begin(); cc < node.x.end(); cc++)
											for (auto xx = geo->x.begin(); xx < geo->x.end(); xx++)
												if (cc->index == xx->index)
													cc->value = xx->value;
										for (auto value = geo->customValues.begin(); value < geo->customValues.end(); value++)
											node.customValues.add_double(value->index, value->value);
										for (auto flag = geo->customFlags.begin(); flag < geo->customFlags.end(); flag++)
											node.customFlags.add_int(flag->index, flag->value);
										for (auto vec = geo->customVec3s.begin(); vec < geo->customVec3s.end(); vec++)
											node.customVec3s.add_vec(vec->index, vec->vec);
										for (auto vec = geo->customVec6s.begin(); vec < geo->customVec6s.end(); vec++)
											node.customVec6s.add_vec(vec->index, vec->vec);
										node.temperature.T = geo->temperature;
									}
									else if (geo->isReverseRegion && !check0) {
										if (geo->isNormalized) {
											double sum_phis = 0.0;
											for (auto p = node.begin(); p < node.end(); p++)
												if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
													sum_phis += p->phi;
											if (sum_phis > SYS_EPSILON) {
												for (auto p = node.begin(); p < node.end(); p++) {
													if (p->phi > SYS_EPSILON && p->index != geo->phaseIndex)
														p->phi *= (1.0 - geo->phi) / sum_phis;
													if (p->phi < Phi_Num_Cut_Off) {
														for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
															comp->value = 0.0;
													}
												}
											}
											else {
												for (auto p = node.begin(); p < node.end(); p++) {
													for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
														comp->value = 0.0;
												}
											}
										}
										node[geo->phaseIndex].phi = geo->phi;
										if (node[geo->phaseIndex].phi > Phi_Num_Cut_Off) {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = geo->x[cc->index].value;
										}
										else {
											for (auto cc = node[geo->phaseIndex].x.begin(); cc < node[geo->phaseIndex].x.end(); cc++)
												cc->value = 0.0;
										}
										for (auto cc = node.x.begin(); cc < node.x.end(); cc++)
											for (auto xx = geo->x.begin(); xx < geo->x.end(); xx++)
												if (cc->index == xx->index)
													cc->value = xx->value;
										for (auto value = geo->customValues.begin(); value < geo->customValues.end(); value++)
											node.customValues.add_double(value->index, value->value);
										for (auto flag = geo->customFlags.begin(); flag < geo->customFlags.end(); flag++)
											node.customFlags.add_int(flag->index, flag->value);
										for (auto vec = geo->customVec3s.begin(); vec < geo->customVec3s.end(); vec++)
											node.customVec3s.add_vec(vec->index, vec->vec);
										for (auto vec = geo->customVec6s.begin(); vec < geo->customVec6s.end(); vec++)
											node.customVec6s.add_vec(vec->index, vec->vec);
										node.temperature.T = geo->temperature;
									}
								}
								stringstream report;
								report << "> A new Polygon for grain : " << to_string(geo->phaseIndex) << " phase: " << Solvers::get_instance()->parameters.Phases[geo->phaseProperty].phi_name << " has been initialized at position ( " 
									<< geo->polyhedron.point_inside_polyhedron.x << ", " << geo->polyhedron.point_inside_polyhedron.y << ", " << geo->polyhedron.point_inside_polyhedron.z << " )." << std::endl;
								Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
					}
					geo = nucleation_box.geometry_box.erase(geo);
				}
				else {
					geo++;
				}
			}
			for (auto point_set = nucleation_box.point_set_box.begin(); point_set < nucleation_box.point_set_box.end();) {
				if (point_set->generate_step == Solvers::get_instance()->current_istep) {
					_creatNewStorage = true;
					for (auto phase = phaseMesh(0, 0, 0).begin(); phase < phaseMesh(0, 0, 0).end(); phase++)
						if (phase->index == point_set->phaseIndex) {
							if (phase->property != point_set->phaseProperty) {
								string error_report = "> init microstructure error, the phase name of new generate structure should be the same as the old one as they are in one grain (phi index) !\n";
								InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
								std::exit(0);
							}
							_creatNewStorage = false;
						}
					if (_creatNewStorage)
						add_new_phi_by_index(point_set->phaseIndex, point_set->phaseProperty, 0.0);

					if (point_set->is_normalized) {
						for (auto phi = point_set->points_phi.begin(); phi < point_set->points_phi.end(); phi++)
						if (*phi > 1.0)
							*phi = 1.0;
						else if (*phi < 0.0)
							*phi = 0.0;
					}

					for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
						PhaseNode& node = phaseMesh(double_to_int(point_set->points[point_index].x), 
							double_to_int(point_set->points[point_index].y), double_to_int(point_set->points[point_index].z));
						if (point_set->is_normalized) {
							double sum_phis = 0.0;
							for (auto p = node.begin(); p < node.end(); p++)
								if (p->phi > SYS_EPSILON && p->index != point_set->phaseIndex)
									sum_phis += p->phi;
							if (sum_phis > SYS_EPSILON) {
								for (auto p = node.begin(); p < node.end(); p++) {
									if (p->phi > SYS_EPSILON && p->index != point_set->phaseIndex)
										p->phi *= (1.0 - point_set->points_phi[point_index]) / sum_phis;
									if (p->phi < Phi_Num_Cut_Off) {
										for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
											comp->value = 0.0;
									}
								}
							}
							else {
								for (auto p = node.begin(); p < node.end(); p++) {
									for (auto comp = p->x.begin(); comp < p->x.end(); comp++)
										comp->value = 0.0;
								}
							}
						}
						node[point_set->phaseIndex].phi = point_set->points_phi[point_index];
						if (node[point_set->phaseIndex].phi > Phi_Num_Cut_Off) {
							for (auto cc = node[point_set->phaseIndex].x.begin(); cc < node[point_set->phaseIndex].x.end(); cc++)
								cc->value = point_set->x[cc->index].value;
						}
						else {
							for (auto cc = node[point_set->phaseIndex].x.begin(); cc < node[point_set->phaseIndex].x.end(); cc++)
								cc->value = 0.0;
						}
						node.temperature.T = point_set->temperature;
					}
					stringstream report;
					report << "> A new PointSet for grain : " << to_string(point_set->phaseIndex) << " phase: " << Solvers::get_instance()->parameters.Phases[point_set->phaseProperty].phi_name << " has been initialized," << std::endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
					point_set = nucleation_box.point_set_box.erase(point_set);
				}
				else {
					point_set++;
				}
			}
		}

		static void generate_voronoi_structure(Vector3 box_position, Vector3 box_size, int grain_0_index, int grain_number, int generate_step,
			vector<int> phases_properties, vector<double> phases_weight, XNode x, double temperature) {
			Dimension dimention = Dimension::Three_Dimension;
			if (box_size[0] == 0 && box_size[1] == 0 && box_size[2] == 0)
				return;
			else if (box_size[0] == 0 || box_size[1] == 0 || box_size[2] == 0)
				dimention = Dimension::Two_Dimension;

			if (phases_properties.size() == 0) {
				string error_report = "> generate voronoi structure error, please set the phases' name of each grains (phi index) !\n";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				return;
			}
			else if (phases_properties.size() != phases_weight.size()) {
				phases_weight.clear();
				for (unsigned int i = 0; i < phases_properties.size(); i++)
					phases_weight.push_back(1.0 / phases_properties.size());
			}
			// > normalize weight
			double sum_weight = 0.0;
			for (auto dd = phases_weight.begin(); dd < phases_weight.end(); dd++)
				sum_weight += *dd;
			for (auto dd = phases_weight.begin(); dd < phases_weight.end(); dd++)
				*dd = *dd / sum_weight;
			// > generate points
			if (is_voronoi_rand)
				RAND_INIT_ALL;
			else
				RAND_INIT(rand_seed);
			vector<Point> points;
			int grain_index = 0;
			double distance2 = voronoi_point_distance * voronoi_point_distance;
			while (grain_index < grain_number) {
				bool is_point_add = false;
				double rand_x = RAND_0_1, rand_y = RAND_0_1, rand_z = RAND_0_1;
				Point p(rand_x * box_size[0] + box_position[0], rand_y * box_size[1] + box_position[1], rand_z * box_size[2] + box_position[2]);
				if (voronoi_point_distance > 1.0) {
					is_point_add = true;
					for (auto ip = points.begin(); ip < points.end(); ip++) {
						double D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
						if (D2 < distance2) {
							is_point_add = false;
						}
					}
				}
				else {
					is_point_add = true;
				}
				if (is_point_add) {
					points.push_back(p);
					grain_index++;
					string str_report = "> Voronoi: generate point index : " + to_string(grain_index) + ", at : (x, y, z) (" + to_string(p.x) + ", " + to_string(p.y) + ", " + to_string(p.x) + ")\n";
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(str_report, LOG_FILE_NAME);
				}
			}
			// > points index and property
			vector<int> grains_property; // = points.size()
			for (auto grain = points.begin(); grain < points.end(); grain++) {
				double rand = RAND_0_1, sum_weight2 = 0.0;
				for (unsigned int ii = 0; ii < phases_weight.size(); ii++) {
					sum_weight2 += phases_weight[ii];
					if (sum_weight2 > rand) {
						grains_property.push_back(phases_properties[ii]);
						break;
					}
				}
			}
			// > periodic boundary condition
			vector<vector<Point>> mirror_points; // = 27 * points.size()
			int region_number = 0;
			if (dimention == Dimension::Three_Dimension)
				region_number = 27;
			else
				region_number = 9;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = points[grain] + Point(0, 0, 0);
				mirror_points[1][grain] = points[grain] + Point(int(box_size[0]), 0, 0);
				mirror_points[2][grain] = points[grain] + Point(int(-box_size[0]), 0, 0);
				mirror_points[3][grain] = points[grain] + Point(0, int(box_size[1]), 0);
				mirror_points[4][grain] = points[grain] + Point(0, int(-box_size[1]), 0);
				mirror_points[5][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), 0);
				mirror_points[6][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), 0);
				mirror_points[7][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), 0);
				mirror_points[8][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), 0);
				if (dimention == Dimension::Three_Dimension) {
					mirror_points[9][grain] = points[grain] + Point(0, 0, int(box_size[2]));
					mirror_points[10][grain] = points[grain] + Point(0, 0, int(-box_size[2]));
					mirror_points[11][grain] = points[grain] + Point(int(box_size[0]), 0, int(box_size[2]));
					mirror_points[12][grain] = points[grain] + Point(int(-box_size[0]), 0, int(box_size[2]));
					mirror_points[13][grain] = points[grain] + Point(int(box_size[0]), 0, int(-box_size[2]));
					mirror_points[14][grain] = points[grain] + Point(int(-box_size[0]), 0, int(-box_size[2]));
					mirror_points[15][grain] = points[grain] + Point(0, int(box_size[1]), int(box_size[2]));
					mirror_points[16][grain] = points[grain] + Point(0, int(-box_size[1]), int(box_size[2]));
					mirror_points[17][grain] = points[grain] + Point(0, int(box_size[1]), int(-box_size[2]));
					mirror_points[18][grain] = points[grain] + Point(0, int(-box_size[1]), int(-box_size[2]));
					mirror_points[19][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), int(box_size[2]));
					mirror_points[20][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), int(box_size[2]));
					mirror_points[23][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), int(box_size[2]));
					mirror_points[24][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), int(-box_size[2]));
					mirror_points[26][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), int(-box_size[2]));
					mirror_points[21][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), int(box_size[2]));
					mirror_points[22][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), int(-box_size[2]));
					mirror_points[25][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), int(-box_size[2]));
				}
			}
#pragma omp parallel for
			for (int region = 0; region < region_number; region++)
				for (int grain = 0; grain < grain_number; grain++) {
					Polyhedron poly(mirror_points[region][grain]);
					vector<point_in_region_index> record_points;
					point_in_region_index rp(region, grain);
					record_points.push_back(rp);
					for (unsigned int region_index = 0; region_index < mirror_points.size(); region_index++)
						for (unsigned int grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
							// Avoid inclusion points
							bool is_point_contained = false;
							for (auto re = record_points.begin(); re < record_points.end(); re++)
								if (re->region == region_index && re->grain_index == grain_index)
									is_point_contained = true;
							if (is_point_contained)
								continue;
							// prepare
							Point norm, mid_point;
							norm = poly.point_inside_polyhedron - mirror_points[region_index][grain_index];
							mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
							// Judged ipsilateral to poly center
							if (poly.check_point_inside_polyhedron(mid_point) == false)
								continue;
							// < add point
							poly.add_surf(norm, mid_point);
							// < Eliminate meaningless points in poly
							for (auto re = record_points.begin(); re < record_points.end();) {
								Point check = (points[grain] + mirror_points[re->region][re->grain_index]) / 2;
								if (poly.check_point_inside_polyhedron(check) == false) {
									re = record_points.erase(re);
								}
								else {
									++re;
								}
							}
						}
					string str_report = "> Voronoi: One polyhedron in region : " + to_string(region) + ", grain : " + to_string(grain) + " has been generated ! \n";
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(str_report, LOG_FILE_NAME);
					GeometricRegion geo;
					geo.geometryProperty = Geometry::Geo_Polyhedron;
					geo.generate_step = generate_step;
					geo.polyhedron = poly;
					geo.phaseIndex = grain_0_index + grain;
					geo.phaseProperty = grains_property[grain];
					geo.temperature = temperature;
					geo.x = x;
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						nucleation_box.geometry_box.push_back(geo);
					}
				}

		}

		static void generate_voronoi_structure_in_phis(Vector3 box_position, Vector3 box_size, int grain_0_index, int grain_number, int generate_step,
			vector<int> phases_properties, vector<double> phases_weight, XNode x, double temperature, vector<int> phis_index, FieldStorage_forPhaseNode& phaseMesh) {
			Dimension dimention = Dimension::Three_Dimension;
			if (box_size[0] == 0 && box_size[1] == 0 && box_size[2] == 0)
				return;
			else if (box_size[0] == 0 || box_size[1] == 0 || box_size[2] == 0)
				dimention = Dimension::Two_Dimension;

			if (phases_properties.size() == 0) {
				string error_report = "> generate voronoi structure error, please set the phases' name of each grains (phi index) !\n";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				return;
			}
			else if (phases_properties.size() != phases_weight.size()) {
				phases_weight.clear();
				for (unsigned int i = 0; i < phases_properties.size(); i++)
					phases_weight.push_back(1.0 / phases_properties.size());
			}
			// > normalize weight
			double sum_weight = 0.0;
			for (auto dd = phases_weight.begin(); dd < phases_weight.end(); dd++)
				sum_weight += *dd;
			for (auto dd = phases_weight.begin(); dd < phases_weight.end(); dd++)
				*dd = *dd / sum_weight;
			// > generate points
			if (is_voronoi_rand)
				RAND_INIT_ALL;
			else
				RAND_INIT(rand_seed);
			vector<Point> points;
			int grain = 0;
			double distance2 = voronoi_point_distance * voronoi_point_distance;
			while (grain < grain_number)
			{
				bool is_point_add = false;
				double rand_x = RAND_0_1, rand_y = RAND_0_1, rand_z = RAND_0_1;
				Point p(rand_x * box_size[0] + box_position[0], rand_y * box_size[1] + box_position[1], rand_z * box_size[2] + box_position[2]);
				PhaseNode& node = phaseMesh(int(p.x), int(p.y), int(p.z));
				double sum_phi = 0.0;
				for (auto phase = node.begin(); phase < node.end(); phase++)
					for (auto index = phis_index.begin(); index < phis_index.end(); index++)
						if (phase->index == *index)
							sum_phi += phase->phi;
				if (voronoi_point_distance > 1.0) {
					is_point_add = true;
					for (auto ip = points.begin(); ip < points.end(); ip++) {
						double D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
						if (D2 < distance2) {
							is_point_add = false;
						}
					}
				}
				else {
					is_point_add = true;
				}
				if (is_point_add && sum_phi > SYS_EPSILON) {
					grain++;
					points.push_back(p);
					string str_report = "> Voronoi: generate point index : " + to_string(grain) + ",  at : (x, y, z) (" + to_string(p.x) + ", " + to_string(p.y) + ", " + to_string(p.x) + ")\n";
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(str_report, LOG_FILE_NAME);
				}
			}
			// > points index and property
			vector<int> grains_property; // = points.size()
			for (auto grain = points.begin(); grain < points.end(); grain++) {
				double rand = RAND_0_1, sum_weight2 = 0.0;
				for (unsigned int ii = 0; ii < phases_weight.size(); ii++) {
					sum_weight2 += phases_weight[ii];
					if (sum_weight2 > rand) {
						grains_property.push_back(phases_properties[ii]);
						break;
					}
				}
			}
			// > periodic boundary condition
			vector<vector<Point>> mirror_points; // = 27 * points.size()
			int region_number = 0;
			if (dimention == Dimension::Three_Dimension)
				region_number = 27;
			else
				region_number = 9;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = points[grain] + Point(0, 0, 0);
				mirror_points[1][grain] = points[grain] + Point(int(box_size[0]), 0, 0);
				mirror_points[2][grain] = points[grain] + Point(int(-box_size[0]), 0, 0);
				mirror_points[3][grain] = points[grain] + Point(0, int(box_size[1]), 0);
				mirror_points[4][grain] = points[grain] + Point(0, int(-box_size[1]), 0);
				mirror_points[5][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), 0);
				mirror_points[6][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), 0);
				mirror_points[7][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), 0);
				mirror_points[8][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), 0);
				if (dimention == Dimension::Three_Dimension) {
					mirror_points[9][grain] = points[grain] + Point(0, 0, int(box_size[2]));
					mirror_points[10][grain] = points[grain] + Point(0, 0, int(-box_size[2]));
					mirror_points[11][grain] = points[grain] + Point(int(box_size[0]), 0, int(box_size[2]));
					mirror_points[12][grain] = points[grain] + Point(int(-box_size[0]), 0, int(box_size[2]));
					mirror_points[13][grain] = points[grain] + Point(int(box_size[0]), 0, int(-box_size[2]));
					mirror_points[14][grain] = points[grain] + Point(int(-box_size[0]), 0, int(-box_size[2]));
					mirror_points[15][grain] = points[grain] + Point(0, int(box_size[1]), int(box_size[2]));
					mirror_points[16][grain] = points[grain] + Point(0, int(-box_size[1]), int(box_size[2]));
					mirror_points[17][grain] = points[grain] + Point(0, int(box_size[1]), int(-box_size[2]));
					mirror_points[18][grain] = points[grain] + Point(0, int(-box_size[1]), int(-box_size[2]));
					mirror_points[19][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), int(box_size[2]));
					mirror_points[20][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), int(box_size[2]));
					mirror_points[23][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), int(box_size[2]));
					mirror_points[24][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), int(-box_size[2]));
					mirror_points[26][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), int(-box_size[2]));
					mirror_points[21][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), int(box_size[2]));
					mirror_points[22][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), int(-box_size[2]));
					mirror_points[25][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), int(-box_size[2]));
				}
			}
#pragma omp parallel for
			for (int region = 0; region < region_number; region++)
				for (int grain = 0; grain < grain_number; grain++) {
					Polyhedron poly(mirror_points[region][grain]);
					vector<point_in_region_index> record_points;
					point_in_region_index rp(region, grain);
					record_points.push_back(rp);
					for (unsigned int region_index = 0; region_index < mirror_points.size(); region_index++)
						for (unsigned int grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
							// Avoid inclusion points
							bool is_point_contained = false;
							for (auto re = record_points.begin(); re < record_points.end(); re++)
								if (re->region == region_index && re->grain_index == grain_index)
									is_point_contained = true;
							if (is_point_contained)
								continue;
							// prepare
							Point norm, mid_point;
							norm = poly.point_inside_polyhedron - mirror_points[region_index][grain_index];
							mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
							// Judged ipsilateral to poly center
							if (poly.check_point_inside_polyhedron(mid_point) == false)
								continue;
							// < add point
							poly.add_surf(norm, mid_point);
							// < Eliminate meaningless points in poly
							for (auto re = record_points.begin(); re < record_points.end();) {
								Point check = (points[grain] + mirror_points[re->region][re->grain_index]) / 2;
								if (poly.check_point_inside_polyhedron(check) == false) {
									re = record_points.erase(re);
								}
								else {
									++re;
								}
							}
						}
					PointSet set;
					double sum_sum_phi = 0.0;
					for (auto node = phaseMesh._mesh.begin(); node < phaseMesh._mesh.end(); node++)
						if (poly.check_point_inside_polyhedron(pf::Point(node->_x, node->_y, node->_z)) == true) {
							double sum_phi = 0.0;
							for (auto phase = node->begin(); phase < node->end(); phase++)
								for (auto index = phis_index.begin(); index < phis_index.end(); index++)
									if (phase->index == *index && phase->phi > SYS_EPSILON) {
										sum_phi += phase->phi;
										phase->phi = 0.0;
									}
							if (sum_phi > SYS_EPSILON) {
								set.add_point(node->_x, node->_y, node->_z, sum_phi);
								sum_sum_phi += sum_phi;
							}
						}
					string str_report = "> Voronoi: One polyhedron in region : " + to_string(region) + ", grain : " + to_string(grain) + " could been generated ! \n";
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(str_report, LOG_FILE_NAME);
					set.generate_step = generate_step;
					set.phaseIndex = grain_0_index + grain;
					set.phaseProperty = grains_property[grain];
					set.temperature = temperature;
					set.x = x;
					set.is_normalized = false;
					if (sum_sum_phi > SYS_EPSILON) {
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							nucleation_box.point_set_box.push_back(set);
						}
					}
				}
		}
#ifdef _WIN32
		static void generate_structure_from_BMP_pic(int phaseIndex, int phase_property, int generate_step, XNode x, double threshold[2], double temperature, string fileName, double phi = 1.0, bool isNormalized = true) {
			FieldStorage_forPhaseNode& phaseMesh = Solvers::get_instance()->phaseMesh;
			if (phaseMesh.limit_z > 1) {
				string error_report = "> generate structure error : Nz > 1, cant init structure by a BMP picture !\n";
				InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
				return;
			}

			BMP24reader bmpReader;
			if (fileName == "") {
				SelectFilePath(fileName);
			}
			else {
				string input_file_folder_path = pf::GetFolderOfPath(Solvers::get_instance()->Working_Folder_Path);
				fileName = input_file_folder_path + dirSeparator + fileName;
			}
			bmpReader.safe(fileName);
			bmpReader.read(fileName);
			PointSet set;
#pragma omp parallel for
			for (int y = 0; y < Solvers::get_instance()->phaseMesh.limit_y; y++)
				for (int x = 0; x < Solvers::get_instance()->phaseMesh.limit_x; x++) {
					int xx = double_to_int(x * bmpReader.bmp_width / Solvers::get_instance()->phaseMesh.limit_x), yy = double_to_int(y * bmpReader.bmp_height / Solvers::get_instance()->phaseMesh.limit_y);
					double graypercent = bmpReader.getGrayPercentage(xx, yy);
					if (graypercent > threshold[0] && graypercent < threshold[1])
						set.add_point(x, y, 0, phi);
				}
			set.generate_step = generate_step;
			set.phaseIndex = phaseIndex;
			set.phaseProperty = phase_property;
			set.temperature = temperature;
			set.x = x;
			set.is_normalized = isNormalized;
			nucleation_box.point_set_box.push_back(set);
		}
#endif

		static void init_mesh_with_datafile(FieldStorage_forPhaseNode& phaseMesh, Data_report& report, string datafile_path, bool infile_debug) {
			if (!Solvers::get_instance()->data_writer.read_dataFile(phaseMesh, datafile_path, datafile_report)) {
				string _report = "> Error : datafile_path can't be opened ! \n";
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt(_report, InputFileReader::get_instance()->debug_file);
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(_report, LOG_FILE_NAME);
				std::exit(0);
			}
			else {
				// check datafile with this input file
				vector<string> bc_type; bc_type.push_back("FIXED"); bc_type.push_back("PERIODIC"); bc_type.push_back("ADIABATIC");
				vector<string> bool_type; bool_type.push_back("IN-VALID"); bool_type.push_back("VALID");
				vector<string> custom_type; custom_type.push_back("DOUBLE"); custom_type.push_back("INT"); custom_type.push_back("VEC3");
				stringstream _report;
				_report << "> " << endl;
				bool line_valid = true, all_valid = true;
				_report << "> | init microstructure with datafile :                    (check in this simulation)" << endl;
				if (datafile_report.Nx == phaseMesh.limit_x && datafile_report.Ny == phaseMesh.limit_y && datafile_report.Nz == phaseMesh.limit_z && datafile_report.dr == phaseMesh.dr)
					line_valid = true;
				else {
					line_valid = false;
					all_valid = false;
				}
				_report << "> | mesh size               : Nx - " << datafile_report.Nx << ", Ny - " << datafile_report.Ny << ", Nz - " << datafile_report.Nz << ", dr - " << datafile_report.dr << " (" << bool_type[line_valid] << ")" << endl;
				if (datafile_report.x_down_bc == phaseMesh._bc_x_down && datafile_report.x_up_bc == phaseMesh._bc_x_up)
					line_valid = true;
				else {
					line_valid = false;
				}
				_report << "> | mesh boundary           : x_down - " << bc_type[datafile_report.x_down_bc] << ", x_up - " << bc_type[datafile_report.x_up_bc] << " (" << bool_type[line_valid] << ")" << endl;
				if (datafile_report.y_down_bc == phaseMesh._bc_y_down && datafile_report.y_up_bc == phaseMesh._bc_y_up)
					line_valid = true;
				else {
					line_valid = false;
				}
				_report << "> |                           y_down - " << bc_type[datafile_report.y_down_bc] << ", y_up - " << bc_type[datafile_report.y_up_bc] << " (" << bool_type[line_valid] << ")" << endl;
				if (datafile_report.z_down_bc == phaseMesh._bc_z_down && datafile_report.z_up_bc == phaseMesh._bc_z_up)
					line_valid = true;
				else {
					line_valid = false;
				}
				_report << "> |                           z_down - " << bc_type[datafile_report.z_down_bc] << ", z_up - " << bc_type[datafile_report.z_up_bc] << " (" << bool_type[line_valid] << ")" << endl;
				for (auto phi = datafile_report.phi_property.begin(); phi < datafile_report.phi_property.end(); phi++) {
					string _name = "";
					for (auto p = Solvers::get_instance()->parameters.Phases.begin(); p < Solvers::get_instance()->parameters.Phases.end(); p++)
						if (p->phi_property == phi->value) {
							line_valid = true;
							_name = p->phi_name;
						}
					if (_name.size() == 0) {
						line_valid = false;
						all_valid = false;
					}
					if (line_valid)
						_report << "> | phi_" << phi->index << "                   : property - " << phi->value << " (" << _name << ")" << endl;
					else
						_report << "> | phi_" << phi->index << "                   : property - " << phi->value << " (" << " " << ")" << endl;
					_report << "> |                           con - ";
					for (auto phi_x = datafile_report.phi_comps[phi->index].begin(); phi_x < datafile_report.phi_comps[phi->index].end(); phi_x++)
						_report << phi_x->index << ", ";
					_report << endl;
				}
				if (datafile_report.comps.size() != 0) {
					string _name = "";
					for (auto c = Solvers::get_instance()->parameters.Components.begin(); c < Solvers::get_instance()->parameters.Components.end(); c++)
						if (c->index == *datafile_report.comps.begin()) {
							line_valid = true;
							_name = c->name;
						}
					if (_name.size() == 0) {
						line_valid = false;
						all_valid = false;
					}
					if (line_valid)
						_report << "> | components in domain        : con - " << *datafile_report.comps.begin() << " (" << _name << ")" << endl;
					else
						_report << "> | components in domain        : con - " << *datafile_report.comps.begin() << " (" << " " << ")" << endl;
					for (auto comp = datafile_report.comps.begin() + 1; comp < datafile_report.comps.end(); comp++) {
						string _name = "";
						for (auto c = Solvers::get_instance()->parameters.Components.begin(); c < Solvers::get_instance()->parameters.Components.end(); c++)
							if (c->index == *comp) {
								line_valid = true;
								_name = c->name;
							}
						if (_name.size() == 0) {
							line_valid = false;
							all_valid = false;
						}
						if (line_valid)
							_report << "> |                             : con - " << *comp << " (" << _name << ")" << endl;
						else
							_report << "> |                             : con - " << *comp << " (" << " " << ")" << endl;
					}
				}
				if (datafile_report.custom_vars.size() != 0) {
					_report << "> | custom value                : index - " << datafile_report.custom_vars.begin()->index << ", type - " << custom_type[datafile_report.custom_vars.begin()->value] << endl;
					for (auto vars = datafile_report.custom_vars.begin() + 1; vars < datafile_report.custom_vars.end(); vars++)
						_report << "> |                             : index - " << vars->index << ", type - " << custom_type[vars->value] << endl;
				}
				_report << ">" << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(_report.str(), LOG_FILE_NAME);
				if (all_valid == false) {
					string error_report = "> Error : mesh data structure from datafile and inputfile mismatch ! \n";
					if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt(error_report, InputFileReader::get_instance()->debug_file);
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(error_report, LOG_FILE_NAME);
					std::exit(0);
				}
			}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			InputFileReader::get_instance()->read_bool_value("Preprocess.Microstructure.is_datafile_init", is_datafile_init, infile_debug);
			if (!is_datafile_init) {
				vector<input_value> init_comps;
				for (auto comp = Solvers::get_instance()->parameters.Components.begin(); comp < Solvers::get_instance()->parameters.Components.end(); comp++) {
					input_value init_comp;
					init_comp.double_value = 0.0;
					init_comps.push_back(init_comp);
				}
				total_component(init_comps);
				string matrix_key = "Preprocess.Microstructure.matrix", matrix_string = "{[()]}";
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .matrix = {[(phi_index),(phi_name),(phi_comp_0_value, phi_comp_1_value, ... )],[(total_comp_0_value, total_comp_1_value, ... )],[(temp_value)]} \n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_string, infile_debug)) {
					vector<vector<InputValueType>> matrix_structure; matrix_structure.resize(3); matrix_structure[0].resize(3); matrix_structure[1].resize(1); matrix_structure[2].resize(1);
					matrix_structure[0][0] = InputValueType::IVType_INT; matrix_structure[0][1] = InputValueType::IVType_STRING; matrix_structure[0][2] = InputValueType::IVType_DOUBLE;
					matrix_structure[1][0] = InputValueType::IVType_DOUBLE; matrix_structure[2][0] = InputValueType::IVType_DOUBLE;
					vector<vector<vector<input_value>>> matrix_value;
					matrix_value = InputFileReader::get_instance()->trans_matrix_3d_array_array_const_to_input_value(matrix_structure, matrix_key, matrix_string, infile_debug);
					double default_phi = 1.0;
					InputFileReader::get_instance()->read_double_value("Preprocess.Microstructure.Matrix.phi", default_phi, infile_debug);
					total_component(matrix_value[1][0]);
					temperature(matrix_value[2][0][0].double_value);
					add_new_phi_by_index(matrix_value[0][0][0].int_value, Solvers::get_instance()->parameters.Phases[matrix_value[0][1][0].string_value].phi_property, default_phi);
					matrix_phase_component(matrix_value[0][0][0].int_value, Solvers::get_instance()->parameters.Phases[matrix_value[0][1][0].string_value].phi_property, matrix_value[0][2]);
				}

				// voronoi structure
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .property = [(phi_index_begin, phi_index_end), (phi_name, ... ),(phi_weight, ...)] \n", InputFileReader::get_instance()->debug_file);
				string voronoi_property_key = "Preprocess.Microstructure.voronoi.property", voronoi_property_input = "[()]";
				vector<InputValueType> voronoi_property_structure; voronoi_property_structure.push_back(InputValueType::IVType_INT); voronoi_property_structure.push_back(InputValueType::IVType_STRING);
				voronoi_property_structure.push_back(InputValueType::IVType_DOUBLE);
				if (InputFileReader::get_instance()->read_string_value(voronoi_property_key, voronoi_property_input, infile_debug)) {
					vector<vector<input_value>> voronoi_value = InputFileReader::get_instance()->trans_matrix_2d_array_const_to_input_value(voronoi_property_structure, voronoi_property_key, voronoi_property_input, infile_debug);
					vector<int> _vor_property; vector<double> _vor_weight;
					for (int vor_index = 0; vor_index < voronoi_value[1].size(); vor_index++)
						_vor_property.push_back(Solvers::get_instance()->parameters.Phases[voronoi_value[1][vor_index].string_value].phi_property);
					for (int vor_index = 0; vor_index < voronoi_value[2].size(); vor_index++)
						_vor_weight.push_back(voronoi_value[2][vor_index].double_value);

					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .box = [(box_origin_point),(box_end_point)] \n", InputFileReader::get_instance()->debug_file);
					string voronoi_box_key = "Preprocess.Microstructure.voronoi.box", voronoi_box_input = "[(0,0,0),(0,0,0)]";
					InputFileReader::get_instance()->read_string_value(voronoi_box_key, voronoi_box_input, infile_debug);
					vector<vector<input_value>> voronoi_box_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, voronoi_box_key, voronoi_box_input, infile_debug);
					Vector3 box_pos, box_size; box_pos[0] = voronoi_box_value[0][0].double_value; box_pos[1] = voronoi_box_value[0][1].double_value; box_pos[2] = voronoi_box_value[0][2].double_value;
					box_size[0] = voronoi_box_value[1][0].double_value; box_size[1] = voronoi_box_value[1][1].double_value; box_size[2] = voronoi_box_value[1][2].double_value;

					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .x = [(comp_0_name,comp_0_value), (comp_1_name,comp_1_value), ...] \n", InputFileReader::get_instance()->debug_file);
					string voronoi_x_key = "Preprocess.Microstructure.voronoi.x", voronoi_x_input = "[()]";
					pf::XNode vor_x;
					if (InputFileReader::get_instance()->read_string_value(voronoi_x_key, voronoi_x_input, infile_debug)) {
						vector<InputValueType> voronoi_x_structure; voronoi_x_structure.push_back(InputValueType::IVType_STRING); voronoi_x_structure.push_back(InputValueType::IVType_DOUBLE);
						vector<vector<input_value>> voronoi_x_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(voronoi_x_structure, voronoi_x_key, voronoi_x_input, infile_debug);
						for (int x_index = 0; x_index < voronoi_x_value.size(); x_index++)
							vor_x.add_x(Solvers::get_instance()->parameters.Components[voronoi_x_value[x_index][0].string_value].index, voronoi_x_value[x_index][1].double_value);
					}
					string voronoi_temperature_key = "Preprocess.Microstructure.voronoi.temperature";
					double vor_temp = 0.0;
					InputFileReader::get_instance()->read_double_value(voronoi_temperature_key, vor_temp, infile_debug);
					if (InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.voronoi.rand_seed", rand_seed, infile_debug))
						is_voronoi_rand = false;

					InputFileReader::get_instance()->read_double_value("Preprocess.Microstructure.voronoi.point_distance", voronoi_point_distance, infile_debug);

					generate_voronoi_structure(box_pos, box_size, voronoi_value[0][0].int_value, voronoi_value[0][1].int_value - voronoi_value[0][0].int_value + 1, 0, _vor_property, _vor_weight, vor_x, vor_temp);
				}

				// read bmp24
#ifdef _WIN32
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .property = (bmp_layer, file_name) \n", InputFileReader::get_instance()->debug_file);
				string bmp24_property_key = "Preprocess.Microstructure.bmp24.property", bmp24_property_input = "()";
				vector<InputValueType> bmp24_property_structure; bmp24_property_structure.push_back(InputValueType::IVType_INT); bmp24_property_structure.push_back(InputValueType::IVType_STRING);
				if (InputFileReader::get_instance()->read_string_value(bmp24_property_key, bmp24_property_input, infile_debug)) {
					vector<input_value> bmp24_property_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(bmp24_property_structure, bmp24_property_key, bmp24_property_input, infile_debug);
					for (int bmp_layer = 0; bmp_layer < bmp24_property_value[0].int_value; bmp_layer++) {
						string bmp24_phi_key = "Preprocess.Microstructure.bmp24_" + to_string(bmp_layer) + ".phi", bmp24_phi_input = "()";
						if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .phi = (phi_index, phi_name, phi_fraction) \n", InputFileReader::get_instance()->debug_file);
						if (InputFileReader::get_instance()->read_string_value(bmp24_phi_key, bmp24_phi_input, infile_debug)) {
							vector<InputValueType> bmp24_phi_structure; bmp24_phi_structure.push_back(InputValueType::IVType_INT); bmp24_phi_structure.push_back(InputValueType::IVType_STRING); 
							bmp24_phi_structure.push_back(InputValueType::IVType_DOUBLE);
							vector<input_value> bmp24_phi_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(bmp24_phi_structure, bmp24_phi_key, bmp24_phi_input, infile_debug);

							if (infile_debug)
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .x = [(comp_0_name,comp_0_value), ...] \n", InputFileReader::get_instance()->debug_file);
							string bmp24_x_key = "Preprocess.Microstructure.bmp24_" + to_string(bmp_layer) + ".x", bmp24_x_input = "[()]";
							pf::XNode bmp24_x;
							if (InputFileReader::get_instance()->read_string_value(bmp24_x_key, bmp24_x_input, infile_debug)) {
								vector<InputValueType> bmp24_x_structure; bmp24_x_structure.push_back(InputValueType::IVType_STRING); bmp24_x_structure.push_back(InputValueType::IVType_DOUBLE);
								vector<vector<input_value>> bmp24_x_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(bmp24_x_structure, bmp24_x_key, bmp24_x_input, infile_debug);
								for (int x_index = 0; x_index < bmp24_x_value.size(); x_index++)
									bmp24_x.add_x(Solvers::get_instance()->parameters.Components[bmp24_x_value[x_index][0].string_value].index, bmp24_x_value[x_index][1].double_value);
							}
							string bmp24_threshold_key = "Preprocess.Microstructure.bmp24_" + to_string(bmp_layer) + ".gray_threshold", bmp24_threshold_input = "(0.0,1.0)";
							InputFileReader::get_instance()->read_string_value(bmp24_threshold_key, bmp24_threshold_input, infile_debug);
							vector<input_value> bmp24_threshold_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, bmp24_threshold_key, bmp24_threshold_input, infile_debug);
							double bmp24_threshold[] = { bmp24_threshold_value[0].double_value, bmp24_threshold_value[1].double_value };
							string bmp24_temp_key = "Preprocess.Microstructure.bmp24_" + to_string(bmp_layer) + ".temperature"; double bmp_temp = 0.0;
							InputFileReader::get_instance()->read_double_value(bmp24_temp_key, bmp_temp, infile_debug);
							string bmp24_norm_key = "Preprocess.Microstructure.bmp24_" + to_string(bmp_layer) + ".is_normalized"; bool is_normalized = true;
							InputFileReader::get_instance()->read_bool_value(bmp24_norm_key, is_normalized, infile_debug);
							generate_structure_from_BMP_pic(bmp24_phi_value[0].int_value, Solvers::get_instance()->parameters.Phases[bmp24_phi_value[1].string_value].phi_property,
								0, bmp24_x, bmp24_threshold, bmp_temp, bmp24_property_value[1].string_value, bmp24_phi_value[2].double_value, is_normalized);
						}
					}
				}
#endif
			}
			else {
				bool is_read_datafile_by_path = false;
				string datafile_path = "DATA.dat";
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Preprocess.Microstructure.datafile_path : relative path from infile folder.\n", InputFileReader::get_instance()->debug_file);
				if (InputFileReader::get_instance()->read_string_value("Preprocess.Microstructure.datafile_path", datafile_path, infile_debug))
					is_read_datafile_by_path = true;
				datafile_path = Solvers::get_instance()->Infile_Folder_Path + dirSeparator + datafile_path;
#ifdef _WIN32
				if (!is_read_datafile_by_path) {
					Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Please select a datafile (in .dat format) to initialize the simulation mesh...", LOG_FILE_NAME);
					SelectFilePath(datafile_path);
				}
#else
				if (!is_read_datafile_by_path) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("> Error : Preprocess.Microstructure.datafile_path should be set here ! \n", InputFileReader::get_instance()->debug_file);
					std::exit(0);
				}
#endif
				string out = "> Open datafile: " + datafile_path + "\n";
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(out, LOG_FILE_NAME);
				init_mesh_with_datafile(phaseMesh, datafile_report, datafile_path, infile_debug);
				
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								for (auto c = phase->x.begin(); c < phase->x.end(); c++) {
									phase->potential.add_con(c->index);
									phaseMesh.info_node[phase->index].potential.add_con(c->index);
								}
								for (auto comp1 = phase->x.begin(); comp1 < phase->x.end(); comp1++)
									for (auto comp2 = phase->x.begin(); comp2 < phase->x.end(); comp2++)
										phase->kinetics_coeff.set(comp1->index, comp2->index, 0.0);
							}
							for (auto comp1 = node.x.begin(); comp1 < node.x.end(); comp1++)
								for (auto comp2 = node.x.begin(); comp2 < node.x.end(); comp2++)
									node.kinetics_coeff.set(comp1->index, comp2->index, 0.0);
						}
			}

			nucleation_box.nucleation_property = NucleationProperty::DefiniteNucleation;

			// geometry structure
			int nulceation_layer = 0;
			InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.geometry_layer_number", nulceation_layer, infile_debug);
			for (int layer_index = 0; layer_index < nulceation_layer; layer_index++) {
				GeometricRegion geo;
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .property = (phi_index, phi_name, geometry_type, rotation_gauge, reverse_region) \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              geometry_type  : 0 - None, 1 - Ellipsoid, 2 - Polyhedron \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#              rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n", InputFileReader::get_instance()->debug_file);
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n", InputFileReader::get_instance()->debug_file);
				}
				string layer_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".property";
				string layer_input_property = "(0,phi_name,0,1,false)";
				if (InputFileReader::get_instance()->read_string_value(layer_key, layer_input_property, infile_debug)) {
					vector<InputValueType> layer_input_property_structure; layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_STRING);
					layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_BOOL);
					vector<input_value> layer_input_property_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(layer_input_property_structure, layer_key, layer_input_property, infile_debug);
					int phi_property = Solvers::get_instance()->parameters.Phases[layer_input_property_value[1].string_value].phi_property;
					geo.init(pf::Geometry(layer_input_property_value[2].int_value), 0, layer_input_property_value[0].int_value, phi_property, layer_input_property_value[4].bool_value);
					pf::RotationGauge rotation_gauge = pf::RotationGauge(layer_input_property_value[3].int_value);
					if (geo.geometryProperty == pf::Geometry::Geo_Ellipsoid) {
						if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .ellipsoid = [(core_x,core_y,core_z),(radius_x,radius_y,radius_z),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n", InputFileReader::get_instance()->debug_file);
						string layer_ellipsoid_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".ellipsoid";
						string layer_input_ellipsoid = "[(0,0,0),(0,0,0),(0,0,0)]";
						InputFileReader::get_instance()->read_string_value(layer_ellipsoid_key, layer_input_ellipsoid, infile_debug);
						vector<vector<input_value>> layer_input_ellipsoid_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, layer_ellipsoid_key, layer_input_ellipsoid, infile_debug);
						geo.ellipSolid.set_core(layer_input_ellipsoid_value[0][0].double_value, layer_input_ellipsoid_value[0][1].double_value, layer_input_ellipsoid_value[0][2].double_value);
						geo.ellipSolid.set_radius(layer_input_ellipsoid_value[1][0].double_value, layer_input_ellipsoid_value[1][1].double_value, layer_input_ellipsoid_value[1][2].double_value);
						double radian[] = { AngleToRadians(layer_input_ellipsoid_value[2][0].double_value), AngleToRadians(layer_input_ellipsoid_value[2][1].double_value), AngleToRadians(layer_input_ellipsoid_value[2][2].double_value) };
						geo.ellipSolid.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
					}
					if (geo.geometryProperty == pf::Geometry::Geo_Polyhedron) {
						if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .polyhedron = {[inside_point],[surf_point,surf_point,surf_point], .... ,[(rotation_angle_1,rotation_angle_2,rotation_angle_3)]} \n", InputFileReader::get_instance()->debug_file);
						if (infile_debug)
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("#                point = (position_x,position_y,position_z) \n", InputFileReader::get_instance()->debug_file);
						string layer_polyhedron_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".polyhedron";
						string layer_input_polyhedron = "{[(-1,0,0)],[(0,0,0),(0,1,0),(0,0,1)],[(0,0,0)]}";
						InputFileReader::get_instance()->read_string_value(layer_polyhedron_key, layer_input_polyhedron, infile_debug);
						vector<vector<vector<input_value>>> layer_input_polyhedron_value = InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_DOUBLE, layer_polyhedron_key, layer_input_polyhedron, infile_debug);
						geo.polyhedron.set_a_point_inside_polyhedron(pf::Point(layer_input_polyhedron_value[0][0][0].double_value, layer_input_polyhedron_value[0][0][1].double_value, layer_input_polyhedron_value[0][0][2].double_value));
						for (int surf_index = 1; surf_index < layer_input_polyhedron_value.size() - 1; surf_index++) {
							geo.polyhedron.add_surf(pf::Point(layer_input_polyhedron_value[surf_index][0][0].double_value, layer_input_polyhedron_value[surf_index][0][1].double_value, layer_input_polyhedron_value[surf_index][0][2].double_value),
								pf::Point(layer_input_polyhedron_value[surf_index][1][0].double_value, layer_input_polyhedron_value[surf_index][1][1].double_value, layer_input_polyhedron_value[surf_index][1][2].double_value),
								pf::Point(layer_input_polyhedron_value[surf_index][2][0].double_value, layer_input_polyhedron_value[surf_index][2][1].double_value, layer_input_polyhedron_value[surf_index][2][2].double_value));
						}
						int rIndex = int(layer_input_polyhedron_value.size() - 1);
						double radian[] = { AngleToRadians(layer_input_polyhedron_value[rIndex][0][0].double_value), AngleToRadians(layer_input_polyhedron_value[rIndex][0][1].double_value), AngleToRadians(layer_input_polyhedron_value[rIndex][0][2].double_value) };
						geo.polyhedron.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
					}
					string layer_temp_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".T";
					InputFileReader::get_instance()->read_double_value(layer_temp_key, geo.temperature, infile_debug);
					string layer_phi_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".phi";
					InputFileReader::get_instance()->read_double_value(layer_phi_key, geo.phi, infile_debug);
					string layer_norm_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".is_normalized";
					InputFileReader::get_instance()->read_bool_value(layer_norm_key, geo.isNormalized, infile_debug);
					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .x = [(comp_0_name,comp_0_value),(comp_1_name,comp_1_value), ...] \n", InputFileReader::get_instance()->debug_file);
					string layer_x_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".x", layer_x_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(layer_x_key, layer_x_input, infile_debug)) {
						vector<InputValueType> layer_input_x_structure; layer_input_x_structure.push_back(InputValueType::IVType_STRING); layer_input_x_structure.push_back(InputValueType::IVType_DOUBLE);
						vector<vector<input_value>> layer_input_x_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(layer_input_x_structure, layer_x_key, layer_x_input, infile_debug);
						for (int x_index = 0; x_index < layer_input_x_value.size(); x_index++)
							geo.x.add_x(Solvers::get_instance()->parameters.Components[layer_input_x_value[x_index][0].string_value].index, layer_input_x_value[x_index][1].double_value);
					}
					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .custom_int = [(custom_0_index, custom_0_value),(custom_1_index, custom_1_value), ...] \n", InputFileReader::get_instance()->debug_file);
					string layer_custom_int_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".custom_int", layer_custom_int_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(layer_custom_int_key, layer_custom_int_input, infile_debug)) {
						vector<InputValueType> layer_input_custom_int_structure; layer_input_custom_int_structure.push_back(InputValueType::IVType_INT); layer_input_custom_int_structure.push_back(InputValueType::IVType_INT);
						vector<vector<input_value>> layer_input_custom_int_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(layer_input_custom_int_structure, layer_custom_int_key, layer_custom_int_input, infile_debug);
						for (int custom_int_index = 0; custom_int_index < layer_input_custom_int_value.size(); custom_int_index++)
							geo.customFlags.add_int(layer_input_custom_int_value[custom_int_index][0].int_value, layer_input_custom_int_value[custom_int_index][1].int_value);
					}
					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .custom_double = [(custom_0_index, custom_0_value),(custom_1_index, custom_1_value), ...] \n", InputFileReader::get_instance()->debug_file);
					string layer_custom_double_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".custom_double", layer_custom_double_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(layer_custom_double_key, layer_custom_double_input, infile_debug)) {
						vector<InputValueType> layer_input_custom_double_structure; layer_input_custom_double_structure.push_back(InputValueType::IVType_INT); layer_input_custom_double_structure.push_back(InputValueType::IVType_DOUBLE);
						vector<vector<input_value>> layer_input_custom_double_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(layer_input_custom_double_structure, layer_custom_double_key, layer_custom_double_input, infile_debug);
						for (int custom_double_index = 0; custom_double_index < layer_input_custom_double_value.size(); custom_double_index++)
							geo.customValues.add_double(layer_input_custom_double_value[custom_double_index][0].int_value, layer_input_custom_double_value[custom_double_index][1].double_value);
					}
					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .custom_vec3 = [(custom_0_index, custom_0_value_0, custom_0_value_1, custom_0_value_2), (custom_1_index, custom_1_value_0, custom_1_value_1, custom_1_value_2), ...] \n", InputFileReader::get_instance()->debug_file);
					string layer_custom_vec3_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".custom_vec3", layer_custom_vec3_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(layer_custom_vec3_key, layer_custom_vec3_input, infile_debug)) {
						vector<InputValueType> layer_input_custom_vec3_structure; layer_input_custom_vec3_structure.push_back(InputValueType::IVType_INT); layer_input_custom_vec3_structure.push_back(InputValueType::IVType_DOUBLE);
						layer_input_custom_vec3_structure.push_back(InputValueType::IVType_DOUBLE); layer_input_custom_vec3_structure.push_back(InputValueType::IVType_DOUBLE);
						vector<vector<input_value>> layer_input_custom_vec3_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(layer_input_custom_vec3_structure, layer_custom_vec3_key, layer_custom_vec3_input, infile_debug);
						for (int custom_vec3_index = 0; custom_vec3_index < layer_input_custom_vec3_value.size(); custom_vec3_index++) {
							double vec[] = { layer_input_custom_vec3_value[custom_vec3_index][1].double_value,  layer_input_custom_vec3_value[custom_vec3_index][2].double_value, layer_input_custom_vec3_value[custom_vec3_index][3].double_value };
							geo.customVec3s.add_vec(layer_input_custom_vec3_value[custom_vec3_index][0].int_value, vec);
						}
					}
					if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .custom_vec6 = [(custom_vec6_index, custom_vec6_value_0, custom_vec6_value_1, custom_vec6_value_2, custom_vec6_value_3, custom_vec6_value_4, custom_vec6_value_5), ...] \n", InputFileReader::get_instance()->debug_file);
					string layer_custom_vec6_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".custom_vec6", layer_custom_vec6_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(layer_custom_vec6_key, layer_custom_vec6_input, infile_debug)) {
						vector<InputValueType> layer_input_custom_vec6_structure; layer_input_custom_vec6_structure.push_back(InputValueType::IVType_INT); layer_input_custom_vec6_structure.push_back(InputValueType::IVType_DOUBLE);
						layer_input_custom_vec6_structure.push_back(InputValueType::IVType_DOUBLE); layer_input_custom_vec6_structure.push_back(InputValueType::IVType_DOUBLE); layer_input_custom_vec6_structure.push_back(InputValueType::IVType_DOUBLE);
						layer_input_custom_vec6_structure.push_back(InputValueType::IVType_DOUBLE); layer_input_custom_vec6_structure.push_back(InputValueType::IVType_DOUBLE);
						vector<vector<input_value>> layer_input_custom_vec6_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(layer_input_custom_vec6_structure, layer_custom_vec6_key, layer_custom_vec6_input, infile_debug);
						for (int custom_vec6_index = 0; custom_vec6_index < layer_input_custom_vec6_value.size(); custom_vec6_index++) {
							double vec[] = { layer_input_custom_vec6_value[custom_vec6_index][1].double_value,  layer_input_custom_vec6_value[custom_vec6_index][2].double_value, layer_input_custom_vec6_value[custom_vec6_index][3].double_value,
							layer_input_custom_vec6_value[custom_vec6_index][4].double_value,  layer_input_custom_vec6_value[custom_vec6_index][5].double_value, layer_input_custom_vec6_value[custom_vec6_index][6].double_value };
							geo.customVec6s.add_vec(layer_input_custom_vec6_value[custom_vec6_index][0].int_value, vec);
						}
					}
					nucleation_box.geometry_box.push_back(geo);
				}
			}

			// generate voronoi in phases
			// voronoi structure
			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .property = [(phi_index_begin, phi_index_end), (phi_name, ... ),(phi_weight, ...),(in phi_index, ... )] \n", InputFileReader::get_instance()->debug_file);
			string voronoi_property_key = "Preprocess.Microstructure.voronoi_inPhis.property", voronoi_property_input = "[()]";
			vector<InputValueType> voronoi_property_structure; voronoi_property_structure.push_back(InputValueType::IVType_INT); voronoi_property_structure.push_back(InputValueType::IVType_STRING);
			voronoi_property_structure.push_back(InputValueType::IVType_DOUBLE); voronoi_property_structure.push_back(InputValueType::IVType_INT);
			if (InputFileReader::get_instance()->read_string_value(voronoi_property_key, voronoi_property_input, infile_debug)) {
				vector<vector<input_value>> voronoi_value = InputFileReader::get_instance()->trans_matrix_2d_array_const_to_input_value(voronoi_property_structure, voronoi_property_key, voronoi_property_input, infile_debug);
				vector<int> _vor_property; vector<double> _vor_weight; vector<int> _in_phi_index;
				for (int vor_index = 0; vor_index < voronoi_value[1].size(); vor_index++)
					_vor_property.push_back(Solvers::get_instance()->parameters.Phases[voronoi_value[1][vor_index].string_value].phi_property);
				for (int vor_index = 0; vor_index < voronoi_value[2].size(); vor_index++)
					_vor_weight.push_back(voronoi_value[2][vor_index].double_value);
				for (int in_phi_index = 0; in_phi_index < voronoi_value[3].size(); in_phi_index++)
					_in_phi_index.push_back(voronoi_value[3][in_phi_index].int_value);

				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .box = [(box_origin_point),(box_end_point)] \n", InputFileReader::get_instance()->debug_file);
				string voronoi_box_key = "Preprocess.Microstructure.voronoi_inPhis.box", voronoi_box_input = "[(0,0,0),(0,0,0)]";
				InputFileReader::get_instance()->read_string_value(voronoi_box_key, voronoi_box_input, infile_debug);
				vector<vector<input_value>> voronoi_box_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, voronoi_box_key, voronoi_box_input, infile_debug);
				Vector3 box_pos, box_size; box_pos[0] = voronoi_box_value[0][0].double_value; box_pos[1] = voronoi_box_value[0][1].double_value; box_pos[2] = voronoi_box_value[0][2].double_value;
				box_size[0] = voronoi_box_value[1][0].double_value; box_size[1] = voronoi_box_value[1][1].double_value; box_size[2] = voronoi_box_value[1][2].double_value;

				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .x = [(comp_0_name,comp_0_value), (comp_1_name,comp_1_value), ...] \n", InputFileReader::get_instance()->debug_file);
				string voronoi_x_key = "Preprocess.Microstructure.voronoi_inPhis.x", voronoi_x_input = "[()]";
				pf::XNode vor_x;
				if (InputFileReader::get_instance()->read_string_value(voronoi_x_key, voronoi_x_input, infile_debug)) {
					vector<InputValueType> voronoi_x_structure; voronoi_x_structure.push_back(InputValueType::IVType_STRING); voronoi_x_structure.push_back(InputValueType::IVType_DOUBLE);
					vector<vector<input_value>> voronoi_x_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(voronoi_x_structure, voronoi_x_key, voronoi_x_input, infile_debug);
					for (int x_index = 0; x_index < voronoi_x_value.size(); x_index++)
						vor_x.add_x(Solvers::get_instance()->parameters.Components[voronoi_x_value[x_index][0].string_value].index, voronoi_x_value[x_index][1].double_value);
				}
				string voronoi_temperature_key = "Preprocess.Microstructure.voronoi_inPhis.temperature";
				double vor_temp = 0.0;
				InputFileReader::get_instance()->read_double_value(voronoi_temperature_key, vor_temp, infile_debug);
				if (InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.voronoi_inPhis.rand_seed", rand_seed, infile_debug))
					is_voronoi_rand = false;

				InputFileReader::get_instance()->read_double_value("Preprocess.Microstructure.voronoi.point_distance", voronoi_point_distance, infile_debug);

				generate_voronoi_structure_in_phis(box_pos, box_size, voronoi_value[0][0].int_value, voronoi_value[0][1].int_value - voronoi_value[0][0].int_value + 1, 0, _vor_property, _vor_weight, vor_x, vor_temp, _in_phi_index, phaseMesh);
			}

			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP || Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				vector<int> phi_indexes = Solvers::get_instance()->C_Solver.phase_indexes;
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phi_indexes));
							node.customValues.add_double(ExternalFields::CON_Smooth_Old_Phi, node.customValues[ExternalFields::CON_Smooth_Phi]);
						}
			}

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Microstructure !\n", LOG_FILE_NAME);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			definiteNucleation();
			if (Solvers::get_instance()->parameters.is_Normalize_Phi)
				phaseMesh.normalize_phi_in_mesh();
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(report, LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			nucleation_box.condition_phi_box.clear();
			nucleation_box.geometry_box.clear();
			nucleation_box.nucleus_box.clear();
			nucleation_box.point_set_box.clear();
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}