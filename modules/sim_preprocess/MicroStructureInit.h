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
					else if (geo->geometryProperty == Geometry::Geo_Cylindricity) {
						for (int z = 0; z < phaseMesh.limit_z; z++)
							for (int y = 0; y < phaseMesh.limit_y; y++)
								for (int x = 0; x < phaseMesh.limit_x; x++) {
									PhaseNode& node = phaseMesh(x, y, z);
									// Vector3 p(x, y, z);
									Vector3 p(x - geo->cylindricity.core.x, y - geo->cylindricity.core.y, z - geo->cylindricity.core.z);
									p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->cylindricity.radian_x, geo->cylindricity.radian_y, geo->cylindricity.radian_z),
										geo->cylindricity.rotationGauge));
									p[0] += geo->cylindricity.core.x;
									p[1] += geo->cylindricity.core.y;
									p[2] += geo->cylindricity.core.z;
									Point po(p[0], p[1], p[2]);
									bool check0 = geo->cylindricity.check_point_inside_cylindricity(po);
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
						report << "> A new Cylindricity for grain : " << to_string(geo->phaseIndex) << " phase: " << Solvers::get_instance()->parameters.Phases[geo->phaseProperty].phi_name << " has been initialized at position ( "
							<< geo->cylindricity.core.x << ", " << geo->cylindricity.core.y << ", " << geo->cylindricity.core.z << " )." << std::endl;
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

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			if (true) {
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
				
				// read bmp24

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
					if (geo.geometryProperty == pf::Geometry::Geo_Cylindricity) {
						if (infile_debug)
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# .cylindricity = [(core_x,core_y,core_z),(radius_x,radius_z,half_height),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n", InputFileReader::get_instance()->debug_file);
						string layer_cylindricity_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".cylindricity";
						string layer_input_cylindricity = "[(0,0,0),(0,0,0),(0,0,0)]";
						InputFileReader::get_instance()->read_string_value(layer_cylindricity_key, layer_input_cylindricity, infile_debug);
						vector<vector<input_value>> layer_input_cylindricity_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, layer_cylindricity_key, layer_input_cylindricity, infile_debug);
						geo.cylindricity.set_core(layer_input_cylindricity_value[0][0].double_value, layer_input_cylindricity_value[0][1].double_value, layer_input_cylindricity_value[0][2].double_value);
						geo.cylindricity.set_radius(layer_input_cylindricity_value[1][0].double_value, layer_input_cylindricity_value[1][1].double_value);
						geo.cylindricity.set_half_height(layer_input_cylindricity_value[1][2].double_value);
						double radian[] = { AngleToRadians(layer_input_cylindricity_value[2][0].double_value), AngleToRadians(layer_input_cylindricity_value[2][1].double_value), AngleToRadians(layer_input_cylindricity_value[2][2].double_value) };
						geo.cylindricity.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
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
			if (Solvers::get_instance()->parameters.is_Normalize_Phi && Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise)
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