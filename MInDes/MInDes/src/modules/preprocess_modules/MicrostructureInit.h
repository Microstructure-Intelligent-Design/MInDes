#pragma once
#include "../input_modules/ioFiles_Params.h"
#include "../input_modules/inputfiles/selectFile.h"
#include "../model_modules/Model_Params.h"
#include "../Module.h"
#include "WriteMeshData.h"
#include "MicrostructureInit/GeometryStructure.h"
#include "MicrostructureInit/Bmp24Structure.h"
#include "MicrostructureInit/PorousStructure.h"
#include "MicrostructureInit/VoronoiStructure.h"
namespace pf {
	namespace microstructure_init {
		inline std::vector<size_t> porous_phis_indexs;
		inline std::vector<size_t> voronoi_phis_indexs;
		namespace functions {
			inline void init_mesh_with_datafile(pf::Data_MeshInfo& report, std::string dat_path) {
				if (!write_mesh_data::read_dataFile(dat_path, report)) {
					string _report = "> Error : datafile_path can't be opened ! \n";
					WriteLog(_report);
					std::exit(0);
				}
				else {
					// check datafile with this input file
					vector<string> bool_type; bool_type.push_back("MISMATCH"); bool_type.push_back("MATCH");
					stringstream _report;
					_report << "> " << std::endl;
					_report << "> | init microstructure with datafile :                    (check in this simulation)" << std::endl;
					bool all_valid = true;
					if (report.is_phi_mesh) {
						bool line_valid = true;
						if (report.PNx == main_field::phase_field.Nx()
							&& report.PNy == main_field::phase_field.Ny()
							&& report.PNz == main_field::phase_field.Nz()
							&& report.phi_number == main_field::phi_number) {
							line_valid = true;
						}
						else {
							line_valid = false;
							all_valid = false;
						}
						_report << "> | phase-field size               : Nx - " << report.PNx - 2
							    << ", Ny - " << report.PNy - 2 << ", Nz - " << report.PNz - 2
							    << ", phi number - " << report.phi_number << " (" << bool_type[line_valid] << ")" << std::endl;
					}
					if (report.is_con_mesh) {
						bool line_valid = true;
						if (report.CNx == main_field::concentration_field.Nx()
							&& report.CNy == main_field::concentration_field.Ny()
							&& report.CNz == main_field::concentration_field.Nz()
							&& report.con_number == main_field::con_number) {
							line_valid = true;
						}
						else {
							line_valid = false;
							all_valid = false;
						}
						_report << "> | concentration-field size               : Nx - " << report.CNx - 2
							<< ", Ny - " << report.CNy - 2 << ", Nz - " << report.CNz - 2
							<< ", con number - " << report.con_number << " (" << bool_type[line_valid] << ")" << std::endl;
					}
					if (report.is_temp_mesh) {
						bool line_valid = true;
						if (report.TNx == main_field::temperature_field.Nx()
							&& report.TNy == main_field::temperature_field.Ny()
							&& report.TNz == main_field::temperature_field.Nz()) {
							line_valid = true;
						}
						else {
							line_valid = false;
							all_valid = false;
						}
						_report << "> | temperature-field size               : Nx - " << report.PNx - 2
							<< ", Ny - " << report.PNy - 2 << ", Nz - " << report.PNz - 2 << " (" << bool_type[line_valid] << ")" << std::endl;
					}
					// if (report.is_fluid_mesh) {
					// 	
					// }
					// _report << ">" << std::endl;
					_report << std::endl;
					WriteLog(_report.str());
					if (all_valid == false) {
						string error_report = "> Wainning : mesh data structure from datafile and inputfile mismatch ! \n";
						WriteLog(error_report);
					}
				}
			}
			void definiteNucleation(size_t Current_ITE_step) {
				bool _creatNewStorage = true;
				for (auto geo = pf::geometry_structure::nucleation_box.geometry_box.begin(); geo < pf::geometry_structure::nucleation_box.geometry_box.end();) {
					if (geo->generate_step == Current_ITE_step) {
						if (geo->isNormalized) {
							if (geo->phi > 1.0)
								geo->phi = 1.0;
							else if (geo->phi < 0.0)
								geo->phi = 0.0;
						}
						if (geo->geometryProperty == pf::geometry_structure::Geometry::Geo_Ellipsoid) {
							if (pf::main_field::is_phi_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::phase_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::phase_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::phase_field.Nx(); x++) {
											std::vector<REAL>& phi = pf::main_field::phase_field(x, y, z);
											// -
											pf::Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
											p.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->ellipSolid.radian_x, 
												geo->ellipSolid.radian_y, geo->ellipSolid.radian_z), pf::RotationGauge(geo->ellipSolid.rotationGauge)));
											p[0] += geo->ellipSolid.core.x;
											p[1] += geo->ellipSolid.core.y;
											p[2] += geo->ellipSolid.core.z;
											pf::geometry_structure::Point po(p[0], p[1], p[2]);
											bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
											if (!geo->isReverseRegion && check0) {
												if (geo->isNormalized) {
													REAL sum_phis = 0.0;
													for (int index = 0; index < pf::main_field::phi_number; index++)
														if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															sum_phis += phi[index];
													if (sum_phis > SYS_EPSILON) {
														for (int index = 0; index < pf::main_field::phi_number; index++) {
															if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
																phi[index] *= (1 - geo->phi) / sum_phis;
														}
													}
												}
												phi[geo->phaseIndex] = geo->phi;
											}
											else if (geo->isReverseRegion && !check0) {
												if (geo->isNormalized) {
													REAL sum_phis = 0.0;
													for (int index = 0; index < pf::main_field::phi_number; index++)
														if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															sum_phis += phi[index];
													if (sum_phis > SYS_EPSILON) {
														for (int index = 0; index < pf::main_field::phi_number; index++) {
															if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
																phi[index] *= (1 - geo->phi) / sum_phis;
														}
													}
												}
												phi[geo->phaseIndex] = geo->phi;
											}
										}
							}
							if (pf::main_field::is_con_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::concentration_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::concentration_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::concentration_field.Nx(); x++) {
											std::vector<REAL>& con = pf::main_field::concentration_field(x, y, z);
											// -
											pf::Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
											p.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->ellipSolid.radian_x, geo->ellipSolid.radian_y, geo->ellipSolid.radian_z),
												pf::RotationGauge(geo->ellipSolid.rotationGauge)));
											p[0] += geo->ellipSolid.core.x;
											p[1] += geo->ellipSolid.core.y;
											p[2] += geo->ellipSolid.core.z;
											pf::geometry_structure::Point po(p[0], p[1], p[2]);
											bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
											if (!geo->isReverseRegion && check0) {
												con = geo->con;
											}
											else if (geo->isReverseRegion && !check0) {
												con = geo->con;
											}
										}
							}
							if (pf::main_field::is_temp_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::temperature_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::temperature_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::temperature_field.Nx(); x++) {
											REAL& temp = pf::main_field::temperature_field(x, y, z);
											// -
											pf::Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
											p.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->ellipSolid.radian_x, geo->ellipSolid.radian_y, geo->ellipSolid.radian_z),
												pf::RotationGauge(geo->ellipSolid.rotationGauge)));
											p[0] += geo->ellipSolid.core.x;
											p[1] += geo->ellipSolid.core.y;
											p[2] += geo->ellipSolid.core.z;
											pf::geometry_structure::Point po(p[0], p[1], p[2]);
											bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
											if (!geo->isReverseRegion && check0) {
												temp = geo->temperature;
											}
											else if (geo->isReverseRegion && !check0) {
												temp = geo->temperature;
											}
										}
							}
							std::stringstream report;
							report << "> A new Ellipsoid for grain : " << std::to_string(geo->phaseIndex) << " has been initialized at position ( "
								<< geo->ellipSolid.core.x << ", " << geo->ellipSolid.core.y << ", " << geo->ellipSolid.core.z << " )." << std::endl;
							WriteLog(report.str());
						}
						else if (geo->geometryProperty == pf::geometry_structure::Geometry::Geo_Polyhedron) {
							if (pf::main_field::is_phi_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::phase_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::phase_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::phase_field.Nx(); x++) {
											std::vector<REAL>& phi = pf::main_field::phase_field(x, y, z);
											// -
											pf::Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
											pv.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
												pf::RotationGauge(geo->polyhedron.rotationGauge)));
											pv[0] += geo->polyhedron.point_inside_polyhedron.x;
											pv[1] += geo->polyhedron.point_inside_polyhedron.y;
											pv[2] += geo->polyhedron.point_inside_polyhedron.z;
											pf::geometry_structure::Point p(pv[0], pv[1], pv[2]);
											//Point p(x, y, z);
											bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
											if (!geo->isReverseRegion && check0) {
												if (geo->isNormalized) {
													REAL sum_phis = 0.0;
													for (int index = 0; index < pf::main_field::phi_number; index++)
														if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															sum_phis += phi[index];
													if (sum_phis > SYS_EPSILON) {
														for (int index = 0; index < pf::main_field::phi_number; index++) {
															if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
																phi[index] *= (1 - geo->phi) / sum_phis;
														}
													}
												}
												phi[geo->phaseIndex] = geo->phi;
											}
											else if (geo->isReverseRegion && !check0) {
												if (geo->isNormalized) {
													REAL sum_phis = 0.0;
													for (int index = 0; index < pf::main_field::phi_number; index++)
														if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															sum_phis += phi[index];
													if (sum_phis > SYS_EPSILON) {
														for (int index = 0; index < pf::main_field::phi_number; index++) {
															if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
																phi[index] *= (1 - geo->phi) / sum_phis;
														}
													}
												}
												phi[geo->phaseIndex] = geo->phi;
											}
										};
							}
							if (pf::main_field::is_con_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::concentration_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::concentration_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::concentration_field.Nx(); x++) {
											std::vector<REAL>& con = pf::main_field::concentration_field(x, y, z);
											pf::Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
											pv.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
												pf::RotationGauge(geo->polyhedron.rotationGauge)));
											pv[0] += geo->polyhedron.point_inside_polyhedron.x;
											pv[1] += geo->polyhedron.point_inside_polyhedron.y;
											pv[2] += geo->polyhedron.point_inside_polyhedron.z;
											pf::geometry_structure::Point p(pv[0], pv[1], pv[2]);
											//Point p(x, y, z);
											bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
											if (!geo->isReverseRegion && check0) {
												con = geo->con;
											}
											else if (geo->isReverseRegion && !check0) {
												con = geo->con;
											}
										};
							}
							if (pf::main_field::is_temp_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::temperature_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::temperature_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::temperature_field.Nx(); x++) {
											REAL& temp = pf::main_field::temperature_field(x, y, z);
											pf::Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
											pv.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
												pf::RotationGauge(geo->polyhedron.rotationGauge)));
											pv[0] += geo->polyhedron.point_inside_polyhedron.x;
											pv[1] += geo->polyhedron.point_inside_polyhedron.y;
											pv[2] += geo->polyhedron.point_inside_polyhedron.z;
											pf::geometry_structure::Point p(pv[0], pv[1], pv[2]);
											//Point p(x, y, z);
											/*p.do_boundary(phaseMesh._bc_x_up, phaseMesh._bc_y_up, phaseMesh._bc_z_up, phaseMesh._bc_x_down, phaseMesh._bc_y_down, phaseMesh._bc_z_down,
												phaseMesh.limit_x, phaseMesh.limit_y, phaseMesh.limit_z);*/
											bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
											if (!geo->isReverseRegion && check0) {
												temp = geo->temperature;
											}
											else if (geo->isReverseRegion && !check0) {
												temp = geo->temperature;
											}
										};
							}
							std::stringstream report;
							report << "> A new Polygon for grain : " << std::to_string(geo->phaseIndex) << " has been initialized at position ( "
								<< geo->polyhedron.point_inside_polyhedron.x << ", " << geo->polyhedron.point_inside_polyhedron.y << ", " << geo->polyhedron.point_inside_polyhedron.z << " )." << std::endl;
							WriteLog(report.str());
						}
						else if (geo->geometryProperty == pf::geometry_structure::Geometry::Geo_SegmentedCylinder) {
							if (pf::main_field::is_phi_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::phase_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::phase_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::phase_field.Nx(); x++) {
											std::vector<REAL>& phi = pf::main_field::phase_field(x, y, z);
											// Vector3 p(x, y, z);
											pf::Vector3 p(x - geo->cylinder.geometric_center.x, y - geo->cylinder.geometric_center.y, z - geo->cylinder.geometric_center.z);
											p.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->cylinder.radian_x, geo->cylinder.radian_y, geo->cylinder.radian_z),
												pf::RotationGauge(geo->cylinder.rotationGauge)));
											p[0] += geo->cylinder.geometric_center.x;
											p[1] += geo->cylinder.geometric_center.y;
											p[2] += geo->cylinder.geometric_center.z;
											pf::geometry_structure::Point po(p[0], p[1], p[2]);
											bool check0 = geo->cylinder.check_point_inside_segmented_cylinder(po);
											if (!geo->isReverseRegion && check0) {
												if (geo->isNormalized) {
													REAL sum_phis = 0.0;
													for (int index = 0; index < pf::main_field::phi_number; index++)
														if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															sum_phis += phi[index];
													if (sum_phis > SYS_EPSILON) {
														for (int index = 0; index < pf::main_field::phi_number; index++) {
															if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
																phi[index] *= (1 - geo->phi) / sum_phis;
														}
													}
												}
												phi[geo->phaseIndex] = geo->phi;
											}
											else if (geo->isReverseRegion && !check0) {
												if (geo->isNormalized) {
													REAL sum_phis = 0.0;
													for (int index = 0; index < pf::main_field::phi_number; index++)
														if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															sum_phis += phi[index];
													if (sum_phis > SYS_EPSILON) {
														for (int index = 0; index < pf::main_field::phi_number; index++) {
															if (phi[index] > SYS_EPSILON && index != geo->phaseIndex)
																phi[index] *= (1 - geo->phi) / sum_phis;
														}
													}
												}
												phi[geo->phaseIndex] = geo->phi;
											}
										};
							}
							if (pf::main_field::is_con_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::concentration_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::concentration_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::concentration_field.Nx(); x++) {
											std::vector<REAL>& con = pf::main_field::concentration_field(x, y, z);
											// Vector3 p(x, y, z);
											pf::Vector3 p(x - geo->cylinder.geometric_center.x, y - geo->cylinder.geometric_center.y, z - geo->cylinder.geometric_center.z);
											p.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->cylinder.radian_x, geo->cylinder.radian_y, geo->cylinder.radian_z),
												pf::RotationGauge(geo->cylinder.rotationGauge)));
											p[0] += geo->cylinder.geometric_center.x;
											p[1] += geo->cylinder.geometric_center.y;
											p[2] += geo->cylinder.geometric_center.z;
											pf::geometry_structure::Point po(p[0], p[1], p[2]);
											bool check0 = geo->cylinder.check_point_inside_segmented_cylinder(po);
											if (!geo->isReverseRegion && check0) {
												con = geo->con;
											}
											else if (geo->isReverseRegion && !check0) {
												con = geo->con;
											}
										};
							}
							if (pf::main_field::is_temp_field_on) {
#pragma omp parallel for
								for (long long z = 0; z < pf::main_field::temperature_field.Nz(); z++)
									for (long long y = 0; y < pf::main_field::temperature_field.Ny(); y++)
										for (long long x = 0; x < pf::main_field::temperature_field.Nx(); x++) {
											REAL& temp = pf::main_field::temperature_field(x, y, z);
											// Vector3 p(x, y, z);
											pf::Vector3 p(x - geo->cylinder.geometric_center.x, y - geo->cylinder.geometric_center.y, z - geo->cylinder.geometric_center.z);
											p.do_rotate(pf::RotationMatrix::rotationMatrix(pf::Vector3(geo->cylinder.radian_x, geo->cylinder.radian_y, geo->cylinder.radian_z),
												pf::RotationGauge(geo->cylinder.rotationGauge)));
											p[0] += geo->cylinder.geometric_center.x;
											p[1] += geo->cylinder.geometric_center.y;
											p[2] += geo->cylinder.geometric_center.z;
											pf::geometry_structure::Point po(p[0], p[1], p[2]);
											bool check0 = geo->cylinder.check_point_inside_segmented_cylinder(po);
											if (!geo->isReverseRegion && check0) {
												temp = geo->temperature;
											}
											else if (geo->isReverseRegion && !check0) {
												temp = geo->temperature;
											}
										}
							}
							std::stringstream report;
							report << "> A new Segmented Cylinder for grain : " << std::to_string(geo->phaseIndex) << " has been initialized at position ( "
								<< geo->cylinder.geometric_center.x << ", " << geo->cylinder.geometric_center.y << ", " << geo->cylinder.geometric_center.z << " )." << std::endl;
							WriteLog(report.str());
						}
						geo = pf::geometry_structure::nucleation_box.geometry_box.erase(geo);
					}
					else {
						geo++;
					}
				}
				for (auto point_set = pf::geometry_structure::nucleation_box.point_set_box.begin(); 
					point_set < pf::geometry_structure::nucleation_box.point_set_box.end();) {
					if (point_set->generate_step == Current_ITE_step) {
						if (point_set->is_normalized) {
#pragma omp parallel for
							for (int index = 0; index < point_set->points_phi.size(); index++) {
								if (point_set->points_phi[index] > 1.0)
									point_set->points_phi[index] = 1.0;
								else if (point_set->points_phi[index] < 0.0)
									point_set->points_phi[index] = 0.0;
							}
						}
						if (pf::main_field::is_phi_field_on) {
#pragma omp parallel for
							for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
								std::vector<REAL>& phi = pf::main_field::phase_field(pf::REAL_to_int(point_set->points[point_index].x),
									pf::REAL_to_int(point_set->points[point_index].y), pf::REAL_to_int(point_set->points[point_index].z));
								if (point_set->is_normalized) {
									REAL sum_phis = 0.0;
									for (int index = 0; index < pf::main_field::phi_number; index++)
										if (phi[index] > SYS_EPSILON && index != point_set->phaseIndex)
											sum_phis += phi[index];
									if (sum_phis > SYS_EPSILON) {
										for (int index = 0; index < pf::main_field::phi_number; index++) {
											if (phi[index] > SYS_EPSILON && index != point_set->phaseIndex)
												phi[index] *= (1 - point_set->points_phi[point_index]) / sum_phis;
										}
									}
								}
								phi[point_set->phaseIndex] = point_set->points_phi[point_index];
							}
						}
						if (pf::main_field::is_con_field_on) {
#pragma omp parallel for
							for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
								std::vector<REAL>& con = pf::main_field::concentration_field(pf::REAL_to_int(point_set->points[point_index].x),
									pf::REAL_to_int(point_set->points[point_index].y), pf::REAL_to_int(point_set->points[point_index].z));
								con = point_set->con;
							}
						}
						if (pf::main_field::is_temp_field_on) {
#pragma omp parallel for
							for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
								REAL& temp = pf::main_field::temperature_field(pf::REAL_to_int(point_set->points[point_index].x),
									pf::REAL_to_int(point_set->points[point_index].y), pf::REAL_to_int(point_set->points[point_index].z));
								temp = point_set->temperature;
							}
						}
						std::stringstream report;
						report << "> A new PointSet for grain : " << std::to_string(point_set->phaseIndex) << " has been initialized;" << std::endl;
						WriteLog(report.str());
						point_set = pf::geometry_structure::nucleation_box.point_set_box.erase(point_set);
					}
					else {
						point_set++;
					}
				}
			}
		}
		inline void check_phi_index(size_t phi_index) {
			if (phi_index >= main_field::phi_number) {
				std::string _report = "> ERROR : Phi index = " + to_string(phi_index) + " should smaller than phi number !\n";
				WriteLog(_report);
				SYS_PROGRAM_STOP;
			}
		}
		inline void check_con_index(size_t con_index) {
			if (con_index >= main_field::con_number) {
				std::string _report = "> ERROR : Con index = " + to_string(con_index) + " should smaller than con number !\n";
				WriteLog(_report);
				SYS_PROGRAM_STOP;
			}
		}
		inline void check_con_size(size_t con_size) {
			if (con_size != main_field::con_number) {
				std::string _report = "> ERROR : Con size = " + to_string(con_size)
					+ " is not equal to con number = " + to_string(main_field::con_number) + " !\n";
				WriteLog(_report);
				SYS_PROGRAM_STOP;
			}
		}
		// -
		inline void init_microstructure_pre_i() {
			if (write_mesh_data::is_datafile_init) {
				WriteLog("> Open datafile : " + write_mesh_data::datafile_path + "\n");
				functions::init_mesh_with_datafile(write_mesh_data::datafile_report, write_mesh_data::datafile_path);
			}
			else {
				if (main_field::is_phi_field_on) {
#pragma omp parallel for
					for (long long x = 0; x < main_field::phase_field.Nx(); x++)
						for (long long y = 0; y < main_field::phase_field.Ny(); y++)
							for (long long z = 0; z < main_field::phase_field.Nz(); z++)
								main_field::phase_field(x, y, z)[geometry_structure::matrix_phi_index] = geometry_structure::matrix_phi_value;
				}
				if (main_field::is_con_field_on) {
#pragma omp parallel for
					for (long long x = 0; x < main_field::concentration_field.Nx(); x++)
						for (long long y = 0; y < main_field::concentration_field.Ny(); y++)
							for (long long z = 0; z < main_field::concentration_field.Nz(); z++)
								main_field::concentration_field(x, y, z) = geometry_structure::matrix_con;
				}
				if (main_field::is_temp_field_on) {
#pragma omp parallel for
					for (long long x = 0; x < main_field::temperature_field.Nx(); x++)
						for (long long y = 0; y < main_field::temperature_field.Ny(); y++)
							for (long long z = 0; z < main_field::temperature_field.Nz(); z++)
								main_field::temperature_field(x, y, z) = geometry_structure::matrix_temperature;
				}
				// - normal init structure
				functions::definiteNucleation(main_iterator::Current_ITE_step);
				if (bmp24_structure::is_read_bmp24file) {
					bmp24_structure::generate_structure_from_BMP_pic(
						mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ,
						main_field::is_phi_field_on, main_field::is_con_field_on, main_field::is_temp_field_on);
					functions::definiteNucleation(main_iterator::Current_ITE_step);
				}
				// - porous
				if (porous_structure::is_porous) {
					if (porous_phis_indexs.size() != 0 && main_field::is_phi_field_on) {
						std::vector<std::vector<std::vector<float>>> aim_phi;
						aim_phi.resize(mesh_parameters::MESH_NX);
						for (size_t x = 0; x < mesh_parameters::MESH_NX; x++) {
							aim_phi[x].resize(mesh_parameters::MESH_NY);
							for (size_t y = 0; y < mesh_parameters::MESH_NY; y++) {
								aim_phi[x][y].resize(mesh_parameters::MESH_NZ);
								for (size_t z = 0; z < mesh_parameters::MESH_NZ; z++) {
									std::vector<REAL>& phi = main_field::phase_field(x + 1, y + 1, z + 1);
									aim_phi[x][y][z] = 0.0;
									for (size_t index = 0; index < main_field::phi_number; index++)
										for (size_t phi_index : porous_phis_indexs)
											if (index == phi_index && phi[index] > 1e-6) {
												aim_phi[x][y][z] += phi[index];
												phi[index] = 0.0;
											}
								}
							}
						}
						porous_structure::quartet_structure_generation_in_phis(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, aim_phi);
					}
					else {
						porous_structure::quartet_structure_generation(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ);
					}
					functions::definiteNucleation(main_iterator::Current_ITE_step);
				}
				// - voronoi
				if (voronoi_structure::is_voronoi) {
					if (voronoi_phis_indexs.size() != 0 && main_field::is_phi_field_on) {
						std::vector<std::vector<std::vector<float>>> aim_phi;
						aim_phi.resize(main_field::phase_field.Nx());
						for (int x = 0; x < main_field::phase_field.Nx(); x++) {
							aim_phi[x].resize(main_field::phase_field.Ny());
							for (int y = 0; y < main_field::phase_field.Ny(); y++) {
								aim_phi[x][y].resize(main_field::phase_field.Nz(), 0);
								for (int z = 0; z < main_field::phase_field.Nz(); z++) {
									std::vector<REAL>& phi = main_field::phase_field(x, y, z);
									aim_phi[x][y][z] = 0.0;
									for (size_t index = 0; index < main_field::phi_number; index++)
										for (size_t phi_index : voronoi_phis_indexs)
											if (index == phi_index && phi[index] > 1e-6) {
												aim_phi[x][y][z] += phi[index];
												phi[index] = 0.0;
											}
								}
							}
						}
						voronoi_structure::generate_voronoi_structure_in_phis(aim_phi);
					}
					else {
						voronoi_structure::generate_voronoi_structure();
					}
					functions::definiteNucleation(main_iterator::Current_ITE_step);
				}
				// - others
			}
		}
		// - 
		void write_data_pre_iii() {
			write_mesh_data::write_dataFile(input_output_files_parameters::WorkingFolder_Path + dirSeparator + write_mesh_data::mainName + "_step0" + write_mesh_data::format);
		}
		void write_data_pos_iii() {
			if (write_mesh_data::output_frequence == 0)
				return;
			if (main_iterator::Current_ITE_step % write_mesh_data::output_frequence != 0)
				return;
			write_mesh_data::write_dataFile(input_output_files_parameters::WorkingFolder_Path + dirSeparator + write_mesh_data::mainName + "_step" + to_string(main_iterator::Current_ITE_step) + write_mesh_data::format);
		}
		inline void dinit() {
			geometry_structure::nucleation_box.geometry_box.clear();
			geometry_structure::nucleation_box.point_set_box.clear();
		}
		// -
		inline void init() {
			if (main_field::is_phi_field_on || main_field::is_con_field_on || main_field::is_temp_field_on) {
				// - datafile module init
				InputFileReader::get_instance()->read_bool_value("Preprocess.Microstructure.is_datafile_init", write_mesh_data::is_datafile_init, true);
				if (write_mesh_data::is_datafile_init) {
					WriteDebugFile("# Preprocess.Microstructure.datafile_path : relative path from infile folder.\n");
					if (InputFileReader::get_instance()->read_string_value("Preprocess.Microstructure.datafile_path", write_mesh_data::datafile_path, true))
						write_mesh_data::is_read_datafile_by_path = true;
					write_mesh_data::datafile_path = input_output_files_parameters::InFileFolder_Path + dirSeparator + write_mesh_data::datafile_path;
#ifdef _WIN32
					if (!write_mesh_data::is_read_datafile_by_path) {
						WriteLog("> Please select a datafile (in .dat format) to initialize the simulation mesh ...");
						selectDataFile(write_mesh_data::datafile_path);
					}
#else
					if (!write_mesh_data::is_read_datafile_by_path) {
						WriteLog("> Error : Preprocess.Microstructure.datafile_path should be defined ! \n");
						std::exit(0);
					}
#endif
				}
				else {
					// - init matrix 
					if (main_field::is_phi_field_on) {
						string matrix_key = "Preprocess.Microstructure.matrix_phi", matrix_string = "(0,1)";
						WriteDebugFile("# .matrix_phi = ( phi_index, phi_value ) \n");
						if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_string, true)) {
							vector<InputValueType> matrix_structure;
							matrix_structure.push_back(InputValueType::IVType_INT);
							matrix_structure.push_back(InputValueType::IVType_REAL);
							vector<input_value> matrix_value =
								InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(matrix_structure, matrix_key, matrix_string, true);
							geometry_structure::matrix_phi_index = matrix_value[0].int_value;
							geometry_structure::matrix_phi_value = matrix_value[1].REAL_value;
							check_phi_index(geometry_structure::matrix_phi_index);
						}
					}
					if (main_field::is_con_field_on) {
						string matrix_key = "Preprocess.Microstructure.matrix_con", matrix_string = "()";
						WriteDebugFile("# .matrix_con = ( con_0_value, con_1_value, ... ) \n");
						geometry_structure::matrix_con.resize(main_field::con_number, 0);
						if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_string, true)) {
							vector<input_value> matrix_value =
								InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, matrix_key, matrix_string, true);
							for (int index = 0; index < matrix_value.size(); index++)
								geometry_structure::matrix_con[index] = matrix_value[index].REAL_value;
							check_con_size(int(geometry_structure::matrix_con.size()));
						}
					}
					if (main_field::is_temp_field_on) {
						string matrix_key = "Preprocess.Microstructure.matrix_temperature";
						InputFileReader::get_instance()->read_REAL_value(matrix_key, geometry_structure::matrix_temperature, true);
					}
				}
				// - init geometry structure
				// geometry structure
				int nucleation_layer = 0;
				vector<int> geometry_layer;
				InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.geometry_layer_number", nucleation_layer, true);
				for (int layer_index = 0; layer_index < nucleation_layer; layer_index++) {
					geometry_structure::GeometricRegion geo;
					if (layer_index == 0) {
						WriteDebugFile("# .property = (phi_index, geometry_type, rotation_gauge, reverse_region) \n");
						WriteDebugFile("#              geometry_type  : 0 - None, 1 - Ellipsoid, 2 - Polyhedron, 3 - Geo_SegmentedCylinder, 4 - RectangularCuboid \n");
						WriteDebugFile("#              rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
						WriteDebugFile("#                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
						WriteDebugFile("#  The compute box of the simulation mesh is [( 1 , 1 , 1 ) , ( " 
							+ to_string(mesh_parameters::MESH_NX) + " , " + to_string(mesh_parameters::MESH_NY) 
							+ " , " + to_string(mesh_parameters::MESH_NZ) + " )] \n");
						WriteDebugFile("#  The boundary of the simulation mesh is x = 0 and "
							+ to_string(mesh_parameters::MESH_NX + 1) + " , y = 0 and " + to_string(mesh_parameters::MESH_NY + 1)
							+ " , z = 0 and " + to_string(mesh_parameters::MESH_NZ + 1) + " \n");
					}
					string layer_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".property";
					string layer_input_property = "(0,0,1,false)";
					if (InputFileReader::get_instance()->read_string_value(layer_key, layer_input_property, true)) {
						vector<InputValueType> layer_input_property_structure; layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_INT);
						layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_BOOL);
						vector<input_value> layer_input_property_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(layer_input_property_structure, layer_key, layer_input_property, true);
						check_phi_index(layer_input_property_value[0].int_value);
						geo.init(geometry_structure::Geometry(layer_input_property_value[1].int_value), 0, layer_input_property_value[0].int_value, layer_input_property_value[3].bool_value);
						RotationGauge rotation_gauge = RotationGauge(layer_input_property_value[2].int_value);
						if (geo.geometryProperty == geometry_structure::Geometry::Geo_Ellipsoid) {
							WriteDebugFile("# .ellipsoid = [(core_x,core_y,core_z),(radius_x,radius_y,radius_z),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n");
							string layer_ellipsoid_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".ellipsoid";
							string layer_input_ellipsoid = "[(0,0,0),(0,0,0),(0,0,0)]";
							InputFileReader::get_instance()->read_string_value(layer_ellipsoid_key, layer_input_ellipsoid, true);
							vector<vector<input_value>> layer_input_ellipsoid_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, layer_ellipsoid_key, layer_input_ellipsoid, true);
							geo.ellipSolid.set_core(layer_input_ellipsoid_value[0][0].REAL_value, layer_input_ellipsoid_value[0][1].REAL_value, layer_input_ellipsoid_value[0][2].REAL_value);
							geo.ellipSolid.set_radius(layer_input_ellipsoid_value[1][0].REAL_value, layer_input_ellipsoid_value[1][1].REAL_value, layer_input_ellipsoid_value[1][2].REAL_value);
							REAL radian[] = { AngleToRadians(layer_input_ellipsoid_value[2][0].REAL_value), AngleToRadians(layer_input_ellipsoid_value[2][1].REAL_value), AngleToRadians(layer_input_ellipsoid_value[2][2].REAL_value) };
							geo.ellipSolid.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
						}
						if (geo.geometryProperty == geometry_structure::Geometry::Geo_Polyhedron) {
							WriteDebugFile("# .polyhedron = {[inside_point],[surf_point,surf_point,surf_point], .... ,[(rotation_angle_1,rotation_angle_2,rotation_angle_3)]} \n");
							WriteDebugFile("#                surf_point = (position_x,position_y,position_z) \n");
							string layer_polyhedron_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".polyhedron";
							string layer_input_polyhedron = "{[(-1,0,0)],[(0,0,0),(0,1,0),(0,0,1)],[(0,0,0)]}";
							InputFileReader::get_instance()->read_string_value(layer_polyhedron_key, layer_input_polyhedron, true);
							vector<vector<vector<input_value>>> layer_input_polyhedron_value = InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_REAL, layer_polyhedron_key, layer_input_polyhedron, true);
							geo.polyhedron.set_a_point_inside_polyhedron(geometry_structure::Point(layer_input_polyhedron_value[0][0][0].REAL_value, layer_input_polyhedron_value[0][0][1].REAL_value, layer_input_polyhedron_value[0][0][2].REAL_value));
							for (int surf_index = 1; surf_index < layer_input_polyhedron_value.size() - 1; surf_index++) {
								geo.polyhedron.add_surf(geometry_structure::Point(layer_input_polyhedron_value[surf_index][0][0].REAL_value, layer_input_polyhedron_value[surf_index][0][1].REAL_value, layer_input_polyhedron_value[surf_index][0][2].REAL_value),
									geometry_structure::Point(layer_input_polyhedron_value[surf_index][1][0].REAL_value, layer_input_polyhedron_value[surf_index][1][1].REAL_value, layer_input_polyhedron_value[surf_index][1][2].REAL_value),
									geometry_structure::Point(layer_input_polyhedron_value[surf_index][2][0].REAL_value, layer_input_polyhedron_value[surf_index][2][1].REAL_value, layer_input_polyhedron_value[surf_index][2][2].REAL_value));
							}
							int rIndex = int(layer_input_polyhedron_value.size() - 1);
							REAL radian[] = { AngleToRadians(layer_input_polyhedron_value[rIndex][0][0].REAL_value), AngleToRadians(layer_input_polyhedron_value[rIndex][0][1].REAL_value), AngleToRadians(layer_input_polyhedron_value[rIndex][0][2].REAL_value) };
							geo.polyhedron.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
						}
						if (geo.geometryProperty == geometry_structure::Geometry::Geo_SegmentedCylinder) {
							WriteDebugFile("# .segmented_cylinder = [(radius),(central_axis_point_x,central_axis_point_y,central_axis_point_z), ... at least two point ... ,(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n");
							string layer_cylindricity_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".segmented_cylinder";
							string layer_input_cylindricity = "[(0),(0,0,0),(0,0,0),(0,0,0)]";
							InputFileReader::get_instance()->read_string_value(layer_cylindricity_key, layer_input_cylindricity, true);
							vector<vector<input_value>> layer_input_cylindricity_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, layer_cylindricity_key, layer_input_cylindricity, true);
							int vec_size = int(layer_input_cylindricity_value.size());
							if (vec_size < 4) {
								WriteDebugFile("> ERROR .segmented_cylinder : the central axis line need at least two point ! \n");
								exit(0);
							}
							geo.cylinder.set_radius(layer_input_cylindricity_value[0][0].REAL_value);
							for (int index = 1; index < vec_size - 1; index++) {
								geo.cylinder.add_point(geometry_structure::Point(layer_input_cylindricity_value[index][0].REAL_value, layer_input_cylindricity_value[index][1].REAL_value, layer_input_cylindricity_value[index][2].REAL_value));
							}
							REAL radian[] = { AngleToRadians(layer_input_cylindricity_value[vec_size - 1][0].REAL_value),
								AngleToRadians(layer_input_cylindricity_value[vec_size - 1][1].REAL_value),
								AngleToRadians(layer_input_cylindricity_value[vec_size - 1][2].REAL_value) };
							geo.cylinder.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
						}
						if (geo.geometryProperty == geometry_structure::Geometry::Geo_RectangularCuboid) {
							geo.geometryProperty = geometry_structure::Geometry::Geo_Polyhedron;
							WriteDebugFile("# .rectangular_cuboid = [(central_point),(x_length,y_length,z_length),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n");
							string layer_polyhedron_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".rectangular_cuboid";
							string layer_input_polyhedron = "[(-1,-1,-1),(1,1,1),(0,0,0)]";
							InputFileReader::get_instance()->read_string_value(layer_polyhedron_key, layer_input_polyhedron, true);
							vector<vector<input_value>> layer_input_polyhedron_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, layer_polyhedron_key, layer_input_polyhedron, true);
							geometry_structure::Point central_point(layer_input_polyhedron_value[0][0].REAL_value, layer_input_polyhedron_value[0][1].REAL_value, layer_input_polyhedron_value[0][2].REAL_value);
							geo.polyhedron.set_a_point_inside_polyhedron(central_point);

							REAL x_length = layer_input_polyhedron_value[1][0].REAL_value + SYS_EPSILON,
								y_length = layer_input_polyhedron_value[1][1].REAL_value + SYS_EPSILON,
								z_length = layer_input_polyhedron_value[1][2].REAL_value + SYS_EPSILON;
							// - 1
							geometry_structure::Point point_on_surf_1(central_point.x - x_length / 2, central_point.y, central_point.z);
							geometry_structure::Vector3 norm_1 = get_vector(central_point, point_on_surf_1);
							geo.polyhedron.add_surf(norm_1, point_on_surf_1);
							// - 2
							geometry_structure::Point point_on_surf_2(central_point.x + x_length / 2, central_point.y, central_point.z);
							geometry_structure::Vector3 norm_2 = get_vector(central_point, point_on_surf_2);
							geo.polyhedron.add_surf(norm_2, point_on_surf_2);
							// - 3
							geometry_structure::Point point_on_surf_3(central_point.x, central_point.y - y_length / 2, central_point.z);
							geometry_structure::Vector3 norm_3 = get_vector(central_point, point_on_surf_3);
							geo.polyhedron.add_surf(norm_3, point_on_surf_3);
							// - 4
							geometry_structure::Point point_on_surf_4(central_point.x, central_point.y + y_length / 2, central_point.z);
							geometry_structure::Vector3 norm_4 = get_vector(central_point, point_on_surf_4);
							geo.polyhedron.add_surf(norm_4, point_on_surf_4);
							// - 5
							geometry_structure::Point point_on_surf_5(central_point.x, central_point.y, central_point.z - z_length / 2);
							geometry_structure::Vector3 norm_5 = get_vector(central_point, point_on_surf_5);
							geo.polyhedron.add_surf(norm_5, point_on_surf_5);
							// - 6
							geometry_structure::Point point_on_surf_6(central_point.x, central_point.y, central_point.z + z_length / 2);
							geometry_structure::Vector3 norm_6 = get_vector(central_point, point_on_surf_6);
							geo.polyhedron.add_surf(norm_6, point_on_surf_6);
							// -
							REAL radian[] = { AngleToRadians(layer_input_polyhedron_value[2][0].REAL_value), AngleToRadians(layer_input_polyhedron_value[2][1].REAL_value), AngleToRadians(layer_input_polyhedron_value[2][2].REAL_value) };
							geo.polyhedron.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
						}
						if (main_field::is_phi_field_on) {
							string layer_phi_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".phi";
							InputFileReader::get_instance()->read_REAL_value(layer_phi_key, geo.phi, true);
							string layer_norm_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".is_normalized";
							InputFileReader::get_instance()->read_bool_value(layer_norm_key, geo.isNormalized, true);
						}
						if (main_field::is_con_field_on) {
							WriteDebugFile("# .con = ( con_0_value, con_1_value, ... ) \n");
							string layer_x_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".con", layer_x_input = "()";
							geo.con.resize(0);
							if (InputFileReader::get_instance()->read_string_value(layer_x_key, layer_x_input, true)) {
								vector<input_value> layer_input_x_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, layer_x_key, layer_x_input, true);
								for (int x_index = 0; x_index < layer_input_x_value.size(); x_index++)
									geo.con.push_back(layer_input_x_value[x_index].REAL_value);
								check_con_size(int(geo.con.size()));
							}
						}
						if (main_field::is_temp_field_on) {
							string layer_temp_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".temp";
							InputFileReader::get_instance()->read_REAL_value(layer_temp_key, geo.temperature, true);
						}
						geometry_layer.push_back(int(geometry_structure::nucleation_box.geometry_box.size()));
						geometry_structure::nucleation_box.geometry_box.push_back(geo);
					}
				}
				// batch geometry
				if (nucleation_layer > 0) {
					int homo_rotation_gauge = 1;
					WriteDebugFile("# .rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
					WriteDebugFile("#                   6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
					InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.GeometryLayerBatch.rotation_gauge", homo_rotation_gauge, true);
					WriteDebugFile("# .transformation = {[(trans_layer_index), (move_x, move_y, move_z), (rotation_angle_1,rotation_angle_2,rotation_angle_3), (rotation_angle_1,rotation_angle_2,rotation_angle_3,rotation_center_x,rotation_center_y,rotation_center_z)], ... } \n");
					string layer_transform_key = "Preprocess.Microstructure.GeometryLayerBatch.transformation", layer_transform_input = "{[(0),(0,0,0),(0,0,0),(0,0,0,0,0,0)]}";
					if (InputFileReader::get_instance()->read_string_value(layer_transform_key, layer_transform_input, true)) {
						vector<InputValueType> layer_transform_structure = { InputValueType::IVType_INT, InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL };
						vector<vector<vector<input_value>>> layer_transform_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(layer_transform_structure, layer_transform_key, layer_transform_input, true);
						for (int index = 0; index < layer_transform_value.size(); index++) {
							int geo_index = geometry_layer[layer_transform_value[index][0][0].int_value];
							REAL move_x = layer_transform_value[index][1][0].REAL_value, move_y = layer_transform_value[index][1][1].REAL_value, move_z = layer_transform_value[index][1][2].REAL_value,
								rotate_1 = AngleToRadians(layer_transform_value[index][2][0].REAL_value), rotate_2 = AngleToRadians(layer_transform_value[index][2][1].REAL_value), rotate_3 = AngleToRadians(layer_transform_value[index][2][2].REAL_value),
								Rotate_1 = AngleToRadians(layer_transform_value[index][3][0].REAL_value), Rotate_2 = AngleToRadians(layer_transform_value[index][3][1].REAL_value), Rotate_3 = AngleToRadians(layer_transform_value[index][3][2].REAL_value),
								Rotate_x = layer_transform_value[index][3][3].REAL_value, Rotate_y = layer_transform_value[index][3][4].REAL_value, Rotate_z = layer_transform_value[index][3][5].REAL_value;
							geometry_structure::GeometricRegion& geo = geometry_structure::nucleation_box.geometry_box[geo_index];
							if (geo.geometryProperty == geometry_structure::Geometry::Geo_Ellipsoid) {
								geo.ellipSolid.move(move_x, move_y, move_z);
								geo.ellipSolid.add_rotation_radian(rotate_1, rotate_2, rotate_3);
								Vector3 pv(geo.ellipSolid.core.x - Rotate_x, geo.ellipSolid.core.y - Rotate_y, geo.ellipSolid.core.z - Rotate_z);
								Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
								geo.ellipSolid.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
							}
							else if (geo.geometryProperty == geometry_structure::Geometry::Geo_Polyhedron) {
								geo.polyhedron.move(move_x, move_y, move_z);
								geo.polyhedron.add_rotation_radian(rotate_1, rotate_2, rotate_3);
								Vector3 pv(geo.polyhedron.point_inside_polyhedron.x - Rotate_x, geo.polyhedron.point_inside_polyhedron.y - Rotate_y, geo.polyhedron.point_inside_polyhedron.z - Rotate_z);
								Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
								geo.polyhedron.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
							}
							else if (geo.geometryProperty == geometry_structure::Geometry::Geo_SegmentedCylinder) {
								geo.cylinder.move(move_x, move_y, move_z);
								geo.cylinder.add_rotation_radian(rotate_1, rotate_2, rotate_3);
								Vector3 pv(geo.cylinder.geometric_center.x - Rotate_x, geo.cylinder.geometric_center.y - Rotate_y, geo.cylinder.geometric_center.z - Rotate_z);
								Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
								geo.cylinder.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
							}
						}
					}
					WriteDebugFile("# .copy = {[(copy_layer_index), (paste_phi_index), (move_x, move_y, move_z), (rotation_angle_1,rotation_angle_2,rotation_angle_3), (rotation_angle_1,rotation_angle_2,rotation_angle_3,rotation_center_x,rotation_center_y,rotation_center_z)], ... } \n");
					string layer_copy_key = "Preprocess.Microstructure.GeometryLayerBatch.copy", layer_copy_input = "{[(0),(0),(0,0,0),(0,0,0),(0,0,0,0,0,0)]}";
					if (InputFileReader::get_instance()->read_string_value(layer_copy_key, layer_copy_input, true)) {
						vector<InputValueType> layer_copy_structure = { InputValueType::IVType_INT, InputValueType::IVType_INT, InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL };
						vector<vector<vector<input_value>>> layer_copy_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(layer_copy_structure, layer_copy_key, layer_copy_input, true);
						for (int index = 0; index < layer_copy_value.size(); index++) {
							int geo_index = geometry_layer[layer_copy_value[index][0][0].int_value], paste_phi_index = layer_copy_value[index][1][0].int_value;
							REAL move_x = layer_copy_value[index][2][0].REAL_value, move_y = layer_copy_value[index][2][1].REAL_value, move_z = layer_copy_value[index][2][2].REAL_value,
								rotate_1 = AngleToRadians(layer_copy_value[index][3][0].REAL_value), rotate_2 = AngleToRadians(layer_copy_value[index][3][1].REAL_value), rotate_3 = AngleToRadians(layer_copy_value[index][3][2].REAL_value),
								Rotate_1 = AngleToRadians(layer_copy_value[index][4][0].REAL_value), Rotate_2 = AngleToRadians(layer_copy_value[index][4][1].REAL_value), Rotate_3 = AngleToRadians(layer_copy_value[index][4][2].REAL_value),
								Rotate_x = layer_copy_value[index][4][3].REAL_value, Rotate_y = layer_copy_value[index][4][4].REAL_value, Rotate_z = layer_copy_value[index][4][5].REAL_value;
							geometry_structure::GeometricRegion geo = geometry_structure::nucleation_box.geometry_box[geo_index];
							geo.phaseIndex = paste_phi_index;
							check_phi_index(paste_phi_index);
							if (geo.geometryProperty == geometry_structure::Geometry::Geo_Ellipsoid) {
								geo.ellipSolid.move(move_x, move_y, move_z);
								geo.ellipSolid.add_rotation_radian(rotate_1, rotate_2, rotate_3);
								Vector3 pv(geo.ellipSolid.core.x - Rotate_x, geo.ellipSolid.core.y - Rotate_y, geo.ellipSolid.core.z - Rotate_z);
								Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
								geo.ellipSolid.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
							}
							else if (geo.geometryProperty == geometry_structure::Geometry::Geo_Polyhedron) {
								geo.polyhedron.move(move_x, move_y, move_z);
								geo.polyhedron.add_rotation_radian(rotate_1, rotate_2, rotate_3);
								Vector3 pv(geo.polyhedron.point_inside_polyhedron.x - Rotate_x, geo.polyhedron.point_inside_polyhedron.y - Rotate_y, geo.polyhedron.point_inside_polyhedron.z - Rotate_z);
								Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
								geo.polyhedron.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
							}
							else if (geo.geometryProperty == geometry_structure::Geometry::Geo_SegmentedCylinder) {
								geo.cylinder.move(move_x, move_y, move_z);
								geo.cylinder.add_rotation_radian(rotate_1, rotate_2, rotate_3);
								Vector3 pv(geo.cylinder.geometric_center.x - Rotate_x, geo.cylinder.geometric_center.y - Rotate_y, geo.cylinder.geometric_center.z - Rotate_z);
								Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
								geo.cylinder.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
							}
							geometry_structure::nucleation_box.geometry_box.push_back(geo);
						}
					}
				}
				// - init bmp24 structure
				if (InputFileReader::get_instance()->read_string_value("Preprocess.Microstructure.bmp24.file_path", bmp24_structure::bmp24file_path, true)) {
					bmp24_structure::is_read_bmp24file = true;
					std::fstream bmp_reader(input_output_files_parameters::InFileFolder_Path + dirSeparator + bmp24_structure::bmp24file_path, std::ios::in);
					if (!bmp_reader.is_open()) {
						WriteLog("> Error : bmp24.file_path can not be opened ! \n");
#ifdef _WIN32
						selectAllFile(bmp24_structure::bmp24file_path);
#else
						SYS_PROGRAM_STOP;
#endif // _WIN32
					}
					else {
						bmp24_structure::bmp24file_path = input_output_files_parameters::InFileFolder_Path + dirSeparator + bmp24_structure::bmp24file_path;
					}
					InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.bmp24.layer_number", bmp24_structure::bmp24_layer, true);
					for (int bmp_layer = 0; bmp_layer < bmp24_structure::bmp24_layer; bmp_layer++) {
						// - gray_threshold
						WriteDebugFile("# .gray_threshold = (range_left, range_right) \n");
						string bmp24_threshold_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".gray_threshold", bmp24_threshold_input = "(0,0)";
						InputFileReader::get_instance()->read_string_value(bmp24_threshold_key, bmp24_threshold_input, true);
						vector<input_value> bmp24_threshold_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, bmp24_threshold_key, bmp24_threshold_input, true);
						bmp24_structure::bmp24_threshold.push_back({ bmp24_threshold_value[0].REAL_value, bmp24_threshold_value[1].REAL_value });
						// - phi
						if (main_field::is_phi_field_on) {
							WriteDebugFile("# .phi = ( phi_index, phi_value, is_normalized ) \n");
							string bmp24_phi_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".phi", bmp24_phi_input = "(0,0,false)";
							InputFileReader::get_instance()->read_string_value(bmp24_phi_key, bmp24_phi_input, true);
							vector<InputValueType> bmp24_phi_structure; bmp24_phi_structure.push_back(InputValueType::IVType_INT);
							bmp24_phi_structure.push_back(InputValueType::IVType_REAL); bmp24_phi_structure.push_back(InputValueType::IVType_BOOL);
							vector<input_value> bmp24_phi_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(bmp24_phi_structure, bmp24_phi_key, bmp24_phi_input, true);
							bmp24_structure::bmp24_phi_index.push_back(bmp24_phi_value[0].int_value);
							check_phi_index(bmp24_phi_value[0].int_value);
							bmp24_structure::bmp24_phi_value.push_back(bmp24_phi_value[1].REAL_value);
							bmp24_structure::bmp24_phi_normalized.push_back(bmp24_phi_value[2].bool_value);
						}
						// - con
						if (main_field::is_con_field_on) {
							vector<REAL> con;
							WriteDebugFile("# .con = ( con_0_value, con_1_value, ... ) \n");
							string bmp24_x_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".con", bmp24_x_input = "()";
							if (InputFileReader::get_instance()->read_string_value(bmp24_x_key, bmp24_x_input, true)) {
								vector<input_value> bmp24_x_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, bmp24_x_key, bmp24_x_input, true);
								for (int x_index = 0; x_index < bmp24_x_value.size(); x_index++)
									con.push_back(bmp24_x_value[x_index].REAL_value);
								check_con_size(int(con.size()));
								bmp24_structure::bmp24_con.push_back(con);
							}
						}
						// - temperature
						if (main_field::is_temp_field_on) {
							string bmp24_temp_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".temperature"; REAL bmp_temp = 0.0;
							InputFileReader::get_instance()->read_REAL_value(bmp24_temp_key, bmp_temp, true);
							bmp24_structure::bmp24_temperature.push_back(bmp_temp);
						}
					}
				}
				// - init porous structure
				WriteDebugFile("# Preprocess.Microstructure.porous = (first_phi_index, second_phi_index, porosity, noise_level) \n");
				string porous_key = "Preprocess.Microstructure.porous", porous_input = "()";
				if (InputFileReader::get_instance()->read_string_value(porous_key, porous_input, true)) {
					porous_structure::is_porous = true;
					vector<InputValueType> porous_structure = { InputValueType::IVType_INT ,
						InputValueType::IVType_INT ,InputValueType::IVType_REAL, InputValueType::IVType_REAL };
					vector<input_value> porous_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(porous_structure, porous_key, porous_input, true);
					porous_structure::porous_first_phi_index = porous_value[0].int_value;
					porous_structure::porous_second_phi_index = porous_value[1].int_value;
					check_phi_index(porous_structure::porous_first_phi_index);
					check_phi_index(porous_structure::porous_second_phi_index);
					porous_structure::porosity = porous_value[2].REAL_value;
					porous_structure::porous_init_noise = porous_value[3].REAL_value;
					if (main_field::is_con_field_on) {
						WriteDebugFile("# .con = [(first_con_0_value, first_con_1_value, ... ), (second_con_0_value, second_con_1_value, ... )] \n");
						string porous_x_key = "Preprocess.Microstructure.Porous.con", porous_x_input = "[()]";
						if (InputFileReader::get_instance()->read_string_value(porous_x_key, porous_x_input, true)) {
							vector<vector<input_value>> porous_x_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, porous_x_key, porous_x_input, true);
							for (int x_index = 0; x_index < porous_x_value[0].size(); x_index++)
								porous_structure::porous_first_con.push_back(porous_x_value[0][x_index].REAL_value);
							for (int x_index = 0; x_index < porous_x_value[1].size(); x_index++)
								porous_structure::porous_second_con.push_back(porous_x_value[1][x_index].REAL_value);
							check_con_size(int(porous_structure::porous_first_con.size()));
							check_con_size(int(porous_structure::porous_second_con.size()));
						}
					}
					if (main_field::is_temp_field_on) {
						WriteDebugFile("# .temperature = (first_temperature, second_temperature) \n");
						string porous_temp_key = "Preprocess.Microstructure.Porous.temperature", porous_temp_input = "(0,0)";
						InputFileReader::get_instance()->read_string_value(porous_temp_key, porous_temp_input, true);
						vector<input_value> porous_temp_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, porous_temp_key, porous_temp_input, true);
						porous_structure::porous_first_temperature = porous_temp_value[0].REAL_value;
						porous_structure::porous_second_temperature = porous_temp_value[1].REAL_value;
					}

					string porous_norm_key = "Preprocess.Microstructure.Porous.is_normalized";
					InputFileReader::get_instance()->read_bool_value(porous_norm_key, porous_structure::is_porous_normalized, true);

					if (InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.Porous.rand_seed", porous_structure::porous_rand_seed, true))
						porous_structure::is_porous_rand = false;

					WriteDebugFile("# .in_phi_indexs = ( phi_index_1, phi_index_2, ... ) \n");
					string porous_in_phis_key = "Preprocess.Microstructure.Porous.in_phi_indexs", porous_in_phis_input = "()";
					InputFileReader::get_instance()->read_string_value(porous_in_phis_key, porous_in_phis_input, true);
					vector<input_value> porous_in_phis_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, porous_in_phis_key, porous_in_phis_input, true);
					for (int in_phi_index = 0; in_phi_index < porous_in_phis_value.size(); in_phi_index++) {
						porous_phis_indexs.push_back(size_t(porous_in_phis_value[in_phi_index].int_value));
						check_phi_index(porous_in_phis_value[in_phi_index].int_value);
					}
				}
				// - init voronoi structure
				WriteDebugFile("# .voronoi = (phi_index_begin, phi_index_end) \n");
				string voronoi_property_key = "Preprocess.Microstructure.voronoi", voronoi_property_input = "(0,0)";
				if (InputFileReader::get_instance()->read_string_value(voronoi_property_key, voronoi_property_input, true)) {
					voronoi_structure::is_voronoi = true;
					vector<input_value> voronoi_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, voronoi_property_key, voronoi_property_input, true);
					voronoi_structure::voronoi_phi_index_range[0] = voronoi_value[0].int_value;
					voronoi_structure::voronoi_phi_index_range[1] = voronoi_value[1].int_value;
					check_phi_index(voronoi_structure::voronoi_phi_index_range[0]);
					check_phi_index(voronoi_structure::voronoi_phi_index_range[1]);
					voronoi_structure::voronoi_box_position[0] = 0;
					voronoi_structure::voronoi_box_position[1] = 0;
					voronoi_structure::voronoi_box_position[2] = 0;
					voronoi_structure::voronoi_box_size[0] = int(main_field::phase_field.Nx() - 1);
					voronoi_structure::voronoi_box_size[1] = int(main_field::phase_field.Ny() - 1);
					voronoi_structure::voronoi_box_size[2] = int(main_field::phase_field.Nz() - 1);

					InputFileReader::get_instance()->read_bool_value("Preprocess.Microstructure.Voronoi.is_periodic", voronoi_structure::is_voronoi_mirror_generation, true);

					WriteDebugFile("# .con = (con_0_value, con_1_value, ... )] \n");
					string voronoi_x_key = "Preprocess.Microstructure.Voronoi.con", voronoi_x_input = "()";
					InputFileReader::get_instance()->read_string_value(voronoi_x_key, voronoi_x_input, true);
					vector<input_value> voronoi_x_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, voronoi_x_key, voronoi_x_input, true);
					for (int x_index = 0; x_index < voronoi_x_value.size(); x_index++)
						voronoi_structure::voronoi_con.push_back(voronoi_x_value[x_index].REAL_value);
					check_con_size(int(voronoi_structure::voronoi_con.size()));

					string voronoi_temperature_key = "Preprocess.Microstructure.Voronoi.temperature";
					InputFileReader::get_instance()->read_REAL_value(voronoi_temperature_key, voronoi_structure::voronoi_temperature, true);
					if (InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.Voronoi.rand_seed", voronoi_structure::voronoi_rand_seed, true))
						voronoi_structure::is_voronoi_rand = false;

					WriteDebugFile("# Preprocess.Microstructure.Voronoi.const_distance = 0.0 \n");
					WriteDebugFile("#                                  .ref_dot = [(REF_DOT), (distance, min_points_distance, max_points_distance)] \n");
					WriteDebugFile("#                                  .ref_surface = [(REF_SURF_POINT_1), (REF_SURF_POINT_2), (REF_SURF_POINT_3), (distance, min_points_distance, max_points_distance)] \n");
					WriteDebugFile("#                                  .dots_matrix = [(DOT1, points_distances1), (DOT2, points_distances2), ... ] \n");
					string voronoi_const_distance_key = "Preprocess.Microstructure.Voronoi.const_distance";
					string voronoi_ref_dot_key = "Preprocess.Microstructure.Voronoi.ref_dot", voronoi_ref_dot_input = "[()]";
					string voronoi_ref_surface_key = "Preprocess.Microstructure.Voronoi.ref_surface", voronoi_ref_surface_input = "[()]";
					string voronoi_dots_matrix_key = "Preprocess.Microstructure.Voronoi.dots_matrix", voronoi_dots_matrix_input = "[()]";
					if (InputFileReader::get_instance()->read_REAL_value(voronoi_const_distance_key, voronoi_structure::voronoi_const_pointsDistance, true)) {
						voronoi_structure::voronoi_type = voronoi_structure::VoronoiType::VDT_CONST;
					}
					else if (InputFileReader::get_instance()->read_string_value(voronoi_ref_dot_key, voronoi_ref_dot_input, true)) {
						voronoi_structure::voronoi_type = voronoi_structure::VoronoiType::VDT_REF_DOT;
						vector<vector<input_value>> voronoi_ref_dot_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, voronoi_ref_dot_key, voronoi_ref_dot_input, true);
						voronoi_structure::voronoi_reference_dot.set(voronoi_ref_dot_value[0][0].REAL_value, voronoi_ref_dot_value[0][1].REAL_value, voronoi_ref_dot_value[0][2].REAL_value);
						voronoi_structure::voronoi_reference_dot_distance = voronoi_ref_dot_value[1][0].REAL_value;
						voronoi_structure::voronoi_reference_dot_min_pointsDistance = voronoi_ref_dot_value[1][1].REAL_value;
						voronoi_structure::voronoi_reference_dot_max_pointsDistance = voronoi_ref_dot_value[1][2].REAL_value;
					}
					else if (InputFileReader::get_instance()->read_string_value(voronoi_ref_surface_key, voronoi_ref_surface_input, true)) {
						voronoi_structure::voronoi_type = voronoi_structure::VoronoiType::VDT_REF_SURFACE;
						vector<vector<input_value>> voronoi_ref_surface_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, voronoi_ref_surface_key, voronoi_ref_surface_input, true);
						voronoi_structure::voronoi_reference_surface.init(voronoi_ref_surface_value[0][0].REAL_value, voronoi_ref_surface_value[0][1].REAL_value, voronoi_ref_surface_value[0][2].REAL_value,
							voronoi_ref_surface_value[1][0].REAL_value, voronoi_ref_surface_value[1][1].REAL_value, voronoi_ref_surface_value[1][2].REAL_value,
							voronoi_ref_surface_value[2][0].REAL_value, voronoi_ref_surface_value[2][1].REAL_value, voronoi_ref_surface_value[2][2].REAL_value);
						voronoi_structure::voronoi_reference_surface_distance = voronoi_ref_surface_value[3][0].REAL_value;
						voronoi_structure::voronoi_reference_surface_min_pointsDistance = voronoi_ref_surface_value[3][1].REAL_value;
						voronoi_structure::voronoi_reference_surface_max_pointsDistance = voronoi_ref_surface_value[3][2].REAL_value;
					}
					else if (InputFileReader::get_instance()->read_string_value(voronoi_dots_matrix_key, voronoi_dots_matrix_input, true)) {
						voronoi_structure::voronoi_type = voronoi_structure::VoronoiType::VDT_DOTS_MATRIX;
						vector<vector<input_value>> voronoi_dots_matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, voronoi_dots_matrix_key, voronoi_dots_matrix_input, true);
						if (voronoi_dots_matrix_value.size() < 2) {
							WriteDebugFile(" ERROR: Preprocess.Microstructure.Voronoi.dots_matrix, the number of DOTS should be larger than 1. \n");
							exit(0);
						}
						for (int index = 0; index < voronoi_dots_matrix_value.size(); index++) {
							voronoi_structure::voronoi_matrix_dots.push_back(geometry_structure::Point(voronoi_dots_matrix_value[index][0].REAL_value, voronoi_dots_matrix_value[index][1].REAL_value, voronoi_dots_matrix_value[index][2].REAL_value));
							voronoi_structure::voronoi_matrix_dots_pointsDistance.push_back(voronoi_dots_matrix_value[index][3].REAL_value);
						}
					}
					WriteDebugFile("# .in_phi_indexs = ( phi_index_1, phi_index_2, ... ) \n");
					string voronoi_in_phis_key = "Preprocess.Microstructure.Voronoi.in_phi_indexs", voronoi_in_phis_input = "()";
					InputFileReader::get_instance()->read_string_value(voronoi_in_phis_key, voronoi_in_phis_input, true);
					vector<input_value> voronoi_in_phis_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, voronoi_in_phis_key, voronoi_in_phis_input, true);
					for (int in_phi_index = 0; in_phi_index < voronoi_in_phis_value.size(); in_phi_index++) {
						voronoi_phis_indexs.push_back(size_t(voronoi_in_phis_value[in_phi_index].int_value));
						check_phi_index(voronoi_in_phis_value[in_phi_index].int_value);
					}
				}
				// - others
				load_a_new_module(init_microstructure_pre_i, default_module_function, default_module_function,
					default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, default_module_function, default_module_function);
				// - 
				if (infile_reader::read_int_value("Preprocess.Microstructure.Output.frequence", write_mesh_data::output_frequence, true)) {
					if (write_mesh_data::output_frequence == 0) {
						load_a_new_module(default_module_function, default_module_function, write_data_pre_iii,
							default_module_function, default_module_function, default_module_function,
							default_module_function, default_module_function, default_module_function, default_module_function);
					}
					else if (write_mesh_data::output_frequence > 0) {
						load_a_new_module(default_module_function, default_module_function, write_data_pre_iii,
							default_module_function, default_module_function, default_module_function,
							default_module_function, default_module_function, write_data_pos_iii, default_module_function);
					}
				}
			}
		}
	}
}