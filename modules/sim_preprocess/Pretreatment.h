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
#include"../sim_models/InterfaceEnergy.h"
#include"../sim_models/BulkEnergy.h"
#include"../sim_models/Mobility.h"
using namespace std;
namespace pf {
	enum MemoryOptimization { MO_NONE, MO_OPT, MO_RESTORE };
	namespace pretreatment {
		// for relax interface
		static bool is_relax_interface_on = false;
		static bool is_interface_movable{};
		static int relaxation_steps = 0;
		static int out_step = 100;
		static bool is_fix_phi_in_loop = false;
		// for merge phases

		// for remove nonexistent phases
		static bool is_remove_inexistent_phis = false;
		// for phi indexs re-ordering
		static bool is_phis_reordering = false;
		static int reordering_index_0 = 0;
		// filling based on phi index
		static bool is_filling_by_phi = false;
		static vector <vector<int>> target_phi;
		static vector<double_box> fill_phi_con;
		static vector<double_box> fill_total_con;
		static vector<double> fill_temperature;
		// re_construct_mesh_by_phi
		static bool is_reconstruct = false;
		static vector<vector<int>> recons_phi;
		static vector<int> recons_property;
		// auto_merge_phis

		// optimization memory
		static MemoryOptimization memory_optimization = MemoryOptimization::MO_NONE;
		vector<int> fixed_phis;

		static void reconstruct_phis(FieldStorage_forPhaseNode& phaseMesh) {
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Do reconstruct phis:\n", LOG_FILE_NAME);
			Info_Phases& phases_info = Solvers::get_instance()->parameters.Phases;
			Info_Node& comps_info = Solvers::get_instance()->parameters.Components;
			for (int index = 0; index < recons_phi.size(); index++) {
				int phi_property = recons_property[index];
				if (phaseMesh.info_node.x.size() != comps_info.size()) {
					phaseMesh.info_node.x.clear();
					for (auto comp = comps_info.begin(); comp < comps_info.end(); comp++)
						phaseMesh.info_node.x.add_con(comp->index, 0.0);
				}
				if (phaseMesh.info_node.potential.size() != comps_info.size()) {
					phaseMesh.info_node.potential.clear();
					for (auto comp = comps_info.begin(); comp < comps_info.end(); comp++)
						phaseMesh.info_node.potential.add_con(comp->index, 0.0);
				}
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					for (auto phi_index = recons_phi[index].begin(); phi_index < recons_phi[index].end(); phi_index++) {
						if (phase->index == *phi_index) {
							phase->property = phi_property;
							if (phase->x.size() != phases_info[phi_property].x.size()) {
								phase->x.clear();
								for (auto comp = phases_info[phi_property].x.begin(); comp < phases_info[phi_property].x.end(); comp++)
									phase->x.add_con(comp->index, 0.0);
							}
							if (phase->potential.size() != phases_info[phi_property].x.size()) {
								phase->potential.clear();
								for (auto comp = phases_info[phi_property].x.begin(); comp < phases_info[phi_property].x.end(); comp++)
									phase->potential.add_con(comp->index, 0.0);
							}
							continue;
						}
					}
				}
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							if (node.x.size() != comps_info.size()) {
								node.x.clear();
								for (auto comp = comps_info.begin(); comp < comps_info.end(); comp++)
									node.x.add_con(comp->index, 0.0);
							}
							if (node.potential.size() != comps_info.size()) {
								node.potential.clear();
								for (auto comp = comps_info.begin(); comp < comps_info.end(); comp++)
									node.potential.add_con(comp->index, 0.0);
							}
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								for (auto phi_index = recons_phi[index].begin(); phi_index < recons_phi[index].end(); phi_index++) {
									if (phase->index == *phi_index) {
										phase->property = phi_property;
										if (phase->x.size() != phases_info[phi_property].x.size()) {
											phase->x.clear();
											for (auto comp = phases_info[phi_property].x.begin(); comp < phases_info[phi_property].x.end(); comp++)
												phase->x.add_con(comp->index, 0.0);
										}
										if (phase->potential.size() != phases_info[phi_property].x.size()) {
											phase->potential.clear();
											for (auto comp = phases_info[phi_property].x.begin(); comp < phases_info[phi_property].x.end(); comp++)
												phase->potential.add_con(comp->index, 0.0);
										}
										continue;
									}
								}
							}
						}
				stringstream report;
				report << "> Already reconstructed phi index    = ";
				for (auto phi_index = recons_phi[index].begin(); phi_index < recons_phi[index].end(); phi_index++)
					report << *phi_index << ", ";
				report << "                        phi property = " << phases_info[phi_property].phi_name << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
			}
		}

		static void relaxation_interface(FieldStorage_forPhaseNode& phaseMesh) {
			int interface_gradient = interface_energy::get_interface_gradient(), interface_potential = interface_energy::get_interface_potential();
			for (auto phase = phaseMesh(0, 0, 0).begin(); phase < phaseMesh(0, 0, 0).end(); phase++) {
				phaseMesh.add_customFlag_to_allnodes(ExternalFields::RELAX_interface_buff + phase->index, 0);
			}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& currentNode = phaseMesh(x, y, z);
						if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
							if (Solvers::get_instance()->parameters.is_Normalize_Phi)
								currentNode.normalized_phi();
							for (auto phase = currentNode.begin(); phase < currentNode.end(); phase++) {
								phase->old_phi = phase->phi;
								phase->_flag = phaseMesh.currentFlag(currentNode, phase->index);
							}
						}
						for (auto phase = currentNode.begin(); phase < currentNode.end(); phase++) {
							if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
								Solvers::get_instance()->Phi_Solver_AC.Boundary_Condition(currentNode, *phase);
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
								Solvers::get_instance()->Phi_Solver_AC.Boundary_Condition(currentNode, *phase);
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
								Solvers::get_instance()->Phi_Solver_CH.Boundary_Condition_Phi(currentNode, *phase);
							}
							if (phase->phi > 0.5)
								currentNode.customFlags[ExternalFields::RELAX_interface_buff + phase->index] = 1; // > 0.5
							else if (phase->phi < 0.5)
								currentNode.customFlags[ExternalFields::RELAX_interface_buff + phase->index] = -1; // < 0.5
						}
					}
			stringstream log;
			log << "> Do interface relaxation:" << endl;
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(log.str(), LOG_FILE_NAME);
			for (int r_istep = 1; r_istep <= relaxation_steps; r_istep++) {
				Solvers::get_instance()->current_istep = r_istep;
				Solvers::get_instance()->real_time += Solvers::get_instance()->parameters.dt;
				double iwidth = interface_energy::get_interface_width(), dt = Solvers::get_instance()->parameters.dt;
				double init_normal[] = { 0.0, 0.0, 0.0 }, dx = phaseMesh.dr;
				DifferenceMethod diff_method = DifferenceMethod(Solvers::get_instance()->parameters.Difference_Method);
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
								for (auto phase = node.Goast_Phase.begin(); phase < node.Goast_Phase.end(); phase++)
									(*phase) = nullptr;
								node.Goast_Phase.clear();
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									//intphase and near intphase
									if (phase->_flag) {
										phase->bulk_increment = 0.0;
										phase->int_increment = 0.0;
										node.Goast_Phase.push_back(&(*phase));
										phase->phi_grad[0] = (node.get_neighbor_node(Direction::x_up)[phase->index].phi -
											node.get_neighbor_node(Direction::x_down)[phase->index].phi) / 2.0 / phaseMesh.dr;
										phase->phi_grad[1] = (node.get_neighbor_node(Direction::y_up)[phase->index].phi -
											node.get_neighbor_node(Direction::y_down)[phase->index].phi) / 2.0 / phaseMesh.dr;
										phase->phi_grad[2] = (node.get_neighbor_node(Direction::z_up)[phase->index].phi -
											node.get_neighbor_node(Direction::z_down)[phase->index].phi) / 2.0 / phaseMesh.dr;
										if (diff_method == DifferenceMethod::FIVE_POINT) {
											phase->laplacian = (node.get_neighbor_node(Direction::x_down)[phase->index].phi
												+ node.get_neighbor_node(Direction::x_up)[phase->index].phi
												+ node.get_neighbor_node(Direction::y_down)[phase->index].phi
												+ node.get_neighbor_node(Direction::y_up)[phase->index].phi
												+ node.get_neighbor_node(Direction::z_down)[phase->index].phi
												+ node.get_neighbor_node(Direction::z_up)[phase->index].phi - 6 * phase->phi) / phaseMesh.dr / phaseMesh.dr;
										}
										else if (diff_method == DifferenceMethod::NINE_POINT) {
											phase->laplacian = (4.0 * node.get_neighbor_node(Direction::x_down)[phase->index].phi
												+ 4.0 * node.get_neighbor_node(Direction::x_up)[phase->index].phi
												+ 4.0 * node.get_neighbor_node(Direction::y_down)[phase->index].phi
												+ 4.0 * node.get_neighbor_node(Direction::y_up)[phase->index].phi
												+ 4.0 * node.get_neighbor_node(Direction::z_down)[phase->index].phi
												+ 4.0 * node.get_neighbor_node(Direction::z_up)[phase->index].phi
												+ node.get_long_range_node(-1, -1, 0)[phase->index].phi
												+ node.get_long_range_node(-1, 1, 0)[phase->index].phi
												+ node.get_long_range_node(1, -1, 0)[phase->index].phi
												+ node.get_long_range_node(1, 1, 0)[phase->index].phi
												+ node.get_long_range_node(-1, 0, -1)[phase->index].phi
												+ node.get_long_range_node(-1, 0, 1)[phase->index].phi
												+ node.get_long_range_node(1, 0, -1)[phase->index].phi
												+ node.get_long_range_node(1, 0, 1)[phase->index].phi
												+ node.get_long_range_node(0, -1, -1)[phase->index].phi
												+ node.get_long_range_node(0, -1, 1)[phase->index].phi
												+ node.get_long_range_node(0, 1, -1)[phase->index].phi
												+ node.get_long_range_node(0, 1, 1)[phase->index].phi - 36 * phase->phi) / 6.0 / phaseMesh.dr / phaseMesh.dr;
										}
									}
									else {
										phase->phi_grad[0] = 0.0;
										phase->phi_grad[1] = 0.0;
										phase->phi_grad[2] = 0.0;
										phase->laplacian = 0.0;
									}
								}
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									//intphase and near intphase
									phase->bulk_increment = 0.0;
									phase->int_increment = 0.0;
									phase->phi_grad[0] = (node.get_neighbor_node(Direction::x_up)[phase->index].phi -
										node.get_neighbor_node(Direction::x_down)[phase->index].phi) / 2.0 / phaseMesh.dr;
									phase->phi_grad[1] = (node.get_neighbor_node(Direction::y_up)[phase->index].phi -
										node.get_neighbor_node(Direction::y_down)[phase->index].phi) / 2.0 / phaseMesh.dr;
									phase->phi_grad[2] = (node.get_neighbor_node(Direction::z_up)[phase->index].phi -
										node.get_neighbor_node(Direction::z_down)[phase->index].phi) / 2.0 / phaseMesh.dr;
									if (diff_method == DifferenceMethod::FIVE_POINT) {
										phase->laplacian = (node.get_neighbor_node(Direction::x_down)[phase->index].phi
											+ node.get_neighbor_node(Direction::x_up)[phase->index].phi
											+ node.get_neighbor_node(Direction::y_down)[phase->index].phi
											+ node.get_neighbor_node(Direction::y_up)[phase->index].phi
											+ node.get_neighbor_node(Direction::z_down)[phase->index].phi
											+ node.get_neighbor_node(Direction::z_up)[phase->index].phi - 6 * phase->phi) / phaseMesh.dr / phaseMesh.dr;
									}
									else if (diff_method == DifferenceMethod::NINE_POINT) {
										phase->laplacian = (4.0 * node.get_neighbor_node(Direction::x_down)[phase->index].phi
											+ 4.0 * node.get_neighbor_node(Direction::x_up)[phase->index].phi
											+ 4.0 * node.get_neighbor_node(Direction::y_down)[phase->index].phi
											+ 4.0 * node.get_neighbor_node(Direction::y_up)[phase->index].phi
											+ 4.0 * node.get_neighbor_node(Direction::z_down)[phase->index].phi
											+ 4.0 * node.get_neighbor_node(Direction::z_up)[phase->index].phi
											+ node.get_long_range_node(-1, -1, 0)[phase->index].phi
											+ node.get_long_range_node(-1, 1, 0)[phase->index].phi
											+ node.get_long_range_node(1, -1, 0)[phase->index].phi
											+ node.get_long_range_node(1, 1, 0)[phase->index].phi
											+ node.get_long_range_node(-1, 0, -1)[phase->index].phi
											+ node.get_long_range_node(-1, 0, 1)[phase->index].phi
											+ node.get_long_range_node(1, 0, -1)[phase->index].phi
											+ node.get_long_range_node(1, 0, 1)[phase->index].phi
											+ node.get_long_range_node(0, -1, -1)[phase->index].phi
											+ node.get_long_range_node(0, -1, 1)[phase->index].phi
											+ node.get_long_range_node(0, 1, -1)[phase->index].phi
											+ node.get_long_range_node(0, 1, 1)[phase->index].phi - 36 * phase->phi) / 6.0 / phaseMesh.dr / phaseMesh.dr;
									}
								}
							}
						}
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
								if (node.Goast_Phase.size() == 0)
									continue;
								//double first_term, second_term;
								double iwidth = interface_energy::get_interface_width();
								auto mobility = Solvers::get_instance()->Phi_Solver_AC.Lij;
								///< simplified version, to avoid issues
								if (interface_gradient == Int_Gradient::Steinbach_G2009 || interface_potential == Int_Potential::Steinbach_P2009) {
									for (auto alpha = node.Goast_Phase.begin(); alpha < node.Goast_Phase.end() - 1; alpha++)
										for (auto beta = alpha + 1; beta < node.Goast_Phase.end(); beta++) {
											double int_incre_b_a = 0.0;
											int_incre_b_a = mobility(node, **alpha, **beta) / iwidth * interface_energy::dfint_dphi_S2009(node, **alpha, **beta, iwidth);
											if (!is_interface_movable) {
												//< fix phi change begin
												if (node.customFlags[ExternalFields::RELAX_interface_buff + (*alpha)->index] == 0 || node.customFlags[ExternalFields::RELAX_interface_buff + (*beta)->index] == 0)
													int_incre_b_a = 0.0;
											}
											//< fix phi change end
											(*alpha)->int_increment += int_incre_b_a;
											(*beta)->int_increment -= int_incre_b_a;
#ifdef _DEBUG
											if (_isnan(int_incre_b_a)) {
												cout << "DEBUG: int_incre_b_a interface error !" << endl;
												SYS_PROGRAM_STOP;
											}
#endif
										}
								}
								else {
									for (auto alpha = node.Goast_Phase.begin(); alpha < node.Goast_Phase.end() - 1; alpha++)
										for (auto beta = alpha + 1; beta < node.Goast_Phase.end(); beta++) {
											double int_incre_b_a = mobility(node, **alpha, **beta) / iwidth
												* (interface_energy::dfint_dphi_pairwise_acc(node, **beta, iwidth) - interface_energy::dfint_dphi_pairwise_acc(node, **alpha, iwidth));
											//< fix phi change begin
											if (!is_interface_movable) {
												if (node.customFlags[ExternalFields::RELAX_interface_buff + (*alpha)->index] == 0 || node.customFlags[ExternalFields::RELAX_interface_buff + (*beta)->index] == 0)
													int_incre_b_a = 0.0;
											}
											//< fix phi change end
											(*alpha)->int_increment += int_incre_b_a;
											(*beta)->int_increment -= int_incre_b_a;
#ifdef _DEBUG
											if (_isnan(int_incre_b_a)) {
												cout << "DEBUG: int_incre_b_a interface error !" << endl;
												SYS_PROGRAM_STOP;
											}
#endif
										}
								}
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
								AllenCahnSolver& solver = Solvers::get_instance()->Phi_Solver_AC;
								solver.dfint_dphi(node, Solvers::get_instance()->parameters.is_Normalize_Phi);
								// equation //
								for (auto alpha = node.begin(); alpha < node.end(); alpha++) {
									alpha->bulk_increment = 0.0;
									if (alpha->laplacian > SYS_EPSILON || alpha->laplacian < -SYS_EPSILON)
										for (auto beta = node.begin(); beta < node.end(); beta++)
											if (beta->laplacian > SYS_EPSILON || beta->laplacian < -SYS_EPSILON) {
												alpha->bulk_increment += -solver.Lij(node, *alpha, *beta) * solver.dfbulk_dphi(node, *beta);
											}
									//< fix phi change begin
									if (!is_interface_movable) {
										if (node.customFlags[ExternalFields::RELAX_interface_buff + alpha->index] == 0) {
											alpha->int_increment = 0.0;
											alpha->bulk_increment = 0.0;
										}
									}
									//< fix phi change end
								}
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
								CahnHilliardSolver& solver = Solvers::get_instance()->Phi_Solver_CH;
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
										phase->int_increment = solver.dfint_dphi(node, *phase);
										phase->bulk_increment = solver.dfbulk_dphi(node, *phase);
										//< fix phi change begin
										if (!is_interface_movable){
											if (node.customFlags[ExternalFields::RELAX_interface_buff + phase->index] == 0) {
												phase->int_increment = 0.0;
												phase->bulk_increment = 0.0;
											}
										}
										//< fix phi change end
									}
								}
							}
						}
				double MAX_PHI_INCRE = 0.0;
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& currentNode = phaseMesh(x, y, z);
							if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
								for (auto phase = currentNode.Goast_Phase.begin(); phase < currentNode.Goast_Phase.end(); phase++) {
									double old_phi = (*phase)->phi;
									(*phase)->phi += dt * ((*phase)->int_increment + (*phase)->bulk_increment);
									Solvers::get_instance()->Phi_Solver_AC.Boundary_Condition(currentNode, **phase);
									if (Solvers::get_instance()->parameters.is_Normalize_Phi) {
										if ((*phase)->phi > (1.0 - SYS_EPSILON))
											(*phase)->phi = 1.0;
										else if ((*phase)->phi < SYS_EPSILON)
											(*phase)->phi = 0.0;
									}
									if ((*phase)->phi > 0.5 && currentNode.customFlags[ExternalFields::RELAX_interface_buff + (*phase)->index] == -1)
										currentNode.customFlags[ExternalFields::RELAX_interface_buff + (*phase)->index] = 0;
									else if ((*phase)->phi < 0.5 && currentNode.customFlags[ExternalFields::RELAX_interface_buff + (*phase)->index] == 1)
										currentNode.customFlags[ExternalFields::RELAX_interface_buff + (*phase)->index] = 0;
#ifdef _OPENMP
#pragma omp critical
#endif
									{
										if (abs(old_phi - (*phase)->phi) > MAX_PHI_INCRE)
											MAX_PHI_INCRE = abs(old_phi - (*phase)->phi);
									}
#ifdef _DEBUG
									if (_isnan((*phase)->phi)) {
										cout << "DEBUG: phase->phi error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
									// change _flag
									if ((*phase)->phi >= Phi_Num_Cut_Off && (*phase)->phi <= (1.0 - Phi_Num_Cut_Off)) {
										if ((*phase)->_flag != pf_INTERFACE) {
											int index = (*phase)->index;
											(*phase)->_flag = pf_INTERFACE;
											if (currentNode.get_neighbor_node(Direction::x_up)[index]._flag == pf_BULK)
												currentNode.get_neighbor_node(Direction::x_up)[index]._flag = pf_NEAR_INTERFACE;
											if (currentNode.get_neighbor_node(Direction::x_down)[index]._flag == pf_BULK)
												currentNode.get_neighbor_node(Direction::x_down)[index]._flag = pf_NEAR_INTERFACE;
											if (currentNode.get_neighbor_node(Direction::y_up)[index]._flag == pf_BULK)
												currentNode.get_neighbor_node(Direction::y_up)[index]._flag = pf_NEAR_INTERFACE;
											if (currentNode.get_neighbor_node(Direction::y_down)[index]._flag == pf_BULK)
												currentNode.get_neighbor_node(Direction::y_down)[index]._flag = pf_NEAR_INTERFACE;
											if (currentNode.get_neighbor_node(Direction::z_up)[index]._flag == pf_BULK)
												currentNode.get_neighbor_node(Direction::z_up)[index]._flag = pf_NEAR_INTERFACE;
											if (currentNode.get_neighbor_node(Direction::z_down)[index]._flag == pf_BULK)
												currentNode.get_neighbor_node(Direction::z_down)[index]._flag = pf_NEAR_INTERFACE;
										}
									}
									if ((*phase)->phi < Phi_Num_Cut_Off) {
										int index = (*phase)->index;
										if (currentNode.get_neighbor_node(Direction::x_up)[index].phi >= Phi_Num_Cut_Off
											|| currentNode.get_neighbor_node(Direction::x_down)[index].phi >= Phi_Num_Cut_Off
											|| currentNode.get_neighbor_node(Direction::y_up)[index].phi >= Phi_Num_Cut_Off
											|| currentNode.get_neighbor_node(Direction::y_down)[index].phi >= Phi_Num_Cut_Off
											|| currentNode.get_neighbor_node(Direction::z_down)[index].phi >= Phi_Num_Cut_Off
											|| currentNode.get_neighbor_node(Direction::z_up)[index].phi >= Phi_Num_Cut_Off)
											(*phase)->_flag = pf_NEAR_INTERFACE;
										else
											(*phase)->_flag = pf_BULK;
									}
									else if ((*phase)->phi > (1.0 - Phi_Num_Cut_Off)) {
										int index = (*phase)->index;
										if (currentNode.get_neighbor_node(Direction::x_up)[index].phi <= (1.0 - Phi_Num_Cut_Off)
											|| currentNode.get_neighbor_node(Direction::x_down)[index].phi <= (1.0 - Phi_Num_Cut_Off)
											|| currentNode.get_neighbor_node(Direction::y_up)[index].phi <= (1.0 - Phi_Num_Cut_Off)
											|| currentNode.get_neighbor_node(Direction::y_down)[index].phi <= (1.0 - Phi_Num_Cut_Off)
											|| currentNode.get_neighbor_node(Direction::z_up)[index].phi <= (1.0 - Phi_Num_Cut_Off)
											|| currentNode.get_neighbor_node(Direction::z_down)[index].phi <= (1.0 - Phi_Num_Cut_Off))
											(*phase)->_flag = pf_NEAR_INTERFACE;
										else
											(*phase)->_flag = pf_BULK;
									}
								}
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
								for (auto phase = currentNode.begin(); phase < currentNode.end(); phase++) {
									if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
										double old_phi = phase->phi;
										phase->phi += dt * (phase->int_increment + phase->bulk_increment);
										Solvers::get_instance()->Phi_Solver_AC.Boundary_Condition(currentNode, *phase);
										if (Solvers::get_instance()->parameters.is_Normalize_Phi) {
											if (phase->phi > (1.0 - SYS_EPSILON))
												phase->phi = 1.0;
											else if (phase->phi < SYS_EPSILON)
												phase->phi = 0.0;
										}
										if (phase->phi > 0.5 && currentNode.customFlags[ExternalFields::RELAX_interface_buff + phase->index] == -1)
											currentNode.customFlags[ExternalFields::RELAX_interface_buff + phase->index] = 0;
										else if (phase->phi < 0.5 && currentNode.customFlags[ExternalFields::RELAX_interface_buff + phase->index] == 1)
											currentNode.customFlags[ExternalFields::RELAX_interface_buff + phase->index] = 0;
#ifdef _OPENMP
#pragma omp critical
#endif
										{
											if (abs(old_phi - phase->phi) > MAX_PHI_INCRE)
												MAX_PHI_INCRE = abs(old_phi - phase->phi);
										}
#ifdef _DEBUG
										if (_isnan(phase->phi)) {
											cout << "DEBUG: phase->phi error !" << endl;
											SYS_PROGRAM_STOP;
										}
#endif
									}
								}
							}
							else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
								CahnHilliardSolver& solver = Solvers::get_instance()->Phi_Solver_CH;
								vec3_box grad_dfdphis;
								double_box lap_dfdphis;
								for (auto phase = currentNode.begin(); phase < currentNode.end(); phase++) {
									if (phase->laplacian > SYS_EPSILON || phase->laplacian < -SYS_EPSILON) {
										//intphase and near intphase
										PhaseEntry& phase_xdown = currentNode.get_neighbor_node(Direction::x_down)[phase->index], phase_xup = currentNode.get_neighbor_node(Direction::x_up)[phase->index],
											phase_ydown = currentNode.get_neighbor_node(Direction::y_down)[phase->index], phase_yup = currentNode.get_neighbor_node(Direction::y_up)[phase->index],
											phase_zdown = currentNode.get_neighbor_node(Direction::z_down)[phase->index], phase_zup = currentNode.get_neighbor_node(Direction::z_up)[phase->index];
										Vector3 vec((phase_xdown.int_increment + phase_xdown.bulk_increment - phase_xup.int_increment - phase_xup.bulk_increment) / 2.0 / phaseMesh.dr,
											(phase_ydown.int_increment + phase_ydown.bulk_increment - phase_yup.int_increment - phase_yup.bulk_increment) / 2.0 / phaseMesh.dr,
											(phase_zdown.int_increment + phase_zdown.bulk_increment - phase_zup.int_increment - phase_zup.bulk_increment) / 2.0 / phaseMesh.dr);
										grad_dfdphis.add_vec(phase->index, vec);
										if (diff_method == DifferenceMethod::FIVE_POINT) {
											lap_dfdphis.add_double(phase->index,
												(phase_xdown.int_increment + phase_xdown.bulk_increment
													+ phase_xup.int_increment + phase_xup.bulk_increment
													+ phase_ydown.int_increment + phase_ydown.bulk_increment
													+ phase_yup.int_increment + phase_yup.bulk_increment
													+ phase_zdown.int_increment + phase_zdown.bulk_increment
													+ phase_zup.int_increment + phase_zup.bulk_increment
													- 6 * (phase->int_increment + phase->bulk_increment)) / phaseMesh.dr / phaseMesh.dr);
										}
										else if (diff_method == DifferenceMethod::NINE_POINT) {
											PhaseEntry& phase_xdown_ydown = currentNode.get_long_range_node(-1, -1, 0)[phase->index],
												phase_xdown_yup = currentNode.get_long_range_node(-1, 1, 0)[phase->index],
												phase_xup_ydown = currentNode.get_long_range_node(1, -1, 0)[phase->index],
												phase_xup_yup = currentNode.get_long_range_node(1, 1, 0)[phase->index],
												phase_xdown_zdown = currentNode.get_long_range_node(-1, 0, -1)[phase->index],
												phase_xdown_zup = currentNode.get_long_range_node(-1, 0, 1)[phase->index],
												phase_xup_zdown = currentNode.get_long_range_node(1, 0, -1)[phase->index],
												phase_xup_zup = currentNode.get_long_range_node(1, 0, 1)[phase->index],
												phase_ydown_zdown = currentNode.get_long_range_node(0, -1, -1)[phase->index],
												phase_ydown_zup = currentNode.get_long_range_node(0, -1, 1)[phase->index],
												phase_yup_zdown = currentNode.get_long_range_node(0, 1, -1)[phase->index],
												phase_yup_zup = currentNode.get_long_range_node(0, 1, 1)[phase->index];
											lap_dfdphis.add_double(phase->index,
												(4.0 * (phase_xdown.int_increment + phase_xdown.bulk_increment)
													+ 4.0 * (phase_xup.int_increment + phase_xup.bulk_increment)
													+ 4.0 * (phase_ydown.int_increment + phase_ydown.bulk_increment)
													+ 4.0 * (phase_yup.int_increment + phase_yup.bulk_increment)
													+ 4.0 * (phase_zdown.int_increment + phase_zdown.bulk_increment)
													+ 4.0 * (phase_zup.int_increment + phase_zup.bulk_increment)
													+ phase_xdown_ydown.int_increment + phase_xdown_ydown.bulk_increment
													+ phase_xdown_yup.int_increment + phase_xdown_yup.bulk_increment
													+ phase_xup_ydown.int_increment + phase_xup_ydown.bulk_increment
													+ phase_xup_yup.int_increment + phase_xup_yup.bulk_increment
													+ phase_xdown_zdown.int_increment + phase_xdown_zdown.bulk_increment
													+ phase_xdown_zup.int_increment + phase_xdown_zup.bulk_increment
													+ phase_xup_zdown.int_increment + phase_xup_zdown.bulk_increment
													+ phase_xup_zup.int_increment + phase_xup_zup.bulk_increment
													+ phase_ydown_zdown.int_increment + phase_ydown_zdown.bulk_increment
													+ phase_ydown_zup.int_increment + phase_ydown_zup.bulk_increment
													+ phase_yup_zdown.int_increment + phase_yup_zdown.bulk_increment
													+ phase_yup_zup.int_increment + phase_yup_zup.bulk_increment
													- 36 * (phase->int_increment + phase->bulk_increment)) / 6.0 / phaseMesh.dr / phaseMesh.dr);
										}
									}
								}
								// assignment
								for (auto phi_a = currentNode.begin(); phi_a < currentNode.end(); phi_a++) {
									if (phi_a->laplacian > SYS_EPSILON || phi_a->laplacian < -SYS_EPSILON) {
										double old_phi = phi_a->phi;
										double phi_a_increment = 0.0;
										for (auto phi_b = currentNode.begin(); phi_b < currentNode.end(); phi_b++) {
											if (phi_b->laplacian > SYS_EPSILON || phi_b->laplacian < -SYS_EPSILON) {
												Vector3 vec_Mab((solver.M_ab(currentNode.get_neighbor_node(Direction::x_down), *phi_a, *phi_b) - solver.M_ab(currentNode.get_neighbor_node(Direction::x_up), *phi_a, *phi_b)) / 2.0 / phaseMesh.dr,
													(solver.M_ab(currentNode.get_neighbor_node(Direction::y_down), *phi_a, *phi_b) - solver.M_ab(currentNode.get_neighbor_node(Direction::y_up), *phi_a, *phi_b)) / 2.0 / phaseMesh.dr,
													(solver.M_ab(currentNode.get_neighbor_node(Direction::z_down), *phi_a, *phi_b) - solver.M_ab(currentNode.get_neighbor_node(Direction::z_up), *phi_a, *phi_b)) / 2.0 / phaseMesh.dr);
												phi_a_increment += vec_Mab * grad_dfdphis[phi_b->index] + solver.M_ab(currentNode, *phi_a, *phi_b) * lap_dfdphis[phi_b->index];
#ifdef _DEBUG
												if (_isnan(phi_a_increment)) {
													cout << "DEBUG: phi_a_increment interaction term error !" << endl;
													SYS_PROGRAM_STOP;
											}
#endif
										}
										}
										// source term
										phi_a_increment += solver.Source_a(currentNode, *phi_a);
#ifdef _DEBUG
										if (_isnan(phi_a_increment)) {
											cout << "DEBUG: phi_a_increment source term error !" << endl;
											SYS_PROGRAM_STOP;
									}
#endif

										phi_a->phi += phi_a_increment * dt;
										solver.Boundary_Condition_Phi(currentNode, *phi_a);
										if (Solvers::get_instance()->parameters.is_Normalize_Phi) {
											if (phi_a->phi > (1.0 - SYS_EPSILON))
												phi_a->phi = 1.0;
											else if (phi_a->phi < SYS_EPSILON)
												phi_a->phi = 0.0;
										}
										if (phi_a->phi > 0.5 && currentNode.customFlags[ExternalFields::RELAX_interface_buff + phi_a->index] == -1)
											currentNode.customFlags[ExternalFields::RELAX_interface_buff + phi_a->index] = 0;
										else if (phi_a->phi < 0.5 && currentNode.customFlags[ExternalFields::RELAX_interface_buff + phi_a->index] == 1)
											currentNode.customFlags[ExternalFields::RELAX_interface_buff + phi_a->index] = 0;
#ifdef _OPENMP
#pragma omp critical
#endif
										{
											if (abs(old_phi - phi_a->phi) > MAX_PHI_INCRE) // MAX_phi_increment
												MAX_PHI_INCRE = abs(old_phi - phi_a->phi);
										}

#ifdef _DEBUG
										if (_isnan(phi_a->phi)) {
											cout << "DEBUG: phi_a->phi error !" << endl;
											SYS_PROGRAM_STOP;
										}
#endif
								}
							}
						}
			}
				if (r_istep % out_step == 0) {
					stringstream log2;
					log2 << "> iterated step:" << to_string(r_istep) << ",	max_phi_variation = " << MAX_PHI_INCRE << endl;
					Solvers::get_instance()->writer.add_string_to_txt_and_screen(log2.str(), LOG_FILE_NAME);
				}
		}
			for (auto phase = phaseMesh(0, 0, 0).begin(); phase < phaseMesh(0, 0, 0).end(); phase++)
				phaseMesh.delete_customFlag_in_allnodes(ExternalFields::RELAX_interface_buff + phase->index);
			stringstream log3;
			log3 << "> Finish interface relaxation !" << endl;
			Solvers::get_instance()->current_istep = 0;
			Solvers::get_instance()->real_time = 0.0;
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(log3.str(), LOG_FILE_NAME);
	}

		static void remove_inexistent_phis(FieldStorage_forPhaseNode& phaseMesh) {
			bool_box phis_need_erase;
			for (auto phase = phaseMesh(0, 0, 0).begin(); phase < phaseMesh(0, 0, 0).end(); phase++)
				phis_need_erase.add_bool(phase->index, true);
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						for (auto phase = node.begin(); phase < node.end(); phase++)
							if (phase->phi > Phi_Num_Cut_Off) {
#ifdef _OPENMP
#pragma omp critical
#endif
								{
									phis_need_erase.erase(phase->index);
								}
							}
					}
			string out = "> Do phis removing: ";
			for (auto need_erase = phis_need_erase.begin(); need_erase < phis_need_erase.end(); need_erase++) {
				out += to_string(need_erase->index) + ", ";
				phaseMesh.info_node.erase(need_erase->index);
			}
			out += "\n";
			Solvers::get_instance()->writer.add_string_to_txt_and_screen(out, LOG_FILE_NAME);
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						for (auto index = phis_need_erase.begin(); index < phis_need_erase.end(); index++)
							node.erase(index->index);
					}
		}

		static void filling_by_phis(FieldStorage_forPhaseNode& phaseMesh) {
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Do phis filling:\n", LOG_FILE_NAME);
			for (int index = 0; index < target_phi.size(); index++) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							for (auto phase = node.begin(); phase < node.end(); phase++) {
								for (auto phi_index = target_phi[index].begin(); phi_index < target_phi[index].end(); phi_index++) {
									if (phase->index == *phi_index) {
										if (phase->phi > SYS_EPSILON) {
											for (auto x = phase->x.begin(); x < phase->x.end(); x++)
												for (auto set_x = fill_phi_con[index].begin(); set_x < fill_phi_con[index].end(); set_x++)
													if (x->index == set_x->index)
														x->value = set_x->value;
											for (auto x = node.x.begin(); x < node.x.end(); x++)
												for (auto set_x = fill_total_con[index].begin(); set_x < fill_total_con[index].end(); set_x++)
													if (x->index == set_x->index)
														x->value = set_x->value;
											node.temperature.T = fill_temperature[index] * phase->phi;
										}
									}
								}
							}
						}
				stringstream report;
				report << "> already fill phi index = ";
				for (auto phi_index = target_phi[index].begin(); phi_index < target_phi[index].end(); phi_index++)
					report << *phi_index << ", ";
				report << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
			}
		}

		static void reordering_phis_index(FieldStorage_forPhaseNode& phaseMesh) {
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Do phis index reordering:\n", LOG_FILE_NAME);
			int reindex = reordering_index_0;
			for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++, reindex++) {
				stringstream report;
				report << "> phi_index: " << phase->index << " ---> " << reindex << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
				phase->index = reindex;
			}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						int rreindex = reordering_index_0;
						for (auto phase = node.begin(); phase < node.end(); phase++, rreindex++)
							phase->index = rreindex;
					}
		}

		static void optimization_memory_pair_wise_grand_potential(FieldStorage_forPhaseNode& phaseMesh) {
			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> Do optimization_memory :\n", LOG_FILE_NAME);
			pf::DifferenceMethod diff_method = Solvers::get_instance()->parameters.Difference_Method;
			vector<int> phase_indexes = Solvers::get_instance()->C_Solver.phase_indexes;
			double threshold = Solvers::get_instance()->C_Solver.threshold;
			ConEquationDomain _domainType = Solvers::get_instance()->parameters.ConEDomain;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
						for (int index = 0; index < fixed_phis.size(); index++)
							for (auto phase = node.begin(); phase < node.end(); phase++)
								if (phase->index == fixed_phis[index])
									phase->_flag = phaseMesh.currentFlag(node, phase->index);
					}
			for (int index = 0; index < fixed_phis.size(); index++) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh.limit_x; x++)
					for (int y = 0; y < phaseMesh.limit_y; y++)
						for (int z = 0; z < phaseMesh.limit_z; z++) {
							PhaseNode& node = phaseMesh(x, y, z);
							bool is_erase = false;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								if (phase->index == fixed_phis[index] && phase->phi > (1.0 - SYS_EPSILON) && phase->_flag == pf_BULK) {
									is_erase = true;
									for (int rel_x = -1; rel_x <= 1; rel_x++)
										for (int rel_y = -1; rel_y <= 1; rel_y++)
											for (int rel_z = -1; rel_z <= 1; rel_z++) {
												if (rel_x == 0 && rel_y == 0 && rel_z == 0)
													continue;
												PhaseNode& check_node = node.get_long_range_node(rel_x, rel_y, rel_z);
												for (auto check_phase = check_node.begin(); check_phase < check_node.end(); check_phase++)
													if (check_phase->index == phase->index && check_phase->_flag != pf_BULK)
														is_erase = false;
											}
								}
							if (is_erase) {
								node.erase_until(fixed_phis[index]);
							}
							if ((node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && _domainType == ConEquationDomain::CEDomain_Standard) ||
								(node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && _domainType == ConEquationDomain::CEDomain_Reverse)) {
								node.x.clear();
								node.potential.clear();
								node.kinetics_coeff.clear();
							}
						}
				stringstream report;
				report << "> already fixed phi index = " << index << endl;
				Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
			}
		}

		static void restore_memory_pair_wise_grand_potential(FieldStorage_forPhaseNode& phaseMesh) {
			vector<int> phase_indexes = Solvers::get_instance()->C_Solver.phase_indexes;
			double threshold = Solvers::get_instance()->C_Solver.threshold;
			ConEquationDomain _domainType = Solvers::get_instance()->parameters.ConEDomain;
			pf::Info_Node comps = Solvers::get_instance()->parameters.Components;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh.limit_x; x++)
				for (int y = 0; y < phaseMesh.limit_y; y++)
					for (int z = 0; z < phaseMesh.limit_z; z++) {
						PhaseNode& node = phaseMesh(x, y, z);
						node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
						for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
							bool is_phi_here = false;
							for (auto phase2 = node.begin(); phase2 < node.end(); phase2++)
								if (phase->index == phase2->index)
									is_phi_here = true;
							if (!is_phi_here)
								node.add_phase(phase->index, phase->property, pf_BULK, 0.0);
						}
						if ((node.customValues[ExternalFields::CON_Smooth_Phi] < threshold && _domainType == ConEquationDomain::CEDomain_Standard) ||
							(node.customValues[ExternalFields::CON_Smooth_Phi] > threshold && _domainType == ConEquationDomain::CEDomain_Reverse)) {
							if (node.potential.size() == 0)
								for (auto comp = comps.begin(); comp < comps.end(); comp++) {
									node.potential.add_con(comp->index, 0.0);
									node.x.add_con(comp->index, 0.0);
									for (auto comp2 = comps.begin(); comp2 < comps.end(); comp2++)
										node.kinetics_coeff.set(comp->index, comp2->index, 0.0);
								}
						}
					}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Preprocess.reconstruct_phis = {[(phi_index_0, phi_index_1, ... ), (phi_name)], .... } \n", InputFileReader::get_instance()->debug_file);
			string reconstruct_phis_key = "Preprocess.reconstruct_phis", reconstruct_phis_input = "{[()]}";
			if (InputFileReader::get_instance()->read_string_value(reconstruct_phis_key, reconstruct_phis_input, infile_debug)) {
				is_reconstruct = true;
				vector<InputValueType> reconstruct_phis_structure; reconstruct_phis_structure.push_back(InputValueType::IVType_INT); reconstruct_phis_structure.push_back(InputValueType::IVType_STRING);
				vector<vector<vector<input_value>>> reconstruct_phis_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(reconstruct_phis_structure, reconstruct_phis_key, reconstruct_phis_input, infile_debug);
				for (int index = 0; index < reconstruct_phis_value.size(); index++) {
					vector<int> phiss;
					for (auto phi_indexx = reconstruct_phis_value[index][0].begin(); phi_indexx < reconstruct_phis_value[index][0].end(); phi_indexx++)
						phiss.push_back(phi_indexx->int_value);
					recons_phi.push_back(phiss);
					recons_property.push_back(Solvers::get_instance()->parameters.Phases[reconstruct_phis_value[index][1][0].string_value].phi_property);
				}
			}

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Preprocess.relax_interface = (relax_steps, output_steps, is_interface_movable_in_relax , fix_phi_after_relax) \n", InputFileReader::get_instance()->debug_file);
			string relax_interface_key = "Preprocess.relax_interface", relax_interface_input = "()";
			if (InputFileReader::get_instance()->read_string_value(relax_interface_key, relax_interface_input, infile_debug)) {
				is_relax_interface_on = true;
				vector<InputValueType> relax_interface_structure; relax_interface_structure.push_back(InputValueType::IVType_INT);
				relax_interface_structure.push_back(InputValueType::IVType_INT); relax_interface_structure.push_back(InputValueType::IVType_BOOL); relax_interface_structure.push_back(InputValueType::IVType_BOOL);
				vector<input_value> relax_interface_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(relax_interface_structure, relax_interface_key, relax_interface_input, infile_debug);
				relaxation_steps = relax_interface_value[0].int_value;
				out_step = relax_interface_value[1].int_value;
				is_interface_movable = relax_interface_value[2].bool_value;
				is_fix_phi_in_loop = relax_interface_value[3].bool_value;
			}

			// should be concerned that matrix phases problem
			InputFileReader::get_instance()->read_bool_value("Preprocess.remove_inexistent_phis", is_remove_inexistent_phis, infile_debug);

			// should be concerned that matrix phases problem
			if (InputFileReader::get_instance()->read_int_value("Preprocess.re_ordering_phis_indexs_from", reordering_index_0, infile_debug))
				is_phis_reordering = true;

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Preprocess.fill_phis = {[(phi_index_0, phi_index_1, ... ), (phi_con_1, phi_con_2, ... ), (total_con_1, total_con_2, ... ), (temperature)], .... } \n", InputFileReader::get_instance()->debug_file);
			string fill_phis_key = "Preprocess.fill_phis", fill_phis_input = "{[()]}";
			if (InputFileReader::get_instance()->read_string_value(fill_phis_key, fill_phis_input, infile_debug)) {
				is_filling_by_phi = true;
				vector<InputValueType> fill_phis_structure; fill_phis_structure.push_back(InputValueType::IVType_INT);
				fill_phis_structure.push_back(InputValueType::IVType_DOUBLE); fill_phis_structure.push_back(InputValueType::IVType_DOUBLE);
				fill_phis_structure.push_back(InputValueType::IVType_DOUBLE);
				vector<vector<vector<input_value>>> fill_phis_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(fill_phis_structure, fill_phis_key, fill_phis_input, infile_debug);
				for (int index = 0; index < fill_phis_value.size(); index++) {
					vector<int> phiss; int this_property = phaseMesh.info_node[fill_phis_value[index][0][0].int_value].property;
					for (auto phi_indexx = fill_phis_value[index][0].begin(); phi_indexx < fill_phis_value[index][0].end(); phi_indexx++) {
						phiss.push_back(phi_indexx->int_value);
						if (this_property != phaseMesh.info_node[phi_indexx->int_value].property) {
							InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Property of phis want to be filled should be same ! \n", InputFileReader::get_instance()->debug_file);
							exit(0);
						}
					}
					target_phi.push_back(phiss);
					double_box phi_con, total_con, total_potential;
					if (Solvers::get_instance()->parameters.Phases[this_property].x.size() != fill_phis_value[index][1].size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("> ERROR Preprocess.fill_phis : define phase_x and their values error, size of phase_x mismatch ! \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					for (int x_index = 0; x_index < fill_phis_value[index][1].size(); x_index++)
						phi_con.add_double((Solvers::get_instance()->parameters.Phases[this_property].x.begin() + x_index)->index, fill_phis_value[index][1][x_index].double_value);
					if (Solvers::get_instance()->parameters.Components.size() != fill_phis_value[index][2].size()) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("> ERROR Preprocess.fill_phis : define total_x and their values error, size of total_x mismatch ! \n", InputFileReader::get_instance()->debug_file);
						exit(0);
					}
					for (int x_index = 0; x_index < fill_phis_value[index][2].size(); x_index++)
						total_con.add_double((Solvers::get_instance()->parameters.Components.begin() + x_index)->index, fill_phis_value[index][2][x_index].double_value);

					fill_phi_con.push_back(phi_con);
					fill_total_con.push_back(total_con);
					fill_temperature.push_back(fill_phis_value[index][3][0].double_value);
				}
			}

			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise && Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				int_box fix_phi_index, no_grand_phi_index;
				for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
					fix_phi_index.add_int(phase->index, 1);
					no_grand_phi_index.add_int(phase->index, 1);
				}
				tensor2_double matrix_Lij = interface_mobility::get_matirx_Lij();
				PairValue block_Lij = interface_mobility::get_block_Lij();
				for (auto Li = matrix_Lij.begin(); Li < matrix_Lij.end(); Li++)
					for (auto Lj = Li->begin(); Lj < Li->end(); Lj++) {
						for (auto fix_phi = fix_phi_index.begin(); fix_phi < fix_phi_index.end(); fix_phi++)
							if (fix_phi->index == Li->index || fix_phi->index == Lj->index)
								fix_phi->value = 0;
					}
				for (auto block = block_Lij.begin(); block < block_Lij.end(); block++)
					for (int index = block->pairIndex_1; index <= block->pairIndex_2; index++)
						for (auto fix_phi = fix_phi_index.begin(); fix_phi < fix_phi_index.end(); fix_phi++)
							if (fix_phi->index == index)
								fix_phi->value = 0;
				if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Standard) {
					for (auto grand_phi = no_grand_phi_index.begin(); grand_phi < no_grand_phi_index.end(); grand_phi++) {
						grand_phi->value = 1;
						for (auto index = Solvers::get_instance()->C_Solver.phase_indexes.begin(); index < Solvers::get_instance()->C_Solver.phase_indexes.end(); index++)
							if (grand_phi->index == *index)
								grand_phi->value = 0;
					}
				}
				else if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse) {
					for (auto grand_phi = no_grand_phi_index.begin(); grand_phi < no_grand_phi_index.end(); grand_phi++) {
						grand_phi->value = 0;
						for (auto index = Solvers::get_instance()->C_Solver.phase_indexes.begin(); index < Solvers::get_instance()->C_Solver.phase_indexes.end(); index++)
							if (grand_phi->index == *index)
								grand_phi->value = 1;
					}
				}
				if (infile_debug) {
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Memory optimization, which used to reduce running memory costs. ( carefully !!!! ) \n", InputFileReader::get_instance()->debug_file);
					stringstream out;
					fixed_phis.clear();
					out << "# No mobility phises : ";
					for (auto fix_phi = fix_phi_index.begin(); fix_phi < fix_phi_index.end(); fix_phi++)
						if (fix_phi->value == 1)
							out << fix_phi->index << ", ";
					out << endl << "# No grand-potential phises : ";
					for (auto grand_phi = no_grand_phi_index.begin(); grand_phi < no_grand_phi_index.end(); grand_phi++) {
						if (grand_phi->value == 1)
							out << grand_phi->index << ", ";
						else if (grand_phi->value == 0 && fix_phi_index[grand_phi->index] == 1)
							fix_phi_index[grand_phi->index] = 0;
					}
					out << endl << "# Phies can be fixed : ";
					for (auto fix_phi = fix_phi_index.begin(); fix_phi < fix_phi_index.end(); fix_phi++)
						if (fix_phi->value == 1) {
							out << fix_phi->index << ", ";
							fixed_phis.push_back(fix_phi->index);
						}
					out << endl << "# Preprocess.memory_optimized = 0 - MO_NONE, 1 - MO_OPT, 2 - MO_RESTORE" << endl;
					InputFileReader::get_instance()->debug_writer->add_string_to_txt(out.str(), InputFileReader::get_instance()->debug_file);
					int opt = 0;
					InputFileReader::get_instance()->read_int_value("Preprocess.is_memory_optimized", opt, infile_debug);
					memory_optimization = MemoryOptimization(opt);
				}

			}

			Solvers::get_instance()->writer.add_string_to_txt_and_screen("> MODULE INIT : Pretreatment !\n", LOG_FILE_NAME);
		}

		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_reconstruct) {
				reconstruct_phis(phaseMesh);
			}
			if (is_relax_interface_on) {
				relaxation_interface(phaseMesh);
				if (is_fix_phi_in_loop) {
					Solvers::get_instance()->parameters.PhiEType = PhiEquationType::PEType_Const;
				}
			}
			if (is_remove_inexistent_phis) {
				remove_inexistent_phis(phaseMesh);
			}
			if (is_phis_reordering) {
				reordering_phis_index(phaseMesh);
			}
			if (is_filling_by_phi) {
				filling_by_phis(phaseMesh);
			}
			if (memory_optimization == MemoryOptimization::MO_OPT)
				optimization_memory_pair_wise_grand_potential(phaseMesh);
			else if (memory_optimization == MemoryOptimization::MO_RESTORE)
				restore_memory_pair_wise_grand_potential(phaseMesh);
		}

		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			return report;
		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}

		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}

		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {

		}
}
}