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

namespace pf {
	namespace concentration_init {
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
				vector<int> phase_indexes = Solvers::get_instance()->C_Solver.phase_indexes;
				if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Standard) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
								node.customValues.add_double(ExternalFields::CON_Smooth_Old_Phi, node.customValues[ExternalFields::CON_Smooth_Phi]);
							}
				}
				else if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
								node.customValues.add_double(ExternalFields::CON_Smooth_Old_Phi, node.customValues[ExternalFields::CON_Smooth_Phi]);
							}
				}
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				vector<int> phase_indexes = Solvers::get_instance()->C_Solver.phase_indexes;
				double threshold = Solvers::get_instance()->C_Solver.threshold;
				auto dfbulk_dx = Solvers::get_instance()->C_Solver.dfbulk_dx;
				if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Standard) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
								node.customValues.add_double(ExternalFields::CON_Smooth_Old_Phi, node.customValues[ExternalFields::CON_Smooth_Phi]);
								if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
									for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
										p->gradient.set_to_zero();
										p->increment = 0.0;
										p->laplacian = 0.0;
										p->value = 0.0;
										for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++)
											node.kinetics_coeff.set(p->index, p2->index, 0.0);
									}
									double sum_phi = 0.0;
									for (auto phase = node.begin(); phase < node.end(); phase++)
										for (auto index = phase_indexes.begin(); index < phase_indexes.end(); index++)
											if (phase->index == *index && phase->phi > SYS_EPSILON) {
												sum_phi += phase->phi;
												dfbulk_dx(node, *phase);
												for (auto p = phase->potential.begin(); p < phase->potential.end(); p++)
													node.potential[p->index].value += p->value * phase->phi;
											}
									if (sum_phi > SYS_EPSILON) {
										for (auto p = node.potential.begin(); p < node.potential.end(); p++)
											p->value /= sum_phi;
									}
									else {
										for (auto p = node.potential.begin(); p < node.potential.end(); p++)
											p->value = 0.0;
									}
								}
								else {
									for (auto p = node.potential.begin(); p < node.potential.end(); p++)
										p->value = 0.0;
								}
							}
				}
				else if (Solvers::get_instance()->parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh.limit_x; x++)
						for (int y = 0; y < phaseMesh.limit_y; y++)
							for (int z = 0; z < phaseMesh.limit_z; z++) {
								PhaseNode& node = phaseMesh(x, y, z);
								node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
								node.customValues.add_double(ExternalFields::CON_Smooth_Old_Phi, node.customValues[ExternalFields::CON_Smooth_Phi]);
								if (node.customValues[ExternalFields::CON_Smooth_Phi] < threshold) {
									for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
										p->gradient.set_to_zero();
										p->increment = 0.0;
										p->laplacian = 0.0;
										p->value = 0.0;
										for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++)
											node.kinetics_coeff.set(p->index, p2->index, 0.0);
									}
									double sum_phi = 0.0;
									for (auto phase = node.begin(); phase < node.end(); phase++) {
										double is_cal = true;
										for (auto index = phase_indexes.begin(); index < phase_indexes.end(); index++)
											if (phase->index == *index)
												is_cal = false;
										if (is_cal && phase->phi > SYS_EPSILON) {
											sum_phi += phase->phi;
											dfbulk_dx(node, *phase);
											for (auto p = phase->potential.begin(); p < phase->potential.end(); p++)
												node.potential[p->index].value += p->value * phase->phi;
										}
									}
									if (sum_phi > SYS_EPSILON) {
										for (auto p = node.potential.begin(); p < node.potential.end(); p++)
											p->value /= sum_phi;
									}
									else {
										for (auto p = node.potential.begin(); p < node.potential.end(); p++)
											p->value = 0.0;
									}
								}
								else {
									for (auto p = node.potential.begin(); p < node.potential.end(); p++)
										p->value = 0.0;
								}
							}
				}
			}
		}
	}
}