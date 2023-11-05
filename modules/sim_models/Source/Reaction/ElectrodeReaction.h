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
#include "../../../Base.h"
#include "../../InterfaceEnergy.h"

namespace pf {
	namespace electrode_reaction {
        static bool is_charge_on = false;
        //-parameters
        static bool is_Vapp_const = false;
        static double V_theo = 0.0;   // [-]
        static int electrolyte_phase_index = 0;
        static int active_component = 0;
        static double k_BV = 0.0;   // k_0 / F
        static double kb_over_ela = 8.617333e-5;
        static double Vapp_low_cutoff = 0.0;
        static double Vapp_high_cutoff = 0.0;
        static double V_ref = 0.0;
        static double T = 300;           // [Kelvin]
        static double Crate = 1.0;       // [1/h]
        static double X_full = 1.0; // [-]
        static double V_app = 1.0; // [-]

        static ConEquationDomain _domain_type = ConEquationDomain::CEDomain_Standard;

        static int get_electrolyte_phase_index() {
            return electrolyte_phase_index;
        }

        static int get_active_component() {
            return active_component;
        }

        static void calc_applied_voltage_smooth(FieldStorage_forPhaseNode& phaseMesh) {
            if (!is_charge_on)
                return;
            double inty[2] = { 0.0, 0.0 }, flux_theo = 0.0;
            double electrode_phi_threshold = Phi_Num_Cut_Off;
            ConEquationType _type = Solvers::get_instance()->parameters.ConEType;
            auto abs_grad_phi_AB = Solvers::get_instance()->C_Solver.abs_grad_phi_AB;
            if (_type == ConEquationType::CEType_GrandP || _type == ConEquationType::CEType_TotalX)
                electrode_phi_threshold = Solvers::get_instance()->C_Solver.threshold;
            if (_type == ConEquationType::CEType_PhaseX) {
#pragma omp parallel for
                for (int x = 0; x < phaseMesh.limit_x; x++)
                    for (int y = 0; y < phaseMesh.limit_y; y++)
                        for (int z = 0; z < phaseMesh.limit_z; z++) {
                            PhaseNode& node = phaseMesh(x, y, z);
                            PhaseEntry& electrolyte = node[electrolyte_phase_index];
                            double sum_phi = 0.0;
                            for (auto phase = node.begin(); phase < node.end(); phase++)
                                if (phase->index != electrolyte.index) {
                                    if (phase->phi > electrode_phi_threshold && electrolyte.phi > electrode_phi_threshold) {
                                        sum_phi += phase->phi;
                                        double phi_ab_norm = abs_grad_phi_AB(node, *phase, electrolyte), active_c = phase->x[active_component].value,
                                            active_u = phase->potential[active_component].value;
                                        double int_0 = phi_ab_norm * sqrt((1.0 - active_c) * active_c) * exp(-active_u / 2.0),
                                            int_1 = phi_ab_norm * sqrt((1.0 - active_c) * active_c) * exp(active_u / 2.0);

#ifdef _OPENMP
#pragma omp critical
#endif
                                        {
                                            inty[0] += int_0;
                                            inty[1] += int_1;
                                        }
                                    }
                                }
                            double int_flux = sum_phi;
#ifdef _OPENMP
#pragma omp critical
#endif
                            {
                                flux_theo += int_flux;
                            }
                        }
            }
            else if (_type == ConEquationType::CEType_GrandP || _type == ConEquationType::CEType_TotalX) {
                if (_domain_type == ConEquationDomain::CEDomain_Standard) {
#pragma omp parallel for
                    for (int x = 0; x < phaseMesh.limit_x; x++)
                        for (int y = 0; y < phaseMesh.limit_y; y++)
                            for (int z = 0; z < phaseMesh.limit_z; z++) {
                                PhaseNode& node = phaseMesh(x, y, z);
                                if (node.customValues[ExternalFields::CON_Smooth_Phi] > electrode_phi_threshold) {
                                    Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh.dr);
                                    double phi_grad_norm = s_phi_grad.abs();
                                    double active_phi = node.customValues[ExternalFields::CON_Smooth_Phi];
                                    if (phi_grad_norm > SYS_EPSILON) {
                                        double active_c = node.x[active_component].value, active_u = node.potential[active_component].value;
                                        double  int_0 = phi_grad_norm * sqrt((1.0 - active_c) * active_c) * exp(-active_u / 2.0),
                                                int_1 = phi_grad_norm * sqrt((1.0 - active_c) * active_c) * exp(active_u / 2.0);
#ifdef _OPENMP
#pragma omp critical
#endif
                                        {
                                            inty[0] += int_0;
                                            inty[1] += int_1;
                                        }
#ifdef _DEBUG
                                        if (_isnan(int_0) || _isnan(int_1)) {
                                            cout << "DEBUG: int error !" << endl;
                                            SYS_PROGRAM_STOP;
                                        }
#endif
                                    }
#ifdef _OPENMP
#pragma omp critical
#endif
                                    {
                                        flux_theo += active_phi;
                                    }
                                }
                            }
                }
                else {
#pragma omp parallel for
                    for (int x = 0; x < phaseMesh.limit_x; x++)
                        for (int y = 0; y < phaseMesh.limit_y; y++)
                            for (int z = 0; z < phaseMesh.limit_z; z++) {
                                PhaseNode& node = phaseMesh(x, y, z);
                                if (node.customValues[ExternalFields::CON_Smooth_Phi] < electrode_phi_threshold) {
                                    Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh.dr);
                                    double phi_grad_norm = s_phi_grad.abs();
                                    double active_phi = 1.0 - node.customValues[ExternalFields::CON_Smooth_Phi];
                                    if (phi_grad_norm > SYS_EPSILON) {
                                        double active_c = node.x[active_component].value, active_u = node.potential[active_component].value;
                                        if (active_c < 0.0)
                                            active_c = 0.0;
                                        else if (active_c > 1.0)
                                            active_c = 1.0;
                                        double  int_0 = phi_grad_norm * sqrt((1.0 - active_c) * active_c) * exp(-active_u / 2.0),
                                                int_1 = phi_grad_norm * sqrt((1.0 - active_c) * active_c) * exp(active_u / 2.0);
#ifdef _OPENMP
#pragma omp critical
#endif
                                        {
                                            inty[0] += int_0;
                                            inty[1] += int_1;
                                        }
#ifdef _DEBUG
                                        if (_isnan(int_0) || _isnan(int_1)) {
                                            cout << "DEBUG: int error !" << endl;
                                            SYS_PROGRAM_STOP;
                                        }
#endif
                                    }
#ifdef _OPENMP
#pragma omp critical
#endif
                                    {
                                        flux_theo += active_phi;
                                    }
                                }
                            }
                }
            }
            inty[0] = k_BV * inty[0] * phaseMesh.dr * phaseMesh.dr * phaseMesh.dr;
            inty[1] = k_BV * inty[1] * phaseMesh.dr * phaseMesh.dr * phaseMesh.dr;
            flux_theo = X_full * Crate * flux_theo / 3600.0 * phaseMesh.dr * phaseMesh.dr * phaseMesh.dr;
            double expV = (flux_theo + sqrt(flux_theo * flux_theo + 4.0 * inty[0] * inty[1])) / (2.0 * inty[0]);
            V_app = V_ref - 2.0 * kb_over_ela * T * log(expV);
            if (Is_Equality(inty[0], 0.0))
                V_app = V_ref;
            V_app = min(max(Vapp_low_cutoff, V_app), Vapp_high_cutoff);
        }

        static void calc_battery_charge_smoothed_boundary(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double abs_dPaPb) {
            if (alpha.index != electrolyte_phase_index && beta.index == electrolyte_phase_index) {
                for (auto comp = alpha.x.begin(); comp < alpha.x.end(); comp++) {
                    if (comp->index == active_component) {
                        double overpotential = V_ref - V_app - kb_over_ela * T * alpha.potential[active_component].value,
                            exp_flux = exp(overpotential / (2.0 * kb_over_ela * T)),
                            electrode_reac = abs_dPaPb * k_BV * sqrt((1.0 - comp->value) * comp->value) * (exp_flux - 1 / exp_flux);
                        if (comp->value < 0.0 || comp->value > 1.0)
                            electrode_reac = 0.0;
#ifdef _DEBUG
                        if (electrode_reac > 100.0 || electrode_reac < -100.0) {
                            cout << "DEBUG: battery_charge error !" << endl;
                            SYS_PROGRAM_STOP;
                        }
#endif
                        if (electrode_reac > 0 && comp->value < (1.0 - SYS_EPSILON)) {
                            comp->ChemicalReactionFlux += electrode_reac;
                        }
                        else if (electrode_reac < 0 && comp->value > SYS_EPSILON) {
                            comp->ChemicalReactionFlux += electrode_reac;
                        }
                    }
                }
            }
        }

        double int_flux(pf::PhaseNode& node, int con_i) {
            if (con_i == active_component) {
                double overpotential = V_ref - V_app - kb_over_ela * T * node.potential[active_component].value,
                    exp_flux = exp(overpotential / (2.0 * kb_over_ela * T)), comp_val = node.x[active_component].value;
                if (comp_val < 0.0)
                    comp_val = 0.0;
                else if (comp_val > 1.0)
                    comp_val = 1.0;
#ifdef _DEBUG
                if (exp_flux > -SYS_EPSILON && exp_flux < SYS_EPSILON) {  // problems here !!!!
                    cout << "DEBUG: battery_charge error !" << endl;
                    //SYS_PROGRAM_STOP;
                }
#endif
                return k_BV * sqrt((1.0 - comp_val) * comp_val) * (exp_flux - 1 / exp_flux);
            }
            else {
                return 0.0;
            }
        }

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
            bool infile_debug = false;
            is_charge_on = true;

            InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

            InputFileReader::get_instance()->read_int_value("ModelsManager.Con.ElectrodeReaction.electrolyte_index", electrolyte_phase_index, false);

            string comp_name = "first con";
            if (InputFileReader::get_instance()->read_string_value("ModelsManager.Con.ElectrodeReaction.active_con", comp_name, infile_debug))
                active_component = Solvers::get_instance()->parameters.Components[comp_name].index;

            InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.k_BV", k_BV, infile_debug);

            InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.reference_voltage", V_ref, infile_debug);

            InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.temperature", T, infile_debug);

            if (InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.applied_voltage", V_theo, infile_debug)) {
                is_Vapp_const = true;
            }
            else {
                is_Vapp_const = false;

                InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.AppliedVoltage.cut_off_low", Vapp_low_cutoff, infile_debug);
                InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.AppliedVoltage.cut_off_high", Vapp_high_cutoff, infile_debug);

                InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.Crate", Crate, infile_debug);

                InputFileReader::get_instance()->read_double_value("ModelsManager.Con.ElectrodeReaction.X_full", X_full, infile_debug);
            }

            _domain_type = Solvers::get_instance()->parameters.ConEDomain;
		}

		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
            _domain_type = Solvers::get_instance()->parameters.ConEDomain;
            if (!is_charge_on)
                return;
            if (!is_Vapp_const) {
                if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP ||
                    Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_TotalX) {
                    vector<int> phase_indexes = Solvers::get_instance()->C_Solver.phase_indexes;
#pragma omp parallel for
                    for (int x = 0; x < phaseMesh.limit_x; x++)
                        for (int y = 0; y < phaseMesh.limit_y; y++)
                            for (int z = 0; z < phaseMesh.limit_z; z++) {
                                PhaseNode& node = phaseMesh(x, y, z);
                                node.customValues.add_double(ExternalFields::CON_Smooth_Phi, node.cal_phases_fraction_by_index(phase_indexes));
                            }
                }
                calc_applied_voltage_smooth(phaseMesh);
            }
            stringstream report;
            report << "> Electrode intercalation : V_ref = " << V_ref << ", V_app = " << V_app << ", reaction_rate = " << k_BV << endl;
            Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			stringstream report;
            if (!is_charge_on)
                return report.str();
            if (!is_Vapp_const)
                calc_applied_voltage_smooth(phaseMesh);
            report << "> Electrode intercalation : V_ref = " << V_ref << ", V_app = " << V_app << ", reaction_rate = " << k_BV << endl;
			return report.str();
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
            if (!is_charge_on)
                return;
            double electrode_phi_threshold = Phi_Num_Cut_Off;
            ConEquationType _type = Solvers::get_instance()->parameters.ConEType;
            auto abs_grad_phi_AB = Solvers::get_instance()->C_Solver.abs_grad_phi_AB;
            if (_type == ConEquationType::CEType_GrandP || _type == ConEquationType::CEType_TotalX)
                electrode_phi_threshold = Solvers::get_instance()->C_Solver.threshold;
            fout << "<DataArray type = \"Float64\" Name = \"" << "elec_flux" <<
                "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
            for (int k = 0; k < phaseMesh.limit_z; k++)
                for (int j = 0; j < phaseMesh.limit_y; j++)
                    for (int i = 0; i < phaseMesh.limit_x; i++) {
                        PhaseNode& node = phaseMesh(i, j, k);
                        double val = 0.0;
                        if (_type == ConEquationType::CEType_PhaseX) {
                            double merge = 0.0;
                            PhaseEntry& electrolyte = node[electrolyte_phase_index];
                            for (auto phase = node.begin(); phase < node.end(); phase++)
                                if (phase->index != electrolyte.index 
                                    && phase->phi > electrode_phi_threshold && electrolyte.phi > electrode_phi_threshold) {
                                    merge += phase->phi;
                                    for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
                                        if (comp->index == active_component) {
                                            double overpotential = V_ref - V_app - kb_over_ela * T * phase->potential[active_component].value,
                                                exp_flux = exp(overpotential / (2.0 * kb_over_ela * T)),
                                                electrode_reac = abs_grad_phi_AB(node, *phase, electrolyte) * k_BV * sqrt((1.0 - comp->value) * comp->value) * (exp_flux - 1 / exp_flux);
                                            if (comp->value < 0.0 || comp->value > 1.0)
                                                electrode_reac = 0.0;
                                            if (electrode_reac > 0 && comp->value < (1.0 - SYS_EPSILON)) {
                                                val += electrode_reac * phase->phi;
                                            }
                                            else if (electrode_reac < 0 && comp->value > SYS_EPSILON) {
                                                val += electrode_reac * phase->phi;
                                            }
                                        }
                                    }
                                }
                            if (merge > electrode_phi_threshold)
                                val /= merge;
                            else
                                val = 0.0;
                        }
                        else if (_type == ConEquationType::CEType_GrandP || _type == ConEquationType::CEType_TotalX) {
                            Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh.dr);
                            double phi_grad_norm = s_phi_grad.abs();
                            if (_domain_type == ConEquationDomain::CEDomain_Standard && node.customValues[ExternalFields::CON_Smooth_Phi] > electrode_phi_threshold && phi_grad_norm > SYS_EPSILON) {
                                double overpotential = V_ref - V_app - kb_over_ela * T * node.potential[active_component].value,
                                    exp_flux = exp(overpotential / (2.0 * kb_over_ela * T)), comp_val = node.x[active_component].value;
                                if (comp_val < 0.0)
                                    comp_val = 0.0;
                                else if (comp_val > 1.0)
                                    comp_val = 1.0;
                                val = phi_grad_norm * k_BV * sqrt((1.0 - comp_val) * comp_val) * (exp_flux - 1 / exp_flux);
                            }
                            else if (_domain_type == ConEquationDomain::CEDomain_Reverse && node.customValues[ExternalFields::CON_Smooth_Phi] < electrode_phi_threshold && phi_grad_norm > SYS_EPSILON) {
                                double overpotential = V_ref - V_app - kb_over_ela * T * node.potential[active_component].value,
                                    exp_flux = exp(overpotential / (2.0 * kb_over_ela * T)), comp_val = node.x[active_component].value;
                                if (comp_val < 0.0)
                                    comp_val = 0.0;
                                else if (comp_val > 1.0)
                                    comp_val = 1.0;
                                val = phi_grad_norm * k_BV * sqrt((1.0 - comp_val) * comp_val) * (exp_flux - 1 / exp_flux);
                            }
                        }
#ifdef _DEBUG
                        if (IS_NAN(val)) {  // problems here !!!!
                            cout << "DEBUG: battery_charge error !" << endl;
                            SYS_PROGRAM_STOP;
                        }
#endif
                        fout << val << endl;
                    }
            fout << "</DataArray>" << endl;
		}
	}
}