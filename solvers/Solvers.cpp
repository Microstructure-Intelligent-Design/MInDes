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


#include"Solvers.h"
namespace pf {
	Solvers* Solvers::solver = nullptr;
	Solvers::Solvers() {
		Working_Folder_Path = "";
		Infile_Folder_Path = "";
		current_istep = 0;
		real_time = 0.0;
		statistics_step = -1000000;
		is_stop_loop = false;
		Phi_Solver_AC.init(phaseMesh);
		Phi_Solver_CH.init(phaseMesh);
		C_Solver.init(phaseMesh);
		T_Solver.init(phaseMesh);
	};
	Solvers::~Solvers() {
		for (auto _module = modules.begin(); _module < modules.end(); _module++)
			_module->deinit(phaseMesh);
		modules.clear();
		phaseMesh.free();
	};
	void Solvers::modules_pre_exec() {
		stringstream pre_exec;
		pre_exec << timer::return_cunrrent_time_by_string();
		pre_exec << ">-------------------------------------------" << "Modules pre-Exec" << "-------------------------------------------" << endl;
		writer.add_string_to_txt_and_screen(pre_exec.str(), LOG_FILE_NAME);
		for (auto _module = modules.begin(); _module < modules.end(); _module++)
			_module->exec_pre(phaseMesh);

		if(parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise)
			Phi_Solver_AC.init_phi_pair_wise(parameters.is_Normalize_Phi);
		else
			Phi_Solver_AC.init_phi(parameters.is_Normalize_Phi);

		pre_exec.str("");
		pre_exec << ">-------------------------------------------" << "--------------" << "-------------------------------------------" << endl << endl;
		writer.add_string_to_txt_and_screen(pre_exec.str(), LOG_FILE_NAME);
	}
	void Solvers::solve_phi_x_T_in_loop() {
		if (parameters.OpenMP_Thread_Counts >= omp_get_num_procs())
			parameters.OpenMP_Thread_Counts = omp_get_num_procs() - 1;
		else if (parameters.OpenMP_Thread_Counts < 1)
			parameters.OpenMP_Thread_Counts = 1;
		omp_set_num_threads(parameters.OpenMP_Thread_Counts);
		int MAX_counts = 1;
		if (parameters.SKIP[CON_TEMP_SKIP::CTS_CON] > parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP])
			MAX_counts = parameters.SKIP[CON_TEMP_SKIP::CTS_CON];
		else
			MAX_counts = parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP];
		// timer
		double output_times = 0.0;
		MAX_PHI_INCREMENT = 0.0;
		MAX_CON_INCREMENT = 0.0;
		MAX_TEMP_INCREMENT = 0.0;
		MAX_DIFF_FLUX_INCREMENT = 0.0;
		MAX_REAC_FLUX_INCREMENT = 0.0;
		MAX_PHI_TRANS_FLUX_INCREMENT = 0.0;
		for (int count = 1; count <= MAX_counts; count++) {
			if (count == 1) {
				timer::time_interval_precision_secs_begin(output_times);
				if (parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
					Phi_Solver_AC.pre_calculation_phi_normal(parameters.dt, parameters.is_Normalize_Phi, parameters.Difference_Method);
					MAX_PHI_INCREMENT = Phi_Solver_AC.solve_phi_normal(parameters.dt, parameters.is_Normalize_Phi);
				}
				else if (parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise) {
					Phi_Solver_AC.pre_calculation_phi_pair_wise_normal(parameters.dt, parameters.is_Normalize_Phi, parameters.Difference_Method);
					MAX_PHI_INCREMENT = Phi_Solver_AC.solve_phi_pair_wise_normal(parameters.dt, parameters.is_Normalize_Phi);
				}
				else if (parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
					Phi_Solver_CH.pre_calculation_phis();
					MAX_PHI_INCREMENT = Phi_Solver_CH.solve_phis(parameters.dt, parameters.is_Normalize_Phi);
				}
				timer::time_interval_precision_secs_end(output_times);
				Solvers_Timer.t_interval_evolve_phi += output_times;
			}
			if (count <= parameters.SKIP[CON_TEMP_SKIP::CTS_CON]) {
				timer::time_interval_precision_secs_begin(output_times);
				if (parameters.ConEType == ConEquationType::CEType_PhaseX) {
					C_Solver.pre_calculation_phase_concentration(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
					vector<double> INCREMENT = C_Solver.solve_phase_concentration(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON], parameters.is_Normalize_Con);
					if (INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX] > MAX_DIFF_FLUX_INCREMENT)
						MAX_DIFF_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX];
					if (INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX] > MAX_REAC_FLUX_INCREMENT)
						MAX_REAC_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX];
					if (INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] > MAX_PHI_TRANS_FLUX_INCREMENT)
						MAX_PHI_TRANS_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX];
					if (INCREMENT[CH_RETURN::CH_MAX_X_INCREMENT_FLUX] > MAX_CON_INCREMENT)
						MAX_CON_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_X_INCREMENT_FLUX];
				}
				else if (parameters.ConEType == ConEquationType::CEType_TotalX) {
					if (parameters.ConEDomain == ConEquationDomain::CEDomain_Standard) {
						C_Solver.pre_calculation_total_concentration_inside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
						vector<double> INCREMENT = C_Solver.solve_total_concentration_inside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON], parameters.is_Normalize_Con);
						if (INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX] > MAX_DIFF_FLUX_INCREMENT)
							MAX_DIFF_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX];
						if (INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX] > MAX_REAC_FLUX_INCREMENT)
							MAX_REAC_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX];
						if (INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] > MAX_PHI_TRANS_FLUX_INCREMENT)
							MAX_PHI_TRANS_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX];
						if (INCREMENT[CH_RETURN::CH_MAX_X_INCREMENT_FLUX] > MAX_CON_INCREMENT)
							MAX_CON_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_X_INCREMENT_FLUX];
					}
					else if (parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse) {
						C_Solver.pre_calculation_total_concentration_outside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
						vector<double> INCREMENT = C_Solver.solve_total_concentration_outside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON], parameters.is_Normalize_Con);
						if (INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX] > MAX_DIFF_FLUX_INCREMENT)
							MAX_DIFF_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX];
						if (INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX] > MAX_REAC_FLUX_INCREMENT)
							MAX_REAC_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX];
						if (INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] > MAX_PHI_TRANS_FLUX_INCREMENT)
							MAX_PHI_TRANS_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX];
						if (INCREMENT[CH_RETURN::CH_MAX_X_INCREMENT_FLUX] > MAX_CON_INCREMENT)
							MAX_CON_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_X_INCREMENT_FLUX];
					}
				}
				else if (parameters.ConEType == ConEquationType::CEType_GrandP) {
					vector<double> INCREMENT;
					if (parameters.ConEDomain == ConEquationDomain::CEDomain_Standard) {
						INCREMENT = C_Solver.pre_calculation_grand_potential_functional_inside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
						double CON_INCREMENT = C_Solver.solve_grand_potential_functional_inside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
						if (CON_INCREMENT > MAX_CON_INCREMENT)
							MAX_CON_INCREMENT = CON_INCREMENT;
					}
					else if (parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse) {
						INCREMENT = C_Solver.pre_calculation_grand_potential_functional_outside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
						double CON_INCREMENT = C_Solver.solve_grand_potential_functional_outside_phis(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_CON]);
						if (CON_INCREMENT > MAX_CON_INCREMENT)
							MAX_CON_INCREMENT = CON_INCREMENT;
					}
					if (INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX] > MAX_DIFF_FLUX_INCREMENT)
						MAX_DIFF_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_DIFFUSION_FLUX];
					if (INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX] > MAX_REAC_FLUX_INCREMENT)
						MAX_REAC_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_REACTION_FLUX];
					if (INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] > MAX_PHI_TRANS_FLUX_INCREMENT)
						MAX_PHI_TRANS_FLUX_INCREMENT = INCREMENT[CH_RETURN::CH_MAX_PHI_TRANS_FLUX];
				}
				timer::time_interval_precision_secs_end(output_times);
				Solvers_Timer.t_interval_evolve_con += output_times;
			}
			if (count <= parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP]) {
				timer::time_interval_precision_secs_begin(output_times);
				if (parameters.TempEType == TemperatureEquationType::TType_Standard) {
					T_Solver.pre_calculation_temperature(parameters.Difference_Method);
					double INCREMENT = T_Solver.solve_temperature(parameters.dt / parameters.SKIP[CON_TEMP_SKIP::CTS_TEMP]);
					if (INCREMENT > MAX_TEMP_INCREMENT)
						MAX_TEMP_INCREMENT = INCREMENT;
				}
				timer::time_interval_precision_secs_end(output_times);
				Solvers_Timer.t_interval_evolve_T += output_times;
			}
		}
	}
	void Solvers::output_with_loop() {
		if (current_istep % parameters.screen_loop_step == 0) {
			stringstream log;
			log << "#====================================================================================================" << endl;
			log << "# CURRENT STEP = " << current_istep << ", REAL TIME = " << real_time << endl;
			log << "# MAX PHI : increment = " << setprecision(5) << MAX_PHI_INCREMENT << endl;
			log << "# MAX CON : increment = " << setprecision(5) << MAX_CON_INCREMENT
				<< ", diffusion = " << setprecision(5) << MAX_DIFF_FLUX_INCREMENT
				<< ", reaction = " << setprecision(5) << MAX_REAC_FLUX_INCREMENT
				<< ", phaseTrans = " << setprecision(5) << MAX_PHI_TRANS_FLUX_INCREMENT << endl;
			log << "# MAX TEMP : increment = " << setprecision(5) << MAX_TEMP_INCREMENT << endl;
			log << "#====================================================================================================" << endl;
			log << endl << endl;
			writer.add_string_to_txt_and_screen(log.str(), LOG_FILE_NAME);
		}
	}
	void Solvers::output_to_screen(string _postprocessing) {
		if (current_istep % parameters.screen_output_step == 0) {
			//////////////////////output to screen/////////////////////////
			stringstream log;
			PhaseNode& inf_node = statistics_information_in_phaseMesh();
			log << "=====================================================================================================" << endl;
			log << timer::return_cunrrent_time_by_string();

			log.setf(ios::fixed);
			log << "# Simulation step_" << current_istep << " has been finished!" << endl;
			log << "# This " << parameters.screen_output_step << " steps used " << setprecision(3) << timer::interval_end(Solvers_Timer.t_interval_begin) << "(secs.), " << endl;
			log << "# Total " << current_istep << " steps used " << setprecision(3) << timer::total_duration_sec(Solvers_Timer.t_total_begin) << "(secs.)." << endl;
			log << "# Simulation time is " << setprecision(3) << real_time << " (secs.)" << endl;
			log << "#============================== The percentage of phase and component ===============================" << endl;

			stringstream log2;
			for (auto p = inf_node.begin(); p < inf_node.end(); p++) {
				log2 << "#  Phase_" << p->index << "(" << parameters.Phases[p->property].phi_name << ")" << " : " << setprecision(5) << p->phi * 100 << " %" << endl;

				//double matter = 0.0;
				for (auto c = p->x.begin(); c < p->x.end(); c++) {
					log2 << "#   Con_" << parameters.Components[c->index].name << " : " << c->value << endl;
					log2 << "    Potential : " << p->potential[c->index].value << endl;
					//matter += c->value;
				}
				//log << "#   Matter  : " << matter << "(mol/m^3)  " <<  endl;
				log2 << "-----------------------------------------------------------------------------------------------------" << endl;
			}
			if (parameters.Components.size() != 0) {
				log2 << "#----------------------------------------- Total  Component -----------------------------------------" << endl;
				for (auto c = parameters.Components.begin(); c < parameters.Components.end(); c++)
				{
					log2 << "#  Total Con_" << c->name << " : " << setprecision(8) << inf_node.x[c->index].value << endl;
				}
				log2 << "#------------------------------------------ Total Potential -----------------------------------------" << endl;
				for (auto c = parameters.Components.begin(); c < parameters.Components.end(); c++)
				{
					log2 << "#  Total Potential_" << c->name << " : " << setprecision(8) << inf_node.potential[c->index].value << endl;
				}
			}
			log2 << "#===================================== Phi & Con & Temp Field =======================================" << endl;
			log2 << "# MAX PHI : increment = " << setprecision(5) << MAX_PHI_INCREMENT << endl;
			log2 << "# MAX CON : increment = " << setprecision(5) << MAX_CON_INCREMENT
				<< ", diffusion = " << setprecision(5) << MAX_DIFF_FLUX_INCREMENT
				<< ", reaction = " << setprecision(5) << MAX_REAC_FLUX_INCREMENT
				<< ", phaseTrans = " << setprecision(5) << MAX_PHI_TRANS_FLUX_INCREMENT << endl;
			log2 << "# MAX TEMP : increment = " << setprecision(5) << MAX_TEMP_INCREMENT << endl;
			log2 << "#========================================= Modules Exec =============================================" << endl;
			log2 << _postprocessing;
			log2 << "#====================================================================================================" << endl;
			log2 << endl << endl;
			timer::interval_begin(Solvers_Timer.t_interval_begin);
			writer.add_string_to_txt_and_screen(log.str(), LOG_FILE_NAME);
			writer.add_string_to_txt_and_screen(log2.str(), LOG_FILE_NAME);
		}
	}
	void Solvers::output_to_vts(WriteVTSType _type) {
		if (current_istep % parameters.vts_output_step == 0 || _type == WriteVTSType::WVTSType_INIT || _type == WriteVTSType::WVTSType_END) {
			PhaseNode& inf_node = statistics_information_in_phaseMesh();
			// write scalar file
			ofstream fout;
			if (_type == WriteVTSType::WVTSType_INIT)
				writer.open_vts_scalar_file(fout, phaseMesh, "step" + to_string(parameters.begin_step));
			else
				writer.open_vts_scalar_file(fout, phaseMesh, "step" + to_string(current_istep));

			writer.write_scalar_grains(fout, phaseMesh);

			if (parameters.ConEType == ConEquationType::CEType_PhaseX || parameters.ConEType == ConEquationType::CEType_GrandP)
				C_Solver.summary_phix_to_x();
			if (parameters.ConEType == ConEquationType::CEType_PhaseX)
				C_Solver.summary_phip_to_p();

			for (auto _module = modules.begin(); _module < modules.end(); _module++)
				_module->write_scalar(fout, phaseMesh);

			writer.close_vts_file(fout, phaseMesh);
			//write vector file
			ofstream fout2;
			if (_type == WriteVTSType::WVTSType_INIT)
				writer.open_vts_vec3_file(fout2, phaseMesh, "step" + to_string(parameters.begin_step));
			else
				writer.open_vts_vec3_file(fout2, phaseMesh, "step" + to_string(current_istep));
			for (auto _module = modules.begin(); _module < modules.end(); _module++)
				_module->write_vec3(fout2, phaseMesh);
			writer.close_vts_file(fout2, phaseMesh);
		}
	}
	void Solvers::write_data_file() {
		if (current_istep % parameters.data_output_step == 0) {
			if (data_writer.write_dataFile(phaseMesh, "_step" + to_string(current_istep)))
				writer.add_string_to_txt_and_screen("> An data file has been saved ! \n", LOG_FILE_NAME);
			else
				writer.add_string_to_txt_and_screen("> Error, write data file failed ! \n", LOG_FILE_NAME);
		}
	}
	PhaseNode& Solvers::statistics_information_in_phaseMesh() {
		if (current_istep == statistics_step && current_istep != parameters.begin_step)
			return info_node;
		statistics_step = current_istep;
		double node_number = double(phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z);
		info_node = phaseMesh.info_node;
		info_node.temperature.T = 0.0;
		if (parameters.ConEType == ConEquationType::CEType_PhaseX) {
			C_Solver.summary_phip_to_p();
			C_Solver.summary_phix_to_x();
		}
		for (auto p = info_node.begin(); p < info_node.end(); p++) {
			p->phi = 0.0;
			for (auto c = p->x.begin(); c < p->x.end(); c++) {
				c->value = 0.0;
				c->PhaseTransitionFlux = 0.0;
				c->ChemicalReactionFlux = 0.0;
				c->DiffusionFlux = 0.0;
			}
			for (auto c = p->potential.begin(); c < p->potential.end(); c++) {
				c->value = 0.0;
			}
		}
		for (auto c = info_node.x.begin(); c < info_node.x.end(); c++) {
			c->value = 0.0;
			c->PhaseTransitionFlux = 0.0;
			c->ChemicalReactionFlux = 0.0;
			c->DiffusionFlux = 0.0;
		}
		for (auto p = info_node.potential.begin(); p < info_node.potential.end(); p++) {
			p->value = 0.0;
		}

		for (int x = 0; x < phaseMesh.limit_x; x++)
			for (int y = 0; y < phaseMesh.limit_y; y++)
				for (int z = 0; z < phaseMesh.limit_z; z++) {
					PhaseNode& node = phaseMesh(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						for (auto p_con = phase->x.begin(); p_con < phase->x.end(); p_con++) {
							info_node[phase->index].x[p_con->index].value += phase->phi * p_con->value;
						}
						for (auto p_chem = phase->potential.begin(); p_chem < phase->potential.end(); p_chem++) {
							info_node[phase->index].potential[p_chem->index].value += phase->phi * p_chem->value;
						}
						info_node[phase->index].phi += phase->phi;
					}
					for (auto comp = node.x.begin(); comp < node.x.end(); comp++) {
						if (parameters.ConEType == ConEquationType::CEType_TotalX) {
							double smooth_phi = node.customValues[ExternalFields::CON_Smooth_Phi];
							info_node.x[comp->index].value += comp->value * smooth_phi;
							info_node.x[comp->index].ChemicalReactionFlux += comp->ChemicalReactionFlux * smooth_phi;
							info_node.x[comp->index].DiffusionFlux += comp->DiffusionFlux * smooth_phi;
							info_node.x[comp->index].PhaseTransitionFlux += comp->PhaseTransitionFlux * smooth_phi;
						}
						else {
							info_node.x[comp->index].value += comp->value;
							info_node.x[comp->index].ChemicalReactionFlux += comp->ChemicalReactionFlux;
							info_node.x[comp->index].DiffusionFlux += comp->DiffusionFlux;
							info_node.x[comp->index].PhaseTransitionFlux += comp->PhaseTransitionFlux;
						}
					}
					for (auto pot = node.potential.begin(); pot < node.potential.end(); pot++) {
						if (parameters.ConEType == ConEquationType::CEType_GrandP || parameters.ConEType == ConEquationType::CEType_TotalX) {
							if (parameters.ConEDomain == ConEquationDomain::CEDomain_Standard)
								info_node.potential[pot->index].value += pot->value * node.customValues[ExternalFields::CON_Smooth_Phi];
							else if (parameters.ConEDomain == ConEquationDomain::CEDomain_Reverse) {
								info_node.potential[pot->index].value += pot->value * (1.0 - node.customValues[ExternalFields::CON_Smooth_Phi]);
							}
						}
						else {
							info_node.potential[pot->index].value += pot->value;
						}
					}
				}

		for (auto p = info_node.begin(); p < info_node.end(); p++) {
			if (p->phi > Phi_Num_Cut_Off) {
				for (auto p_con = p->x.begin(); p_con < p->x.end(); p_con++) {
					p_con->value = p_con->value / p->phi;
				}
				for (auto p_chem = p->potential.begin(); p_chem < p->potential.end(); p_chem++) {
					p_chem->value /= p->phi;
				}
			}
			else {
				for (auto x = p->x.begin(); x < p->x.end(); x++)
					x->value = 0.0;
				for (auto p_chem = p->potential.begin(); p_chem < p->potential.end(); p_chem++) {
					p_chem->value = 0.0;
				}
			}

			p->phi = p->phi / node_number;
		}
		for (auto c = info_node.x.begin(); c < info_node.x.end(); c++) {
			c->value = c->value / node_number;
			c->ChemicalReactionFlux = c->ChemicalReactionFlux / node_number;
			c->DiffusionFlux = c->DiffusionFlux / node_number;
		}
		for (auto pot = info_node.potential.begin(); pot < info_node.potential.end(); pot++) {
			pot->value = pot->value / node_number;
		}
		return info_node;
	}
	string Solvers::exec_modules_in_loop() {
		stringstream exec_modules;
		for (auto _module = modules.begin(); _module < modules.end(); _module++)
			exec_modules << _module->exec_loop(phaseMesh);
		return exec_modules.str();
	}
	void Solvers::mid_selection() {
		printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		cout << endl;
		printf_color_on_control("> Select the one you want (press corresponding number and press enter, or press any other keys to exit):");
		cout << endl;
		cout << endl;
		// 1 activation information
		printf_color_on_control("> 1 MInDes: activation information:");
		printf_color_on_control("     activated", 32);
		cout << endl;
		cout << endl;
		// 2 license information
		printf_color_on_control("> 2 MInDes: license information:");
		printf_color_on_control("     activated", 32);
		cout << endl;
		cout << endl;
		// 3 User permissions
		printf_color_on_control("> 3 MInDes: user permissions:");
		cout << endl;
		cout << endl;

		// 3 cpu ID of this computer
		printf_color_on_control("> 4 MInDes: cpu ID ?");
		cout << endl;
		cout << endl;
		// 4 about copyright, developers & e-mail
		printf_color_on_control("> 5 MInDes: about");
		cout << endl;
		cout << endl;

		char c = getchar(); char enter = getchar(); //get enter
		if (c == '1')
			activation_info();
		else if (c == '2')
			license_info();
		else if (c == '3')
			permissions_info();
		else if (c == '4')
			cpu_id_info();
		else if (c == '5')
			about_info();
		exit(0);
	}
	void Solvers::activation_info() {
		printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		cout << endl;
		printf_color_on_control("> 1 MInDes activation information:");
		printf_color_on_control("     activated", 32);
		cout << endl;
		cout << endl;
		printf_color_on_control("> press \"B/b\" to go back to previous selection;");
		cout << endl;
		printf_color_on_control("> press any other keys to exit;");
		cout << endl;
		printf_color_on_control("> press enter to confirm.");
		cout << endl;
		char c = getchar(); char enter = getchar(); //get enter
		if (c == 'B' || c == 'b')
			mid_selection();
		exit(0);
	}
	void Solvers::license_info() {
		printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		cout << endl;
		printf_color_on_control("> 2 MInDes license information:");
		printf_color_on_control("     activated", 32);
		cout << endl;
		cout << endl;
		printf_color_on_control("> press \"B/b\" to go back to previous selection;");
		cout << endl;
		printf_color_on_control("> press any other keys to exit;");
		cout << endl;
		printf_color_on_control("> press enter to confirm.");
		cout << endl;
		char c = getchar(); char enter = getchar(); //get enter
		if (c == 'B' || c == 'b')
			mid_selection();
		exit(0);
	}
	void Solvers::permissions_info() {
		printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		cout << endl;
		printf_color_on_control("> 3 MInDes user permissions:");
		cout << endl;
		cout << endl;
		printf_color_on_control("> You get ");
		printf_color_on_control("full access", 32);
		printf_color_on_control(" to MInDes.");
		cout << endl;
		cout << endl;
		printf_color_on_control("> press \"B/b\" to go back to previous selection;");
		cout << endl;
		printf_color_on_control("> press any other keys to exit;");
		cout << endl;
		printf_color_on_control("> press enter to confirm.");
		cout << endl;
		char c = getchar(); char enter = getchar(); //get enter
		if (c == 'B' || c == 'b')
			mid_selection();
		exit(0);
	}
	void Solvers::cpu_id_info() {
		printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		cout << endl;
		printf_color_on_control("> 4 MInDes cpu ID ?");
		cout << endl;
		cout << endl;
		string cpu_id;
		printf_color_on_control("> cpu ID of this computer is:     ");
		printf_color_on_control(cpu_id, 37, 42);
		cout << endl;
		cout << endl;
		printf_color_on_control("> press \"B/b\" to go back to previous selection;");
		cout << endl;
		printf_color_on_control("> press any other keys to exit;");
		cout << endl;
		printf_color_on_control("> press enter to confirm.");
		cout << endl;
		char c = getchar(); char enter = getchar(); //get enter
		if (c == 'B' || c == 'b')
			mid_selection();
		exit(0);
	}
	void Solvers::about_info() {
		printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		cout << endl;
		printf_color_on_control("> 5 MInDes about"); cout << endl;
		cout << endl;
		printf_color_on_control("> Developer:");
		printf_color_on_control("Qi Huang", 34); cout << endl;
		cout << endl;
		printf_color_on_control("> Email:");
		printf_color_on_control("qihuang0908@163.com", 34); cout << endl;
		cout << endl;
		printf_color_on_control("> Copyright (c) ");
		printf_color_on_control("2019-2023 Science center for phase diagram, phase transition,", 34); cout << endl;
		printf_color_on_control(">               ");
		printf_color_on_control("material intelligent design and manufacture. Central South University. China", 34); cout << endl;
		cout << endl;
		printf_color_on_control("> This program is free software: you can redistribute it and/or modify it under the terms "); cout << endl;
		printf_color_on_control("> of the GNU General Public License as published by the Free Software Foundation, either "); cout << endl;
		printf_color_on_control("> version 3 of the License, or (at your option) any later version."); cout << endl;
		cout << endl;
		printf_color_on_control("> This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; "); cout << endl;
		printf_color_on_control("> without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "); cout << endl;
		printf_color_on_control("> See the GNU General Public License for more details."); cout << endl;
		cout << endl;
		printf_color_on_control("> You should have received a copy of the GNU General Public License along with this program. "); cout << endl;
		printf_color_on_control("> If not, see <http://www.gnu.org/licenses/>."); cout << endl;
		cout << endl;
		printf_color_on_control("    へ　　　　　／|"); cout << endl;
		printf_color_on_control("　　/＼7　　　 ∠＿/"); cout << endl;
		printf_color_on_control("　 /　│　　 ／　／"); cout << endl;
		printf_color_on_control("　│　Z ＿,＜　／　　 /`ヽ"); cout << endl;
		printf_color_on_control("　│　　　　　ヽ　　 /　　〉"); cout << endl;
		printf_color_on_control("　 Y　　　　　`　 /　　/"); cout << endl;
		printf_color_on_control("　ｲ●　､　●　　⊂⊃〈　　/"); cout << endl;
		printf_color_on_control("　() へ　　　　|　＼〈"); cout << endl;
		printf_color_on_control("　　>ｰ ､_　 ィ　 │ ／／"); cout << endl;
		printf_color_on_control("　 / へ　　 /　ﾉ＜| ＼＼       Pikachu says: let's go ! it's time to start a simulation !"); cout << endl;
		printf_color_on_control("　 ヽ_ﾉ　　(_／　 │／／"); cout << endl;
		printf_color_on_control("　　7　　　　　　　|／"); cout << endl;
		printf_color_on_control("　　＞―r￣￣`ｰ―＿) "); cout << endl;
		cout << endl;
		printf_color_on_control("> press \"B/b\" to go back to previous selection;");
		cout << endl;
		printf_color_on_control("> press any other keys to exit;");
		cout << endl;
		printf_color_on_control("> press enter to confirm.");
		cout << endl;
		char c = getchar(); char enter = getchar(); //get enter
		if (c == 'B' || c == 'b')
			mid_selection();
		exit(0);
	}
}