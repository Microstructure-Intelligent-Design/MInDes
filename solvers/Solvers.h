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
#include"Solver_poisson_equation.h"
#include"Solver_Allen_Cahn.h"
#include"Solver_Cahn_Hilliard.h"
#include"Solver_Temperature.h"
#include"Solver_Fluid.h"
#include"Solver_Mechanics.h"
#include"Solver_custom.h"
namespace pf {
	using namespace std;
	enum WriteVTSType { WVTSType_INIT, WVTSType_LOOP, WVTSType_END };
	enum CON_TEMP_SKIP { CTS_CON, CTS_TEMP };
	struct solvers_time {
		double t_interval_modules_init;
		double t_interval_modules_pre_exec;
		double t_interval_modules_loop_exec;
		double t_interval_output;
		double t_interval_evolve_phi;
		double t_interval_evolve_con;
		double t_interval_evolve_T;
		double t_interval_modules_deinit;

		double t_total_begin;
		double t_interval_begin;
		solvers_time() {
			t_interval_modules_init = 0.0;
			t_interval_modules_pre_exec = 0.0;
			t_interval_modules_loop_exec = 0.0;
			t_interval_output = 0.0;
			t_interval_evolve_phi = 0.0;
			t_interval_evolve_con = 0.0;
			t_interval_evolve_T = 0.0;
			t_interval_modules_deinit = 0.0;

			t_total_begin = 0.0;
			t_interval_begin = 0.0;
		}
		string print_time_interval() {
			stringstream report;
			report << endl;
			report << ">------------------- time interval (secs.) -------------------" << endl;
			report << fixed << setprecision(3);
			report << "> Modules::init()            = " << t_interval_modules_init << endl;
			report << "> Modules::pre_exec()        = " << t_interval_modules_pre_exec << endl;
			report << "> Solvers::evolve_phi()      = " << t_interval_evolve_phi << endl;
			report << "> Solvers::evolve_con()      = " << t_interval_evolve_con << endl;
			report << "> Solvers::evolve_T()        = " << t_interval_evolve_T << endl;
			report << "> Modules::loop_exec()       = " << t_interval_modules_loop_exec << endl;
			report << "> Solvers::output()          = " << t_interval_output << endl;
			report << "> Modules::deinit()          = " << t_interval_modules_deinit << endl;
			report << ">-------------------------------------------------------------" << endl;
			return report.str();
		}
	};
	struct solvers_parameters {
		PhiEquationType PhiEType;
		ConEquationType ConEType;
		ConEquationDomain ConEDomain;
		TemperatureEquationType TempEType;
		int OpenMP_Thread_Counts;
		DifferenceMethod Difference_Method;
		bool is_Normalize_Phi;
		bool is_Normalize_Con;
		double dt;
		vector<int> SKIP;
		int begin_step;
		int end_step;
		int screen_loop_step;
		int screen_output_step;
		int vts_output_step;
		int data_output_step;
		Info_Node Components;
		Info_Phases Phases;
		GrainsOrientations Grains;
		solvers_parameters() {
			PhiEType = PhiEquationType::PEType_Const;
			ConEType = ConEquationType::CEType_Const;
			ConEDomain = ConEquationDomain::CEDomain_Standard;
			TempEType = TemperatureEquationType::TType_Const;
			OpenMP_Thread_Counts = 1;
			Difference_Method = DifferenceMethod::FIVE_POINT;
			is_Normalize_Phi = false;
			is_Normalize_Con = true;
			dt = 1.0;
			SKIP.push_back(1); SKIP.push_back(1); // [CTS_CON, CTS_TEMP]
			begin_step = 0;
			end_step = 0;
			screen_loop_step = -1;
			screen_output_step = -1;
			vts_output_step = -1;
			data_output_step = -1;
			Components.clear();
			Phases.clear();
			Grains.init(RotationGauge::RG_XZX);
		}
		solvers_parameters& operator=(const solvers_parameters& n) {
			PhiEType = n.PhiEType;
			ConEType = n.ConEType;
			ConEDomain = n.ConEDomain;
			TempEType = n.TempEType;
			OpenMP_Thread_Counts = n.OpenMP_Thread_Counts;
			Difference_Method = n.Difference_Method;
			is_Normalize_Phi = n.is_Normalize_Phi;
			is_Normalize_Con = n.is_Normalize_Con;
			dt = n.dt;
			SKIP = n.SKIP;
			begin_step = n.begin_step;
			end_step = n.end_step;
			screen_loop_step = n.screen_loop_step;
			screen_output_step = n.screen_output_step;
			vts_output_step = n.vts_output_step;
			data_output_step = n.data_output_step;
			Components = n.Components;
			Phases = n.Phases;
			Grains = n.Grains;
			return *this;
		}
	};
	struct Solver_Module
	{
		// return 0
		void(*init)(FieldStorage_forPhaseNode&);
		// return 0
		void(*exec_pre)(FieldStorage_forPhaseNode&);
		// return report string
		string (*exec_loop)(FieldStorage_forPhaseNode&);
		// 
		// return report string
		void(*deinit)(FieldStorage_forPhaseNode&);
		void(*write_scalar)(ofstream&, FieldStorage_forPhaseNode&);
		void(*write_vec3)(ofstream&, FieldStorage_forPhaseNode&);
	};
	//< Memory Request and evolution of phi, con, T
	class Solvers
	{
	public:
		static Solvers* get_instance(){
			if(solver == NULL)
				solver = new Solvers();
			return solver;
		}
		void create_a_new_module(void(*init)(FieldStorage_forPhaseNode&), void(*exec_pre)(FieldStorage_forPhaseNode&),
			string (*exec_loop)(FieldStorage_forPhaseNode&), void(*deinit)(FieldStorage_forPhaseNode&),
			void(*write_scalar)(ofstream&, FieldStorage_forPhaseNode&), void(*write_vec3)(ofstream&, FieldStorage_forPhaseNode&)) {
			Solver_Module _module;
			_module.init = init;
			_module.exec_pre = exec_pre;
			_module.exec_loop = exec_loop;
			_module.deinit = deinit;
			_module.write_scalar = write_scalar;
			_module.write_vec3 = write_vec3;
			modules.push_back(_module);
		}
		void init(string input_file_path) {
			string working_folder_name = pf::erase_tail_of_input_file_name(pf::GetFileNameOfPath(input_file_path));
			Infile_Folder_Path = pf::GetFolderOfPath(input_file_path);
			Working_Folder_Path = Infile_Folder_Path + dirSeparator + working_folder_name;
			writer.init(Working_Folder_Path);
			writer.init_txt_file(LOG_FILE_NAME);
		}
		void run() {
			// timer
			timer::time_interval_precision_secs_begin(Solvers_Timer.t_interval_modules_init);
			// init treatment functions
			stringstream modules_init;
			modules_init << endl;
			modules_init << timer::return_cunrrent_time_by_string();
			modules_init << ">--------------------------------------------" << "Modules Init" << "--------------------------------------------" << endl;
			writer.add_string_to_txt_and_screen(modules_init.str(), LOG_FILE_NAME);
			for (auto _module = modules.begin(); _module < modules.end(); _module++)
				_module->init(phaseMesh);
			modules_init.str("");
			modules_init << ">--------------------------------------------" << "------------" << "--------------------------------------------" << endl;
			writer.add_string_to_txt_and_screen(modules_init.str(), LOG_FILE_NAME);
			// timer
			timer::time_interval_precision_secs_end(Solvers_Timer.t_interval_modules_init);

			stringstream log;
			timer::init(Solvers_Timer.t_total_begin);
			timer::interval_begin(Solvers_Timer.t_interval_begin);

			current_istep = parameters.begin_step;

			// timer
			timer::time_interval_precision_secs_begin(Solvers_Timer.t_interval_modules_pre_exec);
			modules_pre_exec();
			timer::time_interval_precision_secs_end(Solvers_Timer.t_interval_modules_pre_exec);

			// timer
			double output_times = 0.0;
			timer::time_interval_precision_secs_begin(output_times);

			if (parameters.screen_output_step > 0)
				output_to_screen("");
			if (parameters.vts_output_step > 0)
				output_to_vts(WriteVTSType::WVTSType_INIT);
			timer::time_interval_precision_secs_end(output_times);
			Solvers_Timer.t_interval_output += output_times;

			real_time = 0.0;
			for (int istep = parameters.begin_step + 1; istep <= parameters.end_step; istep++) {
				current_istep = istep;
				real_time += parameters.dt;
				MAX_PHI_INCREMENT = 0.0;
				MAX_CON_INCREMENT = 0.0;
				MAX_TEMP_INCREMENT = 0.0;
				MAX_DIFF_FLUX_INCREMENT = 0.0;
				MAX_REAC_FLUX_INCREMENT = 0.0;
				MAX_PHI_TRANS_FLUX_INCREMENT = 0.0;

				solve_phi_x_T_in_loop();

				// timer
				timer::time_interval_precision_secs_begin(output_times);
				string execloop_str = exec_modules_in_loop();
				timer::time_interval_precision_secs_end(output_times);
				Solvers_Timer.t_interval_modules_loop_exec += output_times;

				// timer
				timer::time_interval_precision_secs_begin(output_times);

				if (parameters.screen_loop_step > 0)
					output_with_loop();

				if (parameters.screen_output_step > 0)
					output_to_screen(execloop_str);

				if (parameters.vts_output_step > 0)
					output_to_vts(WriteVTSType::WVTSType_LOOP);

				if (parameters.data_output_step > 0)
					write_data_file();

				timer::time_interval_precision_secs_end(output_times);
				Solvers_Timer.t_interval_output += output_times;

				if (is_stop_loop)
					break;
			}
			// timer
			timer::time_interval_precision_secs_begin(output_times);

			if (parameters.data_output_step > 0) {
				writer.add_string_to_txt_and_screen("> Error, write data file failed at the end of Program! \n", LOG_FILE_NAME);
			}
			if (parameters.vts_output_step > 0)
				output_to_vts(WriteVTSType::WVTSType_END);
			timer::time_interval_precision_secs_end(output_times);
			Solvers_Timer.t_interval_output += output_times;


			timer::time_interval_precision_secs_begin(Solvers_Timer.t_interval_modules_deinit);
			// init treatment functions
			stringstream modules_deinit;
			modules_deinit << endl;
			modules_deinit << timer::return_cunrrent_time_by_string();
			modules_deinit << ">--------------------------------------------" << "Modules Deinit" << "--------------------------------------------" << endl;
			writer.add_string_to_txt_and_screen(modules_deinit.str(), LOG_FILE_NAME);
			for (auto _module = modules.begin(); _module < modules.end(); _module++)
				_module->deinit(phaseMesh);
			modules_deinit.str("");
			modules_deinit << ">--------------------------------------------" << "------------" << "--------------------------------------------" << endl;
			writer.add_string_to_txt_and_screen(modules_deinit.str(), LOG_FILE_NAME);
			// timer
			timer::time_interval_precision_secs_end(Solvers_Timer.t_interval_modules_deinit);

			log.str("");
			log << Solvers_Timer.print_time_interval() << endl;
			log << timer::return_cunrrent_time_by_string();
			log << "##### Simulation End! The simulation used time " << setprecision(4) << timer::total_duration_sec(Solvers_Timer.t_total_begin) << "(s)! #####" << endl << endl;
			log << endl;
			writer.add_string_to_txt_and_screen(log.str(), LOG_FILE_NAME);
		}
		void mid_info() {
			printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
			cout << endl;
			printf_color_on_control(">>> >  >    >      >          >          MInDes   INFORMATION          <         <      <    <  < <<<", 34);
			cout << endl;
			mid_selection();
		}
		PhaseNode& statistics_information_in_phaseMesh();

		AllenCahnSolver		Phi_Solver_AC;
		CahnHilliardSolver	Phi_Solver_CH;
		CahnHilliardSolver	C_Solver;
		TemperatureSolver	T_Solver;
		// parameters adjust here
		string Working_Folder_Path;
		string Infile_Folder_Path;
		WriteToFile writer;
		int current_istep;
		double real_time;
		// parameters adjusted elsewhere
		FieldStorage_forPhaseNode phaseMesh;
		solvers_parameters parameters;
		// field variables change rate
		double MAX_PHI_INCREMENT;
		double MAX_CON_INCREMENT;
		double MAX_TEMP_INCREMENT;
		double MAX_DIFF_FLUX_INCREMENT;
		double MAX_REAC_FLUX_INCREMENT;
		double MAX_PHI_TRANS_FLUX_INCREMENT;

		bool is_stop_loop;
	private:
		vector<Solver_Module> modules;
		solvers_time Solvers_Timer;
		Solvers();
		~Solvers();
		// statistics
		PhaseNode info_node;
		int statistics_step;
		void modules_pre_exec();
		void solve_phi_x_T_in_loop();
		string exec_modules_in_loop();
		void output_with_loop();
		void output_to_screen(string _postprocessing);
		void output_to_vts(WriteVTSType _type);
		void write_data_file();
		void mid_selection();
		void activation_info();
		void license_info();
		void permissions_info();
		void cpu_id_info();
		void about_info();
		// module
		static Solvers* solver;
	};
	
}
