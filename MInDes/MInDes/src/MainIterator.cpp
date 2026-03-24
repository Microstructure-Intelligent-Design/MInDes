#include "MainIterator.h"
namespace pf {
	using namespace std;
	namespace iterator_times {
		string print_time_interval() {
			stringstream report;
			report << endl;
			report << timer::return_cunrrent_time_by_string();
			report << ">------------------------------------------- Time Interval ------------------------------------------" << endl;
			report << fixed << setprecision(3);
			report << "> Modules::init()         =  " << main_iterator::t_interval_modules_init << " (sec.)  "
				<< timer::time_min(main_iterator::t_interval_modules_init) << " (min.)  " << timer::time_hour(main_iterator::t_interval_modules_init) << " (h.) " << endl;
			report << "> Modules::pre_exec()     =  " << main_iterator::t_interval_modules_pre_exec << " (sec.)  "
				<< timer::time_min(main_iterator::t_interval_modules_pre_exec) << " (min.)  " << timer::time_hour(main_iterator::t_interval_modules_pre_exec) << " (h.) " << endl;
			report << "> Solvers::exec()         =  " << main_iterator::t_interval_modules_exec << " (sec.)  "
				<< timer::time_min(main_iterator::t_interval_modules_exec) << " (min.)  " << timer::time_hour(main_iterator::t_interval_modules_exec) << " (h.) " << endl;
			report << "> Modules::pos_exec()     =  " << main_iterator::t_interval_modules_pos_exec << " (sec.)  "
				<< timer::time_min(main_iterator::t_interval_modules_pos_exec) << " (min.)  " << timer::time_hour(main_iterator::t_interval_modules_pos_exec) << " (h.) " << endl;
			report << "> Modules::deinit()       =  " << main_iterator::t_interval_modules_deinit << " (sec.)  "
				<< timer::time_min(main_iterator::t_interval_modules_deinit) << " (min.)  " << timer::time_hour(main_iterator::t_interval_modules_deinit) << " (h.) " << endl;
			report << "> Total                   =  " << timer::total_duration_sec(main_iterator::t_total_begin) << " (sec.)  "
				<< timer::total_duration_min(main_iterator::t_total_begin) << " (min.)  " << timer::total_duration_hour(main_iterator::t_total_begin) << " (h.) " << endl;
			report << ">----------------------------------------------------------------------------------------------------" << endl;
			return report.str();
		}
	}

	namespace main_iterator {

		void init_modules(int argc, char* argv[]) {
			// - for UTF-8 output win/linux & screen/file
			std::locale::global(std::locale(""));  // 使用系统本地化
			std::wcout.imbue(std::locale(""));
			// get the infile path
			if (!Quick_StartUp(input_output_files_parameters::InFile_Path, main_iterator::main_solver_on)) {
				SimuInfo simu_info{ User_StartUp(argc,argv) };
				input_output_files_parameters::InFile_Path = simu_info.simu_path;
			}
			// init InputFileReader
			std::filesystem::path cwd{ std::filesystem::current_path() };
			string selected_file_path = infile_path_selector(input_output_files_parameters::InFile_Path);
			std::filesystem::current_path(cwd);
			InputFileReader::get_instance()->init(selected_file_path, false, INT_MAX, ' ');
			// init working folder path and create the folder
			input_output_files_parameters::WorkingFolder_Path = pf::erase_tail_of_infile(input_output_files_parameters::InFile_Path);
			input_output_files_parameters::InFileFolder_Path = pf::GetFolderOfPath(input_output_files_parameters::WorkingFolder_Path);
#if defined(_WIN32)
			(void)_mkdir(input_output_files_parameters::WorkingFolder_Path.c_str());
#elif defined(__linux__)
			mkdir(input_output_files_parameters::WorkingFolder_Path.c_str(), 0777);
#endif
			// init log and debug files
			input_output_files_parameters::LogFile_Path = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "log.txt";
			input_output_files_parameters::DebugFile_Path = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "input_report.txt";
			write_string_to_file("", input_output_files_parameters::LogFile_Path);
			write_string_to_file("", input_output_files_parameters::DebugFile_Path);
			// init modules
			stringstream modules_output;
			timer::time_interval_precision_secs_begin(main_iterator::t_interval_modules_init);
			modules_output << endl;
			modules_output << timer::return_cunrrent_time_by_string();
			modules_output << ">------------------------------------------- " << "Modules Init" << " -------------------------------------------" << endl;
			WriteLog(modules_output.str());
			init_input_modules();
			register_all_modules();
			modules_output.str("");
			modules_output << ">--------------------------------------------" << "------------" << "--------------------------------------------" << endl << endl;
			WriteLog(modules_output.str());
			timer::time_interval_precision_secs_end(main_iterator::t_interval_modules_init);
		}

		void run() {
			stringstream modules_output;
			timer::init(main_iterator::t_total_begin);
			if (main_solver_on) {
				timer::time_interval_precision_secs_begin(main_iterator::t_interval_modules_pre_exec);
				modules_output << endl;
				modules_output << timer::return_cunrrent_time_by_string();
				modules_output << ">----------------------------------------- " << "Modules Pre-Exec" << " -----------------------------------------" << endl;
				WriteLog(modules_output.str());
				for (auto& module : module_list) {
					if (module.exec_pre_i) module.exec_pre_i();
				}
				for (auto& module : module_list) {
					if (module.exec_pre_ii) module.exec_pre_ii();
				}
				for (auto& module : module_list) {
					if (module.exec_pre_iii) module.exec_pre_iii();
				}
				modules_output.str("");
				modules_output << ">-------------------------------------------" << "----------------" << "-----------------------------------------" << endl << endl;
				WriteLog(modules_output.str());
				timer::time_interval_precision_secs_end(main_iterator::t_interval_modules_pre_exec);

				// main loop;
				main_iterator::t_interval_modules_exec = 0.0;
				main_iterator::t_interval_modules_pos_exec = 0.0;
				timer::interval_begin(main_iterator::t_interval_begin);
				for (size_t istep = ITE_Begin_Step + 1; istep <= ITE_End_Step; istep++) {
					main_iterator::Current_ITE_step = istep;
					// - license
					if (main_iterator::OpenMP_Thread_Counts >= omp_get_num_procs())
						main_iterator::OpenMP_Thread_Counts = omp_get_num_procs() - 1;
					else if (main_iterator::OpenMP_Thread_Counts < 1)
						main_iterator::OpenMP_Thread_Counts = 1;
					omp_set_num_threads(main_iterator::OpenMP_Thread_Counts);
					// - license

					double cal_times = 0.0;

					timer::time_interval_precision_secs_begin(cal_times);
					for (auto& module : module_list) {
						if (module.exec_i) module.exec_i();
					}
					for (auto& module : module_list) {
						if (module.exec_ii) module.exec_ii();
					}
					for (auto& module : module_list) {
						if (module.exec_iii) module.exec_iii();
					}
					timer::time_interval_precision_secs_end(cal_times);
					main_iterator::t_interval_modules_exec += cal_times;

					timer::time_interval_precision_secs_begin(cal_times);
					for (auto& module : module_list) {
						if (module.exec_pos_i) module.exec_pos_i();
					}
					for (auto& module : module_list) {
						if (module.exec_pos_ii) module.exec_pos_ii();
					}
					for (auto& module : module_list) {
						if (module.exec_pos_iii) module.exec_pos_iii();
					}
					timer::time_interval_precision_secs_end(cal_times);
					main_iterator::t_interval_modules_pos_exec += cal_times;
				}
			}
			timer::time_interval_precision_secs_begin(main_iterator::t_interval_modules_deinit);
			modules_output.str("");
			modules_output << endl;
			modules_output << timer::return_cunrrent_time_by_string();
			modules_output << ">------------------------------------------ " << "Modules Deinit" << " ------------------------------------------" << endl;
			WriteLog(modules_output.str());
			for (auto& module : module_list) {
				if (module.deinit) module.deinit();
			}
			modules_output.str("");
			modules_output << ">-------------------------------------------" << "----------------" << "-----------------------------------------" << endl << endl;
			WriteLog(modules_output.str());
			timer::time_interval_precision_secs_end(main_iterator::t_interval_modules_deinit);
			WriteLog(iterator_times::print_time_interval());
		}
	};

}