#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include "selectFile.h"
namespace pf {

	void print_help_doc() {
		std::cout << "************************************ HELP DOC ************************************" << std::endl;
		std::cout << "-H      this help" << std::endl;
		std::cout << "-D      Print current working directory" << std::endl;
		std::cout << "-L      Use the last simu's infomation" << std::endl;
		std::cout << "-C      Choose a input file" << std::endl;
		std::cout << "\"path\"  input the path to start simu; if path contains space, please quote it" << std::endl;
		std::cout << "-Q      Exit the program" << std::endl;
		std::cout << "**********************************************************************************" << std::endl;
	};

	enum class StartOption { Default, Sequencial, Parallel, CWD, LastSimu, Choose, Path, Number };

	inline int DefaultMultiNum = 3;

	StartOption arg_translator(std::string& _arg) {
		std::string _lower_arg{ _arg };
		bool is_number = true;
		for (auto c : _lower_arg) {
			if (isdigit(c) == 0) {
				is_number = false;
				break;
			}
		}
		if (is_number) {
			return StartOption::Number;
		}
		else {
			std::transform(_arg.begin(), _arg.end(), _lower_arg.begin(), ::toupper);
			if (_lower_arg == "-L") {
				return StartOption::LastSimu;
			}
			if (_lower_arg == "-C") {
#ifdef _WIN32
				return StartOption::Choose;
#else
				pf::printf_color_on_control("Please enter the input file path\n", 34);
				std::getline(std::cin, _arg);
				return StartOption::Path;
#endif // _WIN32
			}
			if (_lower_arg == "-D") {
				return StartOption::CWD;
			}
			if (_lower_arg == "-H") {
				print_help_doc();
				return StartOption::Default;
			}
			if (_lower_arg == "-Q") {
				exit(0);
			}
			else {
				pf::str_char_delete(_arg, '\"');
				pf::str_char_delete(_arg, '\'');
				return StartOption::Path;
			}
		}
	}

	void arg_process(std::vector<std::string> _arg_list, SimuInfo& _simu_info) {
		for (int i = 0; i < _arg_list.size(); i++) {
			std::string arg = _arg_list.at(i);

			switch (arg_translator(arg)) {
			case StartOption::LastSimu: {
				std::fstream _fs_last_simu(program_path.string() + "path.in", std::ios::in);

				if (!_fs_last_simu.is_open()) {
					_fs_last_simu.close();
#ifdef _WIN32
					pf::printf_color_on_control("No simu record, please select a file\n", 34);
					selectPMFile(_simu_info);
					break;
#else
					pf::printf_color_on_control("No simu record, please enter the input file path (or with multi-instruction)\n", 34);
#endif // _WIN32
				}

				std::string old_arg{};
				std::vector<std::string> temp_arg_list{};
				while (std::getline(_fs_last_simu, old_arg)) {
					temp_arg_list.push_back(old_arg);
				}
				_fs_last_simu.close();
				arg_process(temp_arg_list, _simu_info);
				break;
			}
			case StartOption::Choose: {
#ifdef _WIN32
				selectPMFile(_simu_info);
#endif // _WIN32
				break;
			}
			case StartOption::Path: {
				std::fstream fs_test_path(arg);
				if (fs_test_path.is_open()) {
					fs_test_path.close();
					_simu_info.simu_path = arg;
					_simu_info.is_simu_ready = true;
					break;
				}
				else {
					pf::printf_color_on_control("invalid path!\n", 34);
					pf::printf_color_on_control("is your path contains space ? quote it please\n", 31);
#ifdef _WIN32
					pf::printf_color_on_control("please select a valid file: \n", 34);
					selectPMFile(_simu_info);
#else
					_simu_info.is_simu_ready = false;
#endif // _WIN32
					break;
				}
			}
			case StartOption::CWD: {
				pf::printf_color_on_control("Current working directory:");
				std::cout << std::endl << std::filesystem::current_path().string() << std::endl;
				pf::printf_color_on_control("Program installed directory:");
				std::cout << std::endl << program_path.string() << std::endl;
				break;
			}
			default:
				break;
			}
			if (_simu_info.is_simu_ready = true)
				break;
		}
	}

	std::vector<std::string> path_processor(std::string s_input) {
		std::istringstream token_stream{ s_input };
		std::vector<std::string> token_list{};
		std::string token{};
		while (token_stream >> token) {
			token_list.push_back(token);

			if (token_stream.eof()) {
				break;
			}
			if (*token.begin() == '\'' or *token.begin() == '\"') {

				auto path_start{ token_stream.tellg() };
				path_start -= (*token_list.rbegin()).size();
				auto path_end{ token_stream.tellg() };

				token_list.pop_back();

				while ((*(token.end() - 1) != '\'') and (*(token.end() - 1) != '\"') and !token_stream.eof()) {
					token_stream >> token;
				}
				path_end = token_stream.tellg();

				token_list.push_back(token_stream.str().substr(path_start, path_end - path_start));
			}
		}
		return token_list;
	}

	void user_interact_process(SimuInfo& _simu_info) {
		while (_simu_info.is_simu_ready == false) {
			pf::printf_color_on_control("-----------------------------------------------------\n");
			pf::printf_color_on_control("Please give a infile path or an option\n");
			pf::printf_color_on_control("For help, type \"-H\"\n");
			pf::printf_color_on_control("-----------------------------------------------------");
			std::cout << std::endl;
			std::string s_raw_commands{};
			std::getline(std::cin, s_raw_commands);
//#if !defined(_READLINE_H_) || !defined(_HISTORY_H_)
//			std::getline(std::cin, s_raw_commands);
//#else
//			char* c_str_raw_commands = readline("");
//			if (c_str_raw_commands) {
//				add_history(c_str_raw_commands);
//				s_raw_commands = c_str_raw_commands;
//				free(c_str_raw_commands);
//			}
//#endif

			std::vector<std::string> command_list{ path_processor(s_raw_commands) };
			arg_process(command_list, _simu_info);
		}
	}

	SimuInfo User_StartUp(int _argc, char** _argv) {
		std::ifstream last_simu_info(program_path / "path.in");
		if (last_simu_info) {
			std::string last_path{};
			getline(last_simu_info, last_path);
			std::filesystem::path p_last_path(last_path);
			if (!std::filesystem::exists(p_last_path)) {
				pf::printf_color_on_control("The last input file is not existed.\n");
				pf::printf_color_on_control("It might be moved or deleted.\n");
			}
			else {
				p_last_path = std::filesystem::canonical(p_last_path);
				pf::printf_color_on_control("The last input file is:\n");
				pf::printf_color_on_control(p_last_path.string() + "\n");
			}
		}
		std::vector<std::string> arg_list{};
		SimuInfo _simu_info{};
		for (int i = 1; i < _argc; i++) {
			arg_list.push_back(std::string(_argv[i]));
		}

		arg_process(arg_list, _simu_info);
		user_interact_process(_simu_info);
		_simu_info.write_info();

		return _simu_info;
	}

}