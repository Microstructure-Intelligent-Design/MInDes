#include "UserStartUp.h"
#include "../../base/license.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace pf {
	namespace {
		enum class StartOption {
			Default,
			CWD,
			LastSimu,
			Choose,
			Path,
			Number,
			Information,
			Help,
			Quit
		};

		void mindes_info();

		std::string trim_copy(const std::string& value) {
			const auto first = std::find_if_not(
				value.begin(), value.end(),
				[](unsigned char character) { return std::isspace(character) != 0; });
			const auto last = std::find_if_not(
				value.rbegin(), value.rend(),
				[](unsigned char character) { return std::isspace(character) != 0; }).base();
			if (first >= last) {
				return {};
			}
			return std::string(first, last);
		}

		bool input_single_char(char& value) {
			std::string input;
			while (std::getline(std::cin, input)) {
				input = trim_copy(input);
				if (input.size() == 1) {
					value = input.front();
					return true;
				}
				pf::printf_color_on_control("only one character accepted, please re-enter.");
				std::cout << std::endl;
			}
			return false;
		}

		void wait_for_enter() {
			pf::printf_color_on_control("> press enter to return to MInDes information.");
			std::cout << std::endl;
			std::string ignored;
			std::getline(std::cin, ignored);
		}

		void print_info() {
			std::cout << "" << std::endl;
			std::cout << " @@@@@@@@@   @@@@@@@@@   @@@@@@@   @@@@@    @@@@@@   @@@@@@@@@@@@@     @@@@@@@@@@@          @@@@@   " << std::endl;
			std::cout << " @@@@@@@@@   @@@@@@@@@   @@@@@@@    @@@@@   @@@@@@   @@@@@@@@@@@@@@@   @@@@@@@@@@@          @@@@@@@ " << std::endl;
			std::cout << " @@@@@@@@@   @@@@@@@@@   @@@@@@@    @@@@@@  @@@@@@   @@@@@@  @@@@@@@   @@@@@@@@@@@           @@@@@@ " << std::endl;
			std::cout << " @@@@@@@@@@ @@@@@@@@@@   @@@@@@@     @@@@@  @@@@@@   @@@@@@   @@@@@@                @@@@@@@  @@@@@@ " << std::endl;
			std::cout << " @@@@@@@@@@ @@@@@@@@@@   @@@@@@@      @@@@@ @@@@@@   @@@@@@   @@@@@@                 @@@@@@@        " << std::endl;
			std::cout << " @@@@@@@@@@ @@@@@@@@@@   @@@@@@@      @@@@@@@@@@@@   @@@@@@   @@@@@@    @@@@@@@@@@   @@@@@@@@@@     " << std::endl;
			std::cout << " @@@@@@ @@@@@@@ @@@@@@   @@@@@@@       @@@@@@@@@@@   @@@@@@   @@@@@@    @@@@@@@@@@     @@@@@@@@@@@  " << std::endl;
			std::cout << " @@@@@@ @@@@@@@ @@@ @    @@@@@@@   @    @@@@@@@@@@   @@@@@@   @@@@@@    @@@@@@@@@@        @@@@@@@@@ " << std::endl;
			std::cout << " @@@@@@ @@@@@@@ @@ @     @@@@@@@   @@    @@@@@@@@@   @@@@@@   @@@@@@                        @@@@@@@@" << std::endl;
			std::cout << " @@@@@@ @@@@@@@ @  @ @   @@@@@@@   @@    @@@@@@@@@   @@@@@@   @@@@@@                 @@@@@@  @@@@@@@" << std::endl;
			std::cout << " @@@@@@  @@@@@  @ @ @@   @@@@@@@   @@@    @@@@@@@@   @@@@@@  @@@@@@@                 @@@@@@   @@@@@@" << std::endl;
			std::cout << " @@@@@@  @@@@@   @@@@@   @@@@@@@   @@@@    @@@@@@@   @@@@@@@@@@@@@@@   @@@@@@@@@@@@  @@@@@@@@@@@@@@@" << std::endl;
			std::cout << " @@@@@@  @@@@@   @@@@@   @@@@@@@   @@@@@   @@@@@@@   @@@@@@@@@@@@@@    @@@@@@@@@@@@   @@@@@@@@@@@@  " << std::endl;
			std::cout << "" << std::endl;
		};

		void print_help_doc() {
			std::cout << "************************************ HELP DOC ************************************" << std::endl;
			std::cout << "-H      this help" << std::endl;
			std::cout << "-I      MInDes information and management" << std::endl;
			std::cout << "-D      Print current working directory" << std::endl;
			std::cout << "-L      Use the last simu's infomation" << std::endl;
			std::cout << "-C      Choose a input file" << std::endl;
			std::cout << "\"path\"  input the path to start simu; if path contains space, please quote it" << std::endl;
			std::cout << "-Q      Exit the program" << std::endl;
			std::cout << "**********************************************************************************" << std::endl;
		}

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
			if (_lower_arg == "-I") {
				return StartOption::Information;
			}
			if (_lower_arg == "-D") {
				return StartOption::CWD;
			}
			if (_lower_arg == "-H") {
				return StartOption::Help;
			}
			if (_lower_arg == "-Q") {
				return StartOption::Quit;
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
			case StartOption::Information: {
				print_info();
				mindes_info();
				break;
			}
			case StartOption::Help: {
				print_help_doc();
				break;
			}
			case StartOption::Quit: {
				std::exit(EXIT_SUCCESS);
			}
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
			if (_simu_info.is_simu_ready == true)
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

			std::vector<std::string> command_list{ path_processor(s_raw_commands) };
			arg_process(command_list, _simu_info);
		}
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

	namespace {
	void activation_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 1 MInDes activation information:");
		auto& license = pf::License::instance();
		if (license.is_activated())
			pf::printf_color_on_control("     activated", 32);
		else
			pf::printf_color_on_control("     inactivated", 31);
		std::cout << std::endl;
		std::cout << std::endl;
		license.print_activation_info();
		std::cout << std::endl;
		pf::printf_color_on_control("> press \"A/a\" to register this user and generate mindes.request;");
		std::cout << std::endl;
		pf::printf_color_on_control("> press any other keys to go back to previous selection.");
		std::cout << std::endl;
		char option = '\0';
		if (!input_single_char(option)) {
			return;
		}
		switch (option) {
		case 'a':
		case 'A': {
			if (license.is_license()) {
				pf::printf_color_on_control("> Re-registering creates a new request and invalidates the current license. Continue? [Y/N]", 33);
				std::cout << std::endl;
				char confirmation = '\0';
				if (!input_single_char(confirmation) || (confirmation != 'y' && confirmation != 'Y')) {
					break;
				}
			}
			const auto result = license.activate_this_user();
			pf::printf_color_on_control("> " + result.message, result.success ? 32 : 31);
			std::cout << std::endl;
			if (result.success) {
				pf::printf_color_on_control("> Request file: " + result.request_file.string(), 36);
				std::cout << std::endl;
			}
			wait_for_enter();
			break;
		}
		default:
			break;
		}
	}
	void license_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 2 MInDes license information:");
		auto& license = pf::License::instance();
		if (license.is_license())
			pf::printf_color_on_control("     valid", 32);
		else
			pf::printf_color_on_control("     invalid", 31);
		std::cout << std::endl;
		std::cout << std::endl;
		license.print_license_info();
		std::cout << std::endl;
		pf::printf_color_on_control("> Production MInDes verifies signed licenses but cannot generate them.");
		std::cout << std::endl;
		pf::printf_color_on_control("> Send mindes.request to the developer and place the returned mindes.license beside MInDes.");
		std::cout << std::endl;
		wait_for_enter();
	}
	void permissions_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 3 MInDes user permissions:");
		std::cout << "\n" << std::endl;
		auto& license = pf::License::instance();
		pf::printf_color_on_control("> Activation: ");
		pf::printf_color_on_control(
			license.is_activated() ? "activated" : "inactivated",
			license.is_activated() ? 32 : 31);
		std::cout << std::endl;
		pf::printf_color_on_control("> License:    ");
		pf::printf_color_on_control(
			license.is_license() ? "valid" : "invalid",
			license.is_license() ? 32 : 31);
		std::cout << std::endl << std::endl;
		if (license.is_activated() && license.is_license()) {
			pf::printf_color_on_control("> MInDes has the required permission for normal solver execution.", 32);
		}
		else {
			pf::printf_color_on_control("> Normal solver execution will be restricted until activation and license checks pass.", 31);
		}
		std::cout << "\n" << std::endl;
		wait_for_enter();
	}
	void machine_binding_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 4 MInDes machine binding information");
		std::cout << "\n" << std::endl;
		pf::MachineIdentity identity;
		if (pf::License::instance().get_machine_identity(identity)) {
			pf::printf_color_on_control("> Binding type:     ");
			pf::printf_color_on_control(pf::binding_type_name(identity.binding_type), 36);
			std::cout << std::endl;
			pf::printf_color_on_control("> Machine code:     ");
			pf::printf_color_on_control(identity.display_code, 37, 42);
		}
		else {
			pf::printf_color_on_control("> Can't obtain a reliable TPM/SMBIOS machine identity.", 31);
		}
		std::cout << "\n" << std::endl;
		wait_for_enter();
	}
	void about_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 5 MInDes about"); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> Developers:");
		pf::printf_color_on_control("[1] Qi Huang; [2] Zihang Wang;", 34); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> Emails:");
		pf::printf_color_on_control("[1] qihuang0908@163.com; [2] w.zihang@qq.com;", 34); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> Copyright (c) ");
		pf::printf_color_on_control("2019-2024 Science center for phase diagram, phase transition,", 34); std::cout << std::endl;
		pf::printf_color_on_control(">               ");
		pf::printf_color_on_control("material intelligent design and manufacture. Central South University. China", 34); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> This program is free software: you can redistribute it and/or modify it under the terms "); std::cout << std::endl;
		pf::printf_color_on_control("> of the GNU General Public License as published by the Free Software Foundation, either "); std::cout << std::endl;
		pf::printf_color_on_control("> version 3 of the License, or (at your option) any later version."); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; "); std::cout << std::endl;
		pf::printf_color_on_control("> without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "); std::cout << std::endl;
		pf::printf_color_on_control("> See the GNU General Public License for more details."); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> You should have received a copy of the GNU General Public License along with this program. "); std::cout << std::endl;
		pf::printf_color_on_control("> If not, see <http://www.gnu.org/licenses/>."); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("    へ　　　　　／|", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　　/＼7　　　 ∠＿/", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　 /　│　　 ／　／", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　│　Z ＿,＜　／　　 /`ヽ", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　│　　　　　ヽ　　 /　　〉", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　 Y　　　　　`　 /　　/", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("   ● .　●　　  | 〈　　/", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　() へ　 ()　 /　＼〈", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　　>- ,_　 ィ　 │ ／／", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　 / へ　　 /　/＜| ＼＼       Pikachu says: let's go ! it's time to start a simulation !", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　 ヽ_/　　(_／　 │／／", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　　7　　　　　　　|／", 33, 40); std::cout << std::endl;
		pf::printf_color_on_control("　　＞―r￣￣`-―＿) ", 33, 40); std::cout << std::endl;
		std::cout << std::endl;
		wait_for_enter();
	}

	void mindes_info() {
		auto& license = pf::License::instance();
		while (true) {
			license.check_mid_active(false);
			license.check_license(false);

			pf::printf_color_on_control("------------------------------------------------------------------------------------------------------\n", 34);
			pf::printf_color_on_control(">>> >  >    >      >          >          MInDes   INFORMATION          <         <      <    <  < <<<", 34);
			std::cout << std::endl;
			pf::printf_color_on_control("> Select the one you want (press corresponding number and press enter):");
			std::cout << "\n" << std::endl;

			pf::printf_color_on_control("> 1 MInDes: activation information:");
			pf::printf_color_on_control(
				license.is_activated() ? "     activated" : "     inactivated",
				license.is_activated() ? 32 : 31);
			std::cout << "\n" << std::endl;

			pf::printf_color_on_control("> 2 MInDes: license information:");
			pf::printf_color_on_control(
				license.is_license() ? "     valid" : "     invalid",
				license.is_license() ? 32 : 31);
			std::cout << "\n" << std::endl;

			pf::printf_color_on_control("> 3 MInDes: user permissions:");
			std::cout << "\n" << std::endl;
			pf::printf_color_on_control("> 4 MInDes: machine binding information");
			std::cout << "\n" << std::endl;
			pf::printf_color_on_control("> 5 MInDes: about");
			std::cout << "\n" << std::endl;
			pf::printf_color_on_control("> 6 back to previous menu");
			std::cout << "\n" << std::endl;
			pf::printf_color_on_control("> 7 Quit MInDes");
			std::cout << "\n" << std::endl;

			char option = '\0';
			if (!input_single_char(option)) {
				return;
			}
			switch (option) {
			case '1':
				activation_info();
				break;
			case '2':
				license_info();
				break;
			case '3':
				permissions_info();
				break;
			case '4':
				machine_binding_info();
				break;
			case '5':
				about_info();
				break;
			case '6':
				return;
			case '7':
				std::exit(EXIT_SUCCESS);
			default:
				std::cout << "invalid option, please re-enter an option." << std::endl;
				break;
			}
		}
	}
	}
}
