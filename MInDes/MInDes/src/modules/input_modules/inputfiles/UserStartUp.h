#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include "selectFile.h"
namespace pf {

	inline void print_help_doc() {
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

	StartOption arg_translator(std::string& _arg);

	void arg_process(std::vector<std::string> _arg_list, SimuInfo& _simu_info);

	std::vector<std::string> path_processor(std::string s_input);

	void user_interact_process(SimuInfo& _simu_info);

	SimuInfo User_StartUp(int _argc, char** _argv);

}