#pragma once
#include <iostream>
#include <string>
#include "../ProgramPath.h"
// Stores the path of the execution. Powered by "whereami.c" & whereami.h, by Gregory Pakosz
// Packaged with lambda expression.
namespace pf {
	struct SimuInfo {
		std::string simu_path{};
		std::string multi_mode{ "-S" };
		bool is_simu_ready{ false };
		int parallel_num{ 0 };

		void write_info() {
			std::ofstream fs_last_simu(program_path.string() + "path.in");
			fs_last_simu << simu_path << std::endl << multi_mode << std::endl << parallel_num;
			fs_last_simu.close();
		}
		void show_info() {
			pf::printf_color_on_control("Input File Path:\t" + simu_path + "\n", 34);
			pf::printf_color_on_control("MultiSImu Mode:\t" + multi_mode + "\n", 34);
			pf::printf_color_on_control("MultiSimu Thread Number:\t" + parallel_num, 34);
			printf("\n");
		}
	};
}