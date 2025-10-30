#pragma once
#include <filesystem>
#include "../ProgramPath.h"
#include "SimuInfo.h"
#include <fstream>
namespace pf {
	inline bool Quick_StartUp(std::string& infile_path, bool& solver_run) {
		std::ifstream quick_simu_info(program_path / "start.in");
		SimuInfo simu_info;
		if (quick_simu_info) {
			infile_path = "";
			getline(quick_simu_info, infile_path);
			std::string str1 = "";
			getline(quick_simu_info, str1);
			quick_simu_info.close();
			std::filesystem::path p_last_path(infile_path);
			if (str1.compare("SOLVER_RUN") == 0)
				solver_run = true;
			else
				solver_run = false;
			if (std::filesystem::exists(p_last_path)) {
				std::filesystem::remove(program_path / "start.in");
				simu_info.simu_path = infile_path;
				simu_info.write_info();
				return true;
			}
		}
		return false;
	}
}