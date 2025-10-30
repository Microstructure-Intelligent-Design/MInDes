#pragma once
#include <filesystem>
#include "../ProgramPath.h"
#include <fstream>
namespace pf {
	inline bool Quick_StartUp(std::string& infile_path) {
		std::ifstream quick_simu_info(program_path / "start.in");
		if (quick_simu_info) {
			infile_path = "";
			getline(quick_simu_info, infile_path);
			quick_simu_info.close();
			std::filesystem::path p_last_path(infile_path);
			if (std::filesystem::exists(infile_path)) {
				std::filesystem::remove(program_path / "start.in");
				return true;
			}
		}
		return false;
	}
}