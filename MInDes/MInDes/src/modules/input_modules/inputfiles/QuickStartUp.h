#pragma once
#include <filesystem>
#include "../ProgramPath.h"
#include "SimuInfo.h"
#include <fstream>
namespace pf {
	bool Quick_StartUp(std::string& infile_path, bool& solver_run);
}