#pragma once
#include <iostream>
#include <string>
#include "SimuInfo.h"
#include "../../base/MACRO_DEF.h"
#ifdef _WIN32
#include <windows.h>
namespace pf {
	void selectPMFile(SimuInfo& _simu_info);
	void selectDataFile(std::string& file_path);
	void selectAllFile(std::string& file_path);
}
#endif // _WIN32
