#pragma once
#include "ProgramPath.h"
#include "../base/MACRO_DEF.h"
namespace pf {
	namespace input_output_files_parameters {
		inline std::string Program_Path = program_path.string();
		inline std::string InFile_Path = "";
		inline std::string LogFile_Path = "";
		inline std::string DebugFile_Path = "";
		inline std::string WorkingFolder_Path = "";
		inline std::string InFileFolder_Path = "";
	}
#define WriteDebugFile(str)  add_string_to_file(str, input_output_files_parameters::DebugFile_Path)
#define WriteLog(str)  add_string_to_screen_and_file(str, input_output_files_parameters::LogFile_Path)
}