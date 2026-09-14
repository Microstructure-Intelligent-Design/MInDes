#pragma once
#include "../../base/MACRO_DEF.h"
#include <filesystem>
#include <algorithm>
#include <vector>
namespace pf {

	void infile_line_process(std::string& str);

	void convert_backslash(std::string& _backslash_string);

	void invalid_path_exit(const std::filesystem::path bad_path, const std::filesystem::path parent_path = "");

	void recursive_read_files(std::filesystem::path _p_this_fpath, std::ofstream& _ofs_outfile, std::vector<std::filesystem::path>& _read_file_list);

	std::string infile_path_selector(std::string s_infile_path);
}
