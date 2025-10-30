#pragma once
#include "../../base/MACRO_DEF.h"
#include <filesystem>
#include <algorithm>
namespace pf {

	void infile_line_process(std::string& str) {
		size_t start = 0;
		while (start < str.length() && std::isspace(str[start])) {
			start++;
		}
		str.erase(0, start);
		str_char_delete(str, '\r');
	}

	void convert_backslash(std::string& _backslash_string) {
		for (auto& c : _backslash_string) {
			if (c == '\\') {
				c = '/';
			}
		}
	}

	void invalid_path_exit(const std::filesystem::path bad_path, const std::filesystem::path parent_path = "") {
		pf::printf_color_on_control("file name error when reading\n");
		pf::printf_color_on_control(bad_path.string());
		if (!parent_path.string().empty()) {
			pf::printf_color_on_control("\nfrom\n" + parent_path.string());
		}
		pf::printf_color_on_control("\n,please check \'INCLUDE\' file name");
		std::cout << std::endl;
#ifdef _WIN32
		system("pause");
#else
		(void)getchar();
#endif //_WIN32	
		exit(0);
	}

	void recursive_read_files(std::filesystem::path _p_this_fpath, std::ofstream& _ofs_outfile, std::vector<std::filesystem::path>& _read_file_list) {

		std::ifstream ifs_this_file(_p_this_fpath);
		if (ifs_this_file.is_open()) {
			std::string line{};
			while (getline(ifs_this_file, line)) {
				infile_line_process(line);

				if (line.substr(0, line.find(" ")) == "INCLUDE") {// if it is nested file
					std::filesystem::current_path(_p_this_fpath.parent_path()); //setting pwd to this folder

					std::string s_included_fpath{ line.substr(line.find(" ") + 1) };
					convert_backslash(s_included_fpath);

					std::filesystem::path p_included_fpath(s_included_fpath); //store the next file's name
					if (!std::filesystem::exists(p_included_fpath)) {
						invalid_path_exit(p_included_fpath, _p_this_fpath);
						return;
					}
					p_included_fpath = std::filesystem::canonical(p_included_fpath);

					if (std::find(_read_file_list.begin(), _read_file_list.end(), p_included_fpath) == _read_file_list.end()) {//if not existed
						_read_file_list.push_back(p_included_fpath);
						recursive_read_files(p_included_fpath, _ofs_outfile, _read_file_list);
					}
					continue;
				}
				_ofs_outfile << line << std::endl;
			}
			ifs_this_file.close();
			return;
		}
		else {
			invalid_path_exit(_p_this_fpath);
			return;
		}
	}


	std::string infile_path_selector(std::string s_infile_path) { // read which one?

		convert_backslash(s_infile_path);
		std::filesystem::path p_infile_path(s_infile_path);
		std::filesystem::path p_infile_name_folder(p_infile_path);
		p_infile_name_folder.replace_extension();
		if (!std::filesystem::exists(p_infile_path)) {
			invalid_path_exit(p_infile_path);
		}
		p_infile_path = std::filesystem::canonical(p_infile_path);

		std::ifstream ifs_infile(p_infile_path);
		if (ifs_infile.is_open()) {
			std::string line{};
			while (getline(ifs_infile, line)) {
				infile_line_process(line);

				if (line.substr(0, line.find(" ")) == "INCLUDE") {// if it is nested file
					std::filesystem::path p_combined_fpath(p_infile_name_folder / std::filesystem::path("combined_infile.mindes"));
					std::ofstream ofs_combined(p_combined_fpath);
					if (!ofs_combined) {
						pf::printf_color_on_control("The file folder is invalid. Please check.\n");
						exit(1);
					}
					ofs_combined.close();
					ofs_combined.open(p_combined_fpath, std::ios::app);

					std::vector<std::filesystem::path> read_file_list{ p_infile_path };
					recursive_read_files(p_infile_path, ofs_combined, read_file_list);

					pf::printf_color_on_control("-----------------------------------------------------\n");
					pf::printf_color_on_control("Read from multiple files:\n");
					for (const auto& read_file_path : read_file_list) {
						pf::printf_color_on_control(read_file_path.string() + "\n");
					}
					pf::printf_color_on_control("-----------------------------------------------------");
					std::cout << std::endl;
					ifs_infile.close();
					return p_combined_fpath.string();
				}
			}
			ifs_infile.close();
			return p_infile_path.string();
		}
		pf::printf_color_on_control("Invalid file path, check the input file path!\n");
		exit(1);
		return{};
	}
}