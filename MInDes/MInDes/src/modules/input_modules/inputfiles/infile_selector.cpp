#include "infile_selector.h"
#include "../ProgramPath.h"
#include <cctype>
#include <cerrno>
#include <cstring>
#include <utility>

#ifdef _WIN32
#include <Windows.h>
#elif defined(__linux__)
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#endif

namespace pf {
	namespace {
		enum class PythonPreprocessResult {
			Success,
			NoIncludes,
			Unavailable,
			InputError
		};

#ifdef _WIN32
		struct SystemPythonCandidate {
			std::wstring executable;
			std::vector<std::wstring> prefix_arguments;
			std::string display_name;
		};

		std::wstring quote_windows_argument(const std::wstring& argument) {
			if (argument.empty())
				return L"\"\"";
			if (argument.find_first_of(L" \t\n\v\"") == std::wstring::npos)
				return argument;

			std::wstring quoted(1, L'\"');
			size_t backslashes = 0;
			for (const wchar_t character : argument) {
				if (character == L'\\') {
					++backslashes;
				}
				else if (character == L'\"') {
					quoted.append(backslashes * 2 + 1, L'\\');
					quoted.push_back(L'\"');
					backslashes = 0;
				}
				else {
					quoted.append(backslashes, L'\\');
					quoted.push_back(character);
					backslashes = 0;
				}
			}
			quoted.append(backslashes * 2, L'\\');
			quoted.push_back(L'\"');
			return quoted;
		}

		PythonPreprocessResult run_python_preprocessor(
			const std::wstring& executable,
			const std::vector<std::wstring>& prefix_arguments,
			const std::string& display_name,
			const std::filesystem::path& script,
			const std::filesystem::path& input,
			const std::filesystem::path& output) {
			std::wstring command_line = quote_windows_argument(executable);
			for (const auto& argument : prefix_arguments)
				command_line += L" " + quote_windows_argument(argument);
			command_line += L" " + quote_windows_argument(script.native()) + L" --input "
				+ quote_windows_argument(input.native()) + L" --output "
				+ quote_windows_argument(output.native());
			std::vector<wchar_t> mutable_command(command_line.begin(), command_line.end());
			mutable_command.push_back(L'\0');

			STARTUPINFOW startup_info{};
			startup_info.cb = sizeof(startup_info);
			PROCESS_INFORMATION process_info{};
			const BOOL started = CreateProcessW(
				nullptr, mutable_command.data(), nullptr, nullptr, FALSE, 0,
				nullptr, program_path.c_str(), &startup_info, &process_info);
			if (!started) {
				pf::printf_color_on_control(
					"Failed to start system Python candidate '" + display_name + "' (Windows error "
					+ std::to_string(GetLastError()) + ").\n");
				return PythonPreprocessResult::Unavailable;
			}

			WaitForSingleObject(process_info.hProcess, INFINITE);
			DWORD exit_code = 1;
			GetExitCodeProcess(process_info.hProcess, &exit_code);
			CloseHandle(process_info.hThread);
			CloseHandle(process_info.hProcess);
			if (exit_code == 0)
				return PythonPreprocessResult::Success;
			if (exit_code == 10)
				return PythonPreprocessResult::NoIncludes;
			if (exit_code == 20)
				return PythonPreprocessResult::InputError;
			pf::printf_color_on_control(
				"System Python candidate '" + display_name + "' exited unexpectedly (exit code "
				+ std::to_string(exit_code) + ").\n");
			return PythonPreprocessResult::Unavailable;
		}
#elif defined(__linux__)
		PythonPreprocessResult run_python_preprocessor(
			const std::string& executable,
			const std::vector<std::string>& prefix_arguments,
			const std::filesystem::path& script,
			const std::filesystem::path& input,
			const std::filesystem::path& output) {
			std::vector<std::string> argument_storage{ executable };
			argument_storage.insert(argument_storage.end(), prefix_arguments.begin(), prefix_arguments.end());
			argument_storage.push_back(script.string());
			argument_storage.push_back("--input");
			argument_storage.push_back(input.string());
			argument_storage.push_back("--output");
			argument_storage.push_back(output.string());

			std::vector<char*> arguments;
			for (auto& argument : argument_storage)
				arguments.push_back(argument.data());
			arguments.push_back(nullptr);

			pid_t child = -1;
			const int spawn_error = posix_spawnp(
				&child, executable.c_str(), nullptr, nullptr, arguments.data(), environ);
			if (spawn_error != 0) {
				pf::printf_color_on_control(
					"Failed to start system Python candidate '" + executable
					+ "': " + std::strerror(spawn_error) + "\n");
				return PythonPreprocessResult::Unavailable;
			}

			int status = 0;
			while (waitpid(child, &status, 0) == -1) {
				if (errno == EINTR)
					continue;
				pf::printf_color_on_control(
					"Failed to wait for system Python candidate '" + executable
					+ "': " + std::strerror(errno) + "\n");
				return PythonPreprocessResult::Unavailable;
			}
			if (WIFSIGNALED(status)) {
				pf::printf_color_on_control(
					"System Python candidate '" + executable + "' was terminated by signal "
					+ std::to_string(WTERMSIG(status)) + ".\n");
				return PythonPreprocessResult::Unavailable;
			}
			if (!WIFEXITED(status))
				return PythonPreprocessResult::Unavailable;

			const int exit_code = WEXITSTATUS(status);
			if (exit_code == 0)
				return PythonPreprocessResult::Success;
			if (exit_code == 10)
				return PythonPreprocessResult::NoIncludes;
			if (exit_code == 20)
				return PythonPreprocessResult::InputError;
			pf::printf_color_on_control(
				"System Python candidate '" + executable
				+ "' exited unexpectedly (exit code " + std::to_string(exit_code) + ").\n");
			return PythonPreprocessResult::Unavailable;
		}
#endif

		PythonPreprocessResult preprocess_with_system_python(
			const std::filesystem::path& input,
			const std::filesystem::path& output) {
			const auto script = program_path / "tools" / "mindes_preprocess" / "preprocess_input.py";
			std::error_code ec;
			if (!std::filesystem::is_regular_file(script, ec)) {
				pf::printf_color_on_control("Python preprocessor script was not found: " + script.string() + "\n");
				return PythonPreprocessResult::Unavailable;
			}
#ifdef _WIN32
			const std::vector<SystemPythonCandidate> candidates = {
				{ L"py.exe", { L"-3" }, "py -3" },
				{ L"python.exe", {}, "python" }
			};
#elif defined(__linux__)
			const std::vector<std::pair<std::string, std::vector<std::string>>> candidates = {
				{ "python3", {} },
				{ "python", {} }
			};
#else
			(void)input;
			(void)output;
			return PythonPreprocessResult::Unavailable;
#endif
#if defined(_WIN32) || defined(__linux__)
			for (const auto& candidate : candidates) {
#ifdef _WIN32
				const auto result = run_python_preprocessor(
					candidate.executable, candidate.prefix_arguments, candidate.display_name,
					script, input, output);
#else
				const auto result = run_python_preprocessor(
					candidate.first, candidate.second, script, input, output);
#endif
				if (result != PythonPreprocessResult::Unavailable)
					return result;
			}
			return PythonPreprocessResult::Unavailable;
#endif
		}
	}

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

	void invalid_path_exit(const std::filesystem::path bad_path, const std::filesystem::path parent_path) {
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
		if (!std::filesystem::exists(p_infile_path)) {
			invalid_path_exit(p_infile_path);
		}
		p_infile_path = std::filesystem::canonical(p_infile_path);
		std::filesystem::path p_infile_name_folder(p_infile_path);
		p_infile_name_folder.replace_extension();
		const std::filesystem::path p_combined_fpath = p_infile_name_folder / "combined_infile.mindes";

		const auto python_result = preprocess_with_system_python(p_infile_path, p_combined_fpath);
		if (python_result == PythonPreprocessResult::Success) {
			std::cout << "Infile Parser : tools/mindes_preprocess.\n" << std::endl;
			return p_combined_fpath.string();
		}
		if (python_result == PythonPreprocessResult::NoIncludes) {
			std::cout << "Infile Parser : tools/mindes_preprocess.\n" << std::endl;
			return p_infile_path.string();
		}
		if (python_result == PythonPreprocessResult::InputError) {
			std::cout << "System Python input preprocessing failed.\n" << std::endl;
			exit(1);
		}

		std::cout << "Infile Parser Default : Mindes.exe.\n" << std::endl;

		std::ifstream ifs_infile(p_infile_path);
		if (ifs_infile.is_open()) {
			std::string line{};
			while (getline(ifs_infile, line)) {
				infile_line_process(line);

				if (line.substr(0, line.find(" ")) == "INCLUDE") {// if it is nested file
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
