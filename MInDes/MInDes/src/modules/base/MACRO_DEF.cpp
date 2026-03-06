#include "MACRO_DEF.h"
namespace pf {
	bool isTwoREALEquality(float a, float b) {
		if ((a - b) < SYS_EPSILON && (a - b) > -SYS_EPSILON)
			return true;
		else
			return false;
	}
	bool isTwoREALEquality(double a, double b) {
		if ((a - b) < SYS_EPSILON && (a - b) > -SYS_EPSILON)
			return true;
		else
			return false;
	}
	int REAL_to_int(float a) {
		if ((a - int(a)) > float(0.5))
			return int(a) + 1;
		else
			return int(a);
	}
	int REAL_to_int(double a) {
		if ((a - int(a)) > double(0.5))
			return int(a) + 1;
		else
			return int(a);
	}
	void str_char_delete(std::string& str, char _word) {
		std::size_t start_position{ str.find(_word) };
		while (start_position != std::string::npos) {
			str.erase(start_position, 1);
			start_position = str.find(_word);
		}
	}
	std::string GetFolderOfPath(std::string file_path) {
		while (file_path.size() != 0)
		{
			char _c = *(file_path.end() - 1);
			if (_c == '/' || _c == '\\') {
				file_path.erase(file_path.end() - 1);
				break;
			}
			else {
				file_path.erase(file_path.end() - 1);
			}
		}
		return file_path;
	}
	std::string GetFileNameOfPath(std::string file_path) {
		std::string name = "";
		int index = 1;
		while (index <= file_path.size())
		{
			char _c = *(file_path.end() - index);
			if (_c == '/' || _c == '\\') {
				break;
			}
			else {
				index++;
				name = _c + name;
			}
		}
		return name;
	}
	std::string erase_tail_of_infile(std::string input_file_name) {
		std::string tail = ".mindes", name_without_tail = input_file_name;
		int index = 0;
		bool is_name_correct = true;
		while (name_without_tail.size() != 0 && index < tail.size())
		{
			index++;
			if (*(name_without_tail.end() - 1) != *(tail.end() - index)) {
				is_name_correct = false;
				break;
			}
			name_without_tail.erase(name_without_tail.end() - 1);
		}
		if (name_without_tail.size() == 0)
			is_name_correct = false;
		if (!is_name_correct) {
			std::cout << "> input file name error, file name = " << input_file_name << ", aim tail is " << tail << std::endl;
			SYS_PROGRAM_STOP;
		}
		return name_without_tail;
	}
	char get_char_not_show() {
		char c;
#ifdef _WIN32
		c = _getch();
#else
		system("stty -echo");
		c = getchar();
		system("stty echo");
#endif
		return c;
	}
	std::string get_string_from_consol(bool is_show, char replace_char) {
		std::string str = "";
		int i = 0;
		if (is_show) {
			while (true) {
				char ch = getchar();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
			}
		}
		else {
			while (true) {
				char ch = get_char_not_show();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
				putchar(replace_char);
			}
		}
		putchar('\n');
		return str;
	};
	void printf_color_on_control(std::string str, int front_color, int back_color) {
		printf("\x1b[%d;%dm%s\x1b[0m", back_color, front_color, str.c_str());
	}
	void write_string_to_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path);
		if (!fout) {
			std::cout << "Failed to write the txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content << std::endl;
		fout.close();
	}
	void add_string_to_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path, std::ios::app);
		if (!fout) {
			std::cout << "Failed to add the string to txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
	}
	void add_string_to_screen_and_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path, std::ios::app);
		if (!fout) {
			std::cout << "Failed to add the string to txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
#ifdef _WIN32
		int len = MultiByteToWideChar(CP_UTF8, 0, content.data(), (int)content.size(), nullptr, 0);
		if (len > 0) {
			std::wstring wstr;
			wstr.resize(len);
			MultiByteToWideChar(CP_UTF8, 0, content.data(), (int)content.size(), &wstr[0], len);
			std::wcout << wstr;
		}
		else {
			// UTF-8 failed
			std::cout << content;
		}
#else
		// Linux/macOS
		std::cout << content;
#endif
	}
}