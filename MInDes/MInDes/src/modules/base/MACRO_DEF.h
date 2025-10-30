#pragma once
#include <limits>
#define _CRT_SECURE_NO_WARNINGS
#ifdef _WIN32
#include <conio.h>
#endif
#include <iostream>
#include <fstream>
#include <string>
// control the REAL data
// #define USE_DOUBLE
#undef max // 取消 max 宏定义
#undef min // 取消 min 宏定义
#ifdef USE_DOUBLE
using REAL = double; // use double

#define SYS_EPSILON   (0.000001)
#define SYS_EPSILON_R (0.999999)

#define Phi_Num_Cut_Off   (0.001)
#define Phi_Num_Cut_Off_R (0.999)

#define PhiCon_Num_Cut_Off   (0.01)
#define PhiCon_Num_Cut_Off_R (0.99)

inline double NaN() { return std::numeric_limits<double>::max(); };
inline double REAL_MAX() { return std::numeric_limits<double>::max(); };

#define PI (3.1415926535897932)

#define AngleToRadians(angle) double(angle/180.0*PI)

#else
using REAL = float;  // use float

#define SYS_EPSILON   (0.000001f)
#define SYS_EPSILON_R (0.999999f)

#define Phi_Num_Cut_Off   (0.001f)
#define Phi_Num_Cut_Off_R (0.999f)

#define PhiCon_Num_Cut_Off   (0.01f)
#define PhiCon_Num_Cut_Off_R (0.99f)

inline float NaN() { return std::numeric_limits<float>::max(); };

inline float REAL_MAX() { return std::numeric_limits<float>::max(); };

#define PI (3.1415927f)

#define AngleToRadians(angle) float(angle/180.0f*PI)

#endif


#define SYS_PROGRAM_STOP std::exit(1);

#if defined(_WIN32)
#define dirSeparator std::string("\\")                                     //< Windows style directory separator
#elif defined(__linux__)
#define dirSeparator std::string("/")                                     //< Windows style directory separator
#endif

namespace pf {
	inline bool isTwoREALEquality(float a, float b) {
		if ((a - b) < SYS_EPSILON && (a - b) > -SYS_EPSILON)
			return true;
		else
			return false;
	}
	inline bool isTwoREALEquality(double a, double b) {
		if ((a - b) < SYS_EPSILON && (a - b) > -SYS_EPSILON)
			return true;
		else
			return false;
	}
	inline int REAL_to_int(float a) {
		if ((a - int(a)) > float(0.5))
			return int(a) + 1;
		else
			return int(a);
	}
	inline int REAL_to_int(double a) {
		if ((a - int(a)) > double(0.5))
			return int(a) + 1;
		else
			return int(a);
	}
	inline void str_char_delete(std::string& str, char _word) {
		std::size_t start_position{ str.find(_word) };
		while (start_position != std::string::npos) {
			str.erase(start_position, 1);
			start_position = str.find(_word);
		}
	}
	inline std::string GetFolderOfPath(std::string file_path) {
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
	inline std::string GetFileNameOfPath(std::string file_path) {
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
	inline std::string erase_tail_of_infile(std::string input_file_name) {
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
	inline char get_char_not_show() {
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
	inline std::string get_string_from_consol(bool is_show = true, char replace_char = '*') {
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
	inline void printf_color_on_control(std::string str, int front_color = 30, int back_color = 43) {
		printf("\x1b[%d;%dm%s\x1b[0m", back_color, front_color, str.c_str());
	}
	inline void write_string_to_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path);
		if (!fout) {
			std::cout << "Failed to write the txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content << std::endl;
		fout.close();
	}
	inline void add_string_to_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path, std::ios::app);
		if (!fout) {
			std::cout << "Failed to add the string to txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
	}
	inline void add_string_to_screen_and_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path, std::ios::app);
		if (!fout) {
			std::cout << "Failed to add the string to txt file!" << std::endl;
			fout.close();
			return;
		}
		std::cout << content;
		fout << content;
		fout.close();
	}

}