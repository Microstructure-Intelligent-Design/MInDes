#pragma once
#include <limits>
#include <cstddef>
#define _CRT_SECURE_NO_WARNINGS
#ifdef _WIN32
#include <conio.h>
#include <Windows.h>
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
	bool isTwoREALEquality(float a, float b);
	bool isTwoREALEquality(double a, double b);
	int REAL_to_int(float a);
	int REAL_to_int(double a);
	void str_char_delete(std::string& str, char _word);
	std::string GetFolderOfPath(std::string file_path);
	std::string GetFileNameOfPath(std::string file_path);
	std::string erase_tail_of_infile(std::string input_file_name);
	char get_char_not_show();
	std::string get_string_from_consol(bool is_show = true, char replace_char = '*');
	void printf_color_on_control(std::string str, int front_color = 30, int back_color = 43);
	void write_string_to_file(std::string content, std::string file_path);
	void add_string_to_file(std::string content, std::string file_path);
	void add_string_to_screen_and_file(std::string content, std::string file_path);
}