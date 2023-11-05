/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include "../../solvers/Solvers.h"
#include "InputFileMath.h"
namespace pf {
	using namespace std;
	const char matrix_separator = ',';
	enum InputValueType { IVType_INT, IVType_DOUBLE, IVType_BOOL, IVType_STRING };
	struct input_value {
		int int_value;
		double double_value;
		bool bool_value;
		string string_value;
		input_value() {
			int_value = 0;
			double_value = 0.0;
			bool_value = false;
			string_value = "";
		}
		~input_value() {};
	};
	class InputFileReader
	{
	public:
		static InputFileReader* get_instance() {
			if (infile == NULL)
				infile = new InputFileReader();
			return infile;
		}
		InputFileReader() {
			_split = ' ';
			input_file.clear();
		}
		~InputFileReader() {
			_split = ' ';
			debug_writer = nullptr;
			infileMath.clear();
			input_file.clear();
			_valid_words.clear();
		}
		void init(string input_file_name, WriteToFile& writter, bool debug = false, int file_max_lines = 1000, char split = ' ', string debug_file_name = "input_report") {
			input_file.clear();
			_split = split;
			defined_parameters_number = 0;
			debug_writer = &writter;
			debug_file = debug_file_name;
			debug_writer->init_txt_file(debug_file);
			std::fstream fin(input_file_name, std::ios::in);
			if (debug)
				cout << "-----------------------------------------------------------------------------------------" << endl;
			if (!fin) {
				cout << "> Failed to read the input_file:" << input_file_name << ", use default value 0 (int, double, bool), \"\" (string)" << endl;
			}
			else {
				int line = 1;
				string strline;
				while (getline(fin, strline) && line <= file_max_lines) {
					input_file.add_string(line, strline);
					line++;
				}
			}
			fin.close();
			get_valid_words();
			get_define_func_and_variable();
			if (infileMath.check_variables_funcs_keys() == false) {
				cout << "> Defined key of variables and functions cannot be the same." << endl;
				SYS_PROGRAM_STOP;
			}
			if (debug) {
				debug_infile_and_valid_words();
				debug_custom_variavle_and_funcs();
			}
			// Remove duplicate definitions
			if(debug)
				cout << "> Reading data from input file! Following is [keyword] + [Value], the valid keyword works." << endl;
			for (auto vec1 = _valid_words.begin(); vec1 < _valid_words.end(); vec1++) {
				if (debug)
					cout << "Valid   :		[" << (*vec1)[0] << "] + [" << (*vec1)[1] << "]" << endl;
				if ((*vec1)[0].compare(infileMath.define_func_key) == 0 || (*vec1)[0].compare(infileMath.define_variable_key) == 0)
					continue;
				for (auto vec2 = vec1 + 1; vec2 < _valid_words.end();) {
					if ((*vec1)[0].compare((*vec2)[0]) == 0) {
						if (debug)
							cout << "in-Valid:		[" << (*vec2)[0] << "] + [" << (*vec2)[1] << "]" << endl;
						vec2 = _valid_words.erase(vec2);
					}
					else
						++vec2;
				}
			}
			if (debug)
				cout << "-----------------------------------------------------------------------------------------" << endl;
		}
		
		void auto_read_infile(string input_file_name, int file_max_lines = 1000) {
			input_file.clear();
			std::fstream fin(input_file_name, std::ios::in);
			if (!fin) {
				cout << "> Failed to read the input_file:" << input_file_name << ", use default value 0 (int, double, bool), \"\" (string)" << endl;
			}
			else {
				int line = 1;
				string strline;
				while (getline(fin, strline) && line <= file_max_lines) {
					input_file.add_string(line, strline);
					line++;
				}
			}
			fin.close();
		}
		bool auto_read_int_value(string value_name, int& int_value, bool debug = false) {
			stringstream report;
			for (int line = 1; line <= input_file.size(); line++) {
				string str = input_file[line];
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else {
						if (*(s->end() - 1) == '\r' || *(s->end() - 1) == '\t' || *(s->end() - 1) == '\n' || *(s->end() - 1) == '\0')
							s->erase(s->end() - 1);
						++s;
					}
				}
				if (head == '#' || head == ' ' || head == '\t' || head == '\n' || head == '\r' || head == '\0' || vec3.size() != 4)
					continue;
				else if (vec3[0].compare("<AUTO>") != 0)
					continue;
				else if (vec3[2].compare("=") != 0)
					continue;
				else {
					if (value_name.compare(vec3[1]) == 0) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_int_value(int_value, "<AUTO> " + vec3[1], vec3[3]), debug_file);
						else
							trans_string_to_int_value(int_value, "<AUTO> " + vec3[1], vec3[3]);
						return true;
					}
				}
			}
			if (debug)
				debug_writer->add_string_to_txt(string("> [DEFAULT] <AUTO> " + value_name + " = " + to_string(int_value) + "\n"), debug_file);
			return false;
		}
		bool auto_read_double_value(string value_name, double& double_value, bool debug = false) {
			stringstream report;
			for (int line = 1; line <= input_file.size(); line++) {
				string str = input_file[line];
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else {
						if (*(s->end() - 1) == '\r' || *(s->end() - 1) == '\t' || *(s->end() - 1) == '\n' || *(s->end() - 1) == '\0')
							s->erase(s->end() - 1);
						++s;
					}
				}
				if (head == '#' || head == ' ' || head == '\t' || head == '\n' || head == '\r' || head == '\0' || vec3.size() != 4)
					continue;
				else if (vec3[0].compare("<AUTO>") != 0)
					continue;
				else if (vec3[2].compare("=") != 0)
					continue;
				else {
					if (value_name.compare(vec3[1]) == 0) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_double_value(double_value, "<AUTO> " + vec3[1], vec3[3]), debug_file);
						else
							trans_string_to_double_value(double_value, "<AUTO> " + vec3[1], vec3[3]);
						return true;
					}
				}
			}
			if (debug)
				debug_writer->add_string_to_txt(string("> [DEFAULT] <AUTO> " + value_name + " = " + to_string(double_value) + "\n"), debug_file);
			return false;
		}
		bool auto_read_bool_value(string value_name, bool& bool_value, bool debug = false) {
			stringstream report;
			for (int line = 1; line <= input_file.size(); line++) {
				string str = input_file[line];
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else {
						if (*(s->end() - 1) == '\r' || *(s->end() - 1) == '\t' || *(s->end() - 1) == '\n' || *(s->end() - 1) == '\0')
							s->erase(s->end() - 1);
						++s;
					}
				}
				if (head == '#' || head == ' ' || head == '\t' || head == '\n' || head == '\r' || head == '\0' || vec3.size() != 4)
					continue;
				else if (vec3[0].compare("<AUTO>") != 0)
					continue;
				else if (vec3[2].compare("=") != 0)
					continue;
				else {
					if (value_name.compare(vec3[1]) == 0) {
						bool_value = false;
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_bool_value(bool_value, "<AUTO> " + vec3[1], vec3[3]), debug_file);
						else
							trans_string_to_bool_value(bool_value, "<AUTO> " + vec3[1], vec3[3]);
						return true;
					}
				}
			}
			if (bool_value) {
				if (debug)
					debug_writer->add_string_to_txt(string("> [DEFAULT] <AUTO> " + value_name + " = TRUE" + "\n"), debug_file);
			}
			else {
				if (debug)
					debug_writer->add_string_to_txt(string("> [DEFAULT] <AUTO> " + value_name + " = FALSE" + "\n"), debug_file);
			}
			return false;
		}

		bool read_int_value(string value_name, int& int_value, bool debug = false) {
			stringstream report;
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_int_value(int_value, (*vec)[0], (*vec)[1]), debug_file);
					else
						trans_string_to_int_value(int_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (debug)
				debug_writer->add_string_to_txt(string("> [DEFAULT] " + value_name + " = " + to_string(int_value) + "\n"), debug_file);
			return false;
		}
		bool read_double_value(string value_name, double& double_value, bool debug = false) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_double_value(double_value, (*vec)[0], (*vec)[1]), debug_file);
					else
						trans_string_to_double_value(double_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (debug)
				debug_writer->add_string_to_txt(string("> [DEFAULT] " + value_name + " = " + to_string(double_value) + "\n"), debug_file);
			return false;
		}
		bool read_bool_value(string value_name, bool& bool_value, bool debug = false) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					bool_value = false;
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_bool_value(bool_value, (*vec)[0], (*vec)[1]), debug_file);
					else
						trans_string_to_bool_value(bool_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (bool_value) {
				if (debug)
					debug_writer->add_string_to_txt(string("> [DEFAULT] " + value_name + " = TRUE" + "\n"), debug_file);
			}
			else {
				if (debug)
					debug_writer->add_string_to_txt(string("> [DEFAULT] " + value_name + " = FALSE" + "\n"), debug_file);
			}
			return false;
		}
		bool read_string_value(string value_name, string& string_value, bool debug = false) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_string_value(string_value, (*vec)[0], (*vec)[1]), debug_file);
					else
						trans_string_to_string_value(string_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (debug)
				debug_writer->add_string_to_txt(string("> [DEFAULT] " + value_name + " = " + string_value + "\n"), debug_file);
			return false;
		}
		
		string trans_string_to_int_value(int& val, string str_key, string str_value) {
			stringstream report;
			try
			{
				val = stoi(str_value);
				report << "> [-VALID-] " << str_key << " = " << str_value << endl;
			}
			catch (double)
			{
				cout << "> Input file error! the value of Keyword: " << str_key << " cannot be translate to int value." << endl;
				SYS_PROGRAM_STOP;
			}
			return report.str();
		}
		string trans_string_to_double_value(double& val, string str_key, string str_value) {
			stringstream report;
			int index = -1;
			if (infile_math_default_funcs::is_string_double(str_value, val)) {
				;
			}
			else if (infileMath.search_var(str_value, index)) {
				val = infileMath.infile_vars[index].var;
			}
			else if (infileMath.search_func(str_value, index)) {
				PhaseNode node;
				vector<int> para_int; vector<double> para_double;
				InFileFunc& inFunc = infileMath.infile_funcs[index];
				val = inFunc.func(node, para_int, para_double, inFunc.func_structure, inFunc.operators, inFunc.terms_type, infileMath.infile_vars, infileMath.infile_funcs);
			}
			else {
				cout << "> Input file error! the value of Keyword: " << str_key << " cannot be translate to double value." << endl;
				SYS_PROGRAM_STOP;
			}
			report << "> [-VALID-] " << str_key << " = " << val << endl;
			return report.str();
		}
		string trans_string_to_bool_value(bool& val, string str_key, string str_value) {
			stringstream report;
			if (str_value.compare("true") == 0 || str_value.compare("TRUE") == 0 || str_value.compare("1") == 0)
				val = true;
			if (val) {
				report << "> [-VALID-] " << str_key << " = TRUE" << endl;
			}
			else {
				report << "> [-VALID-] " << str_key << " = FALSE" << endl;
			}
			return report.str();
		}
		string trans_string_to_string_value(string& val, string str_key, string str_value) {
			stringstream report;
			val = str_value;
			report << "> [-VALID-] " << str_key << " = " << str_value << endl;
			return report.str();
		}
		
		vector<input_value> trans_matrix_1d_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '(') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == ')') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<input_value> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator) {
					if (str.size() == 0) {
						cout << "> Input matirx error ! matrix_1d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					input_value in_val;
					string str_key = key + "(" + to_string(index) + ")";
					if (_type[index] == InputValueType::IVType_INT) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_int_value(in_val.int_value, str_key, str), debug_file);
						else
							trans_string_to_int_value(in_val.int_value, str_key, str);
					}
					else if (_type[index] == InputValueType::IVType_DOUBLE) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_double_value(in_val.double_value, str_key, str), debug_file);
						else
							trans_string_to_double_value(in_val.double_value, str_key, str);
					}
					else if (_type[index] == InputValueType::IVType_BOOL) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_bool_value(in_val.bool_value, str_key, str), debug_file);
						else
							trans_string_to_bool_value(in_val.bool_value, str_key, str);
					}
					else if (_type[index] == InputValueType::IVType_STRING) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_string_value(in_val.string_value, str_key, str), debug_file);
						else
							trans_string_to_string_value(in_val.string_value, str_key, str);
					}
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_1d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			input_value in_val;
			string str_key = key + "(" + to_string(index) + ")";
			if (_type[index] == InputValueType::IVType_INT) {
				if (debug)
					debug_writer->add_string_to_txt(trans_string_to_int_value(in_val.int_value, str_key, str), debug_file);
				else
					trans_string_to_int_value(in_val.int_value, str_key, str);
			}
			else if (_type[index] == InputValueType::IVType_DOUBLE) {
				if (debug)
					debug_writer->add_string_to_txt(trans_string_to_double_value(in_val.double_value, str_key, str), debug_file);
				else
					trans_string_to_double_value(in_val.double_value, str_key, str);
			}
			else if (_type[index] == InputValueType::IVType_BOOL) {
				if (debug)
					debug_writer->add_string_to_txt(trans_string_to_bool_value(in_val.bool_value, str_key, str), debug_file);
				else
					trans_string_to_bool_value(in_val.bool_value, str_key, str);
			}
			else if (_type[index] == InputValueType::IVType_STRING) {
				if (debug)
					debug_writer->add_string_to_txt(trans_string_to_string_value(in_val.string_value, str_key, str), debug_file);
				else
					trans_string_to_string_value(in_val.string_value, str_key, str);
			}
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_1d = " << input_string << ", its size = " << index + 1 << ", defined type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<input_value> trans_matrix_1d_const_to_input_value(InputValueType _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '(') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == ')') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<input_value> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator) {
					if (str.size() == 0) {
						cout << "> Input matirx error ! matrix_1d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					input_value in_val;
					string str_key = key + "(" + to_string(index) + ")";
					if (_type == InputValueType::IVType_INT) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_int_value(in_val.int_value, str_key, str), debug_file);
						else
							trans_string_to_int_value(in_val.int_value, str_key, str);
					}
					else if (_type == InputValueType::IVType_DOUBLE) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_double_value(in_val.double_value, str_key, str), debug_file);
						else
							trans_string_to_double_value(in_val.double_value, str_key, str);
					}
					else if (_type == InputValueType::IVType_BOOL) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_bool_value(in_val.bool_value, str_key, str), debug_file);
						else
							trans_string_to_bool_value(in_val.bool_value, str_key, str);
					}
					else if (_type == InputValueType::IVType_STRING) {
						if (debug)
							debug_writer->add_string_to_txt(trans_string_to_string_value(in_val.string_value, str_key, str), debug_file);
						else
							trans_string_to_string_value(in_val.string_value, str_key, str);
					}
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_1d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				input_value in_val;
				string str_key = key + "(" + to_string(index) + ")";
				if (_type == InputValueType::IVType_INT) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_int_value(in_val.int_value, str_key, str), debug_file);
					else
						trans_string_to_int_value(in_val.int_value, str_key, str);
				}
				else if (_type == InputValueType::IVType_DOUBLE) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_double_value(in_val.double_value, str_key, str), debug_file);
					else
						trans_string_to_double_value(in_val.double_value, str_key, str);
				}
				else if (_type == InputValueType::IVType_BOOL) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_bool_value(in_val.bool_value, str_key, str), debug_file);
					else
						trans_string_to_bool_value(in_val.bool_value, str_key, str);
				}
				else if (_type == InputValueType::IVType_STRING) {
					if (debug)
						debug_writer->add_string_to_txt(trans_string_to_string_value(in_val.string_value, str_key, str), debug_file);
					else
						trans_string_to_string_value(in_val.string_value, str_key, str);
				}
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		vector<vector<input_value>> trans_matrix_2d_array_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '[') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == ']') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<input_value>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ')' && *(c + 1) == '(') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<input_value> in_val;
					string str_key = key + "[" + to_string(index) + "]";
					in_val = trans_matrix_1d_array_to_input_value(_type[index], str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			vector<input_value> in_val;
			string str_key = key + "[" + to_string(index) + "]";
			in_val = trans_matrix_1d_array_to_input_value(_type[index], str_key, str, debug);
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_2d = " << input_string << ", its 2D size = " << index + 1 << ", defined 2D type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<vector<input_value>> trans_matrix_2d_array_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '[') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == ']') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<input_value>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ')' && *(c + 1) == '(') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<input_value> in_val;
					string str_key = key + "[" + to_string(index) + "]";
					in_val = trans_matrix_1d_const_to_input_value(_type[index], str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			vector<input_value> in_val;
			string str_key = key + "[" + to_string(index) + "]";
			in_val = trans_matrix_1d_const_to_input_value(_type[index], str_key, str, debug);
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_2d = " << input_string << ", its 2D size = " << index + 1 << ", defined 2D type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<vector<input_value>> trans_matrix_2d_const_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '[') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == ']') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<input_value>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ')' && *(c + 1) == '(') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<input_value> in_val;
					string str_key = key + "[" + to_string(index) + "]";
					in_val = trans_matrix_1d_array_to_input_value(_type, str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				vector<input_value> in_val;
				string str_key = key + "[" + to_string(index) + "]";
				in_val = trans_matrix_1d_array_to_input_value(_type, str_key, str, debug);
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		vector<vector<input_value>> trans_matrix_2d_const_const_to_input_value(InputValueType _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '[') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == ']') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<input_value>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ')' && *(c + 1) == '(') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<input_value> in_val;
					string str_key = key + "[" + to_string(index) + "]";
					in_val = trans_matrix_1d_const_to_input_value(_type, str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_2d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				vector<input_value> in_val;
				string str_key = key + "[" + to_string(index) + "]";
				in_val = trans_matrix_1d_const_to_input_value(_type, str_key, str, debug);
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_array_array_array_to_input_value(vector<vector<vector<InputValueType>>> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_array_array_to_input_value(_type[index], str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			vector<vector<input_value>> in_val;
			string str_key = key + "{" + to_string(index) + "}";
			in_val = trans_matrix_2d_array_array_to_input_value(_type[index], str_key, str, debug);
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_3d = " << input_string << ", its 3D size = " << index + 1 << ", defined 3D type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_array_array_const_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_array_const_to_input_value(_type[index], str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			vector<vector<input_value>> in_val;
			string str_key = key + "{" + to_string(index) + "}";
			in_val = trans_matrix_2d_array_const_to_input_value(_type[index], str_key, str, debug);
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_3d = " << input_string << ", its 3D size = " << index + 1 << ", defined 3D type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_array_const_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_const_array_to_input_value(_type[index], str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			vector<vector<input_value>> in_val;
			string str_key = key + "{" + to_string(index) + "}";
			in_val = trans_matrix_2d_const_array_to_input_value(_type[index], str_key, str, debug);
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_3d = " << input_string << ", its 3D size = " << index + 1 << ", defined 3D type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_const_array_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_array_array_to_input_value(_type, str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				vector<vector<input_value>> in_val;
				string str_key = key + "{" + to_string(index) + "}";
				in_val = trans_matrix_2d_array_array_to_input_value(_type, str_key, str, debug);
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_const_const_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_const_array_to_input_value(_type, str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				vector<vector<input_value>> in_val;
				string str_key = key + "{" + to_string(index) + "}";
				in_val = trans_matrix_2d_const_array_to_input_value(_type, str_key, str, debug);
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_const_array_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_array_const_to_input_value(_type, str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				vector<vector<input_value>> in_val;
				string str_key = key + "{" + to_string(index) + "}";
				in_val = trans_matrix_2d_array_const_to_input_value(_type, str_key, str, debug);
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_array_const_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_const_const_to_input_value(_type[index], str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			vector<vector<input_value>> in_val;
			string str_key = key + "{" + to_string(index) + "}";
			in_val = trans_matrix_2d_const_const_to_input_value(_type[index], str_key, str, debug);
			vec_intput_value.push_back(in_val);
			if (index + 1 != _type.size()) {
				cout << "> Input matirx size mismatch ! matrix_3d = " << input_string << ", its 3D size = " << index + 1 << ", defined 3D type size = " << _type.size() << endl;
				SYS_PROGRAM_STOP;
			}
			return vec_intput_value;
		}
		vector<vector<vector<input_value>>> trans_matrix_3d_const_const_const_to_input_value(InputValueType _type, string key, string input_string, bool debug = false) {
			while (input_string.size() != 0)
			{
				if ((*input_string.begin()) == '{') {
					input_string.erase(input_string.begin());
					break;
				}
				else {
					input_string.erase(input_string.begin());
				}
			}
			while (input_string.size() != 0)
			{
				if (*(input_string.end() - 1) == '}') {
					input_string.erase(input_string.end() - 1);
					break;
				}
				else {
					input_string.erase(input_string.end() - 1);
				}
			}
			int index = 0;
			string str;
			vector<vector<vector<input_value>>> vec_intput_value;
			for (auto c = input_string.begin(); c < input_string.end(); c++) {
				if (*c == matrix_separator && *(c - 1) == ']' && *(c + 1) == '[') {
					if (str.size() == 0) {
						cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
						SYS_PROGRAM_STOP;
					}
					vector<vector<input_value>> in_val;
					string str_key = key + "{" + to_string(index) + "}";
					in_val = trans_matrix_2d_const_const_to_input_value(_type, str_key, str, debug);
					vec_intput_value.push_back(in_val);
					str.clear();
					index++;
				}
				else {
					str.push_back(*c);
				}
			}
			if (str.size() == 0 && index > 0) {
				cout << "> Input matirx format error ! matrix_3d = " << input_string << endl;
				SYS_PROGRAM_STOP;
			}
			if (str.size() != 0) {
				vector<vector<input_value>> in_val;
				string str_key = key + "{" + to_string(index) + "}";
				in_val = trans_matrix_2d_const_const_to_input_value(_type, str_key, str, debug);
				vec_intput_value.push_back(in_val);
			}
			return vec_intput_value;
		}
		
		void debug_infile_and_valid_words() {
			stringstream _cout;
			vector<int> valie_lines;
			string valid_words = "<valid>		", invalid_words = "<in-valid>	", note_words = "<note>		";
			_cout << "======================================= D E B U G =======================================" << endl;
			_cout << "LINE	PROPERTY	|CONTENT" << endl;
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			for (int line = 1; line <= input_file.size(); line++) {
				string out = to_string(line) + "	|	", str = input_file[line], equal = "=";
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else
						++s;
				}
				str = "|" + str;
				if (head == '#')
					out += note_words + str;
				else if (head == ' ' || head == '	' || head == '\n' || vec3.size() != 3)
					out += invalid_words + str;
				else if (equal.compare(vec3[1]) != 0)
					out += invalid_words + str;
				else {
					out += valid_words + str;
					valie_lines.push_back(line);
				}
				_cout << out << endl;
			}
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			debug_writer->add_string_to_txt(_cout.str(), debug_file);
			_cout.str("");
			for (int index = 0; index < valie_lines.size(); index++) {
				if (index >= _valid_words.size()) {
					cout << "Input File error : Some valid key words has been defined multi-times !" << endl;
					SYS_PROGRAM_STOP;
				}
				else {
					_cout.str("");
					string out = to_string(valie_lines[index]) + "	|	<valid>		|{\"" + _valid_words[index][0] + "\"}, {\"=\"}, {\"" + _valid_words[index][1] + "\"}";
					_cout << out << endl;
					debug_writer->add_string_to_txt(_cout.str(), debug_file);
				}
			}
			_cout.str("");
			_cout << "=========================================================================================" << endl;
			debug_writer->add_string_to_txt(_cout.str(), debug_file);
		}
		void debug_custom_variavle_and_funcs() {
			stringstream _cout;
			vector<int> valie_lines;
			string valid_words = "<valid>		", invalid_words = "<in-valid>	", note_words = "<note>		";
			_cout << "======================================= D E B U G =======================================" << endl;
			_cout << "NO.		VARIABLE	|VALUE" << endl;
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			for (int line = 0; line < infileMath.infile_vars.size(); line++) {
				string out = to_string(line) + "	|	";
				_cout << out << infileMath.infile_vars[line].key << "			|" << infileMath.infile_vars[line].var << endl;
			}
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			_cout << "NO.		FUNCTIONS	|CONTENT" << endl;
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			for (int line = 0; line < infileMath.infile_funcs.size(); line++) {
				string out = to_string(line) + "	|	" + infileMath.infile_funcs[line].key + "			|";
				stringstream _equation; 
				if (infileMath.check_default_func(infileMath.infile_funcs[line].key)) {
					_equation << "Default function";
				}
				else if (infileMath.check_field_variable(infileMath.infile_funcs[line].key)) {
					_equation << "Field Variables";
				}
				else {
					int op_diff_i = 0;
					for (int index_i = 0; index_i < infileMath.infile_funcs[line].func_structure.size(); index_i++) {
						switch (infileMath.infile_funcs[line].operators.operators_1[index_i + op_diff_i])
						{
						case pf::CO_PLus:
							_equation << "+";
							break;
						case pf::CO_Minux:
							_equation << "-";
							break;
						case pf::CO_Multiply:
							_equation << "*";
							break;
						case pf::CO_Divide:
							_equation << "/";
							break;
						case pf::CO_ParaSeparator:
							_equation << ",";
							op_diff_i++;
							switch (infileMath.infile_funcs[line].operators.operators_1[index_i + op_diff_i])
							{
							case pf::CO_PLus:
								_equation << "+";
								break;
							case pf::CO_Minux:
								_equation << "-";
								break;
							case pf::CO_Multiply:
								_equation << "*";
								break;
							case pf::CO_Divide:
								_equation << "/";
								break;
							}
							break;
						}
						if (infileMath.infile_funcs[line].terms_type.index_1[index_i] < 0) {
							_equation << "{";
						}
						else {
							_equation << infileMath.infile_funcs[infileMath.infile_funcs[line].terms_type.index_1[index_i]].key << "{";
						}
						int op_diff_j = 0;
						for (int index_j = 0; index_j < infileMath.infile_funcs[line].func_structure[index_i].size(); index_j++) {
							switch (infileMath.infile_funcs[line].operators.operators_2[index_i][index_j + op_diff_j])
							{
							case pf::CO_PLus:
								_equation << "+";
								break;
							case pf::CO_Minux:
								_equation << "-";
								break;
							case pf::CO_Multiply:
								_equation << "*";
								break;
							case pf::CO_Divide:
								_equation << "/";
								break;
							case pf::CO_ParaSeparator:
								_equation << ",";
								op_diff_j++;
								switch (infileMath.infile_funcs[line].operators.operators_2[index_i][index_j + op_diff_j])
								{
								case pf::CO_PLus:
									_equation << "+";
									break;
								case pf::CO_Minux:
									_equation << "-";
									break;
								case pf::CO_Multiply:
									_equation << "*";
									break;
								case pf::CO_Divide:
									_equation << "/";
									break;
								}
								break;
							}
							if (infileMath.infile_funcs[line].terms_type.index_2[index_i][index_j] < 0) {
								_equation << "[";
							}
							else {
								_equation << infileMath.infile_funcs[infileMath.infile_funcs[line].terms_type.index_2[index_i][index_j]].key << "[";
							}
							int op_diff_k = 0;
							for (int index_k = 0; index_k < infileMath.infile_funcs[line].func_structure[index_i][index_j].size(); index_k++) {
								switch (infileMath.infile_funcs[line].operators.operators_3[index_i][index_j][index_k + op_diff_k])
								{
								case pf::CO_PLus:
									_equation << "+";
									break;
								case pf::CO_Minux:
									_equation << "-";
									break;
								case pf::CO_Multiply:
									_equation << "*";
									break;
								case pf::CO_Divide:
									_equation << "/";
									break;
								case pf::CO_ParaSeparator:
									_equation << ",";
									op_diff_k++;
									switch (infileMath.infile_funcs[line].operators.operators_3[index_i][index_j][index_k + op_diff_k])
									{
									case pf::CO_PLus:
										_equation << "+";
										break;
									case pf::CO_Minux:
										_equation << "-";
										break;
									case pf::CO_Multiply:
										_equation << "*";
										break;
									case pf::CO_Divide:
										_equation << "/";
										break;
									}
									break;
								}
								if (infileMath.infile_funcs[line].terms_type.index_3[index_i][index_j][index_k] < 0) {
									_equation << "(";
								}
								else {
									_equation << infileMath.infile_funcs[infileMath.infile_funcs[line].terms_type.index_3[index_i][index_j][index_k]].key << "(";
								}
								int op_diff_l = 0;
								for (int index_l = 0; index_l < infileMath.infile_funcs[line].func_structure[index_i][index_j][index_k].size(); index_l++) {
									switch (infileMath.infile_funcs[line].operators.operators_4[index_i][index_j][index_k][index_l + op_diff_l])
									{
									case pf::CO_PLus:
										_equation << "+";
										break;
									case pf::CO_Minux:
										_equation << "-";
										break;
									case pf::CO_Multiply:
										_equation << "*";
										break;
									case pf::CO_Divide:
										_equation << "/";
										break;
									case pf::CO_ParaSeparator:
										_equation << ",";
										op_diff_l++;
										switch (infileMath.infile_funcs[line].operators.operators_4[index_i][index_j][index_k][index_l + op_diff_l])
										{
										case pf::CO_PLus:
											_equation << "+";
											break;
										case pf::CO_Minux:
											_equation << "-";
											break;
										case pf::CO_Multiply:
											_equation << "*";
											break;
										case pf::CO_Divide:
											_equation << "/";
											break;
										}
										break;
									}
									if (infileMath.infile_funcs[line].terms_type.index_4[index_i][index_j][index_k][index_l] < 0) {
										_equation << infileMath.infile_vars[infileMath.infile_funcs[line].func_structure[index_i][index_j][index_k][index_l][0]].key;
									}
									else {
										vector<int> para_ints = infileMath.infile_funcs[line].func_structure[index_i][index_j][index_k][index_l];
										para_ints.erase(para_ints.begin());
										_equation << infileMath.infile_funcs[infileMath.infile_funcs[line].func_structure[index_i][index_j][index_k][index_l][0]].key
											<< "<";
										for (int int_index = 0; int_index < para_ints.size(); int_index++) {
											_equation << to_string(para_ints[int_index]);
											if (int_index < para_ints.size() - 1)
												_equation << ",";
										}
										_equation << ">";
									}
								}
								_equation << ")";
							}
							_equation << "]";
						}
						_equation << "}";
					}
				}
				_cout << out << _equation.str() << endl;
			}
			_cout << "=========================================================================================" << endl;
			debug_writer->add_string_to_txt(_cout.str(), debug_file);
		}
		string_box get_whole_file_strings() {
			return input_file;
		}

		InfileMath infileMath;
		string debug_file;
		WriteToFile* debug_writer;
	private:
		string_box input_file;
		// [lines]{name, value}
		vector<vector<string>> _valid_words;
		char _split;
		int defined_parameters_number;
		char get_first_character_of_line(string line) {
			if (line.compare("") != 0)
				return line.at(0);
			else
				return '0';
		}
		void get_valid_words() {
			for (int line = 1; line <= input_file.size(); line++) {
				string str = input_file[line], equal = "=";
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else {
						if (*(s->end() - 1) == '\r' || *(s->end() - 1) == '\t' || *(s->end() - 1) == '\n' || *(s->end() - 1) == '\0')
							s->erase(s->end() - 1);
						++s;
					}
				}
				if (head == '#' || head == ' ' || head == '\t' || head == '\n' || head == '\r' || head == '\0' || vec3.size() != 3)
					continue;
				else if (equal.compare(vec3[1]) != 0)
					continue;
				else {
					vector<string> valie_line;
					valie_line.push_back(vec3[0]);
					valie_line.push_back(vec3[2]);
					_valid_words.push_back(valie_line);
				}
			}
		}
		void get_define_func_and_variable() {
			for (auto valid_word = _valid_words.begin(); valid_word < _valid_words.end(); valid_word++) {
				if ((*valid_word)[0].compare(infileMath.define_func_key) == 0) {
					defined_parameters_number++;
					string func_key = "", func_equation = ""; bool equation_begin = false;
					for (int index = 0; index < (*valid_word)[1].size(); index++) {
						char c = (*valid_word)[1].at(index);
						if (c == '@' && equation_begin == false) {
							equation_begin = true;
							func_equation = "";
						}
						else if (c == '@' && equation_begin == true) {
							break;
						}
						if (c != '@') {
							if (equation_begin) {
								func_equation.push_back(c);
							}
							else {
								func_key.push_back(c);
							}
						}
					}
					if (func_key.size() == 0 || func_equation.size() == 0) {
						cout << "> error! defined func: " << func_key << " @" << func_equation << "@ " << " cant be translated !" << endl;
						SYS_PROGRAM_STOP;
					}
					else {
						infileMath.add_infile_funcs(func_key, func_equation/*, true*/); // debug input math
					}
				}
				else if ((*valid_word)[0].compare(infileMath.define_variable_key) == 0) {
					defined_parameters_number++;
					string var_key = "", var_value = ""; bool value_begin = false;
					for (int index = 0; index < (*valid_word)[1].size(); index++) {
						char c = (*valid_word)[1].at(index);
						if (c == ',') {
							value_begin = true;
							var_value = "";
						}
						else
						{
							if (value_begin) {
								var_value.push_back(c);
							}
							else {
								var_key.push_back(c);
							}
						}
					}
					if (var_key.size() == 0 || var_value.size() == 0) {
						cout << "> error! defined func: " << var_key << " |" << var_value << "| " << " cant be translated !" << endl;
						SYS_PROGRAM_STOP;
					}
					double value = 0.0;
					if (infile_math_default_funcs::is_string_double(var_value, value)) {
						infileMath.add_infile_var(var_key, value);
					}
					else {
						cout << "> error! defined variable: " << var_key << " (" << value << ") " << " cant be translate to double value !" << endl;
						SYS_PROGRAM_STOP;
					}
				}
			}
		}
		vector<string> split_string(string str, const char split) {
			vector<string> res;
			istringstream iss(str);
			string token;
			while (getline(iss, token, split)) {
				res.push_back(token);
			}
			return res;
		}
		// module
		static InputFileReader* infile;
	};
	InputFileReader* InputFileReader::infile = nullptr;
}