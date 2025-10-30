#pragma once
#include "InputFileMath.h"
#include "../ioFiles_Params.h"
namespace pf {
	using namespace std;
	const char matrix_separator = ',';
	const char macro_symble = '$';
	enum InputValueType { IVType_INT, IVType_REAL, IVType_BOOL, IVType_STRING };
	enum InputLineType { ILType_ACTIVE, ILType_INACTIVE, ILType_NOTE };

	inline void str_clean(std::string& str) {
		//windows & linux
		str_char_delete(str, ' ');
		str_char_delete(str, '\t');
		str_char_delete(str, '\r');
		str_char_delete(str, '\n');
		str_char_delete(str, '\0');
	}

	inline void str_string_delete(std::string& str, std::string removed_word) {
		str_char_delete(str, '\r');
		str_char_delete(str, '\n');
		str_char_delete(str, '\0');
		std::size_t start_position{ str.find(removed_word) };
		while (start_position != std::string::npos) {
			str.erase(start_position, removed_word.length());
			start_position = str.find(removed_word);
		}
	}

	struct input_value {
		int int_value;
		REAL REAL_value;
		bool bool_value;
		string string_value;
		input_value() {
			int_value = 0;
			REAL_value = 0.0;
			bool_value = false;
			string_value = "";
		}
		~input_value() {};
	};
	struct string_elem {
		int index;
		std::string value;
		string_elem() {
			index = 0;
			value = "";
		}
		string_elem& operator=(const string_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	class string_box {
	public:
		string_box() {
			_string_box.reserve(0);
		}
		~string_box() {
			_string_box.clear();
		}
		std::vector<string_elem> _string_box;
		typedef std::vector<string_elem>::iterator iterator;
		typedef std::vector<string_elem>::const_iterator citerator;
		iterator  begin() { return _string_box.begin(); };
		iterator  end() { return _string_box.end(); };
		std::string& operator[](const int index) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			std::cout << "string_box error, can't find the value index : " << index << std::endl;
			SYS_PROGRAM_STOP;
		}
		string_box& operator=(const string_box& n) {
			_string_box = n._string_box;
			return *this;
		}
		void add_string(int _index, std::string _value) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			_string_box.reserve(_string_box.size() + 1);
			string_elem elem;
			elem.index = _index;
			elem.value = _value;
			_string_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _string_box.begin(); i < _string_box.end();) {
				if (i->index == index) {
					i = _string_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_string_box.clear();
		}
		int size() const {
			return int(_string_box.size());
		}
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
			defined_parameters_number = 0;
			_split = ' ';
			input_file.clear();
		}
		~InputFileReader() {
			_split = ' ';
			infileMath.clear();
			input_file.clear();
			_valid_words.clear();
		}
		void init(string input_file_name, bool debug = false, int file_max_lines = 1000, char split = ' ') {
			input_file.clear();
			_split = split;
			defined_parameters_number = 0;
			std::fstream fin(input_file_name, std::ios::in);
			if (!fin) {
				cout << "> Failed to read the input_file:" << input_file_name << ", use default value 0 (int, REAL, bool), \"\" (string)" << endl;
			}
			else {
				int line_num = 1;
				std::string strline{};
				std::string temp_line{};
				while (std::getline(fin, strline)) {
					while (strline.find("\\\\") != std::string::npos) {
						str_string_delete(strline, "\\\\");
						temp_line += strline;
						std::getline(fin, strline);
					}
					temp_line += strline; //every lines merge into temp_line
					input_file.add_string(line_num, temp_line);
					line_num++;
					temp_line.clear();
				}
			}
			fin.close();
			// for input_file
			compile_macros();
			// for _valid_words
			get_valid_words();
			// for infileMath
			get_define_func_and_variable();
			/*if (infileMath.check_variables_funcs_keys() == false) {
				cout << "> Defined key of variables and functions cannot be the same." << endl;
				SYS_PROGRAM_STOP;
			}*/
			if (debug) {
				debug_infile_and_valid_words();
				debug_custom_variavle_and_funcs();
			}
			// Remove duplicate definitions, last defined key is active.
			vector<vector<string>> buff_words = _valid_words;
			_valid_words.clear();
			for (auto vec1 = buff_words.begin(); vec1 < buff_words.end(); vec1++) {
				vector<string> words = { (*vec1)[0] ,(*vec1)[1] };
				if (words[0].compare(infileMath.define_func_key) == 0 || words[0].compare(infileMath.define_variable_key) == 0) 
					_valid_words.push_back(words);
				bool is_last_defined = true;
				for (auto vec2 = vec1 + 1; vec2 < buff_words.end(); vec2++)
					if (words[0].compare((*vec2)[0]) == 0)
						is_last_defined = false;
				if(is_last_defined)
					_valid_words.push_back(words);
			}
		}

		bool read_int_value(string value_name, int& int_value, bool debug = false) {
			stringstream report;
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					if (debug)
						add_string_to_file(trans_string_to_int_value(int_value, (*vec)[0], (*vec)[1]), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_int_value(int_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (debug)
				add_string_to_file(string("> [DEFAULT] " + value_name + " = " + to_string(int_value) + "\n"), input_output_files_parameters::DebugFile_Path);
			return false;
		}
		bool read_REAL_value(string value_name, REAL& REAL_value, bool debug = false) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					if (debug)
						add_string_to_file(trans_string_to_REAL_value(REAL_value, (*vec)[0], (*vec)[1]), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_REAL_value(REAL_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (debug)
				add_string_to_file(string("> [DEFAULT] " + value_name + " = " + to_string(REAL_value) + "\n"), input_output_files_parameters::DebugFile_Path);
			return false;
		}
		bool read_bool_value(string value_name, bool& bool_value, bool debug = false) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					bool_value = false;
					if (debug)
						add_string_to_file(trans_string_to_bool_value(bool_value, (*vec)[0], (*vec)[1]), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_bool_value(bool_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (bool_value) {
				if (debug)
					add_string_to_file(string("> [DEFAULT] " + value_name + " = TRUE" + "\n"), input_output_files_parameters::DebugFile_Path);
			}
			else {
				if (debug)
					add_string_to_file(string("> [DEFAULT] " + value_name + " = FALSE" + "\n"), input_output_files_parameters::DebugFile_Path);
			}
			return false;
		}
		bool read_string_value(string value_name, string& string_value, bool debug = false) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					if (debug)
						add_string_to_file(trans_string_to_string_value(string_value, (*vec)[0], (*vec)[1]), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_string_value(string_value, (*vec)[0], (*vec)[1]);
					return true;
				}
			if (debug)
				add_string_to_file(string("> [DEFAULT] " + value_name + " = " + string_value + "\n"), input_output_files_parameters::DebugFile_Path);
			return false;
		}
		
		string trans_string_to_int_value(int& val, string str_key, string str_value) {
			stringstream report;
			int index = -1;
			if (infile_math_default_funcs::is_string_int(str_value, val)) {
				;
			}
			else if (infileMath.search_var(str_value, index)) {
				val = int(infileMath.infile_vars[index].var + 0.5);
			}
			else if (infileMath.search_func(str_value, index)) {
				vector<int> para_int; vector<REAL> para_REAL;
				InFileFunc& inFunc = infileMath.infile_funcs[index];
				val = int(inFunc.func(para_int, para_REAL, inFunc.func_structure, inFunc.operators, inFunc.terms_type, infileMath.infile_vars, infileMath.infile_funcs) + 0.5);
			}
			else {
				cout << "> Input file error! the value of Keyword: " << str_key << " cannot be translate to int value." << endl;
				SYS_PROGRAM_STOP;
			}
			report << "> [-VALID-] " << str_key << " = " << val << endl;
			return report.str();
		}
		string trans_string_to_REAL_value(REAL& val, string str_key, string str_value) {
			stringstream report;
			int index = -1;
			if (infile_math_default_funcs::is_string_REAL(str_value, val)) {
				;
			}
			else if (infileMath.search_var(str_value, index)) {
				val = infileMath.infile_vars[index].var;
			}
			else if (infileMath.search_func(str_value, index)) {
				vector<int> para_int; vector<REAL> para_REAL;
				InFileFunc& inFunc = infileMath.infile_funcs[index];
				val = inFunc.func(para_int, para_REAL, inFunc.func_structure, inFunc.operators, inFunc.terms_type, infileMath.infile_vars, infileMath.infile_funcs);
			}
			else {
				cout << "> Input file error! the value of Keyword: " << str_key << " cannot be translate to REAL value." << endl;
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
							add_string_to_file(trans_string_to_int_value(in_val.int_value, str_key, str), input_output_files_parameters::DebugFile_Path);
						else
							trans_string_to_int_value(in_val.int_value, str_key, str);
					}
					else if (_type[index] == InputValueType::IVType_REAL) {
						if (debug)
							add_string_to_file(trans_string_to_REAL_value(in_val.REAL_value, str_key, str), input_output_files_parameters::DebugFile_Path);
						else
							trans_string_to_REAL_value(in_val.REAL_value, str_key, str);
					}
					else if (_type[index] == InputValueType::IVType_BOOL) {
						if (debug)
							add_string_to_file(trans_string_to_bool_value(in_val.bool_value, str_key, str), input_output_files_parameters::DebugFile_Path);
						else
							trans_string_to_bool_value(in_val.bool_value, str_key, str);
					}
					else if (_type[index] == InputValueType::IVType_STRING) {
						if (debug)
							add_string_to_file(trans_string_to_string_value(in_val.string_value, str_key, str), input_output_files_parameters::DebugFile_Path);
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
					add_string_to_file(trans_string_to_int_value(in_val.int_value, str_key, str), input_output_files_parameters::DebugFile_Path);
				else
					trans_string_to_int_value(in_val.int_value, str_key, str);
			}
			else if (_type[index] == InputValueType::IVType_REAL) {
				if (debug)
					add_string_to_file(trans_string_to_REAL_value(in_val.REAL_value, str_key, str), input_output_files_parameters::DebugFile_Path);
				else
					trans_string_to_REAL_value(in_val.REAL_value, str_key, str);
			}
			else if (_type[index] == InputValueType::IVType_BOOL) {
				if (debug)
					add_string_to_file(trans_string_to_bool_value(in_val.bool_value, str_key, str), input_output_files_parameters::DebugFile_Path);
				else
					trans_string_to_bool_value(in_val.bool_value, str_key, str);
			}
			else if (_type[index] == InputValueType::IVType_STRING) {
				if (debug)
					add_string_to_file(trans_string_to_string_value(in_val.string_value, str_key, str), input_output_files_parameters::DebugFile_Path);
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
							add_string_to_file(trans_string_to_int_value(in_val.int_value, str_key, str), input_output_files_parameters::DebugFile_Path);
						else
							trans_string_to_int_value(in_val.int_value, str_key, str);
					}
					else if (_type == InputValueType::IVType_REAL) {
						if (debug)
							add_string_to_file(trans_string_to_REAL_value(in_val.REAL_value, str_key, str), input_output_files_parameters::DebugFile_Path);
						else
							trans_string_to_REAL_value(in_val.REAL_value, str_key, str);
					}
					else if (_type == InputValueType::IVType_BOOL) {
						if (debug)
							add_string_to_file(trans_string_to_bool_value(in_val.bool_value, str_key, str), input_output_files_parameters::DebugFile_Path);
						else
							trans_string_to_bool_value(in_val.bool_value, str_key, str);
					}
					else if (_type == InputValueType::IVType_STRING) {
						if (debug)
							add_string_to_file(trans_string_to_string_value(in_val.string_value, str_key, str), input_output_files_parameters::DebugFile_Path);
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
						add_string_to_file(trans_string_to_int_value(in_val.int_value, str_key, str), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_int_value(in_val.int_value, str_key, str);
				}
				else if (_type == InputValueType::IVType_REAL) {
					if (debug)
						add_string_to_file(trans_string_to_REAL_value(in_val.REAL_value, str_key, str), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_REAL_value(in_val.REAL_value, str_key, str);
				}
				else if (_type == InputValueType::IVType_BOOL) {
					if (debug)
						add_string_to_file(trans_string_to_bool_value(in_val.bool_value, str_key, str), input_output_files_parameters::DebugFile_Path);
					else
						trans_string_to_bool_value(in_val.bool_value, str_key, str);
				}
				else if (_type == InputValueType::IVType_STRING) {
					if (debug)
						add_string_to_file(trans_string_to_string_value(in_val.string_value, str_key, str), input_output_files_parameters::DebugFile_Path);
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
			vector<int> valid_lines_number;
			vector<vector<string>> valid_lines;
			string valid_words = "<read>      ", invalid_words = "<unread>    ", note_words = "<note>      ";
			_cout << "======================================= D E B U G =======================================" << endl;
			_cout << "LINE	PROPERTY	|CONTENT" << endl;
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			for (int line = 1; line <= input_file.size(); line++) {
				string out = to_string(line) + "	|	", str = input_file[line], equal = "=";
				std::vector<std::string> buff = {};
				InputLineType type = get_valid_word_from_string(str, buff);
				if (type == InputLineType::ILType_NOTE) {
					out += note_words + "|" + str;
				}
				else if (type == InputLineType::ILType_INACTIVE) {
					out += invalid_words + "|" + str;
				}
				else if (type == InputLineType::ILType_ACTIVE) {
					out += valid_words + "|" + str;
					valid_lines_number.push_back(line);
					valid_lines.push_back(buff);
				}

				_cout << out << endl;
			}
			_cout << "-----------------------------------------------------------------------------------------" << endl;
			add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
			_cout.str("");
			for (int index = 0; index < valid_lines.size(); index++) {
				if (valid_lines[index][0].compare(infileMath.define_func_key) == 0 || valid_lines[index][0].compare(infileMath.define_variable_key) == 0) {
					_cout.str("");
					string out = to_string(valid_lines_number[index]) + "	|   <Custom>    |{\"" + valid_lines[index][0] + "\"}, {\"=\"}, {\"" + valid_lines[index][1] + "\"}";
					_cout << out << endl;
					add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
					continue;
				}
				bool is_multi_defined_words = false;
				for (int index2 = index + 1; index2 < valid_lines.size(); index2++)
					if (valid_lines[index][0].compare(valid_lines[index2][0]) == 0)
						is_multi_defined_words = true;
				if (is_multi_defined_words) {
					_cout.str("");
					string out = to_string(valid_lines_number[index]) + "	|   <in-valid>  |{\"" + valid_lines[index][0] + "\"}, {\"=\"}, {\"" + valid_lines[index][1] + "\"}";
					_cout << out << endl;
					add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
				}
				else {
					_cout.str("");
					string out = to_string(valid_lines_number[index]) + "	|   <valid>     |{\"" + valid_lines[index][0] + "\"}, {\"=\"}, {\"" + valid_lines[index][1] + "\"}";
					_cout << out << endl;
					add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
				}
			}
			_cout.str("");
			_cout << "=========================================================================================" << endl;
			add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
		}
		void debug_custom_variavle_and_funcs() {
			stringstream _cout;
			vector<int> valid_lines;
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
			add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
		}
		string_box get_whole_file_strings() {
			return input_file;
		}

		InfileMath infileMath;
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

		string macro_translator(string macro_str) {
			// read_int_value()
			// read_REAL_value()
			string trans_words = "";
			str_clean(macro_str);
			size_t str_size = macro_str.size();
			// - random
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// MACRO: TUBE
			if (macro_str.at(0) == 'T' && macro_str.at(1) == 'U' && macro_str.at(2) == 'B' && macro_str.at(3) == 'E') {
				macro_str.erase(macro_str.begin()); // T
				macro_str.erase(macro_str.begin()); // U
				macro_str.erase(macro_str.begin()); // B
				macro_str.erase(macro_str.begin()); // E
				str_size = macro_str.size();
				if (macro_str.at(0) == '[' && macro_str.at(str_size - 1) == ']') {
					macro_str.erase(macro_str.begin()); // [
					macro_str.erase(macro_str.end() - 1); // ]
					std::vector<std::string> buff = split_string(macro_str, ',', false);
					// MACRO: normal tube
					if (buff.size() == 3) {
						int index_begin = stoi(buff[0]), index_end = stoi(buff[1]), index_jump = stoi(buff[2]);
						for (int index = index_begin; index <= index_end; index = index + index_jump) {
							trans_words += to_string(index);
							if (index + index_jump <= index_end)
								trans_words += ",";
						}
					}
				}
				// others
				// -
			}
			// MACRO: RAND
			else if (macro_str.at(0) == 'R' && macro_str.at(1) == 'A' && macro_str.at(2) == 'N' && macro_str.at(3) == 'D') {
				macro_str.erase(macro_str.begin()); // R
				macro_str.erase(macro_str.begin()); // A
				macro_str.erase(macro_str.begin()); // N
				macro_str.erase(macro_str.begin()); // D
				str_size = macro_str.size();
				// MACRO: RAND_INT
				if (macro_str.at(0) == '_' && macro_str.at(1) == 'I' && macro_str.at(2) == 'N' && macro_str.at(3) == 'T') {
					macro_str.erase(macro_str.begin()); // _
					macro_str.erase(macro_str.begin()); // I
					macro_str.erase(macro_str.begin()); // N
					macro_str.erase(macro_str.begin()); // T
					str_size = macro_str.size();
					if (macro_str.at(0) == '[' && macro_str.at(str_size - 1) == ']') {
						macro_str.erase(macro_str.begin()); // [
						macro_str.erase(macro_str.end() - 1); // ]
						std::vector<std::string> buff = split_string(macro_str, ',', false);
						if (buff.size() == 2) {
							int val_min = stoi(buff[0]), val_max = stoi(buff[1]);
							REAL rand = REAL(real_dist(gen));  // random
							REAL result = REAL(val_max - val_min) * rand + REAL(val_min);
							if (result > 0.0)
								trans_words = to_string(int(result + 0.5));
							else
								trans_words = to_string(int(result - 0.5));
						}
					}
				}
				// MACRO: RAND_DOUBLE
				else if (macro_str.at(0) == '_' && macro_str.at(1) == 'D' && macro_str.at(2) == 'O' && macro_str.at(3) == 'U' && macro_str.at(4) == 'B' && macro_str.at(5) == 'L' && macro_str.at(6) == 'E') {
					macro_str.erase(macro_str.begin()); // _
					macro_str.erase(macro_str.begin()); // D
					macro_str.erase(macro_str.begin()); // O
					macro_str.erase(macro_str.begin()); // U
					macro_str.erase(macro_str.begin()); // B
					macro_str.erase(macro_str.begin()); // L
					macro_str.erase(macro_str.begin()); // E
					str_size = macro_str.size();
					if (macro_str.at(0) == '[' && macro_str.at(str_size - 1) == ']') {
						macro_str.erase(macro_str.begin()); // [
						macro_str.erase(macro_str.end() - 1); // ]
						std::vector<std::string> buff = split_string(macro_str, ',', false);
						// MACRO: normal rand REAL
						if (buff.size() == 2) {
							REAL val_min = REAL(stod(buff[0])), val_max = REAL(stod(buff[1]));
							//read_REAL_value(buff[0], val_min, false);
							//read_REAL_value(buff[1], val_max, false);
							REAL rand = REAL(real_dist(gen));  // random
							REAL result = (val_max - val_min) * rand + val_min;
							trans_words = to_string(result);
						}
					}
				}
				// others
				// -
			}
			// others
			// -
			return trans_words;
		}

		void compile_macros(const char keyword = '$') {
			for (auto line : input_file) {
				std::size_t macro_position{ line.value.find(keyword) };
				std::size_t macro_size{ line.value.substr(macro_position + 1).find(keyword)};
				while (macro_position != std::string::npos && macro_size != std::string::npos) {
					string macro = input_file[line.index].substr(macro_position + 1, macro_size);
					input_file[line.index].replace(macro_position, macro_size + 2, macro_translator(macro));
					macro_position = input_file[line.index].find(keyword);
					macro_size = input_file[line.index].substr(macro_position + 1).find(keyword);
				}
				if (macro_position != std::string::npos && macro_size == std::string::npos) {
					string error = "> MACROS COMPILOR ERROR, line number: " + to_string(line.index) + ", only one MACRO key in this line: \n" + line.value + "\n";
					add_string_to_file(error, input_output_files_parameters::DebugFile_Path);
					exit(0);
				}
			}
		}

		void get_valid_words() {
			for (auto line : input_file) {
				std::vector<std::string> valid_word = {};
				InputLineType type = get_valid_word_from_string(line.value, valid_word);
				if (type == InputLineType::ILType_ACTIVE)
					_valid_words.push_back(valid_word);
			}
		}

		InputLineType get_valid_word_from_string(std::string str, std::vector<std::string>& valid_word) {
			str_clean(str);
			if (str.length() == 0) {
				return InputLineType::ILType_INACTIVE;
			}
			if (str.at(0) == '#') {
				return InputLineType::ILType_NOTE;
			}
			std::vector<std::string> buff = split_string(str, '=', false);
			if (buff.size() == 2) {
				valid_word = buff;
				return InputLineType::ILType_ACTIVE;
			}
			else {
				return InputLineType::ILType_INACTIVE;
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
						int index = 0;
						infileMath.search_func(func_key, index);
						vector<int> para_int; vector<REAL> para_REAL;
						InFileFunc& inFunc = infileMath.infile_funcs[index];
						REAL val = inFunc.func(para_int, para_REAL, inFunc.func_structure, inFunc.operators, inFunc.terms_type, infileMath.infile_vars, infileMath.infile_funcs);
						infileMath.add_infile_var(func_key, val);
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
					REAL value = 0.0;
					if (infile_math_default_funcs::is_string_REAL(var_value, value)) {
						infileMath.add_infile_var(var_key, value);
					}
					else {
						cout << "> error! defined variable: " << var_key << " (" << var_value << ") " << " cant be translate to REAL value !" << endl;
						SYS_PROGRAM_STOP;
					}
				}
			}
		}

		std::vector<std::string> split_string(std::string str, const char keyword = '=', bool preserve_keyword = true) {
			std::size_t keyword_pos = str.find(keyword);
			if (keyword_pos == std::string::npos)
				return { str };
			std::vector<std::string> vec = { str.substr(0, keyword_pos) };
			std::size_t keyword_pos2 = str.substr(keyword_pos + 1).find(keyword);
			// have keyword
			while (keyword_pos2 != std::string::npos) {
				if (preserve_keyword)
					vec.push_back("=");
				vec.push_back(str.substr(keyword_pos + 1, keyword_pos2));
				keyword_pos = keyword_pos2 + keyword_pos + 1;
				keyword_pos2 = str.substr(keyword_pos + 1).find(keyword);
			}
			if (preserve_keyword)
				vec.push_back("=");
			vec.push_back(str.substr(keyword_pos + 1, keyword_pos2 - keyword_pos - 1));
			return vec;
		}
		// module
		inline static InputFileReader* infile = nullptr;
	};
	namespace infile_reader {
		inline bool read_int_value(string value_name, int& int_value, bool debug = false) {
			return InputFileReader::get_instance()->read_int_value(value_name, int_value, debug);
		}
		inline bool read_int_value(string value_name, size_t& int_value, bool debug = false) {
			int val = int(int_value);
			bool result = InputFileReader::get_instance()->read_int_value(value_name, val, debug);
			int_value = size_t(val);
			return result;
		}
		inline bool read_int_value(string value_name, long long& int_value, bool debug = false) {
			int val = int(int_value);
			bool result = InputFileReader::get_instance()->read_int_value(value_name, val, debug);
			int_value = long long(val);
			return result;
		}
		inline bool read_real_value(string value_name, double& REAL_value, bool debug = false) {
			REAL val = REAL(REAL_value);
			bool result = InputFileReader::get_instance()->read_REAL_value(value_name, val, debug);
			REAL_value = double(val);
			return result;
		}
		inline bool read_real_value(string value_name, float& REAL_value, bool debug = false) {
			REAL val = REAL(REAL_value);
			bool result = InputFileReader::get_instance()->read_REAL_value(value_name, val, debug);
			REAL_value = float(val);
			return result;
		}
		inline bool read_bool_value(string value_name, bool& bool_value, bool debug = false) {
			return InputFileReader::get_instance()->read_bool_value(value_name, bool_value, debug);
		}
		inline bool read_string_value(string value_name, string& string_value, bool debug = false) {
			return InputFileReader::get_instance()->read_string_value(value_name, string_value, debug);
		}
	}
}