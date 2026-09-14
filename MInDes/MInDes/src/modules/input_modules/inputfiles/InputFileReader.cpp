#include "InputFileReader.h"
namespace pf {
	using namespace std;
	void InputFileReader::init(string input_file_name, bool debug, int file_max_lines, char split) {
		input_file.clear();
		_current_valid_words.clear();
		_current_valid_lines.clear();
		_unread_lines.clear();
		_source_lines.clear();
		_input_file_name = input_file_name;
		_split = split;
		defined_parameters_number = 0;
		std::fstream fin(input_file_name, std::ios::in);
		if (!fin) {
			cout << "> Failed to read the input_file:" << input_file_name << ", use default value 0 (int, REAL, bool), \"\" (string)" << endl;
		}
		else {
			int line_num = 1, physical_line_num = 0, logical_start_line = 0;
			std::string strline{};
			std::string temp_line{};
			bool continuing = false;
			while ((file_max_lines <= 0 || physical_line_num < file_max_lines) && std::getline(fin, strline)) {
				physical_line_num++;
				const size_t comment_pos = strline.find('#');
				if (comment_pos != string::npos)
					strline.erase(comment_pos);
				while (!strline.empty() && (strline.back() == ' ' || strline.back() == '\t' || strline.back() == '\r'))
					strline.pop_back();
				string check_line = strline;
				str_clean(check_line);
				if (check_line.empty())
					continue;
				if (temp_line.empty())
					logical_start_line = physical_line_num;
				const bool has_continuation = strline.size() >= 2 &&
					strline.compare(strline.size() - 2, 2, "\\\\") == 0;
				if (has_continuation)
					strline.erase(strline.size() - 2);
				temp_line += strline;
				continuing = has_continuation;
				if (!continuing) {
					input_file.add_string(line_num++, temp_line);
					_source_lines.push_back(logical_start_line);
					temp_line.clear();
				}
			}
			if (!temp_line.empty()) {
				input_file.add_string(line_num, temp_line);
				_source_lines.push_back(logical_start_line);
				if (continuing)
					mark_unread_line(line_num);
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
		for (const auto& words : _current_valid_words) {
			if (words[0] == infileMath.define_func_key || words[0] == infileMath.define_variable_key)
				continue;
			auto old = std::find_if(_valid_words.begin(), _valid_words.end(), [&](const vector<string>& item) {
				return item[0] == words[0];
			});
			if (old == _valid_words.end())
				_valid_words.push_back(words);
			else
				*old = words;
		}
	}

	bool InputFileReader::read_int_value(string value_name, int& int_value, bool debug) {
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
	bool InputFileReader::read_REAL_value(string value_name, REAL& REAL_value, bool debug) {
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
	bool InputFileReader::read_bool_value(string value_name, bool& bool_value, bool debug) {
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
	bool InputFileReader::read_string_value(string value_name, string& string_value, bool debug) {
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

	string InputFileReader::trans_string_to_int_value(int& val, string str_key, string str_value) {
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
	string InputFileReader::trans_string_to_REAL_value(REAL& val, string str_key, string str_value) {
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
	string InputFileReader::trans_string_to_bool_value(bool& val, string str_key, string str_value) {
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
	string InputFileReader::trans_string_to_string_value(string& val, string str_key, string str_value) {
		stringstream report;
		val = str_value;
		report << "> [-VALID-] " << str_key << " = " << str_value << endl;
		return report.str();
	}

	vector<input_value> InputFileReader::trans_matrix_1d_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug) {
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
	vector<input_value> InputFileReader::trans_matrix_1d_const_to_input_value(InputValueType _type, string key, string input_string, bool debug) {
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
	vector<vector<input_value>> InputFileReader::trans_matrix_2d_array_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug) {
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
	vector<vector<input_value>> InputFileReader::trans_matrix_2d_array_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug) {
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
	vector<vector<input_value>> InputFileReader::trans_matrix_2d_const_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug) {
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
	vector<vector<input_value>> InputFileReader::trans_matrix_2d_const_const_to_input_value(InputValueType _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_array_array_array_to_input_value(vector<vector<vector<InputValueType>>> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_array_array_const_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_array_const_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_const_array_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_const_const_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_const_array_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_array_const_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug) {
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
	vector<vector<vector<input_value>>> InputFileReader::trans_matrix_3d_const_const_const_to_input_value(InputValueType _type, string key, string input_string, bool debug) {
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

	void InputFileReader::debug_infile_and_valid_words() {
		stringstream _cout;
		string valid_words = "<read>      ", invalid_words = "<unread>    ";
		_cout << "======================================= D E B U G =======================================" << endl;
		_cout << "LINE	PROPERTY	|CONTENT" << endl;
		_cout << "-----------------------------------------------------------------------------------------" << endl;
		for (int line = 1; line <= input_file.size(); line++) {
			string out = to_string(line) + "	|	", str = input_file[line], equal = "=";
			std::vector<std::string> buff = {};
			InputLineType type = is_unread_line(line) ? InputLineType::ILType_INACTIVE : get_valid_word_from_string(str, buff);
			if (type == InputLineType::ILType_INACTIVE) {
				out += invalid_words + "|" + str;
			}
			else if (type == InputLineType::ILType_ACTIVE) {
				out += valid_words + "|" + str;
			}

			_cout << out << endl;
		}
		_cout << "-----------------------------------------------------------------------------------------" << endl;
		add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
		_cout.str("");
		_cout << "=========================================================================================" << endl;
		add_string_to_file(_cout.str(), input_output_files_parameters::DebugFile_Path);
	}
	void InputFileReader::debug_custom_variavle_and_funcs() {
		stringstream _cout;
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
	string_box InputFileReader::get_whole_file_strings() {
		return input_file;
	}

	char InputFileReader::get_first_character_of_line(string line) {
		if (line.compare("") != 0)
			return line.at(0);
		else
			return '0';
	}

	string InputFileReader::macro_translator(string macro_str) {
		str_clean(macro_str);
		auto arguments = [&](const string& prefix, size_t count) {
			if (macro_str.size() <= prefix.size() + 2 || macro_str.compare(0, prefix.size(), prefix) != 0 ||
				macro_str[prefix.size()] != '[' || macro_str.back() != ']')
				throw runtime_error("invalid macro syntax");
			vector<string> values = split_string(macro_str.substr(prefix.size() + 1,
				macro_str.size() - prefix.size() - 2), ',', false);
			if (values.size() != count)
				throw runtime_error("invalid macro argument count");
			for (const auto& value : values)
				if (value.empty()) throw runtime_error("empty macro argument");
			return values;
		};
		auto parse_int = [](const string& value) {
			size_t used = 0;
			int result = stoi(value, &used);
			if (used != value.size()) throw runtime_error("invalid integer argument");
			return result;
		};
		auto parse_real = [](const string& value) {
			size_t used = 0;
			REAL result = REAL(stod(value, &used));
			if (used != value.size()) throw runtime_error("invalid real argument");
			return result;
		};

		if (macro_str.rfind("TUBE", 0) == 0) {
			auto values = arguments("TUBE", 3);
			int begin = parse_int(values[0]), end = parse_int(values[1]), step = parse_int(values[2]);
			if (step == 0) throw runtime_error("TUBE step cannot be zero");
			if ((begin < end && step < 0) || (begin > end && step > 0))
				throw runtime_error("TUBE step direction does not reach the end");
			string result;
			for (long long value = begin; (step > 0) ? value <= end : value >= end; value += step) {
				if (!result.empty()) result += ',';
				result += to_string(value);
			}
			return result;
		}
		if (macro_str.rfind("RAND_INT", 0) == 0) {
			auto values = arguments("RAND_INT", 2);
			int minimum = parse_int(values[0]), maximum = parse_int(values[1]);
			if (minimum > maximum) throw runtime_error("RAND_INT minimum is greater than maximum");
			return to_string(uniform_int_distribution<int>(minimum, maximum)(_random_engine));
		}
		if (macro_str.rfind("RAND_REAL", 0) == 0) {
			auto values = arguments("RAND_REAL", 2);
			REAL minimum = parse_real(values[0]), maximum = parse_real(values[1]);
			if (minimum > maximum) throw runtime_error("RAND_REAL minimum is greater than maximum");
			return to_string(uniform_real_distribution<REAL>(minimum, maximum)(_random_engine));
		}
		throw runtime_error("unknown macro");
	}

	void InputFileReader::compile_macros(const char keyword) {
		for (auto line : input_file) {
			if (is_unread_line(line.index))
				continue;
			string& value = input_file[line.index];
			str_clean(value);
			size_t search_position = 0;
			while (true) {
				size_t begin = value.find(keyword, search_position);
				if (begin == string::npos) break;
				size_t end = value.find(keyword, begin + 1);
				if (end == string::npos) {
					const int source_line = size_t(line.index) <= _source_lines.size() ? _source_lines[size_t(line.index - 1)] : line.index;
					string error = "> MACRO ERROR in " + _input_file_name + ", line " + to_string(source_line) +
						": unmatched macro delimiter: " + value + "\n";
					add_string_to_file(error, input_output_files_parameters::DebugFile_Path);
					cout << error;
					SYS_PROGRAM_STOP;
				}
				string macro = value.substr(begin + 1, end - begin - 1);
				try {
					string translated = macro_translator(macro);
					value.replace(begin, end - begin + 1, translated);
					search_position = begin + translated.size();
				}
				catch (const exception& error_detail) {
					const int source_line = size_t(line.index) <= _source_lines.size() ? _source_lines[size_t(line.index - 1)] : line.index;
					string error = "> MACRO ERROR in " + _input_file_name + ", line " + to_string(source_line) +
						", $" + macro + "$: " + error_detail.what() + "\n";
				add_string_to_file(error, input_output_files_parameters::DebugFile_Path);
					cout << error;
					SYS_PROGRAM_STOP;
				}
			}
		}
	}

	void InputFileReader::get_valid_words() {
		for (auto line : input_file) {
			if (is_unread_line(line.index))
				continue;
			std::vector<std::string> valid_word = {};
			InputLineType type = get_valid_word_from_string(line.value, valid_word);
			if (type == InputLineType::ILType_ACTIVE) {
				bool definition_is_valid = true;
				if (valid_word[0] == infileMath.define_variable_key) {
					const string& value = valid_word[1];
					const size_t comma = value.find(',');
					REAL parsed = 0.0;
					definition_is_valid = comma != string::npos && comma > 0 && comma == value.rfind(',') && comma + 1 < value.size();
					if (definition_is_valid) {
						try {
							definition_is_valid = infile_math_default_funcs::is_string_REAL(value.substr(comma + 1), parsed);
						}
						catch (const exception&) {
							definition_is_valid = false;
						}
					}
				}
				else if (valid_word[0] == infileMath.define_func_key) {
					const string& value = valid_word[1];
					const size_t first = value.find('@'), last = value.rfind('@');
					definition_is_valid = first != string::npos && first > 0 && first != last &&
						last == value.size() - 1 && first + 1 < last && value.find('@', first + 1) == last;
				}
				if (definition_is_valid) {
					_current_valid_words.push_back(valid_word);
					_current_valid_lines.push_back(line.index);
				}
				else
					mark_unread_line(line.index);
			}
			else
				mark_unread_line(line.index);
		}
	}

	InputLineType InputFileReader::get_valid_word_from_string(std::string str, std::vector<std::string>& valid_word) {
		str_clean(str);
		if (str.length() == 0) {
			return InputLineType::ILType_INACTIVE;
		}
		std::vector<std::string> buff = split_string(str, '=', false);
		if (buff.size() == 2 && !buff[0].empty() && !buff[1].empty()) {
			valid_word = buff;
			return InputLineType::ILType_ACTIVE;
		}
		else {
			return InputLineType::ILType_INACTIVE;
		}
	}

	void InputFileReader::get_define_func_and_variable() {
		for (auto valid_word = _current_valid_words.begin(); valid_word < _current_valid_words.end(); valid_word++) {
			const int line_index = _current_valid_lines[size_t(valid_word - _current_valid_words.begin())];
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
				if (func_key.empty() || func_equation.empty()) {
					mark_unread_line(line_index);
				}
				else {
					int existing = -1;
					if (infileMath.search_func(func_key, existing) && infileMath.infile_funcs[existing].func_str == func_equation)
						continue;
					InfileMath candidate = infileMath;
					if (!candidate.add_infile_funcs(func_key, func_equation/*, true*/)) {
						mark_unread_line(line_index);
						continue;
					}
					int index = 0;
					if (!candidate.search_func(func_key, index)) {
						mark_unread_line(line_index);
						continue;
					}
					vector<int> para_int; vector<REAL> para_REAL;
					InFileFunc& inFunc = candidate.infile_funcs[index];
					REAL val = inFunc.func(para_int, para_REAL, inFunc.func_structure, inFunc.operators, inFunc.terms_type, candidate.infile_vars, candidate.infile_funcs);
					int cache_index = -1;
					if (candidate.search_var(func_key, cache_index))
						candidate.infile_vars[cache_index].var = val;
					else
						candidate.add_infile_var(func_key, val);
					infileMath = candidate;
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
				if (var_key.empty() || var_value.empty()) {
					mark_unread_line(line_index);
					continue;
				}
				REAL value = 0.0;
				if (infile_math_default_funcs::is_string_REAL(var_value, value)) {
					infileMath.add_infile_var(var_key, value);
				}
				else
					mark_unread_line(line_index);
			}
		}
		for (auto& inFunc : infileMath.infile_funcs) {
			if (infileMath.check_default_func(inFunc.key))
				continue;
			vector<int> para_int;
			vector<REAL> para_REAL;
			REAL value = inFunc.func(para_int, para_REAL, inFunc.func_structure, inFunc.operators,
				inFunc.terms_type, infileMath.infile_vars, infileMath.infile_funcs);
			int cache_index = -1;
			if (infileMath.search_var(inFunc.key, cache_index))
				infileMath.infile_vars[cache_index].var = value;
			else
				infileMath.add_infile_var(inFunc.key, value);
		}
	}

	bool InputFileReader::is_unread_line(int line_index) const {
		return std::find(_unread_lines.begin(), _unread_lines.end(), line_index) != _unread_lines.end();
	}

	void InputFileReader::mark_unread_line(int line_index) {
		if (!is_unread_line(line_index))
			_unread_lines.push_back(line_index);
	}

	std::vector<std::string> InputFileReader::split_string(std::string str, const char keyword, bool preserve_keyword) {
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
}
