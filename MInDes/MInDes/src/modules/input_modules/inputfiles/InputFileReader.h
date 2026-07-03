#pragma once
#include <algorithm>
#include <stdexcept>
#include "InputFileMath.h"
#include "../ioFiles_Params.h"
namespace pf {
	using namespace std;
	const char matrix_separator = ',';
	const char macro_symble = '$';
	enum InputValueType { IVType_INT, IVType_REAL, IVType_BOOL, IVType_STRING };
	enum InputLineType { ILType_ACTIVE, ILType_INACTIVE };

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
			_random_engine.seed(std::random_device{}());
		}
		~InputFileReader() {
			_split = ' ';
			infileMath.clear();
			input_file.clear();
			_valid_words.clear();
		}
		void init(string input_file_name, bool debug = false, int file_max_lines = 1000, char split = ' ');

		bool read_int_value(string value_name, int& int_value, bool debug = false);
		bool read_REAL_value(string value_name, REAL& REAL_value, bool debug = false);
		bool read_bool_value(string value_name, bool& bool_value, bool debug = false);
		bool read_string_value(string value_name, string& string_value, bool debug = false);
		
		string trans_string_to_int_value(int& val, string str_key, string str_value);
		string trans_string_to_REAL_value(REAL& val, string str_key, string str_value);
		string trans_string_to_bool_value(bool& val, string str_key, string str_value);
		string trans_string_to_string_value(string& val, string str_key, string str_value);

		vector<input_value> trans_matrix_1d_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false);
		vector<input_value> trans_matrix_1d_const_to_input_value(InputValueType _type, string key, string input_string, bool debug = false);
		vector<vector<input_value>> trans_matrix_2d_array_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false);
		vector<vector<input_value>> trans_matrix_2d_array_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false);
		vector<vector<input_value>> trans_matrix_2d_const_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false);
		vector<vector<input_value>> trans_matrix_2d_const_const_to_input_value(InputValueType _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_array_array_array_to_input_value(vector<vector<vector<InputValueType>>> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_array_array_const_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_array_const_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_const_array_array_to_input_value(vector<vector<InputValueType>> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_const_const_array_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_const_array_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_array_const_const_to_input_value(vector<InputValueType> _type, string key, string input_string, bool debug = false);
		vector<vector<vector<input_value>>> trans_matrix_3d_const_const_const_to_input_value(InputValueType _type, string key, string input_string, bool debug = false);
		
		void debug_infile_and_valid_words();
		void debug_custom_variavle_and_funcs();
		string_box get_whole_file_strings();

		InfileMath infileMath;
	private:
		string_box input_file;
		// [lines]{name, value}
		vector<vector<string>> _valid_words;
		vector<vector<string>> _current_valid_words;
		vector<int> _current_valid_lines;
		vector<int> _unread_lines;
		vector<int> _source_lines;
		char _split;
		int defined_parameters_number;
		string _input_file_name;
		std::mt19937 _random_engine;
		char get_first_character_of_line(string line);

		string macro_translator(string macro_str);
		bool is_unread_line(int line_index) const;
		void mark_unread_line(int line_index);

		void compile_macros(const char keyword = '$');

		void get_valid_words();

		InputLineType get_valid_word_from_string(std::string str, std::vector<std::string>& valid_word);

		void get_define_func_and_variable();

		std::vector<std::string> split_string(std::string str, const char keyword = '=', bool preserve_keyword = true);
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
			int_value = static_cast<long long>(val);
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
