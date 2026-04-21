#include "ConProperties.h"
namespace pf {

	ConProperties::ConProperties()
		: _is_init(false),
		_con_number(0) {
	};

	void ConProperties::init() {
		if (!_is_init && main_field::is_con_field_on && main_field::con_number > 0) {
			_con_name.clear();
			_con_number = main_field::con_number;
			string con_name = "CON_0";
			for (int index = 1; index < _con_number; index++)
				con_name += ", CON_" + std::to_string(index);
			WriteDebugFile("# Model.ConProperties = ( con_name_0, ... ) \n");
			string cons_key = "Model.ConProperties", cons_input = "(" + con_name + ")";
			infile_reader::read_string_value(cons_key, cons_input, true);
			std::vector<input_value> cons_value = InputFileReader::get_instance()->
				trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, cons_key, cons_input, true);
			if (cons_value.size() != _con_number) {
				WriteDebugFile("# ERROR: Size of con name should be " + std::to_string(_con_number) + " \n");
				SYS_PROGRAM_STOP;
			}
			_con_name.resize(_con_number);
			for (size_t i = 0; i < _con_number; i++)
				_con_name[i] = cons_value[i].string_value;
			_is_init = true;
		}
	}

	ConProperties& ConProperties::instance() {
		static ConProperties inst;
		return inst;
	};

	size_t ConProperties::con_number() {
		return _con_number;
	}

	size_t ConProperties::con_index(std::string con_name) {
		for (size_t index = 0; index < _con_name.size(); index++)
			if (_con_name[index].compare(con_name) == 0)
				return index;
		return 0;
	}

	bool ConProperties::is_con(std::string con_name) {
		for (size_t index = 0; index < _con_name.size(); index++)
			if (_con_name[index].compare(con_name) == 0)
				return true;
		return false;
	}

	std::string ConProperties::con_name(int con_index) {
		if (con_index < _con_name.size())
			return _con_name[con_index];
		return "DEFAULT_CON";
	}

	std::string ConProperties::con_name(size_t con_index) {
		if (con_index < _con_name.size())
			return _con_name[con_index];
		return "DEFAULT_CON";
	}

}