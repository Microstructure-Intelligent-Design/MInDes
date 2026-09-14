#include "ConRegions.h"
namespace pf {

	ConRegions::ConRegions()
		: _is_init(false),
		_region_number(0) {
	};

	void ConRegions::init() {
		if (!_is_init && main_field::is_phi_field_on && main_field::is_con_field_on && main_field::con_number > 0) {
			PhiProperties::instance().init();
			_con_name.clear();
			string con_name = "CON_0";
			for (int index = 1; index < main_field::con_number; index++)
				con_name += ", CON_" + std::to_string(index);
			WriteDebugFile("# Model.ConNames = ( con_name_0, ... ) \n");
			string cons_key = "Model.ConNames", cons_input = "(" + con_name + ")";
			infile_reader::read_string_value(cons_key, cons_input, true);
			std::vector<input_value> cons_value = InputFileReader::get_instance()->
				trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, cons_key, cons_input, true);
			if (cons_value.size() != main_field::con_number) {
				WriteDebugFile("# ERROR: Size of con name should be " + std::to_string(main_field::con_number) + " \n");
				SYS_PROGRAM_STOP;
			}
			_con_name.resize(main_field::con_number);
			for (size_t i = 0; i < main_field::con_number; i++) {
				if (is_con(cons_value[i].string_value)) {
					WriteDebugFile("# ERROR: Same con name " + cons_value[i].string_value + " has been defined \n");
					SYS_PROGRAM_STOP;
				}
				else {
					_con_name[i] = cons_value[i].string_value;
				}
			}
			WriteDebugFile("# Model.ConRegions = {[(region_name_0),(phi_name_0, phi_name_1, ... ),(con_name_0, con_name_1, ... )],  ... } \n");
			WriteDebugFile("#       Note : 1. phi con default in region_name_0, and 2. same phi_name con_name cannot in different region \n");
			string regions_key = "Model.ConRegions", regions_input = "{[(Region_0),(),()]}";
			infile_reader::read_string_value(regions_key, regions_input, true);
			std::vector<std::vector<std::vector<input_value>>> regions_value = InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value
			(InputValueType::IVType_STRING, regions_key, regions_input, true);
			_con_region.resize(main_field::con_number, 0);
			_phi_region.resize(main_field::phi_number, 0);
			_phi_property_region.resize(PhiProperties::instance().phi_property_number(), 0);
			for (size_t i = 0; i < regions_value.size(); i++) {
				std::string region_name = regions_value[i][0][0].string_value;
				size_t region_index = add_region_name(region_name);
				for (size_t j = 0; j < regions_value[i][1].size(); j++) {
					std::string phi_name = regions_value[i][1][j].string_value;
					if (!PhiProperties::instance().is_phi_property(phi_name)) {
						WriteDebugFile("# ERROR: phi name " + phi_name + " set in region " + region_name + " has not been defined \n");
						SYS_PROGRAM_STOP;
					}
					_phi_property_region[PhiProperties::instance().phi_property(phi_name)] = region_index;
					for (size_t k = 0; k < PhiProperties::instance().property_phi(phi_name).size(); k++)
						_phi_region[PhiProperties::instance().property_phi(phi_name)[k]] = region_index;
				}
				for (size_t j = 0; j < regions_value[i][2].size(); j++) {
					std::string con_name = regions_value[i][2][j].string_value;
					if (!is_con(con_name)) {
						WriteDebugFile("# ERROR: con name " + con_name + " set in region " + region_name + " has not been defined \n");
						SYS_PROGRAM_STOP;
					}
					_con_region[con_index(con_name)] = region_index;
				}
			}
			_region_con.resize(_region_number, main_field::con_number);
			_region_phi.resize(_region_number, main_field::phi_number);
			_region_phi_property.resize(_region_number, PhiProperties::instance().phi_property_number());
			_region_con_number.resize(_region_number, 0);
			_region_phi_number.resize(_region_number, 0);
			_region_phi_property_number.resize(_region_number, 0);
			for (size_t i = 0; i < _region_number; i++) {
				for (size_t j = 0; j < main_field::con_number; j++)
					if (_con_region[j] == i) {
						_region_con(i, _region_con_number[i]) = j;
						_region_con_number[i]++;
					}
				for (size_t j = 0; j < main_field::phi_number; j++)
					if (_phi_region[j] == i) {
						_region_phi(i, _region_phi_number[i]) = j;
						_region_phi_number[i]++;
					}
				for (size_t j = 0; j < PhiProperties::instance().phi_property_number(); j++)
					if (_phi_property_region[j] == i) {
						_region_phi_property(i, _region_phi_property_number[i]) = j;
						_region_phi_property_number[i]++;
					}
			}
			_is_init = true;
		}
	}

	ConRegions& ConRegions::instance() {
		static ConRegions inst;
		return inst;
	};

	size_t ConRegions::con_index(std::string con_name) {
		for (size_t index = 0; index < _con_name.size(); index++)
			if (_con_name[index].compare(con_name) == 0)
				return index;
		return 0;
	}

	size_t ConRegions::region_index(std::string region_name) {
		for (size_t index = 0; index < _region_name.size(); index++)
			if (_region_name[index].compare(region_name) == 0)
				return index;
		return 0;
	}

	std::string ConRegions::con_name(size_t con_index) {
		if (con_index < _con_name.size())
			return _con_name[con_index];
		return "DEFAULT_CON";
	}

	std::string ConRegions::region_name(size_t region_index) {
		if (region_index < _region_name.size())
			return _region_name[region_index];
		return "DEFAULT_CON";
	}

	bool ConRegions::is_con(std::string con_name) {
		for (size_t index = 0; index < _con_name.size(); index++)
			if (_con_name[index].compare(con_name) == 0)
				return true;
		return false;
	}

	bool ConRegions::is_region(std::string region_name) {
		for (size_t index = 0; index < _region_name.size(); index++)
			if (_region_name[index].compare(region_name) == 0)
				return true;
		return false;
	}

	bool ConRegions::is_con_in_region(size_t region_index, size_t con_index) {
		for (size_t index = 0; index < _region_con_number[region_index]; index++)
			if (_region_con(region_index, index) == con_index)
				return true;
		return false;
	}

	bool ConRegions::is_phi_in_region(size_t region_index, size_t phi_index) {
		for (size_t index = 0; index < _region_phi_number[region_index]; index++)
			if (_region_phi(region_index, index) == phi_index)
				return true;
		return false;
	}

	size_t ConRegions::add_region_name(std::string region_name) {
		for (size_t index = 0; index < _region_name.size(); index++)
			if (_region_name[index].compare(region_name) == 0)
				return index;
		_region_name.push_back(region_name);
		_region_number = _region_name.size();
		return _region_number - 1;
	}
}