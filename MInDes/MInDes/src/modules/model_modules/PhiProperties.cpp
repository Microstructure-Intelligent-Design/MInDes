#include "PhiProperties.h"
namespace pf {

	PhiProperties::PhiProperties()
		: _is_init(false),
		_phi_property_number(0){
	};

	void PhiProperties::init() {
		if (!_is_init && main_field::is_phi_field_on && main_field::phi_number > 0) {
			_phi_property.clear();
			_phi_property.resize(main_field::phi_number, 0);
			WriteDebugFile("# Model.PhiProperties = {[(phi_name_0),(phi_index_0, phi_index_1, ... )],  ... } \n");
			string grains_key = "Model.PhiProperties", grains_input = "{[(PHI_0),()]}";
			infile_reader::read_string_value(grains_key, grains_input, true);
			std::vector<std::vector<std::vector<input_value>>> grains_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value
			({ InputValueType::IVType_STRING , InputValueType::IVType_INT}, grains_key, grains_input, true);
			for (size_t i = 0; i < grains_value.size(); i++)
				add_property_name(grains_value[i][0][0].string_value);
			for (size_t i = 0; i < grains_value.size(); i++) {
				for (size_t j = 0; j < grains_value[i][1].size(); j++) {
					int phi_index = grains_value[i][1][j].int_value;
					if (phi_index < main_field::phi_number)
						_phi_property[phi_index] = phi_property(grains_value[i][0][0].string_value);
				}
			}
			_property_phi.resize(_phi_property_number);
			for (size_t index = 0; index < main_field::phi_number; index++)
				_property_phi[_phi_property[index]].push_back(index);
			_is_init = true;
		}
	}

	PhiProperties& PhiProperties::instance() {
		static PhiProperties inst;
		return inst;
	};

	size_t PhiProperties::phi_property(size_t phi_index) {
		if (phi_index < _phi_property.size())
			return _phi_property[phi_index];
		return 0;
	}

	std::vector<size_t> PhiProperties::property_phi(size_t phi_property_index) {
		return _property_phi[phi_property_index];
	}

	std::vector<size_t> PhiProperties::property_phi(std::string _property_name) {
		return _property_phi[phi_property(_property_name)];
	}

	size_t PhiProperties::phi_property_number() {
		return _phi_property_number;
	}

	size_t PhiProperties::phi_property(std::string _property_name) {
		for (size_t index = 0; index < _phi_property_name.size(); index++)
			if (_phi_property_name[index].compare(_property_name) == 0)
				return index;
		return 0;
	}

	bool PhiProperties::is_phi_property(std::string _property_name) {
		for (size_t index = 0; index < _phi_property_name.size(); index++)
			if (_phi_property_name[index].compare(_property_name) == 0)
				return true;
		return false;
	}

	std::string PhiProperties::phi_property_name(size_t phi_property_index) {
		if (phi_property_index < _phi_property_name.size())
			return _phi_property_name[phi_property_index];
		return "DEFAULT_PHASE";
	}

	void PhiProperties::add_property_name(std::string _property_name) {
		for (size_t index = 0; index < _phi_property_name.size(); index++)
			if (_phi_property_name[index].compare(_property_name) == 0)
				return;
		_phi_property_name.push_back(_property_name);
		_phi_property_number = _phi_property_name.size();
	}
}