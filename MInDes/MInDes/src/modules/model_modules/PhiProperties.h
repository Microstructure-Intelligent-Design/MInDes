#pragma once
#include "../base/MACRO_DEF.h"
#include <vector>
namespace pf {
	class PhiProperties
	{
	public:
		PhiProperties(const PhiProperties&) = delete;
		PhiProperties& operator=(const PhiProperties&) = delete;

		static PhiProperties& instance();
		size_t phi_property(size_t phi_index);
		int phi_property(int phi_index);
		size_t phi_property(std::string _property_name);
		size_t phi_property_number();
		bool is_phi_property(std::string _property_name);
		std::string phi_property_name(int phi_property_index);
		void init();

	private:
		void add_property_name(std::string _property_name);
		PhiProperties();
		~PhiProperties() = default;
		bool _is_init = false;
		size_t _phi_property_number = 0;
		std::vector<std::string> _phi_property_name;
		std::vector<size_t> _phi_property;
	};
}