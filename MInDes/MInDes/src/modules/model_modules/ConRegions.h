#pragma once
#include "../base/MACRO_DEF.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Modules_Params.h"
#include "PhiProperties.h"
namespace pf {
	class ConRegions
	{
	public:
		ConRegions(const ConRegions&) = delete;
		ConRegions& operator=(const ConRegions&) = delete;

		static ConRegions& instance();
		void init();
		size_t con_index(std::string con_name);
		std::string con_name(size_t con_index);
		size_t region_index(std::string region_name);
		std::string region_name(size_t region_index);
		bool is_con(std::string con_name);
		bool is_region(std::string region_name);
		bool is_con_in_region(size_t region_index, size_t con_index);
		bool is_phi_in_region(size_t region_index, size_t phi_index);

		inline size_t region_number() {
			return _region_number;
		}

		inline size_t region_con(size_t region_index, size_t index) {
			return _region_con(region_index, index);
		}
		inline size_t region_con(std::string region_name, size_t index) {
			return _region_con(region_index(region_name), index);
		}
		inline size_t region_con_number(size_t region_index) {
			return _region_con_number[region_index];
		}
		inline size_t region_con_number(std::string region_name) {
			return _region_con_number[region_index(region_name)];
		}
		inline size_t con_region(size_t con_index) {
			return _con_region[con_index];
		}
		inline size_t con_region(std::string con_name) {
			return _con_region[con_index(con_name)];
		}
		inline size_t region_phi(size_t region_index, size_t index) {
			return _region_phi(region_index, index);
		}
		inline size_t region_phi(std::string region_name, size_t index) {
			return _region_phi(region_index(region_name), index);
		}
		inline size_t region_phi_number(size_t region_index) {
			return _region_phi_number[region_index];
		}
		inline size_t region_phi_number(std::string region_name) {
			return _region_phi_number[region_index(region_name)];
		}
		inline size_t phi_region(size_t phi_index) {
			return _phi_region[phi_index];
		}
		inline size_t phi_property_region(size_t phi_property) {
			return _phi_property_region[phi_property];
		}
		inline size_t region_phi_property(size_t region_index, size_t index) {
			return _region_phi_property(region_index, index);
		}
		inline size_t region_phi_property(std::string region_name, size_t index) {
			return _region_phi_property(region_index(region_name), index);
		}
		inline size_t region_phi_property_number(size_t region_index) {
			return _region_phi_property_number[region_index];
		}
		inline size_t region_phi_property_number(std::string region_name) {
			return _region_phi_property_number[region_index(region_name)];
		}
	private:
		size_t add_region_name(std::string region_name);
		ConRegions();
		~ConRegions() = default;
		bool _is_init = false;
		size_t _region_number = 0;
		std::vector<std::string> _con_name;
		std::vector<std::string> _region_name;
		std::vector<size_t> _con_region;    // - find region index by con index
		std::vector<size_t> _phi_region;    // - find region index by phi index
		std::vector<size_t> _phi_property_region;        // - find region index by phi property
		Matrix2D<size_t> _region_con;                    // - find con index by region index
		std::vector<size_t> _region_con_number;          // - find con index number by region index
		Matrix2D<size_t> _region_phi;                    // - find phi index by region index
		std::vector<size_t> _region_phi_number;          // - find phi index number by region index
		Matrix2D<size_t> _region_phi_property;           // - find phi property by region index
		std::vector<size_t> _region_phi_property_number; // - find phi property number by region index
	};
}