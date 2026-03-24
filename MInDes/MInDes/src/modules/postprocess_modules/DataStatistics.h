#pragma once
#include "../base/MACRO_DEF.h"
namespace pf {
	namespace data_statistics_functions {
		inline string statistics_separator = "     ";
		inline std::vector<std::pair<std::string, REAL>> statistics_data;
        inline bool update_statistics(const std::string key, REAL newValue) {
            for (auto& item : statistics_data) {
                if (item.first == key) {
                    item.second = newValue;
                    return true;
                }
            }
            return false;
        }
	}
}