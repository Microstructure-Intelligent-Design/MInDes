#pragma once
#include "../base/MACRO_DEF.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Modules_Params.h"
namespace pf {
	class ConProperties
	{
	public:
		ConProperties(const ConProperties&) = delete;
		ConProperties& operator=(const ConProperties&) = delete;

		static ConProperties& instance();
		size_t con_index(std::string _con_name);
		size_t con_number();
		bool is_con(std::string _con_name);
		std::string con_name(int con_index);
		std::string con_name(size_t con_index);
		void init();

	private:
		ConProperties();
		~ConProperties() = default;
		bool _is_init = false;
		size_t _con_number = 0;
		std::vector<std::string> _con_name;
	};
}