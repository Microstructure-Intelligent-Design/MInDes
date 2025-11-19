#pragma once
#include "../../MainIterator_Params.h"
#include "../Modules_Params.h"
#include "../input_modules/ioFiles_Params.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
namespace pf {

	namespace write_vts {
		inline size_t output_frequence = 0;
		inline bool is_show_with_boundary = false;
		inline size_t x_begin = 0;
		inline size_t y_begin = 0;
		inline size_t z_begin = 0;
		inline size_t x_end = 0;
		inline size_t y_end = 0;
		inline size_t z_end = 0;
		inline std::vector<void(*)(std::ofstream& fout)> write_vts_scalar_list;
		inline std::vector<void(*)(std::ofstream& fout)> write_vts_vector_list;
		inline void load_vts_scalar_func(void(*buff)(std::ofstream& fout));
		inline void load_vts_vector_func(void(*buff)(std::ofstream& fout));
		namespace default_functions {
			void open_vts_scalar_file(std::ofstream& fout, std::string tail);
			void open_vts_vec3_file(std::ofstream& fout, std::string tail);
			void write_scalar_grains(std::ofstream& fout);
			void write_scalar_phi_index(std::ofstream& fout);
			void write_scalar_phi_all(std::ofstream& fout);
			void write_scalar_con_all(std::ofstream& fout);
			void write_scalar_temperature(std::ofstream& fout);
			void close_vts_file(std::ofstream& fout);
		}
		void write_vts_pre_iii();

		void write_vts_pos_iii();

		void init_write_vts();
	}
}