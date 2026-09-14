#pragma once
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Modules_Params.h"
#include "../Module.h"
namespace pf {
	namespace pretreatment {
		// filling based on phi index
		inline bool is_fill_grains = false;
		enum FILLING_MODE { FM_BULK, FM_HALF_INT, FM_FULL_INT, FM_INTERPOLATION_INT };
		inline FILLING_MODE fill_mode = FILLING_MODE::FM_HALF_INT;
		inline size_t target_number = 0;
		inline std::vector<std::vector<size_t>> target_phi_index;
		inline std::vector<bool> is_fill_con;
		inline std::vector<std::vector<std::pair<size_t, REAL>>> fill_con;
		inline std::vector<bool> is_fill_temp;
		inline std::vector<REAL> fill_temp;
		// for merge phases
		inline bool is_merge_grains = false;
		inline std::vector<std::vector<size_t>> grains_need_merged;
		// auto_merge_phis
		inline bool is_auto_merge_grains = false;
		inline std::vector<size_t> auto_merged_grains;
		// for re-ordering grains
		inline bool is_reorder_grains = false;
		// -
		namespace functions {
			// - filling function, fill con and temp by phi index
			void filling();
			// - merging function, merging multi-grains into one grain
			void merging();
			// - auto-merging function, merging multi-grains into some grains
			void auto_merging();
			// - re-ordering function, move volume 0 grains to high index
			void reordering();
		}
		// -
		void init_pretreatment();
		// -
		void init_pretreatment_pre_i();
		// - 
		void dinit();
	}
}