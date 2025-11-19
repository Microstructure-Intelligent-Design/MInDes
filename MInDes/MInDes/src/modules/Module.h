#pragma once
#include <vector>
namespace pf {
	struct Solver_Module
	{
		// executes in the preprocess
		void(*exec_pre_i)() = nullptr;
		void(*exec_pre_ii)() = nullptr;
		void(*exec_pre_iii)() = nullptr;
		// executes in the model scheme
		void(*exec_i)() = nullptr;
		void(*exec_ii)() = nullptr;
		void(*exec_iii)() = nullptr;
		// executes in the postprocess
		void(*exec_pos_i)() = nullptr;
		void(*exec_pos_ii)() = nullptr;
		void(*exec_pos_iii)() = nullptr;
		// delete module
		void(*deinit)() = nullptr;
	};

	inline std::vector<Solver_Module> module_list;

	// to creat a new module
	inline void load_a_new_module(void(*exec_pre_i)(), void(*exec_pre_ii)(), void(*exec_pre_iii)(),
		void(*exec_i)(), void(*exec_ii)(), void(*exec_iii)(),
		void(*exec_pos_i)(), void(*exec_pos_ii)(), void(*exec_pos_iii)(), void(*deinit)()) {
		Solver_Module _module;
		_module.exec_pre_i = exec_pre_i;
		_module.exec_pre_ii = exec_pre_ii;
		_module.exec_pre_iii = exec_pre_iii;
		_module.exec_i = exec_i;
		_module.exec_ii = exec_ii;
		_module.exec_iii = exec_iii;
		_module.exec_pos_i = exec_pos_i;
		_module.exec_pos_ii = exec_pos_ii;
		_module.exec_pos_iii = exec_pos_iii;
		_module.deinit = deinit;
		module_list.push_back(_module);
	}

	void register_all_modules();

}

