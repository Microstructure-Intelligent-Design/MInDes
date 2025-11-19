#pragma once
#include "../../base/MACRO_DEF.h"
#include <random>
#include <iostream>
#include <string>
#include <sstream>
namespace pf {
	using namespace std;
	enum CalculationOperator { CO_bottom, CO_PLus, CO_Minux, CO_Multiply, CO_Divide, CO_ParaSeparator, CO_top };
	const vector<char> CalculationOperator_Keys = { ' ', '+', '-', '*', '/', ',', ' '};
	inline char get_CalculationOperator_Keys(int co_index) {
		return CalculationOperator_Keys[co_index];
	}
	inline int get_CalculationOperator_index(char key) {
		for (int i = CO_bottom; i < CO_top; i++)
			if (key == CalculationOperator_Keys[i])
				return i;
		cout << "> CalculationOperator_Keys error, cant find aim CalculationOperator_index !" << endl;
		SYS_PROGRAM_STOP;
	}
	enum FuncType { 
		FType_bottom,
		FType_real,  // defined variable
		FType_custom,  // defined funcs
		FType_fieldVal,
		FType_pow,
		FType_sqrt,
		FType_abs,
		FType_exp,
		FType_ln,
		FType_log,
		FType_sin, 
		FType_cos, 
		FType_tan,
		FType_asin,
		FType_acos,
		FType_atan,
		FType_top
	};
	const vector<std::string> FuncType_Keys = { "", 
		"", 
		"", 
		"", 
		"pow",
		"sqrt",
		"abs",
		"exp",
		"ln",
		"log",
		"sin",
		"cos",
		"tan",
		"asin",
		"acos",
		"atan",
		"" 
	};
	inline std::string get_FuncType_Keys(int FT_index) {
		return FuncType_Keys[FT_index];
	}
	inline int get_FuncType_index(string key) {
		for (int i = FType_bottom; i < FType_top; i++)
			if (key.compare(FuncType_Keys[i]) == 0)
				return i;
		cout << "> FuncType_Keys error, cant find aim FuncType_index !" << endl;
		SYS_PROGRAM_STOP;
	}
	struct Operators_4d {
		vector<vector<vector<vector<CalculationOperator>>>> operators_4;
		vector<vector<vector<CalculationOperator>>> operators_3;
		vector<vector<CalculationOperator>> operators_2;
		vector<CalculationOperator> operators_1;
		Operators_4d() {

		}
		~Operators_4d() {
			clear();
		}
		void clear() {
			operators_4.clear();
			operators_3.clear();
			operators_2.clear();
			operators_1.clear();
		}
	};
	struct FuncIndexes_4d {
		vector<vector<vector<vector<int>>>> index_4;
		vector<vector<vector<int>>> index_3;
		vector<vector<int>> index_2;
		vector<int> index_1;
		FuncIndexes_4d() {

		}
		~FuncIndexes_4d() {
			clear();
		}
		void clear() {
			index_4.clear();
			index_3.clear();
			index_2.clear();
			index_1.clear();
		}
	};
	struct InFileVar {
		string key;
		REAL var;
		InFileVar() {
			key = "";
			var = 0.0;
		}
	};
	struct InFileFunc {
		string key;
		string func_str;
		vector<vector<vector<vector<vector<int>>>>> func_structure;
		Operators_4d operators;
		FuncIndexes_4d terms_type;
		REAL(*func)(vector<int>&, vector<REAL>&
			, vector<vector<vector<vector<vector<int>>>>>&
			, Operators_4d&
			, FuncIndexes_4d&
			, vector<InFileVar>&
			, vector<InFileFunc>&);
		InFileFunc() {
			key = "";
			func_str = "";
			func = nullptr;
		}
		~InFileFunc() {
			func_structure.clear();
			operators.clear();
			terms_type.clear();
			key = "";
			func_str = "";
			func = nullptr;
		}
	};

	namespace infile_math_default_funcs {
		bool is_string_int(string str, int& val);

		bool is_string_REAL(string str, REAL& val);

		REAL pow(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL sqrt(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL abs(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL exp(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL ln(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL log(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL sin(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL cos(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL tan(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL asin(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL acos(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		REAL atan(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);
		vector<REAL> formula(vector<REAL> terms, vector<CalculationOperator> operators);
		// if func_type < 0 mean special function REAL(), 
		REAL equation(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type,
			vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs);

	}
	class InfileMath
	{
	public:
		bool search_var(string var_key, int& var_index);
		bool search_func(string func_key, int& func_index);
		// add REAL in (,) in input value: Definition.DefaultValue = value_name<REAL_value>
		int add_infile_var(string key, REAL var);

		bool check_default_func(string key);

		bool check_variables_funcs_keys();

		int add_default_func(string key);

		string split_equation_str_by_separators_and_operators(string equation, vector<char> separators, vector<string>& formulas
			, vector<string>& formula_keys, vector<CalculationOperator>& cal_operators);

		string split_field_variable_by_separators_and_operators(string equation, vector<char> separators, vector<string>& formulas
			, vector<string>& formula_keys, vector<CalculationOperator>& cal_operators);

		// add string in <> in input value: Definition.Function = func_name<equation_str>
		bool add_infile_funcs(string key, string func_str, bool debug = false);

		InfileMath() {
			vector<char> seperator1 = { '{','}' };
			vector<char> seperator2 = { '[',']' };
			vector<char> seperator3 = { '(',')' };
			_seperators.push_back(seperator1);
			_seperators.push_back(seperator2);
			_seperators.push_back(seperator3);
			define_func_key = "Define.Func";
			define_variable_key = "Define.Var";
		}
		~InfileMath() {
			clear();
		}
		void clear() {
			infile_vars.clear();
			infile_funcs.clear();
			_seperators.clear();
		}
		vector<InFileVar> infile_vars;
		vector<InFileFunc> infile_funcs;
		vector<vector<char>> _seperators;
		string define_func_key;
		string define_variable_key;
	};

}