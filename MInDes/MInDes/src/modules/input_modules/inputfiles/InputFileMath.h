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
		inline bool is_string_int(string str, int& val) {
			// 10000000
			for (int index = 0; index < str.size(); index++) {
				char c = str.at(index);
				if (c >= '0' && c <= '9')
					continue;
				else if (c != '-' && c != '+')
					return false;
			}
			val = stoi(str);
			return true;
		}

		inline bool is_string_REAL(string str, REAL& val) {
			// -123.234E-324  , one E/e one . two sign symble -/+
			int index_of_point = -1, index_of_E = -1;
			vector<int> sign_sym_index = { -1, -1 };
			for (int index = 0; index < str.size(); index++) {
				char c = str.at(index);
				if (c >= '0' && c <= '9')
					continue;
				else if (c != '-' && c != '+' && c != '.' && c != 'e' && c != 'E')
					return false;
				if (c == '-' or c == '+') {
					if (index == 0)
						sign_sym_index[0] = index;
					else if (index_of_E != -1 && index == index_of_E + 1)
						sign_sym_index[1] = index;
					else
						return false;
				}
				if (c == '.') {
					if (index_of_point == -1)
						index_of_point = index;
					else
						return false;
				}
				if (c == 'e' or c == 'E') {
					if (index_of_E == -1)
						index_of_E = index;
					else
						return false;
				}
			}
			if (index_of_E == -1 && sign_sym_index[1] != -1)
				return false;
			if (index_of_point == 0 || index_of_point == str.size() - 1)
				return false;
			if (index_of_E == 0)
				return false;
			if (index_of_point != -1 && index_of_E != -1) {
				if (index_of_point > index_of_E || index_of_point == index_of_E - 1)
					return false;
			}
			val = REAL(stod(str));
			return true;
		}

		inline REAL pow(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::pow(para_vals[0], para_vals[1]);
		}
		inline REAL sqrt(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::sqrt(para_vals[0]);
		}
		inline REAL abs(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::fabs(para_vals[0]);
		}
		inline REAL exp(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::exp(para_vals[0]);
		}
		inline REAL ln(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::log(para_vals[0]);
		}
		inline REAL log(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::log10(para_vals[1]) / std::log10(para_vals[0]);
		}
		static REAL sin(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::sin(para_vals[0]);
		}
		inline REAL cos(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::cos(para_vals[0]);
		}
		inline REAL tan(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::tan(para_vals[0]);
		}
		inline REAL asin(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::asin(para_vals[0]);
		}
		inline REAL acos(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::acos(para_vals[0]);
		}
		inline REAL atan(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type
			, vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			return std::atan(para_vals[0]);
		}
		inline vector<REAL> formula(vector<REAL> terms, vector<CalculationOperator> operators) {
			int term = int(terms.size());
			if (term < 1 || operators[0] == CalculationOperator::CO_ParaSeparator || *(operators.end() - 1) == CalculationOperator::CO_ParaSeparator
				|| operators[0] == CalculationOperator::CO_Multiply || operators[0] == CalculationOperator::CO_Divide) {
				cout << "> formula error, terms.size() < 1 || operators[0] = , or * or /" << endl;
				SYS_PROGRAM_STOP;
			}
			vector<REAL> return_vals;
			bool is_ParaSeparator_exist = false;
			vector<REAL> split_terms; vector<CalculationOperator>  split_operators;
			int diff = 0;
			for (int index = 0; index < operators.size(); index++) {
				if (operators[index] == CalculationOperator::CO_ParaSeparator) {
					vector<REAL> cal_values = formula(split_terms, split_operators);
					return_vals.push_back(cal_values[0]);
					split_terms.clear();
					split_operators.clear();
					is_ParaSeparator_exist = true;
					diff--;
				}
				else {
					split_operators.push_back(operators[index]);
					split_terms.push_back(terms[index + diff]);
				}
			}
			if (is_ParaSeparator_exist)
				return return_vals;
			else {
				return_vals.clear(); return_vals.push_back(0.0);
				if (terms.size() != operators.size()) {
					cout << "> formula error, terms.size() != operators.size()" << endl;
					SYS_PROGRAM_STOP;
				}
				for (int index = 1; index < operators.size();) {
					switch (operators[index])
					{
					case pf::CO_Multiply:
						terms[index - 1] *= terms[index];
						operators.erase(operators.begin() + index);
						terms.erase(terms.begin() + index);
						break;
					case pf::CO_Divide:
						terms[index - 1] /= terms[index];
						operators.erase(operators.begin() + index);
						terms.erase(terms.begin() + index);
						break;
					default:
						index++;
						break;
					}
				}
				for (int index = 0; index < operators.size(); index++) {
					switch (operators[index])
					{
					case pf::CO_PLus:
						return_vals[0] += terms[index];
						break;
					case pf::CO_Minux:
						return_vals[0] -= terms[index];
						break;
					default:
						cout << "> formula error, operators cant be explained" << endl;
						SYS_PROGRAM_STOP;
						break;
					}
				}
				return return_vals;
			}
		}
		// if func_type < 0 mean special function REAL(), 
		inline REAL equation(vector<int>& para_index, vector<REAL>& para_vals
			, vector<vector<vector<vector<vector<int>>>>>& _func_structure, Operators_4d& _operators, FuncIndexes_4d& _terms_type,
			vector<InFileVar>& infile_vars, vector<InFileFunc>& infile_funcs) {
			vector<REAL> val_i;
			for (int index_i = 0; index_i < _func_structure.size(); index_i++) {
				vector<REAL> val_j;
				for (int index_j = 0; index_j < _func_structure[index_i].size(); index_j++) {
					vector<REAL> val_k;
					for (int index_k = 0; index_k < _func_structure[index_i][index_j].size(); index_k++) {
						vector<REAL> val_l;
						for (int index_l = 0; index_l < _func_structure[index_i][index_j][index_k].size(); index_l++) {
							if (_func_structure[index_i][index_j][index_k][index_l].size() < 1) {
								cout << "> input func error, func_structure error !" << endl;
								SYS_PROGRAM_STOP;
							}
							if (_terms_type.index_4[index_i][index_j][index_k][index_l] < 0) {
								val_l.push_back(infile_vars[_func_structure[index_i][index_j][index_k][index_l][0]].var);
							}
							else {
								vector<int> para_ints = _func_structure[index_i][index_j][index_k][index_l];
								para_ints.erase(para_ints.begin());
								vector<REAL> para_REALs;
								pf::InFileFunc& infile_func = infile_funcs[_func_structure[index_i][index_j][index_k][index_l][0]];
								val_l.push_back(infile_func.func(para_ints, para_REALs, infile_func.func_structure, infile_func.operators, infile_func.terms_type, infile_vars, infile_funcs));
							}
						}
						val_l = formula(val_l, _operators.operators_4[index_i][index_j][index_k]);
						if (_terms_type.index_3[index_i][index_j][index_k] < 0 && val_l.size() == 1) {
							val_k.push_back(val_l[0]);
						}
						else {
							vector<int> para_ints;
							pf::InFileFunc& infile_func = infile_funcs[_terms_type.index_3[index_i][index_j][index_k]];
							val_k.push_back(infile_func.func(para_ints, val_l, infile_func.func_structure, infile_func.operators, infile_func.terms_type, infile_vars, infile_funcs));
						}
					}
					val_k = formula(val_k, _operators.operators_3[index_i][index_j]);
					if (_terms_type.index_2[index_i][index_j] < 0 && val_k.size() == 1) {
						val_j.push_back(val_k[0]);
					}
					else {
						vector<int> para_ints;
						pf::InFileFunc& infile_func = infile_funcs[_terms_type.index_2[index_i][index_j]];
						val_j.push_back(infile_func.func(para_ints, val_k, infile_func.func_structure, infile_func.operators, infile_func.terms_type, infile_vars, infile_funcs));
					}
				}
				val_j = formula(val_j, _operators.operators_2[index_i]);
				if (_terms_type.index_1[index_i] < 0 && val_j.size() == 1) {
					val_i.push_back(val_j[0]);
				}
				else {
					vector<int> para_ints;
					pf::InFileFunc& infile_func = infile_funcs[_terms_type.index_1[index_i]];
					val_i.push_back(infile_func.func(para_ints, val_j, infile_func.func_structure, infile_func.operators, infile_func.terms_type, infile_vars, infile_funcs));
				}
			}
			val_i = formula(val_i, _operators.operators_1);
			if (val_i.size() != 1) {
				cout << "> input func error, func_structure error !" << endl;
				SYS_PROGRAM_STOP;
			}
			return val_i[0];
		}

	}
	class InfileMath
	{
	public:
		bool search_var(string var_key, int& var_index) {
			for (int index = 0; index < infile_vars.size(); index++)
				if (infile_vars[index].key.compare(var_key) == 0) {
					var_index = index;
					return true;
				}
			return false;
		}
		bool search_func(string func_key, int& func_index) {
			for (int index = 0; index < infile_funcs.size(); index++)
				if (infile_funcs[index].key.compare(func_key) == 0) {
					func_index = index;
					return true;
				}
			return false;
		}
		// add REAL in (,) in input value: Definition.DefaultValue = value_name<REAL_value>
		int add_infile_var(string key, REAL var) {
			for (int in_var_index = 0; in_var_index < infile_vars.size(); in_var_index++)
				if (key.compare(infile_vars[in_var_index].key) == 0) {
					infile_vars[in_var_index].var = var;
					cout << "> Wainning, custom variable: " << key << " has been defined multi-times !" << endl;
					return in_var_index;
				}
			InFileVar new_var;
			new_var.key = key;
			new_var.var = var;
			infile_vars.push_back(new_var);
			return int(infile_vars.size()) - 1;
		}

		bool check_default_func(string key) {
			if (key.compare("") == 0)
				return false;
			for (int index = FuncType::FType_bottom; index < FuncType::FType_top; index++)
				if (FuncType_Keys[index].compare(key) == 0)
					return true;
			return false;
		}

		bool check_variables_funcs_keys() {
			for (int index1 = 0; index1 < infile_vars.size(); index1++)
				for (int index2 = 0; index2 < infile_funcs.size(); index2++)
					if (infile_vars[index1].key.compare(infile_funcs[index2].key) == 0) {
						return false;
					}
			return true;
		}

		int add_default_func(string key) {
			if (key.compare("") == 0)
				return false;
			for (int func_index = 0; func_index < infile_funcs.size(); func_index++)
				if (infile_funcs[func_index].key.compare(key) == 0)
					return func_index;
			InFileFunc newFunc;
			newFunc.key = key;
			switch (get_FuncType_index(key))
			{
			case pf::FType_pow:
				newFunc.func = infile_math_default_funcs::pow;
				break;
			case pf::FType_sqrt:
				newFunc.func = infile_math_default_funcs::sqrt;
				break;
			case pf::FType_abs:
				newFunc.func = infile_math_default_funcs::abs;
				break;
			case pf::FType_exp:
				newFunc.func = infile_math_default_funcs::exp;
				break;
			case pf::FType_ln:
				newFunc.func = infile_math_default_funcs::ln;
				break;
			case pf::FType_log:
				newFunc.func = infile_math_default_funcs::log;
				break;
			case pf::FType_sin:
				newFunc.func = infile_math_default_funcs::sin;
				break;
			case pf::FType_cos:
				newFunc.func = infile_math_default_funcs::cos;
				break;
			case pf::FType_tan:
				newFunc.func = infile_math_default_funcs::tan;
				break;
			case pf::FType_asin:
				newFunc.func = infile_math_default_funcs::asin;
				break;
			case pf::FType_acos:
				newFunc.func = infile_math_default_funcs::acos;
				break;
			case pf::FType_atan:
				newFunc.func = infile_math_default_funcs::atan;
				break;
			default:
				cout << "> Error, add_default_func error !" << endl;
				SYS_PROGRAM_STOP;
				break;
			}
			infile_funcs.push_back(newFunc);
			return int(infile_funcs.size()) - 1;
		}

		string split_equation_str_by_separators_and_operators(string equation, vector<char> separators, vector<string>& formulas
			, vector<string>& formula_keys, vector<CalculationOperator>& cal_operators) {
			bool is_split = false, is_pair_seperators = true;
			// check equation degree
			if (separators.size() != 0) {
				vector<char> check_separator;
				for (auto c = equation.begin(); c < equation.end(); c++)
					if (*c == separators[0] || *c == separators[1])
						check_separator.push_back(*c);
				for (int index = 0; index < check_separator.size(); index++)
					if (check_separator[index % 2] != separators[index % 2])
						is_pair_seperators = false;
				if (check_separator.size() % 2 == 0 && is_pair_seperators) {
					is_split = true;
				}
			}
			else {
				is_split = true;
			}
			// there has separators in equation
			if (is_split) {
				string current_read_str = ""; bool separator_begin = false;
				int separator_type = -1;
				// check every char in equation
				for (int cindex = 0; cindex < equation.size(); cindex++) {
					char c = equation.at(cindex);
					// for the first char Default is a operator '-' or '+', if no, auto set '+'
					if (cindex == 0 || (cal_operators.size() != 0 && *(cal_operators.end() - 1) == CalculationOperator::CO_ParaSeparator)) {
						switch (c)
						{
						case '+':
							cal_operators.push_back(CalculationOperator::CO_PLus);
							break;
						case '-':
							cal_operators.push_back(CalculationOperator::CO_Minux);
							break;
						default:
							cal_operators.push_back(CalculationOperator::CO_PLus);
							break;
						}
						if (c == '+' || c == '-')
							continue;
						for (int sep_index = 0; sep_index < _seperators.size(); sep_index++) {
							if (c == _seperators[sep_index][0]) {
								separator_begin = true;
								separator_type = sep_index;
								formula_keys.push_back("");
								continue;
							}
						}
						if (separator_begin)
							continue;
						current_read_str.push_back(c);
					}
					else {
						// if this char in the separators, if it isnt, equation is seperated for operators and separators[0]
						if (!separator_begin) {
							bool is_operator = false;
							switch (c)
							{
							case '+':
								cal_operators.push_back(CalculationOperator::CO_PLus);
								is_operator = true;
								break;
							case '-':
								cal_operators.push_back(CalculationOperator::CO_Minux);
								is_operator = true;
								break;
							case '*':
								cal_operators.push_back(CalculationOperator::CO_Multiply);
								is_operator = true;
								break;
							case '/':
								cal_operators.push_back(CalculationOperator::CO_Divide);
								is_operator = true;
								break;
							case ',':
								cal_operators.push_back(CalculationOperator::CO_ParaSeparator);
								is_operator = true;
								break;
							default:
								break;
							}
							if (is_operator) {
								// formula
								if (current_read_str.size() != 0) {
									formulas.push_back(current_read_str);
									formula_keys.push_back("");
									current_read_str.clear();
								}
							}
							else {
								for (int sep_index = 0; sep_index < _seperators.size(); sep_index++) {
									if (c == _seperators[sep_index][0]) {
										// functions or formula
										separator_begin = true;
										separator_type = sep_index;
										// functions
										formula_keys.push_back(current_read_str);
										current_read_str.clear();
									}
								}
								if(!separator_begin)
									current_read_str.push_back(c);
							}
						}
						else {
							if (c == _seperators[separator_type][1]) {
								if (current_read_str == "") {
									cout << "> error, no formula between " << separators[0] << separators[1] << " in equation: " << equation << endl;
									SYS_PROGRAM_STOP;
								}
								formulas.push_back(current_read_str);
								current_read_str.clear();
								separator_begin = false;
								separator_type = -1;
							}
							else
								current_read_str.push_back(c);
						}
					}
				}
				if (current_read_str.size() != 0) {
					formulas.push_back(current_read_str);
					formula_keys.push_back("");
				}
			}
			else if (is_pair_seperators) {
				formulas.push_back(equation); // no separators means the whole equation is included in separators, and no formula_keys, no cal_operators
				cal_operators.push_back(CalculationOperator::CO_PLus);
				formula_keys.push_back("");
			}
			else {
				cout << "> error, equation : " << equation << " cant be explained !" << endl;
				SYS_PROGRAM_STOP;
			}
			int comma_size = 0;
			for (int index = 0; index < cal_operators.size(); index++)
				if (cal_operators[index] == CalculationOperator::CO_ParaSeparator)
					comma_size++;
			if ((cal_operators.size() - comma_size) != formula_keys.size() && formula_keys.size() != formula_keys.size()) {
				cout << "> error, equation split funtions size mismatch !" << endl;
				SYS_PROGRAM_STOP;
			}
			stringstream report;
			report << "[FUNREAD] " << equation << " = ";
			int diff = 0;
			if (separators.size() == 2) {
				for (int index = 0; index < cal_operators.size(); index++) {
					if (cal_operators[index] == CalculationOperator::CO_ParaSeparator) {
						report << get_CalculationOperator_Keys(cal_operators[index]);
						diff--;
					}
					else {
						report << get_CalculationOperator_Keys(cal_operators[index]);
						if (formula_keys[index + diff].compare("") == 0)
							report << "REAL" << separators[0];
						else
							report << formula_keys[index + diff] << separators[0];
						report << formulas[index + diff];
						report << separators[1];
					}
				}
			}
			else {
				for (int index = 0; index < cal_operators.size(); index++) {
					if (cal_operators[index] == CalculationOperator::CO_ParaSeparator) {
						report << get_CalculationOperator_Keys(cal_operators[index]);
						diff--;
					}
					else {
						report << get_CalculationOperator_Keys(cal_operators[index]);
						if (formula_keys[index + diff].compare("") == 0)
							report << "REAL(";
						else
							report << formula_keys[index + diff] << "(";
						report << formulas[index + diff];
						report << ")";
					}
				}
			}
			return report.str();
		}

		string split_field_variable_by_separators_and_operators(string equation, vector<char> separators, vector<string>& formulas
			, vector<string>& formula_keys, vector<CalculationOperator>& cal_operators) {
			bool is_split = false, is_pair_seperators = true;
			// check equation degree
			if (separators.size() != 0) {
				vector<char> check_separator;
				for (auto c = equation.begin(); c < equation.end(); c++)
					if (*c == separators[0] || *c == separators[1])
						check_separator.push_back(*c);
				for (int index = 0; index < check_separator.size(); index++)
					if (check_separator[index % 2] != separators[index % 2])
						is_pair_seperators = false;
				if (check_separator.size() % 2 == 0 && is_pair_seperators) {
					is_split = true;
				}
			}
			else {
				is_split = true;
			}
			// there has separators in equation
			if (is_split) {
				string current_read_str = ""; bool separator_begin = false;
				// check every char in equation
				for (int cindex = 0; cindex < equation.size(); cindex++) {
					char c = equation.at(cindex);
					// for the first char Default is a operator '-' or '+', if no, auto set '+'
					if (cindex == 0 || (cal_operators.size() != 0 && *(cal_operators.end() - 1) == CalculationOperator::CO_ParaSeparator)) {
						switch (c)
						{
						case '+':
							cal_operators.push_back(CalculationOperator::CO_PLus);
							break;
						case '-':
							cal_operators.push_back(CalculationOperator::CO_Minux);
							break;
						default:
							cal_operators.push_back(CalculationOperator::CO_PLus);
							break;
						}
						if (c == '+' || c == '-')
							continue;
						if (c == separators[0]) {
							separator_begin = true;
							continue;
						}
						if (separator_begin)
							continue;
						current_read_str.push_back(c);
					}
					else {
						// if this char in the separators, if it isnt, equation is seperated for operators and separators[0]
						if (!separator_begin) {
							bool is_operator = false;
							switch (c)
							{
							case '+':
								cal_operators.push_back(CalculationOperator::CO_PLus);
								is_operator = true;
								break;
							case '-':
								cal_operators.push_back(CalculationOperator::CO_Minux);
								is_operator = true;
								break;
							case '*':
								cal_operators.push_back(CalculationOperator::CO_Multiply);
								is_operator = true;
								break;
							case '/':
								cal_operators.push_back(CalculationOperator::CO_Divide);
								is_operator = true;
								break;
							case ',':
								cal_operators.push_back(CalculationOperator::CO_ParaSeparator);
								is_operator = true;
								break;
							default:
								break;
							}
							if (is_operator) {
								// formula
								if (current_read_str.size() != 0) {
									formulas.push_back(current_read_str);
									formula_keys.push_back("");
									current_read_str.clear();
								}
							}
							else {
								if (c == separators[0]) {
									// functions or formula
									separator_begin = true;
									// functions
									formula_keys.push_back(current_read_str);
									current_read_str.clear();
								}
								if (!separator_begin)
									current_read_str.push_back(c);
							}
						}
						else {
							if (c == separators[1]) {
								if (current_read_str == "") {
									cout << "> error, no formula between " << separators[0] << separators[1] << " in equation: " << equation << endl;
									SYS_PROGRAM_STOP;
								}
								formulas.push_back(current_read_str);
								current_read_str.clear();
								separator_begin = false;
							}
							else
								current_read_str.push_back(c);
						}
					}
				}
				if (current_read_str.size() != 0) {
					formulas.push_back(current_read_str);
					formula_keys.push_back("");
				}
			}
			else if (is_pair_seperators) {
				formulas.push_back(equation); // no separators means the whole equation is included in separators, and no formula_keys, no cal_operators
				cal_operators.push_back(CalculationOperator::CO_PLus);
				formula_keys.push_back("");
			}
			else {
				cout << "> error, equation : " << equation << " cant be explained !" << endl;
				SYS_PROGRAM_STOP;
			}
			int comma_size = 0;
			for (int index = 0; index < cal_operators.size(); index++)
				if (cal_operators[index] == CalculationOperator::CO_ParaSeparator)
					comma_size++;
			if ((cal_operators.size() - comma_size) != formula_keys.size() && formula_keys.size() != formula_keys.size()) {
				cout << "> error, equation split funtions size mismatch !" << endl;
				SYS_PROGRAM_STOP;
			}
			stringstream report;
			report << "[FUNREAD] " << equation << " = ";
			int diff = 0;
			for (int index = 0; index < cal_operators.size(); index++) {
				if (cal_operators[index] == CalculationOperator::CO_ParaSeparator) {
					report << get_CalculationOperator_Keys(cal_operators[index]);
					diff--;
				}
				else {
					report << get_CalculationOperator_Keys(cal_operators[index]);
					if (formula_keys[index + diff].compare("") == 0)
						report << "REAL" << separators[0];
					else
						report << formula_keys[index + diff] << separators[0];
					report << formulas[index + diff];
					report << separators[1];
				}
			}
			return report.str();
		}

		// add string in <> in input value: Definition.Function = func_name<equation_str>
		bool add_infile_funcs(string key, string func_str, bool debug = false) {
			for (int in_func_index = 0; in_func_index < infile_funcs.size(); in_func_index++)
				if (key.compare(infile_funcs[in_func_index].key) == 0) {
					if (func_str.compare(infile_funcs[in_func_index].func_str) == 0) {
						if(debug)
							cout << "> Wainning, custom function: " << key << " has been defined multi-times !" << endl;
						return false;
					}
					else {
						cout << "> Error, custom function: " << key << " has been defined multi-times with different equation !" << endl;
						SYS_PROGRAM_STOP;
					}
				}
			InFileFunc new_func;
			new_func.key = key;
			new_func.func_str = func_str;
			new_func.func = infile_math_default_funcs::equation;
			// explain equations between { }
			vector<char> separators1 = { '{', '}' };
			vector<string> formulas1, formula_keys1; vector<CalculationOperator> cal_operators1;
			string str = split_equation_str_by_separators_and_operators(func_str, separators1, formulas1, formula_keys1, cal_operators1);
			if (debug)
				cout << str << endl;
			new_func.func_structure.resize(formulas1.size());
			new_func.operators.operators_2.resize(formulas1.size());
			new_func.operators.operators_3.resize(formulas1.size());
			new_func.operators.operators_4.resize(formulas1.size());
			new_func.terms_type.index_2.resize(formulas1.size());
			new_func.terms_type.index_3.resize(formulas1.size());
			new_func.terms_type.index_4.resize(formulas1.size());
			for (auto func_key = formula_keys1.begin(); func_key < formula_keys1.end(); func_key++) {
				int func_index = -1;
				if (search_func(*func_key, func_index)) {
					new_func.terms_type.index_1.push_back(func_index);
				}
				else if ((*func_key).compare("") == 0) {
					new_func.terms_type.index_1.push_back(-1);
				}
				else if (check_default_func(*func_key)) {
					func_index = add_default_func(*func_key);
					new_func.terms_type.index_1.push_back(func_index);
				}
				else {
					cout << "> Error, custom function: " << (*func_key) << " hasnt been defined ! or is not defined before function : "
						<< key << endl;
					SYS_PROGRAM_STOP;
				}
			}
			new_func.operators.operators_1 = cal_operators1;

			// explain equations between [ ]
			for (int index1 = 0; index1 < new_func.func_structure.size(); index1++) {
				vector<char> separators2 = { '[', ']' };
				vector<string> formulas2, formula_keys2; vector<CalculationOperator> cal_operators2;
				string str1 = split_equation_str_by_separators_and_operators(formulas1[index1], separators2, formulas2, formula_keys2, cal_operators2);
				if (debug)
					cout << str1 << endl;
				new_func.func_structure[index1].resize(formulas2.size());
				new_func.operators.operators_3[index1].resize(formulas2.size());
				new_func.operators.operators_4[index1].resize(formulas2.size());
				new_func.terms_type.index_3[index1].resize(formulas2.size());
				new_func.terms_type.index_4[index1].resize(formulas2.size());
				for (auto func_key = formula_keys2.begin(); func_key < formula_keys2.end(); func_key++) {
					int func_index = -1;
					if (search_func(*func_key, func_index)) {
						new_func.terms_type.index_2[index1].push_back(func_index);
					}
					else if ((*func_key).compare("") == 0) {
						new_func.terms_type.index_2[index1].push_back(-1);
					}
					else if (check_default_func(*func_key)) {
						func_index = add_default_func(*func_key);
						new_func.terms_type.index_2[index1].push_back(func_index);
					}
					else {
						cout << "> Error, custom function: " << (*func_key) << " hasnt been defined ! or is not defined before function : "
							<< key << endl;
						SYS_PROGRAM_STOP;
					}
				}
				new_func.operators.operators_2[index1] = cal_operators2;
				// explain equations between ( )
				for (int index2 = 0; index2 < new_func.func_structure[index1].size(); index2++) {
					vector<char> separators3 = { '(', ')' };
					vector<string> formulas3, formula_keys3; vector<CalculationOperator> cal_operators3;
					string str2 = split_equation_str_by_separators_and_operators(formulas2[index2], separators3, formulas3, formula_keys3, cal_operators3);
					if (debug)
						cout << str2 << endl;
					new_func.func_structure[index1][index2].resize(formulas3.size());
					new_func.operators.operators_4[index1][index2].resize(formulas3.size());
					new_func.terms_type.index_4[index1][index2].resize(formulas3.size());
					for (auto func_key = formula_keys3.begin(); func_key < formula_keys3.end(); func_key++) {
						int func_index = -1;
						if (search_func(*func_key, func_index)) {
							new_func.terms_type.index_3[index1][index2].push_back(func_index);
						}
						else if ((*func_key).compare("") == 0) {
							new_func.terms_type.index_3[index1][index2].push_back(-1);
						}
						else if (check_default_func(*func_key)) {
							func_index = add_default_func(*func_key);
							new_func.terms_type.index_3[index1][index2].push_back(func_index);
						}
						else {
							cout << "> Error, custom function: " << (*func_key) << " hasnt been defined ! or is not defined before function : "
								<< key << endl;
							SYS_PROGRAM_STOP;
						}
					}
					new_func.operators.operators_3[index1][index2] = cal_operators3;
					// explain equations between field_vriable<> and REAL() in ( )
					for (int index3 = 0; index3 < new_func.func_structure[index1][index2].size(); index3++) {
						vector<char> separators4 = { '<', '>' };
						vector<string> formulas4, formula_keys4; vector<CalculationOperator> cal_operators4;
						string str3 = split_field_variable_by_separators_and_operators(formulas3[index3], separators4, formulas4, formula_keys4, cal_operators4);
						if (debug)
							cout << str3 << endl;
						new_func.func_structure[index1][index2][index3].resize(formulas4.size());
						for (int _key_index = 0; _key_index < formula_keys4.size(); _key_index++) {
							int _var_index = -1; string _key = formula_keys4[_key_index]; REAL dvalue = 0.0;
							if (search_func(formula_keys4[_key_index], _var_index)) {
								new_func.terms_type.index_4[index1][index2][index3].push_back(_var_index);
								new_func.func_structure[index1][index2][index3][_key_index].push_back(_var_index);
							}
							else if (search_var(formulas4[_key_index], _var_index)) {  // REAL
								new_func.terms_type.index_4[index1][index2][index3].push_back(-1);
								new_func.func_structure[index1][index2][index3][_key_index].push_back(_var_index);
							}
							else if (infile_math_default_funcs::is_string_REAL(formulas4[_key_index], dvalue)) {  // REAL
								_var_index = add_infile_var(formulas4[_key_index], dvalue);
								new_func.terms_type.index_4[index1][index2][index3].push_back(-1);
								new_func.func_structure[index1][index2][index3][_key_index].push_back(_var_index);
							}
							else {
								cout << "> Error, custom Variables or default Field Variables: " << formulas4[_key_index] << ", function type : " << _key << " hasnt been defined ! or is not defined before function : "
									<< key << endl;
								SYS_PROGRAM_STOP;
							}
						}
						new_func.operators.operators_4[index1][index2][index3] = cal_operators4;

					}
				}
			}

			infile_funcs.push_back(new_func);
			return int(infile_funcs.size()) - 1;
		}

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