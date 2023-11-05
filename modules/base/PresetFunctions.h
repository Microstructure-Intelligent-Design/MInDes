/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include "../../solvers/base.h"

namespace pf {
	using namespace std;
	class HighorderLinearFunction
	{
	public:
		HighorderLinearFunction() {
			order_num = 0;
			variate_num = 0;
		}
		// _parameters[variates][orders]
		HighorderLinearFunction(vector<vector<double>> _parameters) {
			parameters = _parameters;
			variate_num = (unsigned int)(_parameters.size());
			order_num = (unsigned int)(_parameters[0].size());
		}
		~HighorderLinearFunction() {
			parameters.clear();
		}
		void init(unsigned int _variate_num, unsigned int _order_num) {
			for (auto par = parameters.begin(); par < parameters.end(); par++)
				par->clear();
			parameters.clear();
			variate_num = _variate_num;
			order_num = _order_num;
			for (unsigned int var_i = 0; var_i < variate_num; var_i++) {
				vector<double> par;
				par.resize(order_num);
				for (unsigned int index = 0; index < order_num; index++)
					par[index] = 0.0;
				parameters.push_back(par);
			}
		}
		void parameter_j(unsigned int _order_index, double _par) {
			if (_order_index < order_num) {
				for (auto par = parameters.begin(); par < parameters.end(); par++)
					(*par)[_order_index] = _par;
			}
		}
		void parameter_ij(unsigned int _variate_index, unsigned int _order_index, double _par) {
			if (_variate_index < variate_num && _order_index < order_num)
				parameters[_variate_index][_order_index] = _par;
		}
		double cal_func(vector<double> x) {
			double val = 0.0;
			for (unsigned int x_index = 0; x_index < variate_num; x_index++)
				for (unsigned int ord_index = 0; ord_index < order_num; ord_index++)
					val += parameters[x_index][ord_index] * pow(x[x_index], ord_index);
			return val;
		}
		double cal_dfunc_dxi(vector<double> x, int xi) {
			double val = 0.0;
			for (unsigned int ord_index = 1; ord_index < order_num; ord_index++)
				val += ord_index * parameters[xi][ord_index] * pow(x[xi], ord_index - 1);
			return val;
		}
		double cal_dfunc_dxi_dxi(vector<double> x, int xi) {
			double val = 0.0;
			for (unsigned int ord_index = 2; ord_index < order_num; ord_index++)
				val += ord_index * (ord_index - 1) * parameters[xi][ord_index] * pow(x[xi], ord_index - 2);
			return val;
		}

		string print_model_info() {
			return "High Order Multivariate Linear Equation : f(xi) = sum_i{ sum_j[ Aij * xi^j ] }, i: component, j: order";
		}
		string print_model_parameter_info() {
			stringstream info;
			info << setiosflags(ios::scientific) << setprecision(2);  // 1.23e+04
			info << "Aij  =  ";
			for (unsigned int order_j = 0; order_j < order_num; order_j++)
				info << "Ai" << order_j << "";
			info << " || i: component, j: order" << endl;
			for (unsigned int comp_i = 0; comp_i < variate_num; comp_i++) {
				info << "A" << comp_i << "j    ";
				for (unsigned int order_j = 0; order_j < order_num; order_j++)
					info << setiosflags(ios::scientific) << setprecision(2) << parameters[comp_i][order_j] << "";
				info << endl;
			}
			return info.str();
		}
		HighorderLinearFunction& operator=(const HighorderLinearFunction& n) {
			parameters = n.parameters;
			variate_num = n.variate_num;
			order_num = n.order_num;
			return *this;
		}
	private:
		vector<vector<double>> parameters;
		unsigned int variate_num;
		unsigned int order_num;
	};
	typedef HighorderLinearFunction HLFunc;
}
