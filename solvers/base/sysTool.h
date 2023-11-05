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
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <omp.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <initializer_list>
#include <complex>
#include <stdio.h>
#include <float.h>
#include <thread>
#include <chrono>
#ifdef _WIN32
#define SYS_PROGRAM_STOP while(true){getchar();};
//#define SYS_PROGRAM_STOP std::exit(0);
#define IS_NAN(a) _isnan(a)
#include <tchar.h>
#include <shlobj.h>
#include <direct.h>
#include <Windows.h>
#include <io.h>
#include <conio.h>
#include <psapi.h>
#include <process.h>
const std::string dirSeparator = "\\";                                      //< Windows style directory separator
#else
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <unistd.h>
#include <cpuid.h>
#include <sys/types.h>
#define SYS_PROGRAM_STOP std::exit(0);
#define IS_NAN(a) __isnan(a)
const std::string dirSeparator = "/";                                       //< Unix/Linux style directory separator
#endif;

#define SYS_EPSILON 1e-6
#define Tangent_Max_Number double(1e10)
#define Is_Equality(a, b) pf::isTwoNumEquality(a, b)
#define PI double(3.1415926535897932)
#define ElementaryChargeQuantity double(1.6021892e-19)   // C
#define AvogadroConstant double(6.02214076e23)  // mol -1
#define NA double(6.02214076e23)  // mol -1
#define FaradayConstant 96485.333  // C/mol
#define VacuumPermeability_Magnetismus (PI*4e-7) // 
#define AngleToRadians(angle) double(angle/180.0*PI)
#define NaN double(-1e99)

#define SOLVENT_NONE int(-1)

#define RAND_INIT_ALL (srand((unsigned int)time(NULL)))
#define RAND_INIT(seed) (srand((unsigned int)(seed)))
#define RAND_0_1 (((double)std::rand())/RAND_MAX)

#define INFILE_TAIL string(".mindes")

static double phi_numerical_threshold = 1e-3;
#define Phi_Num_Cut_Off phi_numerical_threshold

static void change_phi_numerical_threshold(double numerical_thresold = 1e-3) {
	phi_numerical_threshold = numerical_thresold;
}

namespace pf {
	using namespace std;
	enum DifferenceMethod { FIVE_POINT, NINE_POINT };
	enum NumericalLimit { LowerLimit, UpperLimit };
	enum Boundary { UP_X, UP_Y, UP_Z, DOWN_X, DOWN_Y, DOWN_Z };
	enum BoundaryCondition { FIXED, PERIODIC, ADIABATIC };
	enum Dimension { One_Dimension, Two_Dimension, Three_Dimension };
	enum Direction { x_up, x_down, y_up, y_down, z_up, z_down, num_of_directons };
	enum Axis { AXIS_X, AXIS_Y, AXIS_Z };

	enum PhiEquationType { PEType_Const, PEType_AC_Standard, PEType_AC_Pairwise, PEType_CH_Standard };
	enum ConEquationType { CEType_Const, CEType_TotalX, CEType_PhaseX, CEType_GrandP };
	enum ConEquationDomain { CEDomain_Standard, CEDomain_Reverse };
	enum TemperatureEquationType { TType_Const, TType_Standard };

	enum ExternalFields {
		RELAX_interface_buff = -2000,
		LBM_Symbols_INDEX_0 = -1000,
		EFF_NONE = -100,
		CON_Smooth_Phi, CON_Smooth_Old_Phi,
		MECH_stress, MECH_strain, MECH_eigen_strain, MECH_stiffness, MECH_plastic_strain, MECH_ave_plastic_strain,
		FLUID_pressure, FLUID_Fluid_Domain, FLUID_velocity, FLUID_volume_force, FLUID_mass, FLUID_viscosity,
	};
	enum LBM_Symbols
	  { LBM_f_0, LBM_f_1, LBM_f_2, LBM_f_3, LBM_f_4, LBM_f_5, LBM_f_6, LBM_f_7, LBM_f_8, LBM_f_9, LBM_f_10, LBM_f_11, LBM_f_12, LBM_f_13, LBM_f_14, LBM_f_15, LBM_f_16, LBM_f_17, LBM_f_18, 
		LBM_m_0, LBM_m_1, LBM_m_2, LBM_m_3, LBM_m_4, LBM_m_5, LBM_m_6, LBM_m_7, LBM_m_8, LBM_m_9, LBM_m_10, LBM_m_11, LBM_m_12, LBM_m_13, LBM_m_14, LBM_m_15, LBM_m_16, LBM_m_17, LBM_m_18,
		LBM_f_macro, LBM_fv_macro, LBM_SIZE };
	enum StraggeredGridsFlag { Minus_0_5, Plus_0_5 };

	static bool isTwoNumEquality(double a, double b) {
		if (std::fabs(a - b) < SYS_EPSILON)
			return true;
		else
			return false;
	}

	static int double_to_int(double a) {
		if ((a - int(a)) > 0.5)
			return int(a) + 1;
		else
			return int(a);
	}

	struct Point {
		double x;
		double y;
		double z;
		Point() {
			x = 0;
			y = 0;
			z = 0;
		}
		Point(double _x, double _y, double _z) {
			x = _x;
			y = _y;
			z = _z;
		}
		void set(double _x, double _y, double _z) {
			x = _x;
			y = _y;
			z = _z;
		}
		void do_boundary(BoundaryCondition x_up_bc, BoundaryCondition y_up_bc, BoundaryCondition z_up_bc, BoundaryCondition x_down_bc, BoundaryCondition y_down_bc, BoundaryCondition z_down_bc, int x_limit, int y_limit, int z_limit) {
			if (x_down_bc == BoundaryCondition::ADIABATIC && x < 0)
				x = 0;
			else if (x_up_bc == BoundaryCondition::ADIABATIC && x >= x_limit)
				x = x_limit - 1;
			else if (x_down_bc == BoundaryCondition::FIXED && x < 0)
				x = 0;
			else if (x_up_bc == BoundaryCondition::FIXED && x >= x_limit)
				x = x_limit - 1;
			else if (x_down_bc == BoundaryCondition::PERIODIC && x < 0)
				x += x_limit;
			else if (x_up_bc == BoundaryCondition::PERIODIC && x >= x_limit)
				x -= x_limit;
			else if (y_down_bc == BoundaryCondition::ADIABATIC && y < 0)
				y = 0;
			else if (y_up_bc == BoundaryCondition::ADIABATIC && y >= y_limit)
				y = y_limit - 1;
			else if (y_down_bc == BoundaryCondition::FIXED && y < 0)
				y = 0;
			else if (y_up_bc == BoundaryCondition::FIXED && y >= y_limit)
				y = y_limit - 1;
			else if (y_down_bc == BoundaryCondition::PERIODIC && y < 0)
				y += y_limit;
			else if (y_up_bc == BoundaryCondition::PERIODIC && y >= y_limit)
				y -= y_limit;
			else if (z_down_bc == BoundaryCondition::ADIABATIC && z < 0)
				z = 0;
			else if (z_up_bc == BoundaryCondition::ADIABATIC && z >= z_limit)
				z = z_limit - 1;
			else if (z_down_bc == BoundaryCondition::FIXED && z < 0)
				z = 0;
			else if (z_up_bc == BoundaryCondition::FIXED && z >= z_limit)
				z = z_limit - 1;
			else if (z_down_bc == BoundaryCondition::PERIODIC && z < 0)
				z += z_limit;
			else if (z_up_bc == BoundaryCondition::PERIODIC && z >= z_limit)
				z -= z_limit;
			else
				return;
			return do_boundary(x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc, x_limit, y_limit, z_limit);
		}
		Point& operator=(const Point& n) {
			this->x = n.x;
			this->y = n.y;
			this->z = n.z;
			return *this;
		}
		Point operator+(const Point& n) {
			Point re;
			re.x = this->x + n.x;
			re.y = this->y + n.y;
			re.z = this->z + n.z;
			return re;
		}
		Point operator-(const Point& n) {
			Point re;
			re.x = this->x - n.x;
			re.y = this->y - n.y;
			re.z = this->z - n.z;
			return re;
		}
		Point operator*(const double n) {
			Point re;
			re.x = this->x * n;
			re.y = this->y * n;
			re.z = this->z * n;
			return re;
		}
		Point operator/(const double n) {
			Point re;
			re.x = this->x / n;
			re.y = this->y / n;
			re.z = this->z / n;
			return re;
		}
	};

	struct Point2D {
		double x;
		double y;
		Point2D() {
			x = 0;
			y = 0;
		}
		Point2D(double _x, double _y) {
			x = _x;
			y = _y;
		}
		void set(double _x, double _y) {
			x = _x;
			y = _y;
		}
		Point2D& operator=(const Point2D& n) {
			this->x = n.x;
			this->y = n.y;
			return *this;
		}
	};

	struct LinearCompoundModel {
		double m;
		double x;
		double scale_max;
		void set_parameters(double _m, double _x) {
			m = _m;
			x = _x;
		}
		LinearCompoundModel() {
			m = 1.0;
			x = 1.0;
			scale_max = 0.0;
		}
	};
	struct bool_elem {
		int index;
		bool value;
		bool_elem() {
			index = 0;
			value = 0;
		}
		bool_elem& operator=(const bool_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct bool_box {
		vector<bool_elem> _bool_box;
		typedef std::vector<bool_elem>::iterator iterator;
		typedef std::vector<bool_elem>::const_iterator citerator;
		iterator  begin() { return _bool_box.begin(); };
		iterator  end() { return _bool_box.end(); };
		bool& operator[](const int index) {
			for (auto i = _bool_box.begin(); i < _bool_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "_bool_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		bool_box& operator=(const bool_box& n) {
			_bool_box = n._bool_box;
			return *this;
		}
		void add_bool(int _index, bool _value) {
			for (auto i = _bool_box.begin(); i < _bool_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			bool_elem elem;
			elem.index = _index;
			elem.value = _value;
			_bool_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _bool_box.begin(); i < _bool_box.end();) {
				if (i->index == index) {
					i = _bool_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_bool_box.clear();
		}
		int size() {
			return int(_bool_box.size());
		}
	};
	struct int_elem {
		int index;
		int value;
		int_elem() {
			index = 0;
			value = 0;
		}
		int_elem& operator=(const int_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct int_box {
		vector<int_elem> _int_box;
		int index;
		typedef std::vector<int_elem>::iterator iterator;
		typedef std::vector<int_elem>::const_iterator citerator;
		iterator  begin() { return _int_box.begin(); };
		iterator  end() { return _int_box.end(); };
		int& operator[](const int index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "int_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		int_box& operator=(const int_box& n) {
			_int_box = n._int_box;
			index = n.index;
			return *this;
		}
		void add_int(int _index, int _value) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			int_elem elem;
			elem.index = _index;
			elem.value = _value;
			_int_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _int_box.begin(); i < _int_box.end();) {
				if (i->index == index) {
					i = _int_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_int_box.clear();
		}
		int size() {
			return int(_int_box.size());
		}
		int_box() {
			index = 0;
		}
	};
	struct int2_box {
		vector<int_box> _int_box;
		typedef std::vector<int_box>::iterator iterator;
		typedef std::vector<int_box>::const_iterator citerator;
		iterator  begin() { return _int_box.begin(); };
		iterator  end() { return _int_box.end(); };
		int_box& operator[](const int index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i) {
				if (i->index == index) return (*i);
			}
			cout << "int2_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		int2_box& operator=(const int2_box& n) {
			_int_box = n._int_box;
			return *this;
		}
		void add_int(int _index1, int _index2, int _value) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index1) {
					i->add_int(_index2, _value);
					return;
				}
			int_box box;
			box.index = _index1;
			box.add_int(_index2, _value);
			_int_box.push_back(box);
		}
		void add_int(int _index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index) {
					return;
				}
			int_box box;
			box.index = _index;
			_int_box.push_back(box);
		}
		void erase(int index) {
			for (auto i = _int_box.begin(); i < _int_box.end();) {
				if (i->index == index) {
					i = _int_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_int_box.clear();
		}
		int size() {
			return int(_int_box.size());
		}
	};
	struct double_elem {
		int index;
		double value;
		double_elem() {
			index = 0;
			value = 0.0;
		}
		double_elem& operator=(const double_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct double_box {
		vector<double_elem> _double_box;
		typedef std::vector<double_elem>::iterator iterator;
		typedef std::vector<double_elem>::const_iterator citerator;
		iterator  begin() { return _double_box.begin(); };
		iterator  end() { return _double_box.end(); };
		double& operator[](const int index) {
			for (auto i = _double_box.begin(); i < _double_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "double_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		double_box& operator=(const double_box& n) {
			_double_box = n._double_box;
			return *this;
		}
		void add_double(int _index, double _value) {
			for (auto i = _double_box.begin(); i < _double_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			double_elem elem;
			elem.index = _index;
			elem.value = _value;
			_double_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _double_box.begin(); i < _double_box.end();) {
				if (i->index == index) {
					i = _double_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_double_box.clear();
		}
		int size() {
			return int(_double_box.size());
		}
	};
	struct pair_flag {
		int index1;
		int index2;
		int value;
		pair_flag() {
			index1 = 0;
			index2 = 0;
			value = 0;
		}
		pair_flag& operator=(const pair_flag& n) {
			index1 = n.index1;
			index2 = n.index2;
			value = n.value;
			return *this;
		}
	};
	struct pair_flag_box {
		vector<pair_flag> _flag_box;
		typedef std::vector<pair_flag>::iterator iterator;
		typedef std::vector<pair_flag>::const_iterator citerator;
		iterator  begin() { return _flag_box.begin(); };
		iterator  end() { return _flag_box.end(); };
		int& operator()(const int index1, const int index2) {
			for (auto i = _flag_box.begin(); i < _flag_box.end(); ++i) {
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) return i->value;
			}
			cout << "int_box error, can't find the value indexs : " << index1 << ", " << index2 << endl;
			SYS_PROGRAM_STOP;
		}
		pair_flag_box& operator=(const pair_flag_box& n) {
			_flag_box = n._flag_box;
		}
		void set_flag(int index1, int index2, int flag) {
			for (auto i = _flag_box.begin(); i < _flag_box.end(); ++i)
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) {
					i->value = flag;
					return;
				}
			pair_flag elem;
			elem.index1 = index1;
			elem.index2 = index2;
			elem.value = flag;
			_flag_box.push_back(elem);
		}
		void erase(int index1, int index2) {
			for (auto i = _flag_box.begin(); i < _flag_box.end();) {
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) {
					i = _flag_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_flag_box.clear();
		}
		int size() {
			return int(_flag_box.size());
		}
	};
	struct string_elem {
		int index;
		string value;
		string_elem() {
			index = 0;
			value = "";
		}
		string_elem& operator=(const string_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct string_box {
		vector<string_elem> _string_box;
		typedef std::vector<string_elem>::iterator iterator;
		typedef std::vector<string_elem>::const_iterator citerator;
		iterator  begin() { return _string_box.begin(); };
		iterator  end() { return _string_box.end(); };
		string& operator[](const int index) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "string_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		string_box& operator=(const string_box& n) {
			_string_box = n._string_box;
			return *this;
		}
		void add_string(int _index, string _value) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			string_elem elem;
			elem.index = _index;
			elem.value = _value;
			_string_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _string_box.begin(); i < _string_box.end();) {
				if (i->index == index) {
					i = _string_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_string_box.clear();
		}
		int size() {
			return int(_string_box.size());
		}
	};

	// 组分 及 相 信息
	struct Info_NodeEntry {
		int index;
		double value;
		string name;
		Info_NodeEntry() {
			index = 0;
			value = 0.0;
			name = "";
		}
		Info_NodeEntry& operator=(const Info_NodeEntry& ne) {
			index = ne.index;
			value = ne.value;
			name = ne.name;
			return *this;
		}
	};
	class Info_Node {
	public:
		Info_Node() {
			siteNum = 0.0;
			index = 0;
			name = "";
		}
		Info_NodeEntry& operator[](int index) {
			for (auto i = node.begin(); i < node.end(); ++i) {
				if (i->index == index) return(*i);
			}
			cout << "Info_Node error, can't find the Info_NodeEntry, index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		Info_NodeEntry& operator[](string name) {
			for (auto i = node.begin(); i < node.end(); ++i) {
				if (i->name.compare(name) == 0) return(*i);
			}
			cout << "Info_Node error, can't find the Info_NodeEntry, name : " << name << endl;
			SYS_PROGRAM_STOP;
		}
		Info_Node& operator=(const Info_Node& n) {
			node = n.node;
			index = n.index;
			name = n.name;
			siteNum = n.siteNum;
			return *this;
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			node.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(node.size());
		};

		void add_nodeEntry(int _index, string _name = "", double _value = 0.0) {
			for (auto i = node.begin(); i < node.end(); ++i)
				if (i->index == _index) {
					i->name = _name;
					i->value = _value;
					return;
				}
			Info_NodeEntry Ne;
			Ne.index = _index;
			Ne.name = _name;
			Ne.value = _value;
			node.push_back(Ne);
		}
		void erase(Info_NodeEntry& n) {
			for (auto i = node.begin(); i < node.end(); ++i) {
				if (i->index == n.index) node.erase(i);
				return;
			}
			cout << "Info_Node error, don't have aim NodeEntry to erase, index : " << n.index << endl;
			SYS_PROGRAM_STOP;
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<Info_NodeEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Info_NodeEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return node.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return node.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return node.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return node.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<Info_NodeEntry> node;
		int index;
		string name;
		double siteNum;
	};
	class Info_Phase {
	public:
		int phi_property;
		string phi_name;
		Info_Node x;
		Info_Phase() {
			phi_property = 0;
			phi_name = "";
		}
		Info_Phase& operator=(const Info_Phase& p) {
			phi_property = p.phi_property;
			phi_name = p.phi_name;
			return *this;
		}
	};
	class Info_Phases {
	public:
		Info_Phases() {
		}
		Info_Phase& operator[](int property) {
			for (auto i = phases.begin(); i < phases.end(); ++i) {
				if (i->phi_property == property) return(*i);
			}
			cout << "Info_Phases error, can't find the Info_Phase, property : " << property << endl;
			SYS_PROGRAM_STOP;
		}
		Info_Phase& operator[](string name) {
			for (auto i = phases.begin(); i < phases.end(); ++i) {
				if (i->phi_name.compare(name) == 0) return(*i);
			}
			cout << "Info_Phases error, can't find the Info_Phase, name : " << name << endl;
			SYS_PROGRAM_STOP;
		}
		Info_Phases& operator=(const Info_Phases& n) {
			phases = n.phases;
			return *this;
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			phases.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(phases.size());
		};

		void add_Phase(int _property, string _name) {
			for (auto i = phases.begin(); i < phases.end(); ++i)
				if (i->phi_property == _property) {
					i->phi_name = _name;
					return;
				}
			Info_Phase phase;
			phase.phi_property = _property;
			phase.phi_name = _name;
			phases.push_back(phase);
		}
		void add_Phase(Info_Phase phase) {
			phases.push_back(phase);
		}
		void erase(Info_Phase& n) {
			for (auto i = phases.begin(); i < phases.end(); ++i) {
				if (i->phi_property == n.phi_property) phases.erase(i);
				return;
			}
			cout << "Info_Phases error, don't have aim phase to erase, property : " << n.phi_property << endl;
			SYS_PROGRAM_STOP;
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<Info_Phase>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Info_Phase>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return phases.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return phases.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return phases.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return phases.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<Info_Phase> phases;
	};

	static string erase_tail_of_input_file_name(string input_file_name) {
		string tail = INFILE_TAIL, name_without_tail = input_file_name;
		int index = 0;
		bool is_name_correct = true;
		while (name_without_tail.size() != 0 && index < tail.size())
		{
			index++;
			if (*(name_without_tail.end() - 1) != *(tail.end() - index)) {
				is_name_correct = false;
				break;
			}
			name_without_tail.erase(name_without_tail.end() - 1);
		}
		if (name_without_tail.size() == 0)
			is_name_correct = false;
		if (!is_name_correct) {
			cout << "> input file name error, file name = " << input_file_name << ", aim tail is " << tail << endl;
			SYS_PROGRAM_STOP;
		}
		return name_without_tail;
	}

#ifdef _WIN32
	static std::string TCHAR2STRING(TCHAR* STR)
	{
		int iLen = WideCharToMultiByte(CP_ACP, 0, STR, -1, NULL, 0, NULL, NULL);
		char* chRtn = new char[iLen * sizeof(char)];
		WideCharToMultiByte(CP_ACP, 0, STR, -1, chRtn, iLen, NULL, NULL);
		std::string str(chRtn);
		delete chRtn;
		return str;
	}
	static bool SelectFilePath(std::string& folderName)
	{
		TCHAR szBuffer[MAX_PATH] = { 0 };
		BROWSEINFO bi;
		ZeroMemory(&bi, sizeof(BROWSEINFO));
		bi.hwndOwner = GetForegroundWindow();
		bi.pszDisplayName = szBuffer;
		bi.pidlRoot = NULL;
		bi.lpszTitle = _T("从下面选文件:");
		bi.ulFlags = BIF_BROWSEINCLUDEFILES;;
		LPITEMIDLIST idl = SHBrowseForFolder(&bi);
		if (NULL == idl){
			return false;
		}
		SHGetPathFromIDList(idl, szBuffer);
		folderName = TCHAR2STRING(szBuffer);
		return true;
		
	}
#endif;
	static string GetFolderOfPath(std::string file_path) {
		while (file_path.size() != 0)
		{
			char _c = *(file_path.end() - 1);
			if (_c == '/' || _c == '\\') {
				file_path.erase(file_path.end() - 1);
				break;
			}
			else {
				file_path.erase(file_path.end() - 1);
			}
		}
		return file_path;
	}
	static char get_char_not_show() {
		char c;
#ifdef _WIN32
		c = _getch();
#else
		system("stty -echo");
		c = getchar();
		system("stty echo");
#endif
		return c;
	}
	static string GetFileNameOfPath(std::string file_path) {
		string name = "";
		int index = 1;
		while (index <= file_path.size())
		{
			char _c = *(file_path.end() - index);
			if (_c == '/' || _c == '\\') {
				break;
			}
			else {
				index++;
				name = _c + name;
			}
		}
		return name;
	}
	static string get_string_from_consol(bool is_show = true, char replace_char = '*') {
		string str = "";
		int i = 0;
		if (is_show) {
			while (true) {
				char ch = getchar();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
			}
		}
		else {
			while (true) {
				char ch = get_char_not_show();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
				putchar(replace_char);
			}
		}
		putchar('\n');
		return str;
	}

	//printf("\x1b[%d;%dm%s\x1b[%dm", backcolor, frountcolor, str, control);
	/*第一个% d:backcolor表示显示字符串的背景颜色, 其值如下：
	40 : 黑  
	41 : 红  
	42 : 绿  
	43 : 黄  
	44 : 蓝  
	45 : 紫  
	46 : 深绿
	47 : 白色
	第二个% d : frountcolor表示字体颜色, 其值如下：
	30 : 黑
	31 : 红
	32 : 绿
	33 : 黄
	34 : 蓝
	35 : 紫
	36 : 深绿
	37 : 白色

	第三个 % s : str 表示需要显示的字符串

	第四个 % d : control表示ANSI控制码, 其值如下表所示：

	ANSI控制码:

	QUOTE:
	\x1b[0m     关闭所有属性
	\x1b[1m     设置高亮度
	\x1b[4m     下划线
	\x1b[5m     闪烁
	\x1b[7m     反显
	\x1b[8m     消隐
	\x1b[30m--  \x1b[37m   设置前景色
	\x1b[40m--  \x1b[47m   设置背景色
	\x1b[nA    光标上移n行
	\x1b[nB    光标下移n行
	\x1b[nC    光标右移n行
	\x1b[nD    光标左移n行
	\x1b[y; xH  设置光标位置
	\x1b[2J    清屏
	\x1b[K     清除从光标到行尾的内容
	\x1b[s     保存光标位置
	\x1b[u     恢复光标位置
	\x1b[? 25l  隐藏光标
	\x1b[? 25h  显示光标
	例：
	printf("\x1b[%d;%dmhello world\n\x1b[0m",i, j);
	*/
	static void printf_color_on_control(string str, int front_color = 30, int back_color = 43) {
		printf("\x1b[%d;%dm%s\x1b[0m", back_color, front_color, str.c_str());
	}
}