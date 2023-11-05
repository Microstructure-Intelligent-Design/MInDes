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
#include "sysTool.h"
using namespace std;
namespace pf{
	enum PairValueProperty{pf_NULL, pf_QUANTITY, pf_VECTORS, pf_STORAGE};
	enum GradientDrection{grad_x, grad_y, grad_z};
	struct valueEntry {
		int pairIndex_1;
		int pairIndex_2;

		double value;

		valueEntry() {
			this->pairIndex_1 = 0;
			this->pairIndex_2 = 0;
			this->value = 0.0;
		};
		valueEntry& operator=(const valueEntry& n) {
			this->pairIndex_1 = n.pairIndex_1;
			this->pairIndex_2 = n.pairIndex_2;
			this->value = n.value;
			return *this;
		}
	};
	struct flagEntry {
		int pairIndex_1;
		int pairIndex_2;

		int flag;

		flagEntry() {
			this->pairIndex_1 = 0;
			this->pairIndex_2 = 0;
			this->flag = 0;
		};
		flagEntry& operator=(const flagEntry& n) {
			this->pairIndex_1 = n.pairIndex_1;
			this->pairIndex_2 = n.pairIndex_2;
			this->flag = n.flag;
			return *this;
		}
	};
	struct kineticEntry {
		int pairIndex_1;
		int pairIndex_2;

		double value; 
		union
		{
			double laplacian;
			double value2;
		};
		Vector3 gradient;

		kineticEntry() {
			this->pairIndex_1 = 0;
			this->pairIndex_2 = 0;
			this->value = 0.0;
			this->value2 = 0.0;
			this->gradient.set_to_zero();
		};
		kineticEntry& operator=(const kineticEntry& n) {
			this->pairIndex_1 = n.pairIndex_1;
			this->pairIndex_2 = n.pairIndex_2;
			this->value = n.value;
			this->value2 = n.value2;
			this->gradient = n.gradient;
			return *this;
		}
	};
	struct vec3Entry {
		int pairIndex_1;
		int pairIndex_2;

		Vector3 value;
		
		vec3Entry() {
			this->pairIndex_1 = 0;
			this->pairIndex_2 = 0;
			this->value[0] = 0.0;
			this->value[1] = 0.0;
			this->value[2] = 0.0;
		};
		vec3Entry& operator=(const vec3Entry& n) {
			this->pairIndex_1 = n.pairIndex_1;
			this->pairIndex_2 = n.pairIndex_2;
			this->value[0] = n.value[0];
			this->value[1] = n.value[1];
			this->value[2] = n.value[2];
			return *this;
		}
		vec3Entry& operator=(double vec[3]) {
			this->value[0] = vec[0];
			this->value[1] = vec[1];
			this->value[2] = vec[2];
			return *this;
		}
	};
	struct doubleBoxEntry {
		int pairIndex_1;
		int pairIndex_2;

		double_box box;

		doubleBoxEntry() {
			this->pairIndex_1 = 0;
			this->pairIndex_2 = 0;
		};
		doubleBoxEntry& operator=(const doubleBoxEntry& n) {
			this->pairIndex_1 = n.pairIndex_1;
			this->pairIndex_2 = n.pairIndex_2;
			this->box = n.box;
			return *this;
		}
		doubleBoxEntry& operator=(double_box box) {
			this->box = box;
			return *this;
		}
	};


	class PairValue
	{
	public:
		PairValue(PairValueProperty n) {
			_valueProperty = n;
		};
		PairValue() {
			_valueProperty = pf_NULL;
		};
		~PairValue() {
			clear();
		};

		typedef std::vector<valueEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<valueEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return pairValue.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return pairValue.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return pairValue.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return pairValue.cend(); };                           //< Constant iterator to the end of storage vector

		PairValue&  operator=(const PairValue& n);
		double operator ()(const int n, const int m);
		void add(const int n, const int m, const double t = 0.0);
		void set(const int n, const int m, const double t = 0.0);
		void set(PairValueProperty _valueProperty);
		std::vector<valueEntry>::iterator erase(const int n, const int m);
		void clear() {
			pairValue.clear();
		}

		PairValueProperty _valueProperty;
		std::vector<valueEntry> pairValue;
	};
	inline PairValue&  PairValue::operator=(const PairValue& n) {
		pairValue = n.pairValue;
		_valueProperty = n._valueProperty;
		return *this;
	}
	inline double PairValue::operator ()(const int n, const int m) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return i->value;
				if(i->pairIndex_1 == m && i->pairIndex_2 == n) return i->value;
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return i->value;
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) return -(i->value);
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return i->value;
			}
		}
		std::cout << "PairValue error, property defined problem." << endl;
		SYS_PROGRAM_STOP;
	}
	inline void PairValue::add(const int n, const int m, const double t) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->value += t;
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					i->value += t;
					return;
				}
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->value += t;
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					i->value -= t;
					return;
				}
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->value += t;
					return;
				}
			}
		}
		valueEntry newPair;
		newPair.pairIndex_1 = n;
		newPair.pairIndex_2 = m;
		newPair.value = t;
		pairValue.push_back(newPair);
	}
	inline void PairValue::set(const int n, const int m, const double t) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->value = t;
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					i->value = t;
					return;
				}
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->value = t;
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					i->value = -t;
					return;
				}
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->value = t;
					return;
				}
			}
		}
		valueEntry newPair;
		newPair.pairIndex_1 = n;
		newPair.pairIndex_2 = m;
		newPair.value = t;
		pairValue.push_back(newPair);
	}
	inline std::vector<valueEntry>::iterator PairValue::erase(const int n, const int m) {
		if(_valueProperty == pf_QUANTITY)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if ((i->pairIndex_1 == n && i->pairIndex_2 == m) || (i->pairIndex_1 == m && i->pairIndex_2 == n)) return pairValue.erase(i);
			}
		else if(_valueProperty == pf_VECTORS)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if ((i->pairIndex_1 == n && i->pairIndex_2 == m) || (i->pairIndex_1 == m && i->pairIndex_2 == n)) return pairValue.erase(i);
			}
		else if (_valueProperty == pf_STORAGE)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return pairValue.erase(i);
			}
		std::cout << "PairValue error, don't find aim value to erase!" << std::endl;
		SYS_PROGRAM_STOP;
	}
	inline void PairValue::set(PairValueProperty _vP) {
		_valueProperty = _vP;
	}

	class PairKinetic
	{
	public:
		PairKinetic() {
			
		};
		~PairKinetic() {
			clear();
		};

		typedef std::vector<kineticEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<kineticEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return pairValue.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return pairValue.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return pairValue.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return pairValue.cend(); };                           //< Constant iterator to the end of storage vector

		PairKinetic& operator=(const PairKinetic& n);
		kineticEntry& operator ()(const int n, const int m);
		void set(const int n, const int m, const double t = 0.0);
		Vector3& get_gradientVec3(const int n, const int m);
		std::vector<kineticEntry>::iterator erase(const int n, const int m);
		void clear() {
			pairValue.clear();
		}

		std::vector<kineticEntry> pairValue;
	};
	inline PairKinetic& PairKinetic::operator=(const PairKinetic& n) {
		pairValue = n.pairValue;
		return *this;
	}
	inline kineticEntry& PairKinetic::operator ()(const int n, const int m) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
			if (i->pairIndex_1 == n && i->pairIndex_2 == m) return (*i);
		}
		std::cout << "PairKinetic error, cont find kineticEntry." << endl;
		SYS_PROGRAM_STOP;
	}
	inline void PairKinetic::set(const int n, const int m, const double t) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
			if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
				i->value = t;
				return;
			}
		}
		kineticEntry newPair;
		newPair.pairIndex_1 = n;
		newPair.pairIndex_2 = m;
		newPair.value = t;
		pairValue.push_back(newPair);
	}
	inline Vector3& PairKinetic::get_gradientVec3(const int n, const int m) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++)
			if (i->pairIndex_1 == n && i->pairIndex_2 == m)
				return i->gradient;
		
		std::cout << "PairKinetic error, cont find kineticEntry." << endl;
		SYS_PROGRAM_STOP;
	}
	inline std::vector<kineticEntry>::iterator PairKinetic::erase(const int n, const int m) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
			if (i->pairIndex_1 == n && i->pairIndex_2 == m) return pairValue.erase(i);
		}
		std::cout << "PairKinetic error, don't find aim kineticEntry to erase!" << std::endl;
		SYS_PROGRAM_STOP;
	}

	class PairVec3
	{
	public:
		PairVec3() {

		};
		~PairVec3() {
			clear();
		};

		typedef std::vector<vec3Entry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<vec3Entry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return pairValue.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return pairValue.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return pairValue.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return pairValue.cend(); };                           //< Constant iterator to the end of storage vector

		PairVec3& operator=(const PairVec3& n);
		Vector3& operator ()(const int n, const int m);
		void set(const int n, const int m, const double t[3]);
		std::vector<vec3Entry>::iterator erase(const int n, const int m);
		void clear() {
			pairValue.clear();
		}

		std::vector<vec3Entry> pairValue;
	};
	inline PairVec3& PairVec3::operator=(const PairVec3& n) {
		pairValue = n.pairValue;
		return *this;
	}
	inline Vector3& PairVec3::operator ()(const int n, const int m) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
			if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
				return i->value;
			}
		}
		std::cout << "PairVec3 error, cont find vec3Entry." << endl;
		SYS_PROGRAM_STOP;
	}
	inline void PairVec3::set(const int n, const int m, const double t[3]) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
			if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
				i->value[0] = t[0];
				i->value[1] = t[1];
				i->value[2] = t[2];
				return;
			}
		}
		vec3Entry newPair;
		newPair.pairIndex_1 = n;
		newPair.pairIndex_2 = m;
		newPair.value[0] = t[0];
		newPair.value[1] = t[1];
		newPair.value[2] = t[2];
		pairValue.push_back(newPair);
	}
	inline std::vector<vec3Entry>::iterator PairVec3::erase(const int n, const int m) {
		for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
			if (i->pairIndex_1 == n && i->pairIndex_2 == m) return pairValue.erase(i);
		}
		std::cout << "PairVec3 error, don't find aim vec3Entry to erase!" << std::endl;
		SYS_PROGRAM_STOP;
	}

	class PairFlag
	{
	public:
		PairFlag() {
			_valueProperty = PairValueProperty::pf_QUANTITY;
		};
		~PairFlag() {
			clear();
		};

		typedef std::vector<flagEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<flagEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return pairValue.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return pairValue.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return pairValue.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return pairValue.cend(); };                           //< Constant iterator to the end of storage vector

		PairFlag& operator=(const PairFlag& n);
		int operator ()(const int n, const int m);
		void set(const int n, const int m, const int flag);
		void set(PairValueProperty _flag) {
			_valueProperty = _flag;
		}
		std::vector<flagEntry>::iterator erase(const int n, const int m);
		void clear() {
			pairValue.clear();
		}

		PairValueProperty _valueProperty;
		std::vector<flagEntry> pairValue;
	};
	inline PairFlag& PairFlag::operator=(const PairFlag& n) {
		pairValue = n.pairValue;
		_valueProperty = n._valueProperty;
		return *this;
	}
	inline int PairFlag::operator ()(const int n, const int m) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return i->flag;
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) return i->flag;
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return i->flag;
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) return -(i->flag);
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return i->flag;
			}
		}
		std::cout << "PairValue error, property defined problem." << endl;
		SYS_PROGRAM_STOP;
	}
	inline void PairFlag::set(const int n, const int m, const int flag) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->flag = flag;
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					i->flag = flag;
					return;
				}
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->flag = flag;
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					i->flag = -flag;
					return;
				}
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					i->flag = flag;
					return;
				}
			}
		}
		flagEntry newPair;
		newPair.pairIndex_1 = n;
		newPair.pairIndex_2 = m;
		newPair.flag = flag;
		pairValue.push_back(newPair);
	}
	inline std::vector<flagEntry>::iterator PairFlag::erase(const int n, const int m) {
		if (_valueProperty == pf_QUANTITY)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if ((i->pairIndex_1 == n && i->pairIndex_2 == m) || (i->pairIndex_1 == m && i->pairIndex_2 == n)) return pairValue.erase(i);
			}
		else if (_valueProperty == pf_VECTORS)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if ((i->pairIndex_1 == n && i->pairIndex_2 == m) || (i->pairIndex_1 == m && i->pairIndex_2 == n)) return pairValue.erase(i);
			}
		else if (_valueProperty == pf_STORAGE)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return pairValue.erase(i);
			}
		std::cout << "PairValue error, don't find aim value to erase!" << std::endl;
		SYS_PROGRAM_STOP;
	}

	class PairBox
	{
	public:
		PairBox() {
			_valueProperty = PairValueProperty::pf_QUANTITY;
		};
		~PairBox() {
			clear();
		};

		typedef std::vector<doubleBoxEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<doubleBoxEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return pairValue.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return pairValue.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return pairValue.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return pairValue.cend(); };                           //< Constant iterator to the end of storage vector

		PairBox& operator=(const PairBox& n);
		double operator ()(const int n, const int m, const int box_index);
		void set(const int n, const int m, const int box_index, const double value = 0.0);
		void set(PairValueProperty _flag) {
			_valueProperty = _flag;
		}
		std::vector<doubleBoxEntry>::iterator erase(const int n, const int m);
		void clear() {
			pairValue.clear();
		}

		PairValueProperty _valueProperty;
		std::vector<doubleBoxEntry> pairValue;
	};
	inline PairBox& PairBox::operator=(const PairBox& n) {
		pairValue = n.pairValue;
		_valueProperty = n._valueProperty;
		return *this;
	}
	inline double PairBox::operator ()(const int n, const int m, const int box_index) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m)
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) return b->value;
				if (i->pairIndex_1 == m && i->pairIndex_2 == n)
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) return b->value;
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m)
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) return b->value;
				if (i->pairIndex_1 == m && i->pairIndex_2 == n)
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) return -b->value;
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m)
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) return b->value;
			}
		}
		std::cout << "PairBox error, property defined problem." << endl;
		SYS_PROGRAM_STOP;
	}
	inline void PairBox::set(const int n, const int m, const int box_index, const double value) {
		if (_valueProperty == pf_QUANTITY) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) {
							b->value = value;
							return;
						}
					i->box.add_double(box_index, value);
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) {
							b->value = value;
							return;
						}
					i->box.add_double(box_index, value);
					return;
				}
			}
		}
		else if (_valueProperty == pf_VECTORS) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) {
							b->value = value;
							return;
						}
					i->box.add_double(box_index, value);
					return;
				}
				if (i->pairIndex_1 == m && i->pairIndex_2 == n) {
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) {
							b->value = -value;
							return;
						}
					i->box.add_double(box_index, -value);
					return;
				}
			}
		}
		else if (_valueProperty == pf_STORAGE) {
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) {
					for (auto b = i->box.begin(); b < i->box.end(); b++)
						if (b->index == box_index) {
							b->value = value;
							return;
						}
					i->box.add_double(box_index, value);
					return;
				}
			}
		}
		doubleBoxEntry newPair;
		newPair.pairIndex_1 = n;
		newPair.pairIndex_2 = m;
		newPair.box.add_double(box_index, value);
		pairValue.push_back(newPair);
		return;
	}
	inline std::vector<doubleBoxEntry>::iterator PairBox::erase(const int n, const int m) {
		if (_valueProperty == pf_QUANTITY)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if ((i->pairIndex_1 == n && i->pairIndex_2 == m) || (i->pairIndex_1 == m && i->pairIndex_2 == n)) return pairValue.erase(i);
			}
		else if (_valueProperty == pf_VECTORS)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if ((i->pairIndex_1 == n && i->pairIndex_2 == m) || (i->pairIndex_1 == m && i->pairIndex_2 == n)) return pairValue.erase(i);
			}
		else if (_valueProperty == pf_STORAGE)
			for (auto i = pairValue.begin(); i < pairValue.end(); i++) {
				if (i->pairIndex_1 == n && i->pairIndex_2 == m) return pairValue.erase(i);
			}
		std::cout << "PairBox error, don't find aim value to erase!" << std::endl;
		SYS_PROGRAM_STOP;
	}

}