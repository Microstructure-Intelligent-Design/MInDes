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

namespace pf {
	struct Constinuent {
		int index;
		double value;
		Constinuent() {
			this->index = 0;
			this->value = 0.0;
		};
		Constinuent& operator=(const Constinuent& n) {
			this->index = n.index;
			this->value = n.value;
			return *this;
		}
	};
	/**********************************************************/
	class Sublattice {
	public:
		Sublattice()
		{
			index = 0;
			siteNum = 0.0;
			_isVacancyExist = false;
		}
		~Sublattice() {
			clear();
		}
		Constinuent& operator[](const int n);                                       //< Index operator for accessing the n's field value
		Sublattice&  operator=(const Sublattice& n);                                            //< Assignement operator
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			cons.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(cons.size());
		};
		void push_back(Constinuent& n) {
			cons.push_back(n);
		}
		void add_y(int _index, double _value = 0.0) {
			for (auto i = cons.begin(); i < cons.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			Constinuent con;
			con.index = _index;
			con.value = _value;
			cons.push_back(con);
		}
		void add_y(Constinuent y) {
			for (auto i = cons.begin(); i < cons.end(); ++i)
				if (i->index == y.index) {
					i->value = y.value;
					return;
				}
			Constinuent con;
			con.index = y.index;
			con.value = y.value;
			cons.push_back(con);
		}
		void erase(int index) {
			for (auto j = cons.begin(); j < cons.end(); j++) {
				if (j->index == index) {
					cons.erase(j);
					return;
				}
			}
		}
		void erase(Constinuent& c) {
			for (auto j = cons.begin(); j < cons.end(); j++) {
				if (j->index == c.index) {
					cons.erase(j);
					return;
				}
			}
		}

		typedef std::vector<Constinuent>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Constinuent>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return cons.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return cons.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return cons.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return cons.cend(); };                           //< Constant iterator to the end of storage vector

		double siteNum;
		int index;
		bool _isVacancyExist;
	private:
		std::vector<Constinuent> cons;
	};

	inline Constinuent& Sublattice::operator[](const int n)
	{
		for (auto i = cons.begin(); i < cons.end(); i++) {
			if (i->index == n) {
				return (*i);
			}
		}
		cout << "Sublattice error, don't have the constituent." << endl;
		SYS_PROGRAM_STOP;
	}

	inline Sublattice& Sublattice::operator=(const Sublattice& n)
	{
		cons = n.cons;
		index = n.index;
		siteNum = n.siteNum;
		_isVacancyExist = n._isVacancyExist;
		return *this;
	}
	/**********************************************************/
	class SublatticeNode {
	public:
		SublatticeNode()
		{
		}
		~SublatticeNode()
		{
			clear();
		}
		Sublattice& operator[](const int n);                                       //< Index operator for accessing the n's field value
		SublatticeNode&  operator=(const SublatticeNode& n);                                            //< Assignement operator
		void erase(int index) {
			for (auto j = _subs.begin(); j < _subs.end(); j++)
				if (j->index == index) {
					_subs.erase(j);
					return;
				}
		}
		void erase(Sublattice n) {
			for (auto j = _subs.begin(); j < _subs.end(); j++)
				if (j->index == n.index) {
					_subs.erase(j);
					return;
				}
		}

		void new_sub(int _index, double _siteNum = 1.0, bool isVacancyExist = false) {
			for (auto i = _subs.begin(); i < _subs.end(); ++i)
				if (i->index == _index) {
					i->siteNum = _siteNum;
					i->_isVacancyExist = isVacancyExist;
					return;
				}
			Sublattice sub;
			sub.index = _index;
			sub.siteNum = _siteNum;
			sub._isVacancyExist = isVacancyExist;
			_subs.push_back(sub);
		}
		void new_sub(Sublattice _sub) {
			for (auto i = _subs.begin(); i < _subs.end(); ++i)
				if (i->index == _sub.index) {
					(*i) = _sub;
					return;
				}
			Sublattice sub;
			sub.index = _sub.index;
			sub.siteNum = _sub.siteNum;
			sub._isVacancyExist = _sub._isVacancyExist;
			_subs.push_back(sub);
		}

		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			_subs.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(_subs.size());
		};

		typedef std::vector<Sublattice>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Sublattice>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return _subs.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return _subs.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return _subs.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return _subs.cend(); };                           //< Constant iterator to the end of storage vector

	private:
		std::vector<Sublattice> _subs;
	};

	inline Sublattice& SublatticeNode::operator[](const int n)
	{
		for (auto i = _subs.begin(); i < _subs.end(); i++) {
			if (i->index == n) {
				return (*i);
			}
		}
		cout << "SublatticeNode error, don't have the sublattice." << endl;
		SYS_PROGRAM_STOP;
	}

	inline SublatticeNode& SublatticeNode::operator=(const SublatticeNode& n)
	{
		_subs = n._subs;
		return *this;
	}
}


