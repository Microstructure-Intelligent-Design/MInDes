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
#include "SublatticNode.h"
#include "vectorMatrix.h"
using namespace std;
namespace pf {

	struct XEntry {
		int index;
		double value;
		double value2;
		XEntry() {
			this->index = 0;
			this->value = 0.0;
			this->value2 = 0.0;
		};
		XEntry& operator=(const XEntry& n) {
			this->index = n.index;
			this->value = n.value;
			this->value2 = n.value2;
			return *this;
		}
	};
	struct ChemEntry {
		int index;
		double value;
		double increment;
		double laplacian;
		Vector3 gradient;
		ChemEntry() {
			this->index = 0;
			this->value = 0.0;
			this->laplacian = 0.0;
			this->increment = 0.0;
			this->gradient.set_to_zero();
		};
		ChemEntry& operator=(const ChemEntry& n) {
			this->index = n.index;
			this->value = n.value;
			this->laplacian = n.laplacian;
			this->increment = n.increment;
			this->gradient = n.gradient;
			return *this;
		}
	};
	struct ConEntry
	{
		int index;
		double value;
		double increment;
		double ChemicalReactionFlux;
		double DiffusionFlux;
		double PhaseTransitionFlux;
		ConEntry() {
			this->index = 0;
			this->value = 0.0;
			this->increment = 0.0;
			this->ChemicalReactionFlux = 0.0;
			this->DiffusionFlux = 0.0;
			this->PhaseTransitionFlux = 0.0;
		};
		ConEntry& operator=(const ConEntry& n) {
			this->index = n.index;
			this->value = n.value;
			this->increment = n.increment;
			this->ChemicalReactionFlux = n.ChemicalReactionFlux;
			this->DiffusionFlux = n.DiffusionFlux;
			this->PhaseTransitionFlux = n.PhaseTransitionFlux;
			return *this;
		}
	};
	/**********************************************************/
	class XNode {
	public:
		XNode() {

		}
		~XNode() {
			clear();
		}
		XEntry& operator[](const int index);
		XNode& operator=(const XNode& n);
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			_Node.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(_Node.size());
		};

		void add_x(int _index, double _value = 0.0) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			XEntry entry;
			entry.index = _index;
			entry.value = _value;
			_Node.push_back(entry);
		}
		void add_x(XEntry _sample) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i)
				if (i->index == _sample.index) {
					i->value = _sample.value;
					return;
				}
			XEntry entry;
			entry.index = _sample.index;
			entry.value = _sample.value;
			_Node.push_back(entry);
		}
		void erase(XEntry& n) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i) {
				if (i->index == n.index) _Node.erase(i);
				return;
			}
		}
		void erase(int index) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i) {
				if (i->index == index) _Node.erase(i);
				return;
			}
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<XEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<XEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return _Node.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return _Node.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return _Node.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return _Node.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<XEntry> _Node;
	};
	inline XEntry& XNode::operator[](const int index) {
		for (auto i = _Node.begin(); i < _Node.end(); ++i) {
			if (i->index == index) return(*i);
		}
		cout << "XNode error, can't find the XEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline XNode& XNode::operator=(const XNode& n)
	{
		_Node = n._Node;
		return *this;
	}

	class MechanicNode {
	public:
		//Matrix3x3     VelocityGradient;                               ///< Storage for deformation gradient rate matrix
		vStress       Stresses;                                       ///< Storage for stresses
		vStrain       Strains;                                        ///< Storage for strains
		vStrain       StrainIncrements;                               ///< Storage for strain increments
		vStrain	      VirtualEigenStrains; ///< Storage for virtual eigenstrains
		vStrain       EffectiveEigenStrains;                          ///< Storage for effective eigenstrains
		Matrix6x6     EffectiveElasticConstants;                      ///< Storage for effective elastic constants
		MechanicNode& operator=(const MechanicNode& n);
		MechanicNode() {};
	};
	inline MechanicNode& MechanicNode::operator=(const MechanicNode& n) {
		//VelocityGradient = n.VelocityGradient;
		Stresses = n.Stresses;
		Strains = n.Strains;
		VirtualEigenStrains = n.VirtualEigenStrains;
		StrainIncrements = n.StrainIncrements;
		EffectiveEigenStrains = n.EffectiveEigenStrains;
		EffectiveElasticConstants = n.EffectiveElasticConstants;
		//VirtualEigenStrains = n.VirtualEigenStrains;
	}

	class ConNode {
	public:
		ConNode(){
			
		}
		~ConNode() {
			clear();
		}
		ConEntry& operator[](const int index);
		ConNode&  operator=(const ConNode& n);
		ConNode& operator=(const XNode& n);
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			_Node.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(_Node.size());
		};

		void add_con(int _index, double _value = 0.0) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			ConEntry entry;
			entry.index = _index;
			entry.value = _value;
			_Node.push_back(entry);
		}
		void add_con(ConEntry _sample) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i)
				if (i->index == _sample.index) {
					*i = _sample;
					return;
				}
			ConEntry entry;
			entry = _sample;
			_Node.push_back(entry);
		}
		void erase(ConEntry& n) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i){
				if (i->index == n.index) _Node.erase(i);
				return;
			}
		}
		void erase(int index) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i) {
				if (i->index == index) _Node.erase(i);
				return;
			}
		}
		
		//std::vector<int> interphaseIndex();      

		typedef std::vector<ConEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<ConEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return _Node.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return _Node.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return _Node.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return _Node.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<ConEntry> _Node;
	};
	inline ConEntry& ConNode::operator[](const int index){
		for (auto i = _Node.begin(); i < _Node.end(); ++i) {
			if (i->index == index) return(*i);
		}
		cout << "Node error, can't find the NodeEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline ConNode& ConNode::operator=(const ConNode& n)
	{
		_Node = n._Node;
		return *this;
	}
	inline ConNode& ConNode::operator=(const XNode& n) {
		this->clear();
		for (auto x = n.cbegin(); x < n.cend(); x++)
			this->add_con(x->index, x->value);
		return *this;
	}

	class ChemNode {
	public:
		ChemNode() {

		}
		~ChemNode() {
			clear();
		}
		ChemEntry& operator[](const int index);
		ChemNode& operator=(const ChemNode& n);
		ChemNode& operator=(const XNode& n);
		ChemNode& operator=(const ConNode& n);
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			_Node.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(_Node.size());
		};

		void add_con(int _index, double _value = 0.0) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			ChemEntry entry;
			entry.index = _index;
			entry.value = _value;
			_Node.push_back(entry);
		}
		void add_con(ChemEntry _sample) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i)
				if (i->index == _sample.index) {
					*i = _sample;
					return;
				}
			ChemEntry entry;
			entry = _sample;
			_Node.push_back(entry);
		}
		void erase(ChemEntry& n) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i) {
				if (i->index == n.index) _Node.erase(i);
				return;
			}
		}
		void erase(int index) {
			for (auto i = _Node.begin(); i < _Node.end(); ++i) {
				if (i->index == index) _Node.erase(i);
				return;
			}
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<ChemEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<ChemEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return _Node.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return _Node.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return _Node.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return _Node.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<ChemEntry> _Node;
	};
	inline ChemEntry& ChemNode::operator[](const int index) {
		for (auto i = _Node.begin(); i < _Node.end(); ++i) {
			if (i->index == index) return(*i);
		}
		cout << "ChemNode error, can't find the ChemEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline ChemNode& ChemNode::operator=(const ChemNode& n)
	{
		_Node = n._Node;
		return *this;
	}
	inline ChemNode& ChemNode::operator=(const ConNode& n) {
		this->clear();
		for (auto x = n.cbegin(); x < n.cend(); x++)
			this->add_con(x->index, x->value);
		return *this;
	}
	inline ChemNode& ChemNode::operator=(const XNode& n) {
		this->clear();
		for (auto x = n.cbegin(); x < n.cend(); x++)
			this->add_con(x->index, x->value);
		return *this;
	}
	/**********************************************************/
	struct tensorEntry_int
	{
		int val;
		int index;
		tensorEntry_int() {
			index = 0;
			val = 0;
		}
		tensorEntry_int& operator=(const tensorEntry_int& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	struct tensorEntry_double
	{
		double val;
		int index;
		tensorEntry_double() {
			index = 0;
			val = 0.0;
		}
		tensorEntry_double& operator=(const tensorEntry_double& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	struct tensorEntry_vec3
	{
		Vector3 val;
		int index;
		tensorEntry_vec3() {
			index = 0;
		}
		tensorEntry_vec3& operator=(const tensorEntry_vec3& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	/*struct tensorEntry_vec6
	{
		Vector6 val;
		int index;
		tensorEntry_vec6() {
			index = 0;
		}
		tensorEntry_vec6& operator=(const tensorEntry_vec6& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};*/
	struct tensorEntry_matrix3
	{
		Matrix3x3 val;
		int index;
		tensorEntry_matrix3() {
			index = 0;
		}
		tensorEntry_matrix3& operator=(const tensorEntry_matrix3& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	struct tensorEntry_matrix6
	{
		Matrix6x6 val;
		int index;
		tensorEntry_matrix6() {
			index = 0;
		}
		tensorEntry_matrix6& operator=(const tensorEntry_matrix6& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	struct tensorEntry_strain
	{
		vStrain val;
		int index;
		tensorEntry_strain() {
			index = 0;
		}
		tensorEntry_strain& operator=(const tensorEntry_strain& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	struct tensorEntry_stress
	{
		vStress val;
		int index;
		tensorEntry_stress() {
			index = 0;
		}
		tensorEntry_stress& operator=(const tensorEntry_stress& n) {
			this->index = n.index;
			this->val = n.val;
			return *this;
		}
	};
	/**********************************************************/
	class tensor1_int
	{
	public:
		tensor1_int() {
			index = 0;
		}
		~tensor1_int() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_int(int _index, int _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_int entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		int& operator()(const int index);
		tensor1_int& operator=(const tensor1_int& n);

		typedef std::vector<tensorEntry_int>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_int>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_int> val;
	};
	inline int& tensor1_int::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_int error, can't find the tensorEntry_int, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_int& tensor1_int::operator=(const tensor1_int& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor1_double
	{
	public:
		tensor1_double() {
			index = 0;
		}
		~tensor1_double() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_double(int _index, double _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_double entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		double& operator()(const int index);
		tensor1_double& operator=(const tensor1_double& n);

		typedef std::vector<tensorEntry_double>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_double>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_double> val;
	};
	inline double& tensor1_double::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_double& tensor1_double::operator=(const tensor1_double& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor1_vec3
	{
	public:
		tensor1_vec3() {
			index = 0;
		}
		~tensor1_vec3() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_vec3(int _index, Vector3 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_vec3 entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		Vector3& operator()(const int index);
		tensor1_vec3& operator=(const tensor1_vec3& n);

		typedef std::vector<tensorEntry_vec3>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_vec3>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_vec3> val;
	};
	inline Vector3& tensor1_vec3::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_vec3& tensor1_vec3::operator=(const tensor1_vec3& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor1_matrix3
	{
	public:
		tensor1_matrix3() {
			index = 0;
		}
		~tensor1_matrix3() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_matrix3(int _index, Matrix3x3 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_matrix3 entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		Matrix3x3& operator()(const int index);
		tensor1_matrix3& operator=(const tensor1_matrix3& n);

		typedef std::vector<tensorEntry_matrix3>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_matrix3>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_matrix3> val;
	};
	inline Matrix3x3& tensor1_matrix3::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_matrix3 error, can't find the tensorEntry_matrix3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_matrix3& tensor1_matrix3::operator=(const tensor1_matrix3& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor1_matrix6
	{
	public:
		tensor1_matrix6() {
			index = 0;
		}
		~tensor1_matrix6() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_matrix6(int _index, Matrix6x6 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_matrix6 entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		Matrix6x6& operator()(const int index);
		tensor1_matrix6& operator=(const tensor1_matrix6& n);

		typedef std::vector<tensorEntry_matrix6>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_matrix6>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_matrix6> val;
	};
	inline Matrix6x6& tensor1_matrix6::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_matrix6 error, can't find the tensorEntry_matrix6, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_matrix6& tensor1_matrix6::operator=(const tensor1_matrix6& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor1_strain
	{
	public:
		tensor1_strain() {
			index = 0;
		}
		~tensor1_strain() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_strain(int _index, vStrain _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_strain entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		vStrain& operator()(const int index);
		tensor1_strain& operator=(const tensor1_strain& n);

		typedef std::vector<tensorEntry_strain>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_strain>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_strain> val;
	};
	inline vStrain& tensor1_strain::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_strain error, can't find the tensorEntry_strain, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_strain& tensor1_strain::operator=(const tensor1_strain& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor1_stress
	{
	public:
		tensor1_stress() {
			index = 0;
		}
		~tensor1_stress() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_stress(int _index, vStress _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index) {
					i->val = _value;
					return;
				}
			tensorEntry_stress entry;
			entry.index = _index;
			entry.val = _value;
			val.push_back(entry);
		}
		void erase(int index) {
			for (auto i = val.begin(); i < val.end(); ++i) {
				if (i->index == index) val.erase(i);
				return;
			}
		}
		vStress& operator()(const int index);
		tensor1_stress& operator=(const tensor1_stress& n);

		typedef std::vector<tensorEntry_stress>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensorEntry_stress>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensorEntry_stress> val;
	};
	inline vStress& tensor1_stress::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i) {
			if (i->index == index) return i->val;
		}
		cout << "tensor1_stress error, can't find the tensorEntry_stress, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_stress& tensor1_stress::operator=(const tensor1_stress& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}
	/**********************************************************/
	class tensor2_int
	{
	public:
		tensor2_int() {
			index = 0;
		}
		~tensor2_int() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_int(int _index1, int _index2, int _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_int(_index2, _value);
					return;
				}
			tensor1_int t;
			t.index = _index1;
			t.add_int(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_int t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		int& operator()(const int index1, const int index2);
		tensor1_int& operator()(const int index);
		tensor2_int& operator=(const tensor2_int& n);

		typedef std::vector<tensor1_int>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_int>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_int> val;
	};
	inline int& tensor2_int::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_int& tensor2_int::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor2_int error, can't find the tensor1_int, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_int& tensor2_int::operator=(const tensor2_int& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor2_double
	{
	public:
		tensor2_double() {
			index = 0;
		}
		~tensor2_double() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_double(int _index1, int _index2, double _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_double(_index2, _value);
					return;
				}
			tensor1_double t;
			t.index = _index1;
			t.add_double(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_double t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		double& operator()(const int index1, const int index2);
		tensor1_double& operator()(const int index);
		tensor2_double& operator=(const tensor2_double& n);

		typedef std::vector<tensor1_double>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_double>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_double> val;
	};
	inline double& tensor2_double::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_double& tensor2_double::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor2_vec3 error, can't find the tensor1_double, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_double& tensor2_double::operator=(const tensor2_double& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor2_vec3
	{
	public:
		tensor2_vec3() {
			index = 0;
		}
		~tensor2_vec3() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_vec3(int _index1, int _index2, Vector3 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_vec3(_index2, _value);
					return;
				}
			tensor1_vec3 t;
			t.index = _index1;
			t.add_vec3(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_vec3 t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		Vector3& operator()(const int index1, const int index2);
		tensor1_vec3& operator()(const int index);
		tensor2_vec3& operator=(const tensor2_vec3& n);

		typedef std::vector<tensor1_vec3>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_vec3>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_vec3> val;
	};
	inline Vector3& tensor2_vec3::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_vec3& tensor2_vec3::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor2_vec3 error, can't find the tensor1_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_vec3& tensor2_vec3::operator=(const tensor2_vec3& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor2_matrix3
	{
	public:
		tensor2_matrix3() {
			index = 0;
		}
		~tensor2_matrix3() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_matrix3(int _index1, int _index2, Matrix3x3 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_matrix3(_index2, _value);
					return;
				}
			tensor1_matrix3 t;
			t.index = _index1;
			t.add_matrix3(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_matrix3 t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		Matrix3x3& operator()(const int index1, const int index2);
		tensor1_matrix3& operator()(const int index);
		tensor2_matrix3& operator=(const tensor2_matrix3& n);

		typedef std::vector<tensor1_matrix3>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_matrix3>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_matrix3> val;
	};
	inline Matrix3x3& tensor2_matrix3::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_matrix3 error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_matrix3& tensor2_matrix3::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor2_matrix3 error, can't find the tensor1_matrix3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_matrix3& tensor2_matrix3::operator=(const tensor2_matrix3& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor2_matrix6
	{
	public:
		tensor2_matrix6() {
			index = 0;
		}
		~tensor2_matrix6() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_matrix6(int _index1, int _index2, Matrix6x6 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_matrix6(_index2, _value);
					return;
				}
			tensor1_matrix6 t;
			t.index = _index1;
			t.add_matrix6(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_matrix6 t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		Matrix6x6& operator()(const int index1, const int index2);
		tensor1_matrix6& operator()(const int index);
		tensor2_matrix6& operator=(const tensor2_matrix6& n);

		typedef std::vector<tensor1_matrix6>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_matrix6>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_matrix6> val;
	};
	inline Matrix6x6& tensor2_matrix6::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_matrix6& tensor2_matrix6::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == index) return *(i);
		cout << "tensor2_vec3 error, can't find the tensor1_matrix6, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_matrix6& tensor2_matrix6::operator=(const tensor2_matrix6& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor2_strain
	{
	public:
		tensor2_strain() {
			index = 0;
		}
		~tensor2_strain() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_strain(int _index1, int _index2, vStrain _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_strain(_index2, _value);
					return;
				}
			tensor1_strain t;
			t.index = _index1;
			t.add_strain(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_strain t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		vStrain& operator()(const int index1, const int index2);
		tensor1_strain& operator()(const int index);
		tensor2_strain& operator=(const tensor2_strain& n);

		typedef std::vector<tensor1_strain>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_strain>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_strain> val;
	};
	inline vStrain& tensor2_strain::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_strain& tensor2_strain::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == index) return *(i);
		cout << "tensor2_vec3 error, can't find the tensor1_strain, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_strain& tensor2_strain::operator=(const tensor2_strain& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor2_stress
	{
	public:
		tensor2_stress() {
			index = 0;
		}
		~tensor2_stress() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_stress(int _index1, int _index2, vStress _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_stress(_index2, _value);
					return;
				}
			tensor1_stress t;
			t.index = _index1;
			t.add_stress(_index2, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor1_stress t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) {
					if (i->index == _index1 && j->index == _index2) i->val.erase(j);
					return;
				}
		}
		vStress& operator()(const int index1, const int index2);
		tensor1_stress& operator()(const int index);
		tensor2_stress& operator=(const tensor2_stress& n);

		typedef std::vector<tensor1_stress>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor1_stress>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor1_stress> val;
	};
	inline vStress& tensor2_stress::operator()(const int index1, const int index2) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) {
				if (i->index == index1 && j->index == index2) return j->val;
			}
		cout << "tensor2_vec3 error, can't find the tensorEntry_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor1_stress& tensor2_stress::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == index) return *(i);
		cout << "tensor2_vec3 error, can't find the tensor1_stress, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_stress& tensor2_stress::operator=(const tensor2_stress& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}
	/**********************************************************/
	class tensor3_int
	{
	public:
		tensor3_int() {
			index = 0;
		}
		~tensor3_int() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_int(int _index1, int _index2, int _index3, int _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_int(_index2, _index3, _value);
					return;
				}
			tensor2_int t;
			t.index = _index1;
			t.add_int(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_int t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		int& operator()(const int index1, const int index2, const int index3);
		tensor2_int& operator()(const int index);
		tensor3_int& operator=(const tensor3_int& n);

		typedef std::vector<tensor2_int>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_int>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_int> val;
	};
	inline int& tensor3_int::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_int error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_int& tensor3_int::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_int error, can't find the tensor2_int, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_int& tensor3_int::operator=(const tensor3_int& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor3_double
	{
	public:
		tensor3_double() {
			index = 0;
		}
		~tensor3_double() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_double(int _index1, int _index2, int _index3, double _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_double(_index2, _index3, _value);
					return;
				}
			tensor2_double t;
			t.index = _index1;
			t.add_double(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_double t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		double& operator()(const int index1, const int index2, const int index3);
		tensor2_double& operator()(const int index);
		tensor3_double& operator=(const tensor3_double& n);

		typedef std::vector<tensor2_double>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_double>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_double> val;
	};
	inline double& tensor3_double::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_vec3 error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_double& tensor3_double::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_vec3 error, can't find the tensor2_double, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_double& tensor3_double::operator=(const tensor3_double& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor3_vec3
	{
	public:
		tensor3_vec3() {
			index = 0;
		}
		~tensor3_vec3() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_vec3(int _index1, int _index2, int _index3, Vector3 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_vec3(_index2, _index3, _value);
					return;
				}
			tensor2_vec3 t;
			t.index = _index1;
			t.add_vec3(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_vec3 t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j) 
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		Vector3& operator()(const int index1, const int index2, const int index3);
		tensor2_vec3& operator()(const int index);
		tensor3_vec3& operator=(const tensor3_vec3& n);

		typedef std::vector<tensor2_vec3>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_vec3>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_vec3> val;
	};
	inline Vector3& tensor3_vec3::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j) 
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_vec3 error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_vec3& tensor3_vec3::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_vec3 error, can't find the tensor2_vec3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_vec3& tensor3_vec3::operator=(const tensor3_vec3& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor3_matrix3
	{
	public:
		tensor3_matrix3() {
			index = 0;
		}
		~tensor3_matrix3() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_matrix3(int _index1, int _index2, int _index3, Matrix3x3 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_matrix3(_index2, _index3, _value);
					return;
				}
			tensor2_matrix3 t;
			t.index = _index1;
			t.add_matrix3(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_matrix3 t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		Matrix3x3& operator()(const int index1, const int index2, const int index3);
		tensor2_matrix3& operator()(const int index);
		tensor3_matrix3& operator=(const tensor3_matrix3& n);

		typedef std::vector<tensor2_matrix3>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_matrix3>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_matrix3> val;
	};
	inline Matrix3x3& tensor3_matrix3::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_vec3 error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_matrix3& tensor3_matrix3::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_vec3 error, can't find the tensor2_matrix3, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_matrix3& tensor3_matrix3::operator=(const tensor3_matrix3& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor3_matrix6
	{
	public:
		tensor3_matrix6() {
			index = 0;
		}
		~tensor3_matrix6() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_matrix6(int _index1, int _index2, int _index3, Matrix6x6 _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_matrix6(_index2, _index3, _value);
					return;
				}
			tensor2_matrix6 t;
			t.index = _index1;
			t.add_matrix6(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_matrix6 t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		Matrix6x6& operator()(const int index1, const int index2, const int index3);
		tensor2_matrix6& operator()(const int index);
		tensor3_matrix6& operator=(const tensor3_matrix6& n);

		typedef std::vector<tensor2_matrix6>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_matrix6>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_matrix6> val;
	};
	inline Matrix6x6& tensor3_matrix6::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_vec3 error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_matrix6& tensor3_matrix6::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_vec3 error, can't find the tensor2_matrix6, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_matrix6& tensor3_matrix6::operator=(const tensor3_matrix6& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor3_strain
	{
	public:
		tensor3_strain() {
			index = 0;
		}
		~tensor3_strain() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_strain(int _index1, int _index2, int _index3, vStrain _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_strain(_index2, _index3, _value);
					return;
				}
			tensor2_strain t;
			t.index = _index1;
			t.add_strain(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_strain t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		vStrain& operator()(const int index1, const int index2, const int index3);
		tensor2_strain& operator()(const int index);
		tensor3_strain& operator=(const tensor3_strain& n);

		typedef std::vector<tensor2_strain>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_strain>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_strain> val;
	};
	inline vStrain& tensor3_strain::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_strain error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_strain& tensor3_strain::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_strain error, can't find the tensor2_strain, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_strain& tensor3_strain::operator=(const tensor3_strain& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}

	class tensor3_stress
	{
	public:
		tensor3_stress() {
			index = 0;
		}
		~tensor3_stress() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_stress(int _index1, int _index2, int _index3, vStress _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_stress(_index2, _index3, _value);
					return;
				}
			tensor2_stress t;
			t.index = _index1;
			t.add_stress(_index2, _index3, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor2_stress t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) {
						if (i->index == _index1 && j->index == _index2 && k->index == _index3) j->val.erase(k);
						return;
					}
		}
		vStress& operator()(const int index1, const int index2, const int index3);
		tensor2_stress& operator()(const int index);
		tensor3_stress& operator=(const tensor3_stress& n);

		typedef std::vector<tensor2_stress>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor2_stress>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor2_stress> val;
	};
	inline vStress& tensor3_stress::operator()(const int index1, const int index2, const int index3) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k) {
					if (i->index == index1 && j->index == index2 && k->index == index3) return k->val;
				}
		cout << "tensor3_stress error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor2_stress& tensor3_stress::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_stress error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_stress& tensor3_stress::operator=(const tensor3_stress& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}
	//*********************************************************************
	class tensor4_double
	{
	public:
		tensor4_double() {
			index = 0;
		}
		~tensor4_double() {
			val.clear();
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			val.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(val.size());
		};

		void add_double(int _index1, int _index2, int _index3, int _index4, double _value) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					i->add_double(_index2, _index3, _index4, _value);
					return;
				}
			tensor3_double t;
			t.index = _index1;
			t.add_double(_index2, _index3, _index4, _value);
			val.push_back(t);
		}
		void add_tensor(int _index1) {
			for (auto i = val.begin(); i < val.end(); ++i)
				if (i->index == _index1) {
					return;
				}
			tensor3_double t;
			t.index = _index1;
			val.push_back(t);
		}
		void erase(int _index1, int _index2, int _index3, int _index4) {
			for (auto i = val.begin(); i < val.end(); ++i)
				for (auto j = i->val.begin(); j < i->val.end(); ++j)
					for (auto k = j->val.begin(); k < j->val.end(); ++k) 
						for (auto l = k->val.begin(); l < k->val.end(); ++l) {
							if (i->index == _index1 && j->index == _index2 && k->index == _index3 && l->index == _index4) k->val.erase(l);
							return;
						}
		}
		double& operator()(const int index1, const int index2, const int index3, const int index4);
		tensor3_double& operator()(const int index);
		tensor4_double& operator=(const tensor4_double& n);

		typedef std::vector<tensor3_double>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<tensor3_double>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return val.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return val.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return val.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return val.cend(); };                           //< Constant iterator to the end of storage vector

		int index;
		vector<tensor3_double> val;
	};
	inline double& tensor4_double::operator()(const int index1, const int index2, const int index3, const int index4) {
		for (auto i = val.begin(); i < val.end(); ++i)
			for (auto j = i->val.begin(); j < i->val.end(); ++j)
				for (auto k = j->val.begin(); k < j->val.end(); ++k)
					for (auto l = k->val.begin(); l < k->val.end(); ++l) {
						if (i->index == index1 && j->index == index2 && k->index == index3 && l->index == index4) return l->val;
					}
		cout << "tensor3_vec3 error, can't find the tensorEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor3_double& tensor4_double::operator()(const int index) {
		for (auto i = val.begin(); i < val.end(); ++i)
			if (i->index == index) return *(i);
		cout << "tensor3_vec3 error, can't find the tensor3_double, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline tensor4_double& tensor4_double::operator=(const tensor4_double& n)
	{
		index = n.index;
		val = n.val;
		return *this;
	}
}