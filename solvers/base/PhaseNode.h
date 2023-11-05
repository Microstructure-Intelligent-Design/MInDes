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
#include"NormalNode.h"
#include"PairValueNode.h"
#include"SublatticNode.h"

using namespace std;
namespace pf {
	enum Flag { pf_BULK, pf_NEAR_INTERFACE, pf_INTERFACE }; //check in calculation
	struct tempNode {
		tempNode() {
			set_zero();
		}
		void set_zero() {
			T = 0;
			D = 0;
			increment = 0;
			laplace = 0;
		}
		tempNode& operator=(const tempNode& n) {
			T = n.T;
			D = n.D;
			increment = n.increment;
			laplace = n.laplace;
			return *this;
		}
		double T;
		double D;
		double increment;
		double laplace;
	};

	class PhaseEntry
	{
	public:
		PhaseEntry() {
			_flag = 0;
			index = 0;
			property = 0;
			phi = 0.0;
			old_phi = 0.0;
			int_increment = 0.0;
			bulk_increment = 0.0;
			laplacian = 0.0;
			phi_grad.set_to_zero();
		};
		~PhaseEntry() {
			x.clear();
			potential.clear();
			kinetics_coeff.clear();
		};
		
		PhaseEntry& operator=(const PhaseEntry& n);
		void set(int _index, int _property, int __flag = 0, double _phi = 0.0) {
			index = _index;
			property = _property;
			_flag = __flag;
			phi = _phi;
		}
		//< flags
		int _flag; //Node_Phase_Flag
		int index; //search for this entry
		int property;
		//< flags
		//< data
		double phi;  // main variable
		double old_phi;
		Vector3 phi_grad;
		double laplacian;
		double int_increment;
		double bulk_increment;

		ConNode x;  // main variable
		SublatticeNode y;
		ChemNode potential;
		PairKinetic kinetics_coeff;

		//< data
	};
	inline PhaseEntry& PhaseEntry::operator=(const PhaseEntry& node) {
		_flag = node._flag;
		index = node.index;
		property = node.property;

		phi = node.phi;
		old_phi = node.old_phi;
		laplacian = node.laplacian;
		phi_grad = node.phi_grad;
		int_increment = node.int_increment;
		bulk_increment = node.bulk_increment;

		x = node.x;
		y = node.y;
		potential = node.potential;
		kinetics_coeff = node.kinetics_coeff;

		return *this;
	}

	class PhaseNode
	{
	public:
		PhaseNode() {
			_Phase.reserve(2);
			Goast_Phase.reserve(2);
		};
		~PhaseNode() {
			clear();
		};
		PhaseEntry& operator[](const int index);                                       //< Index operator for accessing the n's field value
		PhaseNode&  operator=(const PhaseNode& n);                                            //< Assignement operator
		typedef std::vector<PhaseEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<PhaseEntry*>::iterator giterator;                        //< Iterator over storage vector
		typedef std::vector<PhaseEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return _Phase.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return _Phase.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return _Phase.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return _Phase.cend(); };                           //< Constant iterator to the end of storage vector
		void add_phase(int _index, int _property = 0, int flag = 0, double _phi = 0.0) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i)
				if (i->index == _index) {
					i->_flag = flag;
					i->property = _property;
					i->phi = _phi;
					return;
				}
			PhaseEntry entry;
			entry.index = _index;
			entry._flag = flag;
			entry.property = _property;
			entry.phi = _phi;
			_Phase.push_back(entry);
		}
		void add_phase(PhaseEntry _sample) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i)
				if (i->index == _sample.index) {
					i->property = _sample.property;
					i->phi = _sample.phi;
					i->x = _sample.x;
					i->_flag = _sample._flag;
					i->potential = _sample.potential;
					i->kinetics_coeff = _sample.kinetics_coeff;
					return;
				}
			PhaseEntry entry;
			entry.index = _sample.index;
			entry.property = _sample.property;
			entry.phi = _sample.phi;
			entry.x = _sample.x;
			entry._flag = _sample._flag;
			entry.potential = _sample.potential;
			entry.kinetics_coeff = _sample.kinetics_coeff;
			_Phase.push_back(entry);
		}
		bool is_phase_exist(int phaseIndex) {
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				if (phase->index == phaseIndex && phase->phi > Phi_Num_Cut_Off)
					return true;
			return false;
		}
		bool is_phase_inexist(int phaseIndex) {
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				if (phase->index == phaseIndex && phase->phi > Phi_Num_Cut_Off)
					return false;
			return true;
		}
		double cal_phases_fraction_by_index(vector<int> pIndex) {
			double sum_fraction = 0.0;
			for (auto index = pIndex.begin(); index < pIndex.end(); index++)
				sum_fraction += (*this)[*index].phi;
			return sum_fraction;
		}
		double cal_phases_fraction_by_property(vector<int> phasesProperty) {
			double sum_fraction = 0.0;
			for (auto index = phasesProperty.begin(); index < phasesProperty.end(); index++) {
				for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
					if (phase->property == *index)
						sum_fraction += phase->phi;
			}
			return sum_fraction;
		}
		void cal_x_from_phase_x() {
			for (auto con = x.begin(); con < x.end(); con++)
				con->value = 0.0;
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				for (auto con = phase->x.begin(); con < phase->x.end(); con++)
					x[con->index].value += phase->phi * con->value;
		}
		void cal_potential_from_phase_potential() {
			for (auto p = potential.begin(); p < potential.end(); p++) {
				p->value = 0.0;
			}
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				for (auto p = phase->potential.begin(); p < phase->potential.end(); p++) {
					potential[p->index].value += phase->phi * p->value;
				}
		}
		double get_averege_phi(int phaseIndex, int average_range = 1) {
			int num = 0;
			double p_f = 0.0;
			for (int i = -average_range; i <= +average_range; i++)
				for (int j = -average_range; j <= +average_range; j++)
					for (int k = -average_range; k <= +average_range; k++) {
						p_f += get_long_range_node(i, j, k)[phaseIndex].phi;
						num += 1;
					}
			if (num > 0)
				p_f = p_f / num;
			else
				p_f = this->operator[](phaseIndex).phi;
			return p_f;
		}
		void normalized_phi() {
			double sum = 0.0;
			for (auto phase = _Phase.begin(); phase < _Phase.end(); phase++)
				sum += phase->phi;
			for (auto phase = _Phase.begin(); phase < _Phase.end(); phase++)
				phase->phi = phase->phi / sum;
		}
		double cal_customValues_laplace(int customValues_index, double dx = 1.0, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT) {
			double laplace = 0.0;
			if (diff_method == DifferenceMethod::FIVE_POINT) {
				laplace = (get_neighbor_node(Direction::x_down).customValues[customValues_index]
					+ get_neighbor_node(Direction::x_up).customValues[customValues_index]
					+ get_neighbor_node(Direction::y_down).customValues[customValues_index]
					+ get_neighbor_node(Direction::y_up).customValues[customValues_index]
					+ get_neighbor_node(Direction::z_down).customValues[customValues_index]
					+ get_neighbor_node(Direction::z_up).customValues[customValues_index]
					- 6.0 * customValues[customValues_index]) / dx / dx;
			}
			else if (diff_method == DifferenceMethod::NINE_POINT) {
				laplace = (4.0 * get_neighbor_node(Direction::x_down).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::x_up).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::y_down).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::y_up).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::z_down).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::z_up).customValues[customValues_index]
					+ get_long_range_node(-1, -1, 0).customValues[customValues_index]
					+ get_long_range_node(-1, 1, 0).customValues[customValues_index]
					+ get_long_range_node(1, -1, 0).customValues[customValues_index]
					+ get_long_range_node(1, 1, 0).customValues[customValues_index]
					+ get_long_range_node(-1, 0, -1).customValues[customValues_index]
					+ get_long_range_node(-1, 0, 1).customValues[customValues_index]
					+ get_long_range_node(1, 0, -1).customValues[customValues_index]
					+ get_long_range_node(1, 0, 1).customValues[customValues_index]
					+ get_long_range_node(0, -1, -1).customValues[customValues_index]
					+ get_long_range_node(0, -1, 1).customValues[customValues_index]
					+ get_long_range_node(0, 1, -1).customValues[customValues_index]
					+ get_long_range_node(0, 1, 1).customValues[customValues_index]
					- 36.0 * customValues[customValues_index]) / 6.0 / dx / dx;
			}
			else {
				cout << "Difference method define error!";
				SYS_PROGRAM_STOP;
			}
			return laplace;
		}
		Vector3 cal_customValues_gradient(int customValues_index, double dx = 1.0) {
			Vector3 grad;
			grad[0] = (get_neighbor_node(Direction::x_up).customValues[customValues_index] - get_neighbor_node(Direction::x_down).customValues[customValues_index]) / 2.0 / dx;
			grad[1] = (get_neighbor_node(Direction::y_up).customValues[customValues_index] - get_neighbor_node(Direction::y_down).customValues[customValues_index]) / 2.0 / dx;
			grad[2] = (get_neighbor_node(Direction::z_up).customValues[customValues_index] - get_neighbor_node(Direction::z_down).customValues[customValues_index]) / 2.0 / dx;
			return grad;
		}
		Vector3 cal_temperature_gradient(double dx = 1.0) {
			Vector3 grad;
			grad[0] = (get_neighbor_node(Direction::x_up).temperature.T - get_neighbor_node(Direction::x_down).temperature.T) / 2.0 / dx;
			grad[1] = (get_neighbor_node(Direction::y_up).temperature.T - get_neighbor_node(Direction::y_down).temperature.T) / 2.0 / dx;
			grad[2] = (get_neighbor_node(Direction::z_up).temperature.T - get_neighbor_node(Direction::z_down).temperature.T) / 2.0 / dx;
			return grad;
		}
		double cal_laplace_temperature(double dx = 1.0, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT) {
			double laplace = 0.0;
			if (diff_method == DifferenceMethod::FIVE_POINT) {
				laplace = (get_neighbor_node(Direction::x_down).temperature.T
					+ get_neighbor_node(Direction::x_up).temperature.T
					+ get_neighbor_node(Direction::y_down).temperature.T
					+ get_neighbor_node(Direction::y_up).temperature.T
					+ get_neighbor_node(Direction::z_down).temperature.T
					+ get_neighbor_node(Direction::z_up).temperature.T
					- 6.0 * temperature.T) / dx / dx;
			}
			else if (diff_method == DifferenceMethod::NINE_POINT) {
				laplace = (4.0 * get_neighbor_node(Direction::x_down).temperature.T
					+ 4.0 * get_neighbor_node(Direction::x_up).temperature.T
					+ 4.0 * get_neighbor_node(Direction::y_down).temperature.T
					+ 4.0 * get_neighbor_node(Direction::y_up).temperature.T
					+ 4.0 * get_neighbor_node(Direction::z_down).temperature.T
					+ 4.0 * get_neighbor_node(Direction::z_up).temperature.T
					+ get_long_range_node(-1, -1, 0).temperature.T
					+ get_long_range_node(-1, 1, 0).temperature.T
					+ get_long_range_node(1, -1, 0).temperature.T
					+ get_long_range_node(1, 1, 0).temperature.T
					+ get_long_range_node(-1, 0, -1).temperature.T
					+ get_long_range_node(-1, 0, 1).temperature.T
					+ get_long_range_node(1, 0, -1).temperature.T
					+ get_long_range_node(1, 0, 1).temperature.T
					+ get_long_range_node(0, -1, -1).temperature.T
					+ get_long_range_node(0, -1, 1).temperature.T
					+ get_long_range_node(0, 1, -1).temperature.T
					+ get_long_range_node(0, 1, 1).temperature.T
					- 36.0 * temperature.T) / 6.0 / dx / dx;
			}
			else {
				cout << "Difference method define error!";
				SYS_PROGRAM_STOP;
			}
			return laplace;
		}
		void automatic_set_flag();
		//< Data
		std::vector<PhaseEntry> _Phase;
		std::vector<PhaseEntry*> Goast_Phase;
		ConNode x;
		ChemNode potential;
		PairKinetic kinetics_coeff;
		tempNode temperature;

		double_box customValues;
		int_box customFlags;
		vec3_box customVec3s;
		vec6_box customVec6s;
		matrix6x6_box customMatrix6x6s;

		//< Data

		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			for (auto phase = Goast_Phase.begin(); phase < Goast_Phase.end(); phase++)
				(*phase) = nullptr;
			Goast_Phase.clear();
			_Phase.clear(); 
			customValues.clear();
			customVec3s.clear();
			customVec6s.clear();
			customFlags.clear();
			customMatrix6x6s.clear();
			_up_x = nullptr;
			_down_x = nullptr;
			_up_y = nullptr;
			_down_y = nullptr;
			_up_z = nullptr;
			_down_z = nullptr;
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(_Phase.size());
		};
		int phaseNum();
		void push_back(PhaseEntry& n) {
			_Phase.push_back(n);
		}
		void erase(PhaseEntry& n) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
				if (i->index == n.index) {
					_Phase.erase(i);
					return;
				}
			}
		}
		void erase(int phi_index) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
				if (i->index == phi_index) {
					_Phase.erase(i);
					return;
				}
			}
		}
		void erase_until(int phi_index) {
			for (auto i = _Phase.begin(); i < _Phase.end();) {
				if (i->index != phi_index) {
					i = _Phase.erase(i);
				}
				else {
					i++;
				}
			}
		}
		void connect_mesh(int x, int y, int z, PhaseNode& up_x, PhaseNode& down_x, PhaseNode& up_y, PhaseNode& down_y, PhaseNode& up_z, PhaseNode& down_z) {
			_x = x;
			_y = y;
			_z = z;
			_up_x = &up_x;
			_down_x = &down_x;
			_up_y = &up_y;
			_down_y = &down_y;
			_up_z = &up_z;
			_down_z = &down_z;
		}
		PhaseNode& get_neighbor_node(Direction _d);
		PhaseNode& get_long_range_node(int relative_x, int relative_y, int relative_z);
		int _x; // 0 - (Nx-1)
		int _y;
		int _z;
	private:
		PhaseNode* _up_x;
		PhaseNode* _down_x;
		PhaseNode* _up_y;
		PhaseNode* _down_y;
		PhaseNode* _up_z;
		PhaseNode* _down_z;
	};
	inline PhaseEntry& PhaseNode::operator[](const int index) {
		for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
			if (i->index == index) return(*i);
		}
		cout << "PhaseNode error, can't find the PhaseEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline int PhaseNode::phaseNum() {
		int result = 0;
		for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
			if (i->_flag != pf_BULK) result++;
		}
		if (result == 0)  ///< both pf_bulk, in the phase's bulk
			result = 1;
		return result;
	}
	inline PhaseNode& PhaseNode::operator=(const PhaseNode& n) {
		_Phase = n._Phase;
		Goast_Phase = n.Goast_Phase;
		x = n.x;
		potential = n.potential;
		kinetics_coeff = n.kinetics_coeff;
		temperature = n.temperature;
		customValues = n.customValues;
		customFlags = n.customFlags;
		customVec3s = n.customVec3s;
		customVec6s = n.customVec6s;
		customMatrix6x6s = n.customMatrix6x6s;
		return *this;
	}
	inline PhaseNode& PhaseNode::get_neighbor_node(Direction _d) {
		switch (_d)
		{
		case pf::x_up:
			return *_up_x;
			break;
		case pf::x_down:
			return *_down_x;
			break;
		case pf::y_up:
			return *_up_y;
			break;
		case pf::y_down:
			return *_down_y;
			break;
		case pf::z_up:
			return *_up_z;
			break;
		case pf::z_down:
			return *_down_z;
			break;
		default:
			cout << "Class get_neighbor_node parameter error !" << endl;
			SYS_PROGRAM_STOP;
			break;
		}
	}
	inline PhaseNode& PhaseNode::get_long_range_node(int relative_x, int relative_y, int relative_z) {
		if (relative_x > 0)
			return (*this).get_neighbor_node(Direction::x_up).get_long_range_node(relative_x - 1, relative_y, relative_z);
		else if (relative_x < 0)
			return (*this).get_neighbor_node(Direction::x_down).get_long_range_node(relative_x + 1, relative_y, relative_z);
		if (relative_y > 0)
			return (*this).get_neighbor_node(Direction::y_up).get_long_range_node(relative_x, relative_y - 1, relative_z);
		else if (relative_y < 0)
			return (*this).get_neighbor_node(Direction::y_down).get_long_range_node(relative_x, relative_y + 1, relative_z);
		if (relative_z > 0)
			return (*this).get_neighbor_node(Direction::z_up).get_long_range_node(relative_x, relative_y, relative_z - 1);
		else if (relative_z < 0)
			return (*this).get_neighbor_node(Direction::z_down).get_long_range_node(relative_x, relative_y, relative_z + 1);
		return (*this);
	}
	inline void PhaseNode::automatic_set_flag() {
		for (auto phase = _Phase.begin(); phase < _Phase.end(); phase++) {
			phase->_flag = pf_BULK;
			if (phase->phi >= Phi_Num_Cut_Off && phase->phi <= (1.0 - Phi_Num_Cut_Off)) {
				phase->_flag = pf_INTERFACE;
				if (get_neighbor_node(Direction::x_up)[phase->index].phi > (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::x_up)[phase->index].phi < Phi_Num_Cut_Off)
					get_neighbor_node(Direction::x_up)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::x_down)[phase->index].phi > (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::x_down)[phase->index].phi < Phi_Num_Cut_Off)
					get_neighbor_node(Direction::x_down)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::y_up)[phase->index].phi > (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::y_up)[phase->index].phi < Phi_Num_Cut_Off)
					get_neighbor_node(Direction::y_up)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::y_down)[phase->index].phi > (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::y_down)[phase->index].phi < Phi_Num_Cut_Off)
					get_neighbor_node(Direction::y_down)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::z_up)[phase->index].phi > (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::z_up)[phase->index].phi < Phi_Num_Cut_Off)
					get_neighbor_node(Direction::z_up)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::z_down)[phase->index].phi > (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::z_down)[phase->index].phi < Phi_Num_Cut_Off)
					get_neighbor_node(Direction::z_down)[phase->index]._flag = pf_NEAR_INTERFACE;
			}
			if (phase->phi < Phi_Num_Cut_Off) {
				if (get_neighbor_node(Direction::x_up)[phase->index].phi >= Phi_Num_Cut_Off
					|| get_neighbor_node(Direction::x_down)[phase->index].phi >= Phi_Num_Cut_Off
					|| get_neighbor_node(Direction::y_up)[phase->index].phi >= Phi_Num_Cut_Off
					|| get_neighbor_node(Direction::y_down)[phase->index].phi >= Phi_Num_Cut_Off
					|| get_neighbor_node(Direction::z_down)[phase->index].phi >= Phi_Num_Cut_Off
					|| get_neighbor_node(Direction::z_up)[phase->index].phi >= Phi_Num_Cut_Off) {
					phase->_flag = pf_NEAR_INTERFACE;
				}
			}
			else if (phase->phi > (1.0 - Phi_Num_Cut_Off)) {
				if (get_neighbor_node(Direction::x_up)[phase->index].phi <= (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::x_down)[phase->index].phi <= (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::y_up)[phase->index].phi <= (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::y_down)[phase->index].phi <= (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::z_up)[phase->index].phi <= (1.0 - Phi_Num_Cut_Off)
					|| get_neighbor_node(Direction::z_down)[phase->index].phi <= (1.0 - Phi_Num_Cut_Off)) {
					phase->_flag = pf_NEAR_INTERFACE;
				}
			}
		}
	}

}