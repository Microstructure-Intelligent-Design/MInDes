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
#include "NormalNode.h"
#include "PairValueNode.h"
#include "PhaseNode.h"
#include "SublatticNode.h"
#define X_MIN 0
#define Y_MIN 0
#define Z_MIN 0

using namespace std;
namespace pf {
	static int Index(int x, int y, int z, int limit_x, int limit_y, int limit_z) {
		return x + y * limit_x + z * limit_x * limit_y;
	}
	class FieldStorage_forPhaseNode
	{
	public:
		FieldStorage_forPhaseNode() {};
		~FieldStorage_forPhaseNode() {
			try {
				_mesh.clear();
			}
			catch(...){
				std::cout << "Destroy FieldStorage_forPhaseNode failed, problems occur." << endl;
			}
		};

		PhaseNode& operator()(const int x, const int y, const int z) ;                                       //< Index operator for accessing the n's field value
		FieldStorage_forPhaseNode&  operator=(const FieldStorage_forPhaseNode& n);

		Flag currentFlag(PhaseNode& node, int phaseIndex);      //<  seek where the node is (set in storage)
		bool is_interface(PhaseNode& node, int phaseIndex);
		void upgrade(PhaseNode& node, int phaseIndex);
		void downgrade(PhaseNode& node, int phaseIndex);
		void init(int Nx, int Ny, int Nz, double _dr, BoundaryCondition bc_x_up, BoundaryCondition bc_y_up, 
			BoundaryCondition bc_z_up, BoundaryCondition bc_x_down, BoundaryCondition bc_y_down, BoundaryCondition bc_z_down) {
			//<maxNum 
			_mesh.resize(Nx*Ny*Nz);
			limit_x = Nx;
			limit_y = Ny;
			limit_z = Nz;
			dr = _dr;
			_bc_x_up = bc_x_up;
			_bc_y_up = bc_y_up;
			_bc_z_up = bc_z_up;
			_bc_x_down = bc_x_down;
			_bc_y_down = bc_y_down;
			_bc_z_down = bc_z_down;
			_dimention = Three_Dimension;
			if (Nx == 1 || Ny == 1 || Nz == 1)
				_dimention = Two_Dimension;
			if ((Nx == 1 && Ny == 1) || (Nx == 1 && Nz == 1) || (Ny == 1 && Nz == 1))
				_dimention = One_Dimension;
			for(int z = 0; z < Nz; z++)
				for (int y = 0; y < Ny; y++)
					for (int x = 0; x < Nx; x++) {
						PhaseNode* _up_x;
						PhaseNode* _down_x;
						PhaseNode* _up_y;
						PhaseNode* _down_y;
						PhaseNode* _up_z;
						PhaseNode* _down_z;
						if (bc_x_up == BoundaryCondition::ADIABATIC || bc_x_up == BoundaryCondition::FIXED) {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
						}
						else {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = 0 + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
						}
						if (bc_x_down == BoundaryCondition::ADIABATIC || bc_x_down == BoundaryCondition::FIXED) {
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							_down_x = &_mesh[down];
						}
						else {
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = (Nx - 1) + y * limit_x + z * limit_x * limit_y;
							_down_x = &_mesh[down];
						}

						if (bc_y_up == BoundaryCondition::ADIABATIC || bc_y_up == BoundaryCondition::FIXED) {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
						}
						else {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + 0 * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
						}
						if (bc_y_down == BoundaryCondition::ADIABATIC || bc_y_down == BoundaryCondition::FIXED) {
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							_down_y = &_mesh[down];
						}
						else {
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + (Ny - 1) * limit_x + z * limit_x * limit_y;
							_down_y = &_mesh[down];
						}

						if (bc_z_up == BoundaryCondition::ADIABATIC || bc_z_up == BoundaryCondition::FIXED) {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_z = &_mesh[up];
						}
						else {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + 0 * limit_x * limit_y;
							_up_z = &_mesh[up];
						}
						if (bc_z_down == BoundaryCondition::ADIABATIC || bc_z_down == BoundaryCondition::FIXED) {
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							_down_z = &_mesh[down];
						}
						else {
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + (Nz - 1) * limit_x * limit_y;
							_down_z = &_mesh[down];
						}
						_mesh[x + y * limit_x + z * limit_x * limit_y].connect_mesh(x, y, z, *_up_x, *_down_x, *_up_y, *_down_y, *_up_z, *_down_z);
					}
		}
		void normalize_phi_in_mesh() {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				double sum_phi = 0.0;
				for (auto phase = node->begin(); phase < node->end(); phase++)
					sum_phi += phase->phi;
				for (auto phase = node->begin(); phase < node->end(); phase++)
					phase->phi /= sum_phi;
			}
		}

		void add_customFlag_to_allnodes(int index, int flag) {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				node->customFlags.add_int(index, flag);
			}
			info_node.customFlags.add_int(index, flag);
		}
		void add_customValue_to_allnodes(int index, double value) {
			for (auto node = _mesh.begin();			node < _mesh.end(); node++) {
				node->customValues.add_double(index, value);
			}
			info_node.customValues.add_double(index, value);
		}
		void add_customVec3_to_allnodes(int index, Vector3 value) {
			for (auto node = _mesh.begin();			node < _mesh.end(); node++) {
				node->customVec3s.add_vec(index, value);
			}
			info_node.customVec3s.add_vec(index, value);
		}
		void add_customVec6_to_allnodes(int index, Vector6 value) {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				node->customVec6s.add_vec(index, value);
			}
			info_node.customVec6s.add_vec(index, value);
		}
		void write_customValue_in_node(ofstream& fout, int index, string _name) {
			string name = "\"cval_" + _name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < limit_z; ++k)
				for (int j = 0; j < limit_y; ++j)
					for (int i = 0; i < limit_x; ++i) {
						PhaseNode& node = operator()(i, j, k);
						for (auto val = node.customValues.begin(); val < node.customValues.end(); val++)
							if (val->index == index)
								fout << val->value << endl;
					}
			fout << "</DataArray>" << endl;
		}
		void write_customFlag_in_node(ofstream& fout, int index, string _name) {
			string name = "\"cflag_" + _name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < limit_z; ++k)
				for (int j = 0; j < limit_y; ++j)
					for (int i = 0; i < limit_x; ++i) {
						PhaseNode& node = operator()(i, j, k);
						for (auto val = node.customFlags.begin(); val < node.customFlags.end(); val++)
							if (val->index == index)
								fout << val->value << endl;
					}
			fout << "</DataArray>" << endl;
		}
		void write_customVec3_in_node(ofstream& fout, int index, string _name) {
			string name = "\"cvec3_" + _name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < limit_z; ++k)
				for (int j = 0; j < limit_y; ++j)
					for (int i = 0; i < limit_x; ++i) {
						PhaseNode& node = operator()(i, j, k);
						for (auto val = node.customVec3s.begin(); val < node.customVec3s.end(); val++)
							if (val->index == index) {
								fout << val->vec[0] << " "
									<< val->vec[1] << " "
									<< val->vec[2] << endl;
							}
					}
			fout << "</DataArray>" << endl;
		}
		void write_customVec6_in_node(ofstream& fout, int index, string _name) {
			for (int elem = 0; elem < 6; elem++) {
				string name = "\"cvec6_" + _name + "_" + to_string(elem) + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
				for (int k = 0; k < limit_z; ++k)
					for (int j = 0; j < limit_y; ++j)
						for (int i = 0; i < limit_x; ++i) {
							PhaseNode& node = operator()(i, j, k);
							for (auto val = node.customVec6s.begin(); val < node.customVec6s.end(); val++)
								if (val->index == index) {
									fout << val->vec[elem] << endl;
								}
						}
				fout << "</DataArray>" << endl;
			}
		}
		void delete_customFlag_in_allnodes(int index) {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				node->customFlags.erase(index);
			}
			info_node.customFlags.erase(index);
		}
		void delete_customValue_in_allnodes(int index) {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				node->customValues.erase(index);
			}
			info_node.customValues.erase(index);
		}
		void delete_customVec3_in_allnodes(int index) {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				node->customVec3s.erase(index);
			}
			info_node.customVec3s.erase(index);
		}
		void delete_customVec6_in_allnodes(int index) {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++) {
				node->customVec6s.erase(index);
			}
			info_node.customVec6s.erase(index);
		}

		void free() {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++)
				node->clear();
			_mesh.clear();
			info_node.clear();
			limit_x = 0;
			limit_y = 0;
			limit_z = 0;
			dr = 0.0;
		}
		int size() {
			return int(_mesh.size());
		}

		// ******   Settings   ****** //
		int limit_x;
		int limit_y;
		int limit_z;
		double dr;
		Dimension _dimention;
		BoundaryCondition _bc_x_up;
		BoundaryCondition _bc_y_up;
		BoundaryCondition _bc_z_up;
		BoundaryCondition _bc_x_down;
		BoundaryCondition _bc_y_down;
		BoundaryCondition _bc_z_down;
		PhaseNode info_node;
		std::vector<PhaseNode> _mesh;
	};
	inline Flag FieldStorage_forPhaseNode::currentFlag(PhaseNode& node, int phaseIndex) {
		if (node[phaseIndex].phi >= Phi_Num_Cut_Off && node[phaseIndex].phi <= (1.0 - Phi_Num_Cut_Off))
			return pf_INTERFACE;
		if (node[phaseIndex].phi < Phi_Num_Cut_Off) {
			for (int direction = 0; direction < 6; direction++) {
				PhaseNode& check_node = node.get_neighbor_node(Direction(direction));
				for (auto check_phase = check_node.begin(); check_phase < check_node.end(); check_phase++)
					if (check_phase->index == phaseIndex && check_phase->phi >= Phi_Num_Cut_Off)
						return pf_NEAR_INTERFACE;
			}
		}
		else if (node[phaseIndex].phi > (1.0 - Phi_Num_Cut_Off)) {
			for (int direction = 0; direction < 6; direction++) {
				PhaseNode& check_node = node.get_neighbor_node(Direction(direction));
				for (auto check_phase = check_node.begin(); check_phase < check_node.end(); check_phase++)
					if (check_phase->index == phaseIndex && check_phase->phi <= (1.0 - Phi_Num_Cut_Off))
						return pf_NEAR_INTERFACE;
			}
		}
		return pf_BULK;
	}
	inline void FieldStorage_forPhaseNode::upgrade(PhaseNode& node, int phaseIndex) {
		node[phaseIndex]._flag = pf_INTERFACE;
		if ((node._x - 1) >= 0 || _bc_x_down == PERIODIC)
			if (node.get_neighbor_node(Direction::x_down)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::x_down)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._x + 1) < limit_x || _bc_x_up == PERIODIC)
			if (node.get_neighbor_node(Direction::x_up)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::x_up)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._y - 1) >= 0 || _bc_y_down == PERIODIC)
			if (node.get_neighbor_node(Direction::y_down)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::y_down)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._y + 1) < limit_y || _bc_y_up == PERIODIC)
			if (node.get_neighbor_node(Direction::y_up)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::y_up)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._z - 1) >= 0 || _bc_z_down == PERIODIC)
			if (node.get_neighbor_node(Direction::z_down)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::z_down)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._z + 1) < limit_z || _bc_z_up == PERIODIC)
			if (node.get_neighbor_node(Direction::z_up)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::z_up)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		
	}
	inline void FieldStorage_forPhaseNode::downgrade(PhaseNode& node, int phaseIndex) {
		if (node[phaseIndex]._flag == pf_INTERFACE) {
			node[phaseIndex]._flag = pf_NEAR_INTERFACE;
			return;
		}
		bool _xu = false, _xd = false, _yu = false, _yd = false, _zu = false, _zd = false;
		if ((node._x + 1) < limit_x || _bc_x_up == PERIODIC)
			_xu = (*this).is_interface(node.get_neighbor_node(Direction::x_up), phaseIndex);
		if ((node._x - 1) >= 0 || _bc_x_down == PERIODIC)
			_xd = (*this).is_interface(node.get_neighbor_node(Direction::x_down), phaseIndex);
		if ((node._y + 1) < limit_y || _bc_y_up == PERIODIC)
			_yu = (*this).is_interface(node.get_neighbor_node(Direction::y_up), phaseIndex);
		if ((node._y - 1) >= 0 || _bc_y_down == PERIODIC)
			_yd = (*this).is_interface(node.get_neighbor_node(Direction::y_down), phaseIndex);
		if ((node._z + 1) < limit_z || _bc_z_up == PERIODIC)
			_zu = (*this).is_interface(node.get_neighbor_node(Direction::z_up), phaseIndex);
		if ((node._z - 1) >= 0 || _bc_z_down == PERIODIC)
			_zd = (*this).is_interface(node.get_neighbor_node(Direction::z_down), phaseIndex);
		if (_xu || _xd || _yu || _yd || _zu || _zd) {
			node[phaseIndex]._flag = pf_NEAR_INTERFACE;
		}
		else {
			node[phaseIndex]._flag = pf_BULK;
			if (node[phaseIndex].phi < Phi_Num_Cut_Off)
				node[phaseIndex].phi = 0.0;
			if (node[phaseIndex].phi > (1 - Phi_Num_Cut_Off))
				node[phaseIndex].phi = 1.0;
		}
	}
	inline bool FieldStorage_forPhaseNode::is_interface(PhaseNode& node, int phaseIndex) {
		if (node[phaseIndex].phi >= Phi_Num_Cut_Off && node[phaseIndex].phi <= (1.0 - Phi_Num_Cut_Off))
			return true;
		else
			return false;
	}
	inline PhaseNode& FieldStorage_forPhaseNode::operator()(const int x, const int y, const int z) {
#ifdef _DEBUG
		if (limit_x <= 0 || limit_y <= 0 || limit_z <= 0) {
			cout << "Settings error ! limit value error !" << endl;
			SYS_PROGRAM_STOP;
		}
#endif
		int _x = x, _y = y, _z = z;
		if (x >= 0 && x < limit_x && y >= 0 && y < limit_y && z >= 0 && z < limit_z) {
			return _mesh[x + y * limit_x + z * limit_x * limit_y];
		}
		if (x < 0) {
			if (_bc_x_down == PERIODIC)
				_x = _x + limit_x;
			else
				_x = 0;
		}
		if (x >= limit_x) {
			if (_bc_x_up == PERIODIC)
				_x = _x - limit_x;
			else
				_x = limit_x - 1;
		}
		if (y < 0) {
			if (_bc_y_down == PERIODIC)
				_y = _y + limit_y;
			else
				_y = 0;
		}
		if (y >= limit_y) {
			if (_bc_y_up == PERIODIC)
				_y = _y - limit_y;
			else
				_y = limit_y - 1;
		}
		if (z < 0) {
			if (_bc_z_down == PERIODIC)
				_z = _z + limit_z;
			else
				_z = 0;
		}
		if (z >= limit_z) {
			if (_bc_z_up == PERIODIC)
				_z = _z - limit_z;
			else
				_z = limit_z - 1;
		}
		return (*this)(_x, _y, _z);
	}
	inline FieldStorage_forPhaseNode& FieldStorage_forPhaseNode::operator=(const FieldStorage_forPhaseNode& n) {
		_mesh = n._mesh;
		info_node = n.info_node;
		limit_x = n.limit_x;
		limit_y = n.limit_y;
		limit_z = n.limit_z;
		_bc_x_down = n._bc_x_down;
		_bc_x_up = n._bc_x_up;
		_bc_y_down = n._bc_y_down;
		_bc_y_up = n._bc_y_up;
		_bc_z_down = n._bc_z_down;
		_bc_z_up = n._bc_z_up;
		dr = n.dr;
		_dimention = n._dimention;
		return *this;
	}

	class FieldStorage_forMechanicNode {
	public:
		FieldStorage_forMechanicNode() {
		};
		~FieldStorage_forMechanicNode() {
			free();
		};

		void init(int _Nx, int _Ny, int _Nz, double _dx, BoundaryCondition bc_x, BoundaryCondition bc_y, BoundaryCondition bc_z) {
			//<maxNum 
			Nx = _Nx;
			Ny = _Ny;
			Nz = _Nz;
			if (bc_x != PERIODIC) {
				Nx = _Nx * 2;
			}
			if (bc_y != PERIODIC) {
				Ny = _Ny * 2;
			}
			if (bc_z != PERIODIC) {
				Nz = _Nz * 2;
			}
			_mesh.resize(Nx * Ny * Nz);
			_virtual_strain_buff.resize(Nx * Ny * Nz);
			dx = _dx;
		}
		void free() {
			_mesh.clear();
			_virtual_strain_buff.clear();
			Nx = 0;
			Ny = 0;
			Nz = 0;
			dx = 0.0;
		}
		int size() {
			return int(_mesh.size());
		}

		MechanicNode& operator()(const int x, const int y, const int z) {
			if (x >= 0 && x < Nx && y >= 0 && y < Ny && z >= 0 && z < Nz)
				return _mesh[x + y * Nx + z * Nx * Ny];
			else {
				cout << "FieldStorage_forMechanicNode error ! limit value error !" << endl;
				SYS_PROGRAM_STOP;
			}
		}

		std::vector<MechanicNode> _mesh;
		std::vector<vStrain> _virtual_strain_buff;
		int Nx;                                                                ///< System size along X direction
		int Ny;                                                                ///< System size along y direction
		int Nz;                                                                ///< System size along z direction
		double dx;
	};

	struct VectorNode {
		VectorNode() {
		}
		~VectorNode() {
			clear();
		}
		void clear() {
			_up_x = nullptr;
			_down_x = nullptr;
			_up_y = nullptr;
			_down_y = nullptr;
			_up_z = nullptr;
			_down_z = nullptr;
			vals.clear();
		}
		void connect_mesh(int x, int y, int z, VectorNode& up_x, VectorNode& down_x
			, VectorNode& up_y, VectorNode& down_y, VectorNode& up_z, VectorNode& down_z) {
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
		VectorNode& get_neighbor_node(Direction _d);
		VectorNode& get_long_range_node(int relative_x, int relative_y, int relative_z);
		vector<double> vals;
		int _x;
		int _y;
		int _z;
	private:
		VectorNode* _up_x;
		VectorNode* _down_x;
		VectorNode* _up_y;
		VectorNode* _down_y;
		VectorNode* _up_z;
		VectorNode* _down_z;
	};
	inline VectorNode& VectorNode::get_neighbor_node(Direction _d) {
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
	inline VectorNode& VectorNode::get_long_range_node(int relative_x, int relative_y, int relative_z) {
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
	class FieldStorage_forVector
	{
	public:
		FieldStorage_forVector() {};
		~FieldStorage_forVector() {
			try {
				_mesh.clear();
			}
			catch (...) {
				std::cout << "Destroy FieldStorage_forVector failed, problems occur." << endl;
			}
		};

		VectorNode& operator()(const int x, const int y, const int z);                                       //< Index operator for accessing the n's field value
		FieldStorage_forVector& operator=(const FieldStorage_forVector& n);
		void init(int Nx, int Ny, int Nz, double _dr, BoundaryCondition bc_x_up, BoundaryCondition bc_y_up, BoundaryCondition bc_z_up, BoundaryCondition bc_x_down, BoundaryCondition bc_y_down, BoundaryCondition bc_z_down) {
			//<maxNum 
			_mesh.resize(Nx * Ny * Nz);
			limit_x = Nx;
			limit_y = Ny;
			limit_z = Nz;
			x_up_bc = bc_x_up;
			y_up_bc = bc_y_up;
			z_up_bc = bc_z_up;
			x_down_bc = bc_x_down;
			y_down_bc = bc_y_down;
			z_down_bc = bc_z_down;
			dr = _dr;
			_dimention = Three_Dimension;
			if (Nx == 1 || Ny == 1 || Nz == 1)
				_dimention = Two_Dimension;
			if ((Nx == 1 && Ny == 1) || (Nx == 1 && Nz == 1) || (Ny == 1 && Nz == 1))
				_dimention = One_Dimension;
			for (int z = 0; z < Nz; z++)
				for (int y = 0; y < Ny; y++)
					for (int x = 0; x < Nx; x++) {
						VectorNode* _up_x;
						VectorNode* _down_x;
						VectorNode* _up_y;
						VectorNode* _down_y;
						VectorNode* _up_z;
						VectorNode* _down_z;
						if (bc_x_up == BoundaryCondition::ADIABATIC || bc_x_up == BoundaryCondition::FIXED) {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
						}
						else {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = 0 + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
						}
						if (bc_x_down == BoundaryCondition::ADIABATIC || bc_x_down == BoundaryCondition::FIXED) {
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							_down_x = &_mesh[down];
						}
						else {
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = (Nx - 1) + y * limit_x + z * limit_x * limit_y;
							_down_x = &_mesh[down];
						}
						if (bc_y_up == BoundaryCondition::ADIABATIC || bc_y_up == BoundaryCondition::FIXED) {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
						}
						else {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + 0 * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
						}
						if (bc_y_down == BoundaryCondition::ADIABATIC || bc_y_down == BoundaryCondition::FIXED) {
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							_down_y = &_mesh[down];
						}
						else {
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + (Ny - 1) * limit_x + z * limit_x * limit_y;
							_down_y = &_mesh[down];
						}
						if (bc_z_up == BoundaryCondition::ADIABATIC || bc_z_up == BoundaryCondition::FIXED) {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_z = &_mesh[up];
						}
						else {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + 0 * limit_x * limit_y;
							_up_z = &_mesh[up];
						}
						if (bc_z_down == BoundaryCondition::ADIABATIC || bc_z_down == BoundaryCondition::FIXED) {
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							_down_z = &_mesh[down];
						}
						else {
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + (Nz - 1) * limit_x * limit_y;
							_down_z = &_mesh[down];
						}
						_mesh[x + y * limit_x + z * limit_x * limit_y].connect_mesh(x, y, z, *_up_x, *_down_x, *_up_y, *_down_y, *_up_z, *_down_z);
					}
		}
		void free() {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++)
				node->clear();
			_mesh.clear();
			limit_x = 0;
			limit_y = 0;
			limit_z = 0;
			dr = 0.0;
		}
		int size() {
			return int(_mesh.size());
		}

		// ******   Settings   ****** //
		BoundaryCondition x_up_bc;
		BoundaryCondition y_up_bc;
		BoundaryCondition z_up_bc;
		BoundaryCondition x_down_bc;
		BoundaryCondition y_down_bc;
		BoundaryCondition z_down_bc;
		int limit_x;
		int limit_y;
		int limit_z;
		double dr;
		Dimension _dimention;
		std::vector<VectorNode> _mesh;
	};
	inline VectorNode& FieldStorage_forVector::operator()(const int x, const int y, const int z) {
#ifdef _DEBUG
		if (limit_x <= 0 || limit_y <= 0 || limit_z <= 0) {
			cout << "Settings error ! mesh size error !" << endl;
			SYS_PROGRAM_STOP;
		}
#endif
		int _x = x, _y = y, _z = z;
		if (x >= 0 && x < limit_x && y >= 0 && y < limit_y && z >= 0 && z < limit_z) {
			return _mesh[x + y * limit_x + z * limit_x * limit_y];
		}
		if (x < 0) {
			if (x_down_bc == PERIODIC)
				_x = _x + limit_x;
			else
				_x = 0;
		}
		if (x >= limit_x) {
			if (x_up_bc == PERIODIC)
				_x = _x - limit_x;
			else
				_x = limit_x - 1;
		}
		if (y < 0) {
			if (y_down_bc == PERIODIC)
				_y = _y + limit_y;
			else
				_y = 0;
		}
		if (y >= limit_y) {
			if (y_up_bc == PERIODIC)
				_y = _y - limit_y;
			else
				_y = limit_y - 1;
		}
		if (z < 0) {
			if (z_down_bc == PERIODIC)
				_z = _z + limit_z;
			else
				_z = 0;
		}
		if (z >= limit_z) {
			if (z_up_bc == PERIODIC)
				_z = _z - limit_z;
			else
				_z = limit_z - 1;
		}
		return (*this)(_x, _y, _z);
	}
	inline FieldStorage_forVector& FieldStorage_forVector::operator=(const FieldStorage_forVector& n) {
		_mesh = n._mesh;
		dr = n.dr;
		_dimention = n._dimention;
		limit_x = n.limit_x;
		limit_y = n.limit_y;
		limit_z = n.limit_z;
		x_up_bc = n.x_up_bc;
		y_up_bc = n.y_up_bc;
		z_up_bc = n.z_up_bc;
		x_down_bc = n.x_down_bc;
		y_down_bc = n.y_down_bc;
		z_down_bc = n.z_down_bc;
		return *this;
	}
}
