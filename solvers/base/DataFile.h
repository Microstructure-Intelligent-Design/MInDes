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
#include "FieldStorage.h"
using namespace std;

#define DataFile_NoneStep -1

namespace pf {
	enum Data_Custom_Flag { DCW_Value, DCW_Flag, DCW_Vec3 };
	struct Data_fieldStorage {
		BoundaryCondition x_down_bc;
		BoundaryCondition y_down_bc;
		BoundaryCondition z_down_bc;
		BoundaryCondition x_up_bc;
		BoundaryCondition y_up_bc;
		BoundaryCondition z_up_bc;
		int Nx;
		int Ny;
		int Nz;
		double dr;
	};
	struct Data_phaseNode {
		double temperature;
		int phaseNum;
		int conNum;
		int potentialNum;
		int custom_value_num;
		int custom_flag_num;
		int custom_vec_num;
	};
	struct Data_custom_info {
		Data_Custom_Flag info;
		int index;
	};
	struct Data_custom_value {
		int index;
		int flag;
		double value;
		double value2;
		double value3;
	};
	struct Data_phaseEntry {
		int _flag;
		int index;
		int phaseProperty;
		double phaseFraction;
		int conNum;
	};
	struct Data_nodeEntry {
		int index;
		double value;
	};
	struct Data_report {
		BoundaryCondition x_down_bc;
		BoundaryCondition y_down_bc;
		BoundaryCondition z_down_bc;
		BoundaryCondition x_up_bc;
		BoundaryCondition y_up_bc;
		BoundaryCondition z_up_bc;
		int Nx;
		int Ny;
		int Nz;
		double dr;
		int_box custom_vars;
		int_box phi_property;
		int2_box phi_comps;
		vector<int> comps;
		Data_report& operator=(const Data_report& n) {
			x_down_bc = n.x_down_bc;
			y_down_bc = n.y_down_bc;
			z_down_bc = n.z_down_bc;
			x_up_bc = n.x_up_bc;
			y_up_bc = n.y_up_bc;
			z_up_bc = n.z_up_bc;
			Nx = n.Nx;
			Ny = n.Ny;
			Nz = n.Nz;
			dr = n.dr;
			custom_vars = n.custom_vars;
			phi_property = n.phi_property;
			phi_comps = n.phi_comps;
			comps = n.comps;
		}
	};
	class DataFile
	{
	public:
		DataFile() {
			mainName = "dataFile";
			format = ".dat";
			_path = "";
			DataFile_header = new char[50];
			string head = "## MInDes PROGRAM FILE VERSION 1.0 ##";
			memcpy(DataFile_header, head.c_str(), 50);
			int i = 2;
		};
		~DataFile() {};
		void set_mainName(string name) {
			mainName = name;
		}
		void set_path(string path) {
			_path = path;
#if defined(_WIN32)
			_mkdir(_path.c_str());
#elif defined(__linux__)
			mkdir(_path.c_str(), 0777);
			// #else more?
#endif
		}
		bool write_dataFile(FieldStorage_forPhaseNode& phaseMesh, string mark);
		bool read_dataFile(FieldStorage_forPhaseNode& phaseMesh, string name, Data_report& report);
		DataFile& operator=(const DataFile& n);

		void write_PhaseNode_to_dataFile(ofstream &fout, PhaseNode &n);
		void read_PhaseNode_from_dataFile(fstream& fout, PhaseNode& n);
		void register_custom_value(int custom_index, Data_Custom_Flag custom_type);
		string _path;
		string mainName;
		string format;
		vector<Data_custom_info> custom_value_register;
		char* DataFile_header;
	};
	inline DataFile& DataFile::operator=(const DataFile& n) {
		_path = n._path;
		mainName = n.mainName;
		format = n.format;
		memcpy(DataFile_header, n.DataFile_header, 50);
	}
	inline void DataFile::write_PhaseNode_to_dataFile(ofstream& fout, PhaseNode& n) {
		Data_phaseNode phaseNode;
		Data_phaseEntry phaseEntry;
		Data_nodeEntry nodeEntry;
		phaseNode.temperature = n.temperature.T;
		phaseNode.phaseNum = n.size();
		phaseNode.conNum = n.x.size();
		phaseNode.potentialNum = n.potential.size();
		phaseNode.custom_flag_num = 0;
		phaseNode.custom_value_num = 0;
		phaseNode.custom_vec_num = 0;
		for (auto _register = custom_value_register.begin(); _register < custom_value_register.end(); _register++) {
			if (_register->info == Data_Custom_Flag::DCW_Flag)
				phaseNode.custom_flag_num++;
			else if (_register->info == Data_Custom_Flag::DCW_Value)
				phaseNode.custom_value_num++;
			else if (_register->info == Data_Custom_Flag::DCW_Vec3)
				phaseNode.custom_vec_num++;
		}
		fout.write((const char*)&phaseNode, sizeof(Data_phaseNode));
		for (auto c = n.x.begin(); c < n.x.end(); c++) {
			nodeEntry.index = c->index;
			nodeEntry.value = c->value;
			fout.write((const char*)&nodeEntry, sizeof(Data_nodeEntry));
		}
		for (auto c = n.potential.begin(); c < n.potential.end(); c++) {
			nodeEntry.index = c->index;
			nodeEntry.value = c->value;
			fout.write((const char*)&nodeEntry, sizeof(Data_nodeEntry));
		}
		for (auto _register = custom_value_register.begin(); _register < custom_value_register.end(); _register++) {
			if (_register->info == Data_Custom_Flag::DCW_Flag) {
				Data_custom_value value;
				value.index = _register->index;
				value.flag = n.customFlags[_register->index];
				value.value = 0.0;
				value.value2 = 0.0;
				value.value3 = 0.0;
				fout.write((const char*)&value, sizeof(Data_custom_value));
			}
		}
		for (auto _register = custom_value_register.begin(); _register < custom_value_register.end(); _register++) {
			if (_register->info == Data_Custom_Flag::DCW_Value) {
				Data_custom_value value;
				value.index = _register->index;
				value.value = n.customValues[_register->index];
				value.value2 = 0.0;
				value.value3 = 0.0;
				value.flag = 0;
				fout.write((const char*)&value, sizeof(Data_custom_value));
			}
		}
		for (auto _register = custom_value_register.begin(); _register < custom_value_register.end(); _register++) {
			if (_register->info == Data_Custom_Flag::DCW_Vec3) {
				Data_custom_value value;
				value.index = _register->index;
				value.flag = 0;
				value.value = n.customVec3s[_register->index][0];
				value.value2 = n.customVec3s[_register->index][1];
				value.value3 = n.customVec3s[_register->index][2];
				fout.write((const char*)&value, sizeof(Data_custom_value));
			}
		}
		///< storage for _Phase
		for (auto p = n._Phase.begin(); p < n._Phase.end(); p++) {
			phaseEntry._flag = p->_flag;
			phaseEntry.index = p->index;
			phaseEntry.phaseProperty = p->property;
			phaseEntry.phaseFraction = p->phi;
			phaseEntry.conNum = p->x.size();
			fout.write((const char*)&phaseEntry, sizeof(Data_phaseEntry));
			///< con
			for (auto c = p->x.begin(); c < p->x.end(); c++) {
				nodeEntry.index = c->index;
				nodeEntry.value = c->value;
				fout.write((const char*)&nodeEntry, sizeof(Data_nodeEntry));
			}
		}
	}
	inline void DataFile::read_PhaseNode_from_dataFile(fstream& fin, PhaseNode& n) {
		Data_phaseNode phaseNode;
		Data_phaseEntry phaseEntry;
		Data_nodeEntry nodeEntry;
		fin.read((char*)&phaseNode, sizeof(Data_phaseNode));
		n.temperature.T = phaseNode.temperature;
		n._Phase.resize(phaseNode.phaseNum);
		for (int c_num = 0; c_num < phaseNode.conNum; c_num++) {
			fin.read((char*)&nodeEntry, sizeof(Data_nodeEntry));
			n.x.add_con(nodeEntry.index, nodeEntry.value);
		}
		for (int p_num = 0; p_num < phaseNode.potentialNum; p_num++) {
			fin.read((char*)&nodeEntry, sizeof(Data_nodeEntry));
			n.potential.add_con(nodeEntry.index, nodeEntry.value);
		}
		for (int i = 0; i < phaseNode.custom_flag_num; i++) {
			Data_custom_value cus_val;
			fin.read((char*)&cus_val, sizeof(Data_custom_value));
			n.customFlags.add_int(cus_val.index, cus_val.flag);
		}
		for (int i = 0; i < phaseNode.custom_value_num; i++) {
			Data_custom_value cus_val;
			fin.read((char*)&cus_val, sizeof(Data_custom_value));
			n.customValues.add_double(cus_val.index, cus_val.value);
		}
		for (int i = 0; i < phaseNode.custom_vec_num; i++) {
			Data_custom_value cus_val;
			fin.read((char*)&cus_val, sizeof(Data_custom_value));
			n.customVec3s.add_vec(cus_val.index, Vector3(cus_val.value, cus_val.value2, cus_val.value3));
		}
		///< storage for _Phase
		for (auto p = n._Phase.begin(); p < n._Phase.end(); p++) {
			fin.read((char*)&phaseEntry, sizeof(Data_phaseEntry));
			p->_flag = phaseEntry._flag;
			p->index = phaseEntry.index;
			p->property = phaseEntry.phaseProperty;
			p->phi = phaseEntry.phaseFraction;
			///< con
			for (int c_num = 0; c_num < phaseEntry.conNum; c_num++) {
				fin.read((char*)&nodeEntry, sizeof(Data_nodeEntry));
				p->x.add_con(nodeEntry.index, nodeEntry.value);
			}
		}
	}
	inline bool DataFile::write_dataFile(FieldStorage_forPhaseNode& phaseMesh, string mark) {
		string fname;
		if (_path == "")
			fname = mainName + mark + format;
		else
			fname = _path + dirSeparator + mainName + mark + format;
		ofstream fout(fname, std::ios::binary);
		if (!fout) {
			cout << "Failed to write the data file!" << endl;
			fout.close();
			return false;
		}
		fout.write((const char*)DataFile_header, 50);
		{ ///< defined in sequence(same with read)
			Data_fieldStorage fieldStorage;
			fieldStorage.dr = phaseMesh.dr;
			fieldStorage.Nx = phaseMesh.limit_x;
			fieldStorage.Ny = phaseMesh.limit_y;
			fieldStorage.Nz = phaseMesh.limit_z;
			fieldStorage.x_down_bc = phaseMesh._bc_x_down;
			fieldStorage.x_up_bc = phaseMesh._bc_x_up;
			fieldStorage.y_down_bc = phaseMesh._bc_y_down;
			fieldStorage.y_up_bc = phaseMesh._bc_y_up;
			fieldStorage.z_down_bc = phaseMesh._bc_z_down;
			fieldStorage.z_up_bc = phaseMesh._bc_z_up;

			fout.write((const char*)&fieldStorage, sizeof(Data_fieldStorage));
			///< storage for _mesh
			for (auto n = phaseMesh._mesh.begin(); n < phaseMesh._mesh.end(); n++)
				write_PhaseNode_to_dataFile(fout, (*n));
		}
		fout.close();
		return true;
	}
	inline bool DataFile::read_dataFile(FieldStorage_forPhaseNode& phaseMesh, string name, Data_report& report) {
		Data_fieldStorage fieldStorage;
		char header[50];
		std::fstream fin(name, std::ios::binary | ios::in);
		if (!fin) {
			cout << "Failed to read the aim file!" << endl;
			fin.close();
			return false;
		}
		fin.read((char*)header, 50);
		if (strcmp(header, DataFile_header)) {
			cout << "File format error, check content or version!" << endl;
			fin.close();
			return false;
		}
		{  ///< defined in sequence(same with write)
			fin.read((char*)&fieldStorage, sizeof(Data_fieldStorage));
			phaseMesh.init(fieldStorage.Nx, fieldStorage.Ny, fieldStorage.Nz, fieldStorage.dr, 
				fieldStorage.x_up_bc, fieldStorage.y_up_bc, fieldStorage.z_up_bc, fieldStorage.x_down_bc, fieldStorage.y_down_bc, fieldStorage.z_down_bc);

			for (auto n = phaseMesh._mesh.begin(); n < phaseMesh._mesh.end(); n++)
				read_PhaseNode_from_dataFile(fin, (*n));
		}
		{
			report.x_down_bc = fieldStorage.x_down_bc;
			report.x_up_bc = fieldStorage.x_up_bc;
			report.y_down_bc = fieldStorage.y_down_bc;
			report.y_up_bc = fieldStorage.y_up_bc;
			report.z_down_bc = fieldStorage.z_down_bc;
			report.z_up_bc = fieldStorage.z_up_bc;
			report.Nx = fieldStorage.Nx;
			report.Ny = fieldStorage.Ny;
			report.Nz = fieldStorage.Nz;
			report.dr = fieldStorage.dr;
			report.phi_property.clear();
			report.phi_comps.clear();
			report.custom_vars.clear();
			report.comps.clear();

			for (auto node = phaseMesh._mesh.begin(); node < phaseMesh._mesh.end(); node++) {
				for (auto phase = node->begin(); phase < node->end(); phase++) {
					phaseMesh.info_node.add_phase(phase->index, phase->property, 0, 0.0);
					for (auto phase_con = phase->x.begin(); phase_con < phase->x.end(); phase_con++)
						phaseMesh.info_node[phase->index].x.add_con(phase_con->index, 0.0);
				}
				for (auto con = node->x.begin(); con < node->x.end(); con++)
					phaseMesh.info_node.x.add_con(con->index, 0.0);
				for (auto pot = node->potential.begin(); pot < node->potential.end(); pot++)
					phaseMesh.info_node.potential.add_con(pot->index, 0.0);
				for (auto flag = node->customFlags.begin(); flag < node->customFlags.end(); flag++)
					phaseMesh.info_node.customFlags.add_int(flag->index, 0);
				for (auto var = node->customValues.begin(); var < node->customValues.end(); var++)
					phaseMesh.info_node.customValues.add_double(var->index, 0.0);
				for (auto vec3 = node->customVec3s.begin(); vec3 < node->customVec3s.end(); vec3++)
					phaseMesh.info_node.customVec3s.add_vec(vec3->index, Vector3(0.0, 0.0, 0.0));
			}

			for (auto phase = phaseMesh.info_node.begin(); phase < phaseMesh.info_node.end(); phase++) {
				report.phi_property.add_int(phase->index, phase->property);
				report.phi_comps.add_int(phase->index);
				for (auto con = phase->x.begin(); con < phase->x.end(); con++)
					report.phi_comps.add_int(phase->index, con->index, 0);
			}
			for (auto con = phaseMesh.info_node.x.begin(); con < phaseMesh.info_node.x.end(); con++)
				report.comps.push_back(con->index);
			for (auto _flag = phaseMesh.info_node.customFlags.begin(); _flag < phaseMesh.info_node.customFlags.end(); _flag++)
				report.custom_vars.add_int(_flag->index, Data_Custom_Flag::DCW_Flag);
			for (auto _value = phaseMesh.info_node.customValues.begin(); _value < phaseMesh.info_node.customValues.end(); _value++)
				report.custom_vars.add_int(_value->index, Data_Custom_Flag::DCW_Value);
			for (auto _vec = phaseMesh.info_node.customVec3s.begin(); _vec < phaseMesh.info_node.customVec3s.end(); _vec++)
				report.custom_vars.add_int(_vec->index, Data_Custom_Flag::DCW_Vec3);
		}
		fin.close();
		return true;
	}
	inline void DataFile::register_custom_value(int custom_index, Data_Custom_Flag custom_type) {
		for (auto i = custom_value_register.begin(); i < custom_value_register.end(); ++i)
			if (i->index == custom_index) {
				if (i->info != custom_type) {
					cout << "> ERROR, Datafile class, register_custom_value, register a value with same index but different type !" << endl;
					exit(0);
				}
				return;
			}
		Data_custom_info _new_info;
		_new_info.index = custom_index;
		_new_info.info = custom_type;
		custom_value_register.push_back(_new_info);
	}
}