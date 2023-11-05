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
	class WriteToFile
	{
	public:
		WriteToFile(string file_path) {
			init(file_path);
		};
		WriteToFile() {};
		~WriteToFile() {};
		void init(string file_path) {
			_path = file_path;
#if defined(_WIN32)
			_mkdir(_path.c_str());
#elif defined(__linux__)
			mkdir(_path.c_str(), 0777);
			// #else more?
#endif
		}

		// open vts file
		void open_vts_scalar_file(ofstream& fout, FieldStorage_forPhaseNode& fs, string tail);
		void open_vts_vec3_file(ofstream& fout, FieldStorage_forPhaseNode& fs, string tail);
		// close vts file
		void close_vts_file(ofstream& fout, FieldStorage_forPhaseNode& fs);
		// write scalar data
		void write_scalar_grains(ofstream& fout, FieldStorage_forPhaseNode& fs);

		void write_string_to_txt(string content, string _fname);
		void write_string_to_txt_and_screen(string content, string _fname);
		void add_string_to_txt(string content, string _fname);
		void add_string_to_txt_and_screen(string content, string _fname);
		void init_txt_file(string _fname);

		void write_string_to_m(string content, string _fname);
		void add_string_to_m(string content, string _fname);
		void init_m_file(string _fname);

		string _path;
	};
	inline void WriteToFile::open_vts_scalar_file(ofstream& fout, FieldStorage_forPhaseNode& fs, string tail) {
		string fname;
		fname = _path + dirSeparator + "scalar_variables_" + tail + ".vts";
		
		fout.open(fname);
		if (!fout) {
			cout << "Failed to write the vtk file..." << endl;
			fout.close();
			return;
		}
		fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
		fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
		fout << "<StructuredGrid WholeExtent=\""
			<< 0 << " " << fs.limit_x - 1 << " "
			<< 0 << " " << fs.limit_y - 1 << " "
			<< 0 << " " << fs.limit_z - 1 << "\"> " << endl;
		fout << "<PointData Scalars= \"ScalarData\">" << endl;
	}
	inline void WriteToFile::open_vts_vec3_file(ofstream& fout, FieldStorage_forPhaseNode& fs, string tail) {
		string fname;
		fname = _path + dirSeparator + "vec3_variables_" + tail + ".vts";
		fout.open(fname);
		if (!fout) {
			cout << "Failed to write the vtk file..." << endl;
			fout.close();
			return;
		}
		fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
		fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
		fout << "<StructuredGrid WholeExtent=\""
			<< 0 << " " << fs.limit_x - 1 << " "
			<< 0 << " " << fs.limit_y - 1 << " "
			<< 0 << " " << fs.limit_z - 1 << "\"> " << endl;
		fout << "<PointData  Vectors= \"VectorData\">" << endl;
	}
	inline void WriteToFile::close_vts_file(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "</PointData>" << endl;
		fout << "<Points>" << endl;
		fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					fout << node._x << " " << node._y << " " << node._z << "\n";
				}
		fout << "</DataArray>" << endl;
		fout << "</Points>" << endl;
		fout << "</StructuredGrid>" << endl;
		fout << "</VTKFile>" << endl;
		fout.close();
	}

	inline void WriteToFile::write_scalar_grains(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "grains" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					double fix = 0.0;
					for (auto phase = fs(i, j, k).begin(); phase < fs(i, j, k).end(); phase++)
						fix += phase->phi * phase->phi;
					fout << fix << endl;
				}
		fout << "</DataArray>" << endl;
	}

	inline void WriteToFile::write_string_to_txt(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to write the txt file!" << endl;
			fout.close();
			return;
		}
		fout << content << endl;
		fout.close();
	}

	inline void WriteToFile::write_string_to_txt_and_screen(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to write the txt file!" << endl;
			fout.close();
			return;
		}
		fout << content;
		cout << content;
		fout.close();
	}

	inline void WriteToFile::add_string_to_txt(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname, ios::app);
		if (!fout) {
			cout << "Failed to add the string to txt file!" << endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
	}

	inline void WriteToFile::add_string_to_txt_and_screen(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname, ios::app);
		if (!fout) {
			cout << "Failed to add the string to txt file!" << endl;
			fout.close();
			return;
		}
		fout << content;
		cout << content;
		fout.close();
	}

	inline void WriteToFile::init_txt_file(string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to init the txt file!" << endl;
			fout.close();
			return;
		}
		fout << "";
		fout.close();
	}

	inline void WriteToFile::write_string_to_m(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".m";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to write the m file!" << endl;
			fout.close();
			return;
		}
		fout << content << endl;
		fout.close();
	}

	inline void WriteToFile::add_string_to_m(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".m";
		ofstream fout(fname, ios::app);
		if (!fout) {
			cout << "Failed to add the string to m file!" << endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
	}

	inline void WriteToFile::init_m_file(string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".m";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to init the m file!" << endl;
			fout.close();
			return;
		}
		fout << "";
		fout.close();
	}

}