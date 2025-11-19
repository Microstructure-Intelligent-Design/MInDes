#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "../../Modules_Params.h"
namespace pf {
	struct Data_MeshInfo {
		bool is_phi_mesh;
		size_t phi_number;
		size_t PNx;
		size_t PNy;
		size_t PNz;
		bool is_con_mesh;
		size_t con_number;
		size_t CNx;
		size_t CNy;
		size_t CNz;
		bool is_temp_mesh;
		size_t TNx;
		size_t TNy;
		size_t TNz;
		bool is_fluid_mesh;
		size_t dis_func_number;
		size_t FNx;
		size_t FNy;
		size_t FNz;
		void operator=(const Data_MeshInfo& n) {
			PNx = n.PNx;
			PNy = n.PNy;
			PNz = n.PNz;
			CNx = n.CNx;
			CNy = n.CNy;
			CNz = n.CNz;
			TNx = n.TNx;
			TNy = n.TNy;
			TNz = n.TNz;
			FNx = n.FNx;
			FNy = n.FNy;
			FNz = n.FNz;
			is_phi_mesh = n.is_phi_mesh;
			phi_number = n.phi_number;
			is_con_mesh = n.is_con_mesh;
			con_number = n.con_number;
			is_temp_mesh = n.is_temp_mesh;
			is_fluid_mesh = n.is_fluid_mesh;
			dis_func_number = n.dis_func_number;
		}
	};
	namespace write_mesh_data {
		// - data file
		inline bool is_datafile_init = false;
		inline bool is_read_datafile_by_path = false;
		const std::string mainName = "MeshData";
		const std::string format = ".dat";
		inline std::string datafile_path = mainName + format;
		inline size_t output_frequence = 0;
		inline Data_MeshInfo datafile_report;
		inline size_t MESH_INDEX(size_t x, size_t y, size_t z, size_t Nx, size_t Ny) { return x + y * Nx + z * Nx * Ny; };
		bool write_dataFile(std::string fname);
		bool read_dataFile(std::string fname, Data_MeshInfo& mesh_info);
	}
}