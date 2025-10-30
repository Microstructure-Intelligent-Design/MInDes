#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "../model_modules/Model_Params.h"
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
		bool write_dataFile(std::string fname) {
			std::ofstream fout(fname, std::ios::binary);
			if (!fout) {
				std::cout << "Failed to write the data file!" << std::endl;
				fout.close();
				return false;
			}
			{ ///< defined in sequence(same with read)
				Data_MeshInfo mesh_info;
				mesh_info.is_phi_mesh = main_field::is_phi_field_on;
				mesh_info.phi_number = main_field::phi_number;
				mesh_info.PNx = main_field::phase_field.Nx();
				mesh_info.PNy = main_field::phase_field.Ny();
				mesh_info.PNz = main_field::phase_field.Nz();
				mesh_info.is_con_mesh = main_field::is_con_field_on;
				mesh_info.con_number = main_field::con_number;
				mesh_info.CNx = main_field::concentration_field.Nx();
				mesh_info.CNy = main_field::concentration_field.Ny();
				mesh_info.CNz = main_field::concentration_field.Nz();
				mesh_info.is_temp_mesh = main_field::is_temp_field_on;
				mesh_info.TNx = main_field::temperature_field.Nx();
				mesh_info.TNy = main_field::temperature_field.Ny();
				mesh_info.TNz = main_field::temperature_field.Nz();
				mesh_info.is_fluid_mesh = false;
				mesh_info.dis_func_number = 0;
				mesh_info.FNx = 0;
				mesh_info.FNy = 0;
				mesh_info.FNz = 0;
				fout.write((const char*)&mesh_info, sizeof(Data_MeshInfo));
				///< storage for mesh
				if (mesh_info.is_phi_mesh) {
					for (size_t x = 0; x < mesh_info.PNx; x++)
						for (size_t y = 0; y < mesh_info.PNy; y++)
							for (size_t z = 0; z < mesh_info.PNz; z++) {
								std::vector<REAL>& point = main_field::phase_field(x, y, z);
								for (size_t index = 0; index < mesh_info.phi_number; index++) {
									float phi = point[index];
									fout.write((const char*)&phi, sizeof(float));
								}
							}
				}
				if (mesh_info.is_con_mesh) {
					for (size_t x = 0; x < mesh_info.CNx; x++)
						for (size_t y = 0; y < mesh_info.CNy; y++)
							for (size_t z = 0; z < mesh_info.CNz; z++) {
								std::vector<REAL>& point = main_field::concentration_field(x, y, z);
								for (size_t index = 0; index < mesh_info.con_number; index++) {
									float con = point[index];
									fout.write((const char*)&con, sizeof(float));
								}
							}
				}
				if (mesh_info.is_temp_mesh) {
					for (size_t x = 0; x < mesh_info.TNx; x++)
						for (size_t y = 0; y < mesh_info.TNy; y++)
							for (size_t z = 0; z < mesh_info.TNz; z++) {
								float temp = main_field::temperature_field(x, y, z);
								fout.write((const char*)&temp, sizeof(float));
							}
				}
				if (mesh_info.is_fluid_mesh) {
					for (size_t x = 0; x < mesh_info.FNx; x++)
						for (size_t y = 0; y < mesh_info.FNy; y++)
							for (size_t z = 0; z < mesh_info.FNz; z++) {
								// std::vector<float>& point = fluid_field[MESH_INDEX(x, y, z, mesh_info.FNx, mesh_info.FNy)];
								// for (size_t index = 0; index < mesh_info.dis_func_number; index++)
								// 	fout.write((const char*)&point[index], sizeof(float));
							}
				}
			}
			fout.close();
			return true;
		}
		bool read_dataFile(std::string fname, Data_MeshInfo& mesh_info) {
			std::fstream fin(fname, std::ios::binary | std::ios::in);
			if (!fin) {
				std::cout << "Failed to read the aim file!" << std::endl;
				fin.close();
				return false;
			}
			fin.read((char*)&mesh_info, sizeof(Data_MeshInfo));
			///< read for mesh
			float data_buff = 0;
			if (mesh_info.is_phi_mesh) {
				size_t pNx = main_field::phase_field.Nx(), pNy = main_field::phase_field.Ny(), pNz = main_field::phase_field.Nz();
				for (size_t x = 0; x < mesh_info.PNx; x++)
					for (size_t y = 0; y < mesh_info.PNy; y++)
						for (size_t z = 0; z < mesh_info.PNz; z++) {
							if (x < pNx && y < pNy && z < pNz) {
								std::vector<REAL>& point = main_field::phase_field(x, y, z);
								for (size_t index = 0; index < mesh_info.phi_number; index++) {
									fin.read((char*)&data_buff, sizeof(float));
									if (index < main_field::phi_number)
										point[index] = data_buff;
								}
							}
							else {
								for (size_t index = 0; index < mesh_info.phi_number; index++)
									fin.read((char*)&data_buff, sizeof(float));
							}
						}
			}
			if (mesh_info.is_con_mesh) {
				size_t cNx = main_field::concentration_field.Nx(), cNy = main_field::concentration_field.Ny(),
					cNz = main_field::concentration_field.Nz();
				for (size_t x = 0; x < mesh_info.CNx; x++)
					for (size_t y = 0; y < mesh_info.CNy; y++)
						for (size_t z = 0; z < mesh_info.CNz; z++) {
							if (x < cNx && y < cNy && z < cNz) {
								std::vector<REAL>& point = main_field::concentration_field(x, y, z);
								for (size_t index = 0; index < mesh_info.con_number; index++) {
									fin.read((char*)&data_buff, sizeof(float));
									if (index < main_field::con_number)
										point[index] = data_buff;
								}
							}
							else {
								for (size_t index = 0; index < mesh_info.con_number; index++)
									fin.read((char*)&data_buff, sizeof(float));
							}
						}
			}
			if (mesh_info.is_temp_mesh) {
				size_t tNx = main_field::temperature_field.Nx(), tNy = main_field::temperature_field.Ny(),
					tNz = main_field::temperature_field.Nz();
				for (size_t x = 0; x < mesh_info.TNx; x++)
					for (size_t y = 0; y < mesh_info.TNy; y++)
						for (size_t z = 0; z < mesh_info.TNz; z++) {
							if (x < tNx && y < tNy && z < tNz) {
								fin.read((char*)&data_buff, sizeof(float));
								main_field::temperature_field(x, y, z) = data_buff;
							}
							else {
								fin.read((char*)&data_buff, sizeof(float));
							}
						}
			}
			if (mesh_info.is_fluid_mesh) {
				//for (int x = 0; x < mesh_info.FNx; x++)
				//	for (int y = 0; y < mesh_info.FNy; y++)
				//		for (int z = 0; z < mesh_info.FNz; z++) {
				//			if (x < fNx && y < fNy && z < fNz && write_mesh_data::is_fluid_mesh) {
				//				std::vector<float>& point = fluid_field[MESH_INDEX(x, y, z, fNx, fNy)];
				//				for (size_t index = 0; index < mesh_info.dis_func_number; index++) {
				//					fin.read((char*)&data_buff, sizeof(float));
				//					if (index < dis_func_number)
				//						point[index] = data_buff;
				//				}
				//			}
				//			else {
				//				for (size_t index = 0; index < mesh_info.dis_func_number; index++)
				//					fin.read((char*)&data_buff, sizeof(float));
				//			}
				//		}
			}
			fin.close();
			return true;
		}
	}
}