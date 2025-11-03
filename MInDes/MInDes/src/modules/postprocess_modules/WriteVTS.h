#pragma once
#include "../model_modules/Model_Params.h"
#include "../input_modules/ioFiles_Params.h"
#include "../../MainIterator_Params.h"
#include "../Module.h"
namespace pf {

	namespace write_vts {
		inline size_t output_frequence = 0;
		inline bool is_show_with_boundary = false;
		inline size_t x_begin = 0;
		inline size_t y_begin = 0;
		inline size_t z_begin = 0;
		inline size_t x_end = 0;
		inline size_t y_end = 0;
		inline size_t z_end = 0;
		inline std::vector<void(*)(std::ofstream& fout)> write_vts_scalar_list;
		inline std::vector<void(*)(std::ofstream& fout)> write_vts_vector_list;
		inline void load_vts_scalar_func(void(*buff)(std::ofstream& fout)) {
			write_vts::write_vts_scalar_list.push_back(buff);
		}
		inline void load_vts_vector_func(void(*buff)(std::ofstream& fout)) {
			write_vts::write_vts_vector_list.push_back(buff);
		}
		namespace default_functions {
			inline void open_vts_scalar_file(std::ofstream& fout, std::string tail) {
				std::string fname;
				fname = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "scalar_variables_" + tail + ".vts";
				fout.open(fname);
				if (!fout) {
					std::cout << "Failed to write the vtk file..." << std::endl;
					fout.close();
					return;
				}
				fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << std::endl;
				fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
				fout << "<StructuredGrid WholeExtent=\""
					<< x_begin << " " << x_end << " "
					<< y_begin << " " << y_end << " "
					<< z_begin << " " << z_end << "\"> " << std::endl;
				fout << "<PointData Scalars= \"ScalarData\">" << std::endl;
			}
			inline void open_vts_vec3_file(std::ofstream& fout, std::string tail) {
				std::string fname;
				fname = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "vec3_variables_" + tail + ".vts";
				fout.open(fname);
				if (!fout) {
					std::cout << "Failed to write the vtk file..." << std::endl;
					fout.close();
					return;
				}
				fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << std::endl;
				fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
				fout << "<StructuredGrid WholeExtent=\""
					<< 0 << " " << x_end - x_begin + 1 << " "
					<< 0 << " " << y_end - y_begin + 1 << " "
					<< 0 << " " << z_end - z_begin + 1 << "\"> " << std::endl;
				fout << "<PointData  Vectors= \"VectorData\">" << std::endl;
			}
			inline void write_scalar_grains(ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "grains" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							std::vector<REAL>& point = main_field::phase_field(i, j, k);
							REAL fix = 0.0;
							for (size_t index = 0; index < main_field::phi_number; index++)
								fix += point[index] * point[index];
							if (std::isnan(fix))
								fout << NaN() << endl;
							else
								fout << 1.0 - fix << endl;
						}
				fout << "</DataArray>" << endl;
			}
			inline void write_scalar_phi_index(ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "phi_index" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							std::vector<REAL>& point = main_field::phase_field(i, j, k);
							REAL fix = 0.0;
							for (size_t index = 0; index < main_field::phi_number; index++)
								fix += point[index] * index;
							fout << fix << endl;
						}
				fout << "</DataArray>" << endl;
			}
			inline void write_scalar_phi_all(ofstream& fout) {
				for (size_t pindex = 0; pindex < main_field::phi_number; pindex++) {
					string phi_name = "phi_" + to_string(pindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << phi_name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i)
								fout << main_field::phase_field(i, j, k)[pindex] << endl;
					fout << "</DataArray>" << endl;
				}
			}
			inline void write_scalar_con_all(ofstream& fout) {
				for (size_t cindex = 0; cindex < main_field::con_number; cindex++) {
					string con_name = "con_" + to_string(cindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << con_name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i)
								fout << main_field::concentration_field(i, j, k)[cindex] << endl;
					fout << "</DataArray>" << endl;
				}
			}
			inline void write_scalar_temperature(ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "temp" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << main_field::temperature_field(i, j, k) << endl;
						}
				fout << "</DataArray>" << endl;
			}
			inline void close_vts_file(std::ofstream& fout) {
				fout << "</PointData>" << std::endl;
				fout << "<Points>" << std::endl;
				fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
				for (size_t k = z_begin; k <= z_end; ++k)
					for (size_t j = y_begin; j <= y_end; ++j)
						for (size_t i = x_begin; i <= x_end; ++i) {
							fout << i * mesh_parameters::delt_r << " " << j * mesh_parameters::delt_r << " " << k * mesh_parameters::delt_r << "\n";
						}
				fout << "</DataArray>" << std::endl;
				fout << "</Points>" << std::endl;
				fout << "</StructuredGrid>" << std::endl;
				fout << "</VTKFile>" << std::endl;
				fout.close();
			}
		}
		inline void write_vts_pre_iii() {
			ofstream fout;
			// - 
			default_functions::open_vts_scalar_file(fout, "step0");
			for (auto writer = write_vts_scalar_list.begin(); writer < write_vts_scalar_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
			// - 
			default_functions::open_vts_vec3_file(fout, "step0");
			for (auto writer = write_vts_vector_list.begin(); writer < write_vts_vector_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
		}

		inline void write_vts_pos_iii() {
			ofstream fout;
			// - 
			default_functions::open_vts_scalar_file(fout, "step" + to_string(main_iterator::Current_ITE_step));
			for (auto writer = write_vts_scalar_list.begin(); writer < write_vts_scalar_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
			// - 
			default_functions::open_vts_vec3_file(fout, "step" + to_string(main_iterator::Current_ITE_step));
			for (auto writer = write_vts_vector_list.begin(); writer < write_vts_vector_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
		}

		inline void init() {
			if (infile_reader::read_int_value("Solver.Output.VTS.frequence", output_frequence, true)) {
				if (output_frequence == 0) {
					load_a_new_module(default_module_function, default_module_function, write_vts_pre_iii,
						default_module_function, default_module_function, default_module_function,
						default_module_function, default_module_function, default_module_function, default_module_function);
				}
				else if (output_frequence > 0) {
					load_a_new_module(default_module_function, default_module_function, write_vts_pre_iii,
						default_module_function, default_module_function, default_module_function,
						default_module_function, default_module_function, write_vts_pos_iii, default_module_function);
				}
			}
			InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.with_boundary", is_show_with_boundary, true);
			if (is_show_with_boundary) {
				x_begin = 0;
				y_begin = 0;
				z_begin = 0;
				if (main_field::is_phi_field_on) {
					x_end = main_field::phase_field.Nx() - 1;
					y_end = main_field::phase_field.Ny() - 1;
					z_end = main_field::phase_field.Nz() - 1;
				}
				else if (main_field::is_con_field_on) {
					x_end = main_field::concentration_field.Nx() - 1;
					y_end = main_field::concentration_field.Ny() - 1;
					z_end = main_field::concentration_field.Nz() - 1;
				}
				else if (main_field::is_temp_field_on) {
					x_end = main_field::temperature_field.Nx() - 1;
					y_end = main_field::temperature_field.Ny() - 1;
					z_end = main_field::temperature_field.Nz() - 1;
				}
			}
			else {
				if (main_field::is_phi_field_on) {
					x_begin = main_field::phase_field.COMP_X_BGN();
					y_begin = main_field::phase_field.COMP_Y_BGN();
					z_begin = main_field::phase_field.COMP_Z_BGN();
					x_end = main_field::phase_field.COMP_X_END();
					y_end = main_field::phase_field.COMP_Y_END();
					z_end = main_field::phase_field.COMP_Z_END();
				}
				else if (main_field::is_con_field_on) {
					x_begin = main_field::concentration_field.COMP_X_BGN();
					y_begin = main_field::concentration_field.COMP_Y_BGN();
					z_begin = main_field::concentration_field.COMP_Z_BGN();
					x_end = main_field::concentration_field.COMP_X_END();
					y_end = main_field::concentration_field.COMP_Y_END();
					z_end = main_field::concentration_field.COMP_Z_END();
				}
				else if (main_field::is_temp_field_on) {
					x_begin = main_field::temperature_field.COMP_X_BGN();
					y_begin = main_field::temperature_field.COMP_Y_BGN();
					z_begin = main_field::temperature_field.COMP_Z_BGN();
					x_end = main_field::temperature_field.COMP_X_END();
					y_end = main_field::temperature_field.COMP_Y_END();
					z_end = main_field::temperature_field.COMP_Z_END();
				}
			}
			bool buff = false;
			if (main_field::is_phi_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.phi_all", buff, true);
				if (buff)
					load_vts_scalar_func(default_functions::write_scalar_phi_all);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.phi_index", buff, true);
				if (buff)
					load_vts_scalar_func(default_functions::write_scalar_phi_index);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.grains", buff, true);
				if (buff)
					load_vts_scalar_func(default_functions::write_scalar_grains);
			}
			if (main_field::is_con_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.con_all", buff, true);
				if (buff)
					load_vts_scalar_func(default_functions::write_scalar_con_all);
			}
			if (main_field::is_temp_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.temp", buff, true);
				if (buff)
					load_vts_scalar_func(default_functions::write_scalar_temperature);
			}
		}
	}

}