#include "WriteVTS.h"
#include <locale>
namespace pf {

	namespace write_vts {
		void load_vts_func(void(*buff)(std::ofstream& fout)) {
			write_vts::write_vts_list.push_back(buff);
		}
		namespace default_functions {
			void open_vts_file(std::ofstream& fout, std::string tail) {
				std::string fname;
				fname = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "SimData_" + tail + ".vts";
				fout.imbue(std::locale::classic());
				fout.open(fname);
				if (!fout) {
					std::cout << "Failed to write the vtk file..." << std::endl;
					fout.close();
					return;
				}
				fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << '\n';
				fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
				fout << "<StructuredGrid WholeExtent=\""
					<< x_begin << " " << x_end << " "
					<< y_begin << " " << y_end << " "
					<< z_begin << " " << z_end << "\"> " << '\n';
				fout << "<PointData Scalars= \"ScalarData\"  Vectors= \"VectorData\">" << '\n';
			}
			void write_scalar_grains(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "phi2_summary" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							Matrix1D<REAL>& point = main_field::phase_field(i, j, k);
							REAL fix = 0.0;
							for (size_t index = 0; index < main_field::phi_number; index++)
								fix += point[index] * point[index];
							if (std::isnan(fix))
								fout << NaN() << '\n';
							else
								fout << 1.0 - fix << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_scalar_phi_index(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "phi_index" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							Matrix1D<REAL>& point = main_field::phase_field(i, j, k);
							REAL fix = 0.0;
							for (size_t index = 0; index < main_field::phi_number; index++)
								fix += point[index] * index;
							fout << fix << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_scalar_phi_all(std::ofstream& fout) {
				for (size_t pindex = 0; pindex < main_field::phi_number; pindex++) {
					std::string phi_name = "phi_" + std::to_string(pindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << phi_name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i)
								fout << main_field::phase_field(i, j, k)[pindex] << '\n';
					fout << "</DataArray>" << '\n';
				}
			}
			void write_scalar_grad_phi_all(std::ofstream& fout) {
				for (size_t pindex = 0; pindex < main_field::phi_number; pindex++) {
					std::string phi_name = "phi_grad_" + std::to_string(pindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << phi_name <<
						"\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								if (i >= size_t(main_field::phase_field.COMP_X_BGN()) && i <= size_t(main_field::phase_field.COMP_X_END()) ||
									j >= size_t(main_field::phase_field.COMP_Y_BGN()) && j <= size_t(main_field::phase_field.COMP_Y_END()) ||
									k >= size_t(main_field::phase_field.COMP_Z_BGN()) && k <= size_t(main_field::phase_field.COMP_Z_END()))
									fout << (main_field::phase_field(i + 1, j, k)[pindex] - main_field::phase_field(i - 1, j, k)[pindex]) / 2 / mesh_parameters::delt_r << " "
									<< (main_field::phase_field(i, j + 1, k)[pindex] - main_field::phase_field(i, j - 1, k)[pindex]) / 2 / mesh_parameters::delt_r << " "
									<< (main_field::phase_field(i, j, k + 1)[pindex] - main_field::phase_field(i, j, k - 1)[pindex]) / 2 / mesh_parameters::delt_r << '\n';
								else
									fout << 0 << " " << 0 << " " << 0 << '\n';
							}
					fout << "</DataArray>" << '\n';
				}
			}
			void write_scalar_con_all(std::ofstream& fout) {
				for (size_t cindex = 0; cindex < main_field::con_number; cindex++) {
					std::string con_name = "con_" + std::to_string(cindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << con_name <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i)
								fout << main_field::concentration_field(i, j, k)[cindex] << '\n';
					fout << "</DataArray>" << '\n';
				}
			}
			void write_scalar_grad_con_all(std::ofstream& fout) {
				for (size_t cindex = 0; cindex < main_field::con_number; cindex++) {
					std::string con_name = "con_grad_" + std::to_string(cindex);
					fout << "<DataArray type = \"Float64\" Name = \"" << con_name <<
						"\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								if (i >= size_t(main_field::concentration_field.COMP_X_BGN()) && i <= size_t(main_field::concentration_field.COMP_X_END()) ||
									j >= size_t(main_field::concentration_field.COMP_Y_BGN()) && j <= size_t(main_field::concentration_field.COMP_Y_END()) ||
									k >= size_t(main_field::concentration_field.COMP_Z_BGN()) && k <= size_t(main_field::concentration_field.COMP_Z_END()))
									fout << (main_field::concentration_field(i + 1, j, k)[cindex] - main_field::concentration_field(i - 1, j, k)[cindex]) / 2 / mesh_parameters::delt_r << " "
									<< (main_field::concentration_field(i, j + 1, k)[cindex] - main_field::concentration_field(i, j - 1, k)[cindex]) / 2 / mesh_parameters::delt_r << " "
									<< (main_field::concentration_field(i, j, k + 1)[cindex] - main_field::concentration_field(i, j, k - 1)[cindex]) / 2 / mesh_parameters::delt_r << '\n';
								else
									fout << 0 << " " << 0 << " " << 0 << '\n';
							}
					fout << "</DataArray>" << '\n';
				}
			}
			void write_scalar_temperature(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "temp" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << main_field::temperature_field(i, j, k) << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_scalar_grad_temperature(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "temp_grad" <<
					"\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							if (i >= size_t(main_field::temperature_field.COMP_X_BGN()) && i <= size_t(main_field::temperature_field.COMP_X_END()) ||
								j >= size_t(main_field::temperature_field.COMP_Y_BGN()) && j <= size_t(main_field::temperature_field.COMP_Y_END()) ||
								k >= size_t(main_field::temperature_field.COMP_Z_BGN()) && k <= size_t(main_field::temperature_field.COMP_Z_END()))
								fout << (main_field::temperature_field(i + 1, j, k) - main_field::temperature_field(i - 1, j, k)) / 2 / mesh_parameters::delt_r << " "
								<< (main_field::temperature_field(i, j + 1, k) - main_field::temperature_field(i, j - 1, k)) / 2 / mesh_parameters::delt_r << " "
								<< (main_field::temperature_field(i, j, k + 1) - main_field::temperature_field(i, j, k - 1)) / 2 / mesh_parameters::delt_r << '\n';
							else
								fout << 0 << " " << 0 << " " << 0 << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_velocity(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"fluid_velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							LBMPoint& point = external_physical_field::lbm_field(i, j, k);
							fout << point.velocity[0] << " "
								<< point.velocity[1] << " "
								<< point.velocity[2] << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_abs_velocity(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_abs_velocity" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << external_physical_field::lbm_field(i, j, k).velocity.abs() << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_pressure(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_pressure" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << external_physical_field::lbm_field(i, j, k).pressure << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_density(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "fluid_density" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << external_physical_field::lbm_field(i, j, k).mass << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_stress(std::ofstream& fout) {
				vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
				for (int ele_index = 0; ele_index < 6; ele_index++)
				{
					string compname = "\"stress_" + compNameV[ele_index] + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << compname <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								fout << external_physical_field::elastic_field(i, j, k).Stress[ele_index] << '\n';
							}
					fout << "</DataArray>" << '\n';
				}
			}
			void write_strain(std::ofstream& fout) {
				vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
				for (int ele_index = 0; ele_index < 6; ele_index++)
				{
					string compname = "\"strain_" + compNameV[ele_index] + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << compname <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								fout << external_physical_field::elastic_field(i, j, k).Strain[ele_index] << '\n';
							}
					fout << "</DataArray>" << '\n';
				}
			}
			void write_nonElasticStrain(std::ofstream& fout) {
				vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
				for (int ele_index = 0; ele_index < 6; ele_index++)
				{
					string compname = "\"nonElasticStrain_" + compNameV[ele_index] + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << compname <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								fout << external_physical_field::elastic_field(i, j, k).EffectiveEigenStrain[ele_index] << '\n';
							}
					fout << "</DataArray>" << '\n';
				}
			}
			void write_J1(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "stress_J1" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << external_physical_field::elastic_field(i, j, k).Stress.J1() << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_vMises(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "stress_vMises" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << external_physical_field::elastic_field(i, j, k).Stress.Mises() << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void write_plastic_strain(std::ofstream& fout) {
				vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
				for (int ele_index = 0; ele_index < 6; ele_index++)
				{
					string compname = "\"plastic_strain_" + compNameV[ele_index] + "\" ";
					fout << "<DataArray type = \"Float64\" Name = " << compname <<
						"NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
					for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
						for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
							for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
								fout << external_physical_field::plastic_field(i, j, k).PlasticStrain[ele_index] << '\n';
							}
					fout << "</DataArray>" << '\n';
				}
			}
			void write_ave_plastic_strain(std::ofstream& fout) {
				string compname = "\"ave_plastic_strain\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << '\n';
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << external_physical_field::plastic_field(i, j, k).AvePlasticStrain << '\n';
						}
				fout << "</DataArray>" << '\n';
			}
			void close_vts_file(std::ofstream& fout) {
				fout << "</PointData>" << '\n';
				fout << "<Points>" << '\n';
				fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
				for (size_t k = z_begin; k <= z_end; ++k)
					for (size_t j = y_begin; j <= y_end; ++j)
						for (size_t i = x_begin; i <= x_end; ++i) {
							fout << i * mesh_parameters::delt_r << " " << j * mesh_parameters::delt_r << " " << k * mesh_parameters::delt_r << "\n";
						}
				fout << "</DataArray>" << '\n';
				fout << "</Points>" << '\n';
				fout << "</StructuredGrid>" << '\n';
				fout << "</VTKFile>" << '\n';
				fout.close();
			}
		}
		void write_vts_pre_iii() {
			std::ofstream fout;
			// - 
			default_functions::open_vts_file(fout, "step0");
			for (auto writer = write_vts_list.begin(); writer < write_vts_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
		}

		void write_vts_pos_iii() {
			if (output_frequence == 0)
				return;
			if (main_iterator::Current_ITE_step % output_frequence == 0) {
				std::ofstream fout;
				// - 
				default_functions::open_vts_file(fout, "step" + std::to_string(main_iterator::Current_ITE_step));
				for (auto writer = write_vts_list.begin(); writer < write_vts_list.end(); writer++)
					(*writer)(fout);
				default_functions::close_vts_file(fout);
			}
		}

		void init_write_vts() {
			if (infile_reader::read_int_value("Solver.Output.VTS.frequence", output_frequence, true)) {
				if (output_frequence == 0) {
					load_a_new_module(nullptr, nullptr, write_vts_pre_iii,
						nullptr, nullptr, nullptr,
						nullptr, nullptr, nullptr, nullptr);
				}
				else if (output_frequence > 0) {
					load_a_new_module(nullptr, nullptr, write_vts_pre_iii,
						nullptr, nullptr, nullptr,
						nullptr, nullptr, write_vts_pos_iii, nullptr);
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
					load_vts_func(default_functions::write_scalar_phi_all);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.phi_grad_all", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_grad_phi_all);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.phi_index", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_phi_index);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.phi2_summary", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_grains);
			}
			if (main_field::is_con_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.con_all", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_con_all);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.con_grad_all", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_grad_con_all);
			}
			if (main_field::is_temp_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.temp", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_temperature);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.temp_grad", buff, true);
				if (buff)
					load_vts_func(default_functions::write_scalar_grad_temperature);
			}
			if (external_physical_field::is_fluid_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.fluid_velocity", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_velocity);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.fluid_abs_velocity", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_abs_velocity);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.fluid_pressure", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_pressure);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.fluid_density", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_density);
			}
			if (external_physical_field::is_mech_field_on) {
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_stress", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_stress);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_strain", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_strain);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_non_elastic_strain", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_nonElasticStrain);
				buff = false;
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_J1", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_J1);
				InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_vMises", buff, true);
				if (buff)
					write_vts::load_vts_func(default_functions::write_vMises);
				if (external_physical_field::is_mech_plastic_field_on) {
					buff = false;
					InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_plastic_strain", buff, true);
					if (buff)
						write_vts::load_vts_func(default_functions::write_plastic_strain);
					buff = false;
					InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.mech_ave_plastic_strain", buff, true);
					if (buff)
						write_vts::load_vts_func(default_functions::write_ave_plastic_strain);
				}
			}
		}
	}
}
