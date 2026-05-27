#include "Pretreatment.h"
namespace pf {
	namespace pretreatment {
		// -
		namespace functions {
			// -
			void filling() {
				WriteLog("> filling Start:\n");
				for (size_t tar_index = 0; tar_index < target_number; tar_index++) {
					WriteLog("> filling phase index : ");
					for (size_t index = 0; index < target_phi_index[tar_index].size(); index++)
						WriteLog(to_string(target_phi_index[tar_index][index]) + ", ");
					WriteLog("\n");
#pragma omp parallel for
					for (long long x = 0; x < main_field::phase_field.Nx(); x++)
						for (long long y = 0; y < main_field::phase_field.Ny(); y++)
							for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
								Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
								Matrix1D<REAL>& con = main_field::concentration_field(x, y, z);
								REAL& temp = main_field::temperature_field(x, y, z);
								REAL target = 0, anti_target = 0;
								for (size_t index = 0; index < main_field::phi_number; index++) {
									bool is_target = false;
									for (size_t index2 = 0; index2 < target_phi_index[tar_index].size(); index2++)
										if (index == target_phi_index[tar_index][index2])
											is_target = true;
									if (is_target)
										target += phi[index];
									else
										anti_target += phi[index];
								}
								switch (fill_mode) {
								case FILLING_MODE::FM_BULK:
									if (target > SYS_EPSILON_R) {
										if (is_fill_con[tar_index])
											for (size_t index = 0; index < fill_con[tar_index].size(); index++)
												con[fill_con[tar_index][index].first] = fill_con[tar_index][index].second;
										if (is_fill_temp[tar_index])
											temp = fill_temp[tar_index];
									}
									break;
								case FILLING_MODE::FM_FULL_INT:
									if (target > SYS_EPSILON) {
										if (is_fill_con[tar_index])
											for (size_t index = 0; index < fill_con[tar_index].size(); index++)
												con[fill_con[tar_index][index].first] = fill_con[tar_index][index].second;
										if (is_fill_temp[tar_index])
											temp = fill_temp[tar_index];
									}
									break;
								case FILLING_MODE::FM_HALF_INT:
									if (target > 0.5) {
										if (is_fill_con[tar_index])
											for (size_t index = 0; index < fill_con[tar_index].size(); index++)
												con[fill_con[tar_index][index].first] = fill_con[tar_index][index].second;
										if (is_fill_temp[tar_index])
											temp = fill_temp[tar_index];
									}
									break;
								case FILLING_MODE::FM_INTERPOLATION_INT:
									if (target < SYS_EPSILON)
										continue;
									if (anti_target < SYS_EPSILON)
										continue;
									if (is_fill_con[tar_index])
										for (size_t index = 0; index < fill_con[tar_index].size(); index++)
											con[fill_con[tar_index][index].first] = (fill_con[tar_index][index].second * target
												+ con[fill_con[tar_index][index].first] * anti_target) / (target + anti_target);
									if (is_fill_temp[tar_index])
										temp = (fill_temp[tar_index] * target + temp * anti_target) / (target + anti_target);
									break;
								default:
									break;
								}
							}
				}
				WriteLog("> Filling End ! \n");
			}
			// -
			void merging() {
				WriteLog("> Merging Start: \n");
				for (size_t merge_index = 0; merge_index < grains_need_merged.size(); merge_index++) {
					size_t tar_phi_index = grains_need_merged[merge_index][0]; 
					WriteLog("> merging step " + to_string(merge_index) + " : target phi index is " + to_string(tar_phi_index) + " , aim grains are ");
					for (size_t index = 1; index < grains_need_merged[merge_index].size(); index++) {
						size_t phi_index = grains_need_merged[merge_index][index]; 
						WriteLog(to_string(phi_index) + ", ");
					}
					WriteLog("\n");
#pragma omp parallel for
					for (long long x = 0; x < main_field::phase_field.Nx(); x++)
						for (long long y = 0; y < main_field::phase_field.Ny(); y++)
							for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
								Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
								for (size_t index = 1; index < grains_need_merged[merge_index].size(); index++) {
									size_t phi_index = grains_need_merged[merge_index][index];
									phi[tar_phi_index] += phi[phi_index];
									phi[phi_index] = 0;
								}
							}
				}
				WriteLog("> Merging End! \n");
			}
			// - 
			void auto_merging() {
				WriteLog("> Auto-Merging Start: \n");
				Matrix2D<int> grain_contacts(main_field::phi_number, main_field::phi_number, 0);
				std::vector<int> grains_need_merged(main_field::phi_number, 0);
				WriteLog("> checking grains distribution ...\n");
				for (long long x = 0; x < main_field::phase_field.Nx(); x++)
					for (long long y = 0; y < main_field::phase_field.Ny(); y++)
						for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							for (size_t index = 0; index < main_field::phi_number; index++)
								if (phi[index] > SYS_EPSILON) {
									grains_need_merged[index] = 1;
									// -
									for (size_t index2 = index + 1; index2 < main_field::phi_number; index2++)
										if (phi[index2] > SYS_EPSILON) {
											grain_contacts(index, index2) = 1;
											grain_contacts(index2, index) = 1;
										}
								}
						}
				for (size_t index = 0; index < main_field::phi_number; index++) {
					bool is_merge = false;
					for (size_t index2 = 0; index2 < auto_merged_grains.size(); index2++)
						if (index == auto_merged_grains[index2])
							is_merge = true;
					if (!is_merge)
						grains_need_merged[index] = 0;
				}
				for (size_t map_index = 0; map_index < main_field::phi_number - 1; map_index++) {
					if (grains_need_merged[map_index] == 0)
						continue;
					std::vector<size_t> merge_grains_buff = { map_index };
					std::string out = "> merging : " + to_string(map_index) + " <-- ";
					for (size_t merge_index = map_index + 1; merge_index < main_field::phi_number; merge_index++) {
						if (grains_need_merged[merge_index] == 0)
							continue;
						// check contact
						bool is_contact = false;
						for (size_t index = 0; index < merge_grains_buff.size(); index++)
							if (grain_contacts(merge_index, merge_grains_buff[index]) == 1)
								is_contact = true;
						if (!is_contact) {
							merge_grains_buff.push_back(merge_index);
							grains_need_merged[map_index] = 0;
							out += to_string(merge_index) + ", ";
						}
					}
					out += "\n";
					WriteLog(out);
					if (merge_grains_buff.size() > 1) {
						size_t map_index = merge_grains_buff[0], merge_number = merge_grains_buff.size();
#pragma omp parallel for
						for (long long x = 0; x < main_field::phase_field.Nx(); x++)
							for (long long y = 0; y < main_field::phase_field.Ny(); y++)
								for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
									Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
									for (size_t index = 1; index < merge_number; index++) {
										phi[map_index] += phi[merge_grains_buff[index]];
										phi[merge_grains_buff[index]] = 0.0;
									}
								}
					}
				}
				WriteLog("> Auto-Merging End! \n");
			}
			// - 
			void reordering() {
				WriteLog("> Re-Ordering Grains Start: \n");
				std::vector<int> grains_empty(main_field::phi_number, 1);
				WriteLog("> checking grains distribution ...\n");
				for (long long x = 0; x < main_field::phase_field.Nx(); x++)
					for (long long y = 0; y < main_field::phase_field.Ny(); y++)
						for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							for (size_t index = 0; index < main_field::phi_number; index++)
								if (phi[index] > SYS_EPSILON) {
									grains_empty[index] = 0;
								}
						}
				for (size_t index = 0; index < main_field::phi_number; index++) {
					if (grains_empty[index] == 0) // jump unempty grain
						continue;
					for (size_t index2 = index + 1; index2 < main_field::phi_number; index2++) {
						if (grains_empty[index2] == 1) // jump empty grain
							continue;
						// index is empty grain, and index2 is unempty grain, move index2 to index
#pragma omp parallel for
						for (long long x = 0; x < main_field::phase_field.Nx(); x++)
							for (long long y = 0; y < main_field::phase_field.Ny(); y++)
								for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
									Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
									phi[index] = phi[index2];
									phi[index2] = 0.0;
								}
						// index2 is empty grain, and index is unempty grain
						grains_empty[index] = 0;
						grains_empty[index2] = 1;
						WriteLog("> reordering : " + to_string(index) + " <-- " + to_string(index2) + "\n");
					}
				}
				WriteLog("> Re-Ordering Grains End! \n");
			}
		}
		// -
		void init_pretreatment_pre_i() {
			if (is_merge_grains)
				functions::merging();
			if (is_auto_merge_grains)
				functions::auto_merging();
			if (is_reorder_grains)
				functions::reordering();
			if (is_fill_grains)
				functions::filling();
		}
		// -
		void init_pretreatment() {
			WriteDebugFile("# Preprocess.merge_grains = [(phi_index_0, phi_index_1, ... ), ... ] \n");
			std::string merge_phis_key = "Preprocess.merge_grains", merge_phis_input = "[()]";
			if (infile_reader::read_string_value(merge_phis_key, merge_phis_input, true)) {
				is_merge_grains = true;
				std::vector<std::vector<input_value>> merge_phis_value = InputFileReader::get_instance()->
					trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, merge_phis_key, merge_phis_input, true);
				for (size_t index = 0; index < merge_phis_value.size(); index++) {
					std::vector<size_t> merge_grains;
					for (size_t index2 = 0; index2 < merge_phis_value[index].size(); index2++)
						merge_grains.push_back(size_t(merge_phis_value[index][index2].int_value));
					grains_need_merged.push_back(merge_grains);
				}
			}

			// -
			WriteDebugFile("# Preprocess.auto_merge_grains = (phi_index_0, phi_index_1, ... ) \n");
			std::string auto_merge_key = "Preprocess.auto_merge_grains", auto_merge_input = "()";
			if (infile_reader::read_string_value(auto_merge_key, auto_merge_input, true)) {
				is_auto_merge_grains = true;
				auto_merged_grains.clear();
				std::vector<input_value> auto_merge_value = InputFileReader::get_instance()->
					trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, auto_merge_key, auto_merge_input, true);
				for (size_t index = 0; index < auto_merge_value.size(); index++) {
					bool is_already_set = false;
					for (size_t index2 = 0; index2 < auto_merged_grains.size(); index2++)
						if (auto_merge_value[index].int_value == auto_merged_grains[index2])
							is_already_set = true;
					if (!is_already_set)
						auto_merged_grains.push_back(auto_merge_value[index].int_value);
				}
			}

			// should be concerned that matrix phases problem
			infile_reader::read_bool_value("Preprocess.reordering_grains", is_reorder_grains, true);

			target_number = 0;
			WriteDebugFile("# Preprocess.fill_times = filling times \n");
			if (infile_reader::read_int_value("Preprocess.fill_times", target_number, true)) {
				is_fill_grains = true;
				target_phi_index.resize(target_number);
				is_fill_con.resize(target_number, false);
				fill_con.resize(target_number);
				is_fill_temp.resize(target_number, false);
				fill_temp.resize(target_number, 0);
				for (size_t fill_index = 0; fill_index < target_number; fill_index++) {
					std::string keys = "Preprocess.fill_time_" + to_string(fill_index);
					// - 
					WriteDebugFile("# " + keys + ".grains_index = (phi_index_0, phi_index_1, ... ) \n");
					std::string fill_phis_key = keys + ".grains_index", fill_phis_input = "()";
					if (infile_reader::read_string_value(fill_phis_key, fill_phis_input, true)) {
						std::vector<input_value> fill_phis_value = InputFileReader::get_instance()->
							trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, fill_phis_key, fill_phis_input, true);
						for (auto phi_index = fill_phis_value.begin(); phi_index < fill_phis_value.end(); phi_index++) {
							if (phi_index->int_value >= main_field::phi_number || phi_index->int_value < 0) {
								WriteDebugFile("> ERROR " + keys + ".grains_index : phi_index should >= 0 and < " 
									+ to_string(main_field::phi_number) + " ! \n");
								exit(0);
							}
							target_phi_index[fill_index].push_back(phi_index->int_value);
						}
					}
					// - 
					WriteDebugFile("# " + keys + ".grains_con = [(con_index_0, con_value_0), (con_index_1, con_value_1), ... ] \n");
					std::string fill_con_key = keys + ".grains_con", fill_con_input = "[()]";
					if (infile_reader::read_string_value(fill_con_key, fill_con_input, true)) {
						is_fill_con[fill_index] = true;
						std::vector<std::vector<input_value>> fill_con_value = InputFileReader::get_instance()->
							trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_INT, InputValueType::IVType_REAL }, fill_con_key, fill_con_input, true);
						for (size_t index = 0; index < fill_con_value.size(); index++)
							fill_con[fill_index].push_back({ fill_con_value[index][0].int_value, fill_con_value[index][1].REAL_value });
					}
					// - 
					WriteDebugFile("# " + keys + ".grains_temp = temp_value \n");
					if (infile_reader::read_real_value(keys + ".grains_temp", fill_temp[fill_index], true))
						is_fill_temp[fill_index] = true;
				}
			}
			load_a_new_module(init_pretreatment_pre_i, nullptr, nullptr,
				nullptr, nullptr, nullptr,
				nullptr, nullptr, nullptr, dinit);
		}
		// - 
		void dinit() {

		}
	}
}