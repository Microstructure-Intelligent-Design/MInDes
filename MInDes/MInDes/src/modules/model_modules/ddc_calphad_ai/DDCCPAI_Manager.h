#pragma once
#include "../../Module.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "DDCCPAI_Params.h"
#include "DDCCPAI_Functions.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		/*
			模型来自文献：
			  原创
		*/
		inline void init_model_modules() {
			// - check
			if (main_field::phi_number == 0) {
				string error_report = "> ERROR, the PCT command line settings do not meet the requirements : PCT = (N>0,K,true) !\n";
				WriteLog(error_report);
				std::exit(0);
			}
			if (main_field::is_temp_field_on == false) {
				string error_report = "> ERROR, the PCT command line settings do not meet the requirements : PCT = (N>0,K,true) !\n";
				WriteLog(error_report);
				std::exit(0);
			}
			WriteLog("> \n");
			WriteLog("> Simulation Model - Data Driven Complex - couple CALPHAD database - accelerated by AI technology - is Activated ! \n");
			WriteLog(u8"> Reference : 原创 \n");
			WriteLog(u8"> DOI: 暂无 \n");
			WriteLog("> \n");
			WriteDebugFile("# Funcs.phi  : d phi_i / d t = XXX \n");
			WriteDebugFile("# Funcs.con  : d c_i / d t   = XXX \n");
			WriteDebugFile("# Funcs.temp : d T / d t     = XXX \n");
			// - 
			PhiProperties::instance().init();
			GrainsOrientations::instance().init();
			ConRegions::instance().init();
			// - 
			Matrix3x3 matrix;
			matrix.set_to_unity();
			parameters::grain_rotation_matrix.resize(main_field::phi_number, matrix);
			for (size_t index = 0; index < main_field::phi_number; index++)
				parameters::grain_rotation_matrix[index] = GrainsOrientations::instance().RotationMatrix(index);
			// - 
			int difference_method = DifferenceMethod::SEVEN_POINT;
			WriteDebugFile("# Model.DDCCPAI.difference_method : 0 - SEVEN_POINT , 1 - NINETEEN_POINT \n");
			infile_reader::read_int_value("Model.DDCCPAI.difference_method", difference_method, true);
			parameters::diff_method = DifferenceMethod(difference_method);
			if (parameters::diff_method == DifferenceMethod::SEVEN_POINT) {
				phase_field_functions::interphase_gradient_lapace_calculation 
					= phase_field_functions::interphase_gradient_lapace_calculation_7P;
				phase_field_functions::currentFlag = phase_field_functions::currentFlag_7P;
				phase_field_functions::upgradeFlag = phase_field_functions::upgradeFlag_7P;
			}
			else if (parameters::diff_method == DifferenceMethod::NINETEEN_POINT) {
				phase_field_functions::interphase_gradient_lapace_calculation
					= phase_field_functions::interphase_gradient_lapace_calculation_19P;
				phase_field_functions::currentFlag = phase_field_functions::currentFlag_19P;
				phase_field_functions::upgradeFlag = phase_field_functions::upgradeFlag_19P;
			}
			infile_reader::read_bool_value("Model.DDCCPAI.Phi.is_normalize", parameters::is_phi_normalized, true);
			infile_reader::read_bool_value("Model.DDCCPAI.Con.is_normalize", parameters::is_con_normalized, true);
			infile_reader::read_bool_value("Model.DDCCPAI.Temp.is_normalize", parameters::is_temp_normalized, true);
			// - statistic
			WriteDebugFile("# Model.DDCCPAI.PhiConTemp.statistic  = ( STA_KEY , ... ) \n");
			WriteDebugFile("#       STA_KEY = \n");
			WriteDebugFile("#                 " + parameters::statistic_phase_volume.first + "  : " + parameters::statistic_phase_volume.second + "\n");
			WriteDebugFile("#                 " + parameters::statistic_region_volume.first + " : " + parameters::statistic_region_volume.second + "\n");
			WriteDebugFile("#                 " + parameters::statistic_average_con.first + "    : " + parameters::statistic_average_con.second + "\n");
			WriteDebugFile("#                 " + parameters::statistic_total_con.first + "    : " + parameters::statistic_total_con.second + "\n");
			std::string statistic_key = "Model.DDCCPAI.PhiConTemp.statistic", statistic_input = "()";
			if (infile_reader::read_string_value(statistic_key, statistic_input, true)) {
				std::vector<input_value> statistic_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value
					(InputValueType::IVType_STRING, statistic_key, statistic_input, true);
				for (size_t index = 0; index < statistic_value.size(); index++) {
					if (statistic_value[index].string_value == parameters::statistic_phase_volume.first)
						parameters::is_phase_volume_statistic = true;
					else if (statistic_value[index].string_value == parameters::statistic_region_volume.first)
						parameters::is_region_volume_statistic = true;
					else if (statistic_value[index].string_value == parameters::statistic_average_con.first)
						parameters::is_average_con_statistic = true;
					else if (statistic_value[index].string_value == parameters::statistic_total_con.first)
						parameters::is_total_con_statistic = true;
				}
			}
			// ==========================================================================================================================
			// - init phi equation
			parameters::PHI_ACC_NUMBER = main_field::phi_number;
			parameters::PAIRWISE_ACC_STOP = main_field::phi_number;
			infile_reader::read_int_value("Model.DDCCPAI.Phi.PairWiseAcc.container_size", parameters::PHI_ACC_NUMBER, true);
			// - 
			WriteDebugFile("# Model.DDCCPAI.Phi.IntMobility.const  = Lij_value \n");
			WriteDebugFile("#                              .matrix = [(phi_index_i, phi_index_j, Lij_value), ... ] \n");
			WriteDebugFile("#                              .block  = [(phi_index_begin, phi_index_end, Lij_value), ... ] \n");
			std::string matrix_key = "Model.DDCCPAI.Phi.IntMobility.matrix", matrix_input = "[()]", 
				block_key = "Model.DDCCPAI.Phi.IntMobility.block", block_input = "[()]";
			REAL const_input = 0;
			parameters::Lij.resize(main_field::phi_number, main_field::phi_number, 0);
			phase_field_functions::_Lij = phase_field_functions::Lij;
			if (infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.const", const_input, true)) {
				WriteDebugFile("# Model.DDCCPAI.Phi.IntMobility.const.block       = [(phi_index_1, phi_index_2, ... ), ... ] \n");
				WriteDebugFile("#                                    .cross_block = {[(phi_alpha_index_1, ... ), (phi_beta_index_1, ... )], ... } \n");
				std::string const_block_key = "Model.DDCCPAI.Phi.IntMobility.const.block", 
					const_cross_block_key = "Model.DDCCPAI.Phi.IntMobility.const.cross_block", const_block_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(const_block_key, const_block_input, true)) {
					std::vector<std::vector<input_value>> const_block_value = 
						InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
					for (int block_index = 0; block_index < const_block_value.size(); block_index++)
						for (auto phi_a_index = const_block_value[block_index].begin(); phi_a_index < const_block_value[block_index].end(); phi_a_index++)
							for (auto phi_b_index = phi_a_index + 1; phi_b_index < const_block_value[block_index].end(); phi_b_index++) {
								parameters::Lij(phi_a_index->int_value, phi_b_index->int_value) = const_input;
								parameters::Lij(phi_b_index->int_value, phi_a_index->int_value) = const_input;
							}
				}
				else if (InputFileReader::get_instance()->read_string_value(const_cross_block_key, const_block_input, true)) {
					std::vector<std::vector<std::vector<input_value>>> const_block_value = 
						InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
					for (int block_index = 0; block_index < const_block_value.size(); block_index++)
						for (auto phi_a_index = const_block_value[block_index][0].begin(); phi_a_index < const_block_value[block_index][0].end(); phi_a_index++)
							for (auto phi_b_index = const_block_value[block_index][1].begin(); phi_b_index < const_block_value[block_index][1].end(); phi_b_index++) {
								parameters::Lij(phi_a_index->int_value, phi_b_index->int_value) = const_input;
								parameters::Lij(phi_b_index->int_value, phi_a_index->int_value) = const_input;
							}
				}
				else {
					for (size_t phi_a_index = 0; phi_a_index < main_field::phi_number; phi_a_index++)
						for (size_t phi_b_index = 0; phi_b_index < main_field::phi_number; phi_b_index++)
							parameters::Lij(phi_a_index, phi_b_index) = const_input;
				}
			}
			if (infile_reader::read_string_value(matrix_key, matrix_input, true)) {
				std::vector<std::vector<input_value>> matrix_value = 
					InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_INT , InputValueType::IVType_INT , InputValueType::IVType_REAL }, matrix_key, matrix_input, true);
				for (int index = 0; index < matrix_value.size(); index++) {
					parameters::Lij(matrix_value[index][0].int_value, matrix_value[index][1].int_value) = matrix_value[index][2].REAL_value;
					parameters::Lij(matrix_value[index][1].int_value, matrix_value[index][0].int_value) = matrix_value[index][2].REAL_value;
				}
			}
			if (infile_reader::read_string_value(block_key, block_input, true)) {
				std::vector<std::vector<input_value>> block_value = 
					InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_INT , InputValueType::IVType_INT , InputValueType::IVType_REAL }, block_key, block_input, true);
				for (int index = 0; index < block_value.size(); index++)
					for (int phi_a_index = block_value[index][0].int_value; phi_a_index <= block_value[index][1].int_value; phi_a_index++)
						for (int phi_b_index = phi_a_index + 1; phi_b_index <= block_value[index][1].int_value; phi_b_index++) {
							parameters::Lij(phi_a_index, phi_b_index) = block_value[index][2].REAL_value;
							parameters::Lij(phi_b_index, phi_a_index) = block_value[index][2].REAL_value;
						}
			}
			// - Temp
			bool is_intMob_temp = false;
			if (main_field::is_temp_field_on) {
				WriteDebugFile("# IntMobility(T) = IntMobility_0 * exp( - IntMobilityQ / R / T) \n");
				WriteDebugFile("# Model.DDCCPAI.Phi.IntMobilityQ.const  = Qij_value \n");
				WriteDebugFile("#                               .matrix = [(phi_index_i, phi_index_j, Qij_value), ... ] \n");
				WriteDebugFile("#                               .block  = [(phi_index_begin, phi_index_end, Qij_value), ... ] \n");
				std::string matrix_key = "Model.DDCCPAI.Phi.IntMobilityQ.matrix", matrix_input = "[()]", 
					block_key = "Model.DDCCPAI.Phi.IntMobilityQ.block", block_input = "[()]";
				const_input = 0;
				if (infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobilityQ.const", const_input, true)) {
					WriteDebugFile("# Model.DDCCPAI.Phi.IntMobilityQ.const.block       = [(phi_index_1, phi_index_2, ... ), ... ] \n");
					WriteDebugFile("#                                     .cross_block = {[(phi_alpha_index_1, ... ), (phi_beta_index_1, ... )], ... } \n");
					std::string const_block_key = "Model.DDCCPAI.Phi.IntMobilityQ.const.block",
						const_cross_block_key = "Model.DDCCPAI.Phi.IntMobilityQ.const.cross_block", const_block_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(const_block_key, const_block_input, true)) {
						std::vector<std::vector<input_value>> const_block_value =
							InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
						for (int block_index = 0; block_index < const_block_value.size(); block_index++)
							for (auto phi_a_index = const_block_value[block_index].begin(); phi_a_index < const_block_value[block_index].end(); phi_a_index++)
								for (auto phi_b_index = phi_a_index + 1; phi_b_index < const_block_value[block_index].end(); phi_b_index++) {
									parameters::Qij(phi_a_index->int_value, phi_b_index->int_value) = const_input;
									parameters::Qij(phi_b_index->int_value, phi_a_index->int_value) = const_input;
								}
					}
					else if (InputFileReader::get_instance()->read_string_value(const_cross_block_key, const_block_input, true)) {
						std::vector<std::vector<std::vector<input_value>>> const_block_value =
							InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
						for (int block_index = 0; block_index < const_block_value.size(); block_index++)
							for (auto phi_a_index = const_block_value[block_index][0].begin(); phi_a_index < const_block_value[block_index][0].end(); phi_a_index++)
								for (auto phi_b_index = const_block_value[block_index][1].begin(); phi_b_index < const_block_value[block_index][1].end(); phi_b_index++) {
									parameters::Qij(phi_a_index->int_value, phi_b_index->int_value) = const_input;
									parameters::Qij(phi_b_index->int_value, phi_a_index->int_value) = const_input;
								}
					}
					else {
						for (size_t phi_a_index = 0; phi_a_index < main_field::phi_number; phi_a_index++)
							for (size_t phi_b_index = 0; phi_b_index < main_field::phi_number; phi_b_index++)
								parameters::Qij(phi_a_index, phi_b_index) = const_input;
					}
					is_intMob_temp = true;
				}
				if (infile_reader::read_string_value(matrix_key, matrix_input, true)) {
					std::vector<std::vector<input_value>> matrix_value =
						InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_INT , InputValueType::IVType_INT , InputValueType::IVType_REAL }, matrix_key, matrix_input, true);
					for (size_t index = 0; index < matrix_value.size(); index++) {
						parameters::Qij(matrix_value[index][0].int_value, matrix_value[index][1].int_value) = matrix_value[index][2].REAL_value;
						parameters::Qij(matrix_value[index][1].int_value, matrix_value[index][0].int_value) = matrix_value[index][2].REAL_value;
					}
					is_intMob_temp = true;
				}
				if (infile_reader::read_string_value(block_key, block_input, true)) {
					std::vector<std::vector<input_value>> block_value =
						InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_INT , InputValueType::IVType_INT , InputValueType::IVType_REAL }, block_key, block_input, true);
					for (int index = 0; index < block_value.size(); index++)
						for (int phi_a_index = block_value[index][0].int_value; phi_a_index <= block_value[index][1].int_value; phi_a_index++)
							for (int phi_b_index = phi_a_index + 1; phi_b_index <= block_value[index][1].int_value; phi_b_index++) {
								parameters::Qij(phi_a_index, phi_b_index) = block_value[index][2].REAL_value;
								parameters::Qij(phi_b_index, phi_a_index) = block_value[index][2].REAL_value;
							}
					is_intMob_temp = true;
				}
				if (is_intMob_temp) {
					parameters::Qij.resize(main_field::phi_number, main_field::phi_number, 0);
					phase_field_functions::_Lij = phase_field_functions::Lij_temp;
				}
			}
			// - anisotropic
			WriteDebugFile("# Model.DDCCPAI.Phi.IntMobility.Anisotropic.property = AnisotropicModelIndex \n");
			WriteDebugFile("#          AnisotropicModelIndex : 0 - ISO, 1 - CUBIC, 2 - HEX_BOETTGER, 3 - HEX_SUN, 4 - HEX_YANG, 5 - DENDRITE_YANG \n");
			int intMobAniso_model = 0;
			if (infile_reader::read_int_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.property", intMobAniso_model, true)) {
				parameters::intMobAniso_model = parameters::Int_Mobility_Anisotropic(intMobAniso_model);
				if (parameters::intMobAniso_model == parameters::Int_Mobility_Anisotropic::IMA_CUBIC) {
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_1", parameters::intMobAniso_param1, true);
					if (is_intMob_temp)
						phase_field_functions::_Lij = phase_field_functions::Lij_temp_cubic;
					else
						phase_field_functions::_Lij = phase_field_functions::Lij_cubic;
				}
				else if (parameters::intMobAniso_model == parameters::Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_1", parameters::intMobAniso_param1, true);
					if (is_intMob_temp)
						phase_field_functions::_Lij = phase_field_functions::Lij_temp_hex_boettger;
					else
						phase_field_functions::_Lij = phase_field_functions::Lij_hex_boettger;
				}
				else if (parameters::intMobAniso_model == parameters::Int_Mobility_Anisotropic::IMA_HEX_SUN) {
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_1", parameters::intMobAniso_param1, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_2", parameters::intMobAniso_param2, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_3", parameters::intMobAniso_param3, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_4", parameters::intMobAniso_param4, true);
					if (is_intMob_temp)
						phase_field_functions::_Lij = phase_field_functions::Lij_temp_hex_sun;
					else
						phase_field_functions::_Lij = phase_field_functions::Lij_hex_sun;
				}
				else if (parameters::intMobAniso_model == parameters::Int_Mobility_Anisotropic::IMA_HEX_YANG) {
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_1", parameters::intMobAniso_param1, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_2", parameters::intMobAniso_param2, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_3", parameters::intMobAniso_param3, true);
					if (is_intMob_temp)
						phase_field_functions::_Lij = phase_field_functions::Lij_temp_hex_yang;
					else
						phase_field_functions::_Lij = phase_field_functions::Lij_hex_yang;
				}
				else if (parameters::intMobAniso_model == parameters::Int_Mobility_Anisotropic::IMA_DENDRITE_YANG) {
					infile_reader::read_real_value("Model.DDCCPAI.Phi.IntMobility.Anisotropic.parameter_1", parameters::intMobAniso_param1, true);
					if (is_intMob_temp)
						phase_field_functions::_Lij = phase_field_functions::Lij_temp_dendrite_yang;
					else
						phase_field_functions::_Lij = phase_field_functions::Lij_dendrite_yang;
				}
			}
			// - interface energy
			WriteDebugFile("# Model.DDCCPAI.Phi.InterfaceEnergy.int_gradient : 0 - Steinbach_1996 , 1 - Steinbach_1999 , 2 - Steinbach_G2009 \n");
			int int_gradient = parameters::interface_gradient;
			infile_reader::read_int_value("Model.DDCCPAI.Phi.InterfaceEnergy.int_gradient", int_gradient, true);
			int int_potential = parameters::interface_potential;
			WriteDebugFile("# Model.DDCCPAI.Phi.InterfaceEnergy.int_potential : 0 - Nestler_Well , 1 - Nestler_Obstacle , 2 - Steinbach_P2009 \n");
			infile_reader::read_int_value("Model.DDCCPAI.Phi.InterfaceEnergy.int_potential", int_potential, true);
			parameters::interface_gradient = Int_Gradient(int_gradient);
			parameters::interface_potential = Int_Potential(int_potential);

			infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.int_width", parameters::interface_width, true);

			if (parameters::interface_gradient == Int_Gradient::Steinbach_1996 && parameters::interface_potential == Int_Potential::Nestler_Obstacle)
				phase_field_functions::dfint_dphi_pairwise_acc = phase_field_functions::dfint_dphi_pairwise_S1996_NO;
			else if (parameters::interface_gradient == Int_Gradient::Steinbach_1996 && parameters::interface_potential == Int_Potential::Nestler_Well)
				phase_field_functions::dfint_dphi_pairwise_acc = phase_field_functions::dfint_dphi_pairwise_S1996_NW;
			else if (parameters::interface_gradient == Int_Gradient::Steinbach_1999 && parameters::interface_potential == Int_Potential::Nestler_Obstacle)
				phase_field_functions::dfint_dphi_pairwise_acc = phase_field_functions::dfint_dphi_pairwise_S1999_NO;
			else if (parameters::interface_gradient == Int_Gradient::Steinbach_1999 && parameters::interface_potential == Int_Potential::Nestler_Well)
				phase_field_functions::dfint_dphi_pairwise_acc = phase_field_functions::dfint_dphi_pairwise_S1999_NW;
			else
				phase_field_functions::dfint_dphi_pairwise_acc = phase_field_functions::dfint_dphi_pairwise_S2009;
			WriteDebugFile("# Model.DDCCPAI.Phi.InterfaceEnergy.const  = xi_ab \n");
			WriteDebugFile("#                                  .matrix = [(phi_a_name, phi_b_name, xi_ab), ...] \n");
			REAL const_xi_ab = 0; size_t phi_property_number = PhiProperties::instance().phi_property_number();
			std::string matrix_key1 = "Model.DDCCPAI.Phi.InterfaceEnergy.matrix", matrix_input1 = "[()]";
			parameters::xi_ab.resize(phi_property_number, phi_property_number, 0);
			phase_field_functions::_xi_ab = phase_field_functions::xi_ab;
			if (infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.const", const_xi_ab, true)) {
				for (size_t alpha_property = 0; alpha_property < phi_property_number; alpha_property++)
					for (size_t beta_property = 0; beta_property < phi_property_number; beta_property++)
						parameters::xi_ab(alpha_property, beta_property) = const_xi_ab;
			}
			if (infile_reader::read_string_value(matrix_key1, matrix_input1, true)) {
				std::vector<std::vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value
				({ InputValueType::IVType_STRING , InputValueType::IVType_STRING , InputValueType ::IVType_REAL}, matrix_key1, matrix_input1, true);
				for (size_t index = 0; index < matrix_value.size(); index++) {
					size_t alpha_property = PhiProperties::instance().phi_property(matrix_value[index][0].string_value), 
						beta_property = PhiProperties::instance().phi_property(matrix_value[index][1].string_value);
					parameters::xi_ab(alpha_property, beta_property) = matrix_value[index][2].REAL_value;
					parameters::xi_ab(beta_property, alpha_property) = matrix_value[index][2].REAL_value;
				}
			}

			WriteDebugFile("# Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.property = AnisotropicModelIndex \n");
			WriteDebugFile("#         AnisotropicModelIndex : 0 - ISO, 1 - CUBIC, 2 - HEX_BOETTGER, 3 - HEX_SUN, 4 - HEX_YANG, 5 - DENDRITE_YANG \n"); 
			int intEnAniso_model = 0;
			if (infile_reader::read_int_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.property", intEnAniso_model, true)) {
				parameters::intEnAniso_model = parameters::Int_Energy_Anisotropic(intEnAniso_model);
				if (parameters::intEnAniso_model == parameters::Int_Energy_Anisotropic::IEA_CUBIC) {
					phase_field_functions::_xi_ab = phase_field_functions::xi_ab_cubic;
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_1", parameters::intEnAniso_param1, true);
				}
				else if (parameters::intEnAniso_model == parameters::Int_Energy_Anisotropic::IEA_HEX_BOETTGER) {
					phase_field_functions::_xi_ab = phase_field_functions::xi_ab_hex_boettger;
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_1", parameters::intEnAniso_param1, true);
				}
				else if (parameters::intEnAniso_model == parameters::Int_Energy_Anisotropic::IEA_HEX_SUN) {
					phase_field_functions::_xi_ab = phase_field_functions::xi_ab_hex_sun;
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_1", parameters::intEnAniso_param1, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_2", parameters::intEnAniso_param2, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_3", parameters::intEnAniso_param3, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_4", parameters::intEnAniso_param4, true);
				}
				else if (parameters::intEnAniso_model == parameters::Int_Energy_Anisotropic::IEA_HEX_YANG) {
					phase_field_functions::_xi_ab = phase_field_functions::xi_ab_hex_yang;
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_1", parameters::intEnAniso_param1, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_2", parameters::intEnAniso_param2, true);
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_3", parameters::intEnAniso_param3, true);
				}
				else if (parameters::intEnAniso_model == parameters::Int_Energy_Anisotropic::IEA_DENDRITE_YANG) {
					phase_field_functions::_xi_ab = phase_field_functions::xi_ab_dendrite_yang;
					infile_reader::read_real_value("Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.parameter_1", parameters::intEnAniso_param1, true);
				}
			}
			WriteDebugFile("# Model.DDCCPAI.Phi.TripleJunctionEnergy.const  = xi_abc \n");
			WriteDebugFile("#                                       .matrix = [(phi_a_name, phi_b_name, phi_c_name, xi_abc), ...] \n");
			std::string matrix_key2 = "Model.DDCCPAI.Phi.TripleJunctionEnergy.matrix", matrix_input2 = "[()]";
			REAL const_xi_abc = 0;
			parameters::xi_abc.resize(phi_property_number, phi_property_number, phi_property_number, 0);
			phase_field_functions::_xi_abc = phase_field_functions::xi_abc;
			if (infile_reader::read_real_value("Model.DDCCPAI.Phi.TripleJunctionEnergy.const", const_xi_abc, true)) {
				for (size_t alpha_property = 0; alpha_property < phi_property_number; alpha_property++)
					for (size_t beta_property = 0; beta_property < phi_property_number; beta_property++)
						for (size_t gamma_property = 0; gamma_property < phi_property_number; gamma_property++)
							parameters::xi_abc(alpha_property, beta_property, gamma_property) = const_xi_abc;
			}
			if (infile_reader::read_string_value(matrix_key2, matrix_input2, true)) {
				std::vector<std::vector<input_value>> matrix_value = InputFileReader::get_instance()->
					trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_STRING , InputValueType::IVType_STRING , 
						InputValueType::IVType_STRING , InputValueType::IVType_REAL }, matrix_key2, matrix_input2, true);
				for (size_t index = 0; index < matrix_value.size(); index++)
					parameters::xi_abc(PhiProperties::instance().phi_property(matrix_value[index][0].string_value),
						PhiProperties::instance().phi_property(matrix_value[index][1].string_value),
						PhiProperties::instance().phi_property(matrix_value[index][2].string_value)) = matrix_value[index][3].REAL_value;
			}
			// - bulk energy const f_0
			WriteDebugFile("# Model.DDCCPAI.Phi.BulkEnergy.const = [(phi_name, bulk_energy), ... ], default energy = 0 \n");
			std::string bulk_energy_const_key = "Model.DDCCPAI.Phi.BulkEnergy.const", bulk_energy_const_input = "[()]";
			parameters::f_bulk_0.resize(phi_property_number, 0);
			if (infile_reader::read_string_value(bulk_energy_const_key, bulk_energy_const_input, true)) {
				std::vector<std::vector<input_value>> bulk_energy_const_value = InputFileReader::get_instance()->
					trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_STRING , InputValueType ::IVType_REAL}, 
						bulk_energy_const_key, bulk_energy_const_input, true);
				for (size_t index = 0; index < bulk_energy_const_value.size(); index++)
					parameters::f_bulk_0[PhiProperties::instance().phi_property(bulk_energy_const_value[index][0].string_value)] = 
						bulk_energy_const_value[index][1].REAL_value;
				parameters::delt_Fbulk_delt_phi.push_back(phase_field_functions::dfbulk_dphi_0);
			}
			// - phi noise
			parameters::is_phi_noise.resize(main_field::phi_number, false);
			WriteDebugFile("# Model.DDCCPAI.Phi.noise = ( phi_index, ... ) , noise will generate between two phi_index\n");
			std::string pn_key = "Model.DDCCPAI.Phi.noise", pn_input = "()";
			infile_reader::read_string_value(pn_key, pn_input, true);
			bool is_noise = false;
			std::vector<input_value> pn_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, pn_key, pn_input, true);
			for (size_t index = 0; index < pn_value.size(); index++)
				if (pn_value[index].int_value < main_field::phi_number) {
					parameters::is_phi_noise[pn_value[index].int_value] = true;
					is_noise = true;
				}
			if (is_noise) {
				infile_reader::read_int_value("Model.DDCCPAI.Phi.noise.begin_step", parameters::phi_noise_begin, true);
				infile_reader::read_int_value("Model.DDCCPAI.Phi.noise.end_step", parameters::phi_noise_end, true);
				infile_reader::read_int_value("Model.DDCCPAI.Phi.noise.frequency", parameters::phi_noise_frequency, true);
				infile_reader::read_real_value("Model.DDCCPAI.Phi.noise.amplitude", parameters::phi_noise_amplitude, true);
				if (infile_reader::read_int_value("Model.DDCCPAI.Phi.noise.seed", parameters::phi_noise_seed, true))
					parameters::is_phi_noise_rand = false;
				parameters::source_alpha_beta.push_back(phase_field_functions::noise_pairwise_acc);
			}
			// ==========================================================================================================================
			// - init con equation
			if (main_field::is_con_field_on) {
				parameters::init_con_in_moving_region = concentration_field_functions::init_con_in_moving_region_default;
				parameters::deinit_con_in_moving_region = concentration_field_functions::deinit_con_in_moving_region_default;
				if (parameters::diff_method == DifferenceMethod::SEVEN_POINT) {
					parameters::cal_mob_miu_grad_lap
						= concentration_field_functions::cal_mob_miu_grad_lap_7P;
				}
				else if (parameters::diff_method == DifferenceMethod::NINETEEN_POINT) {
					parameters::cal_mob_miu_grad_lap
						= concentration_field_functions::cal_mob_miu_grad_lap_19P;
				}
				// - 
				bool mob_with_temp = false;
				{
					REAL Mii_const_value = 0;
					WriteDebugFile("# Model.DDCCPAI.Con.BulkMobility.const  = Mii_value \n");
					WriteDebugFile("#                               .PhaseName.matrix = ( M11_value, ... ) \n");
					infile_reader::read_real_value("Model.DDCCPAI.Con.BulkMobility.const", Mii_const_value, true);
					parameters::Mii.resize(PhiProperties::instance().phi_property_number(), main_field::con_number, Mii_const_value);
					for (size_t phi_property = 0; phi_property < PhiProperties::instance().phi_property_number(); phi_property++) {
						std::string Mii_matrix_key = "Model.DDCCPAI.Con.BulkMobility." + PhiProperties::instance().phi_property_name(phi_property)
							+ ".matrix", Mii_matrix_input = "(";
						size_t region_index = ConRegions::instance().phi_property_region(phi_property);
						for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
							if (index < ConRegions::instance().region_con_number(region_index) - 1)
								Mii_matrix_input += "Mbulk_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ",";
							else
								Mii_matrix_input += "Mbulk_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ")";
						if (infile_reader::read_string_value(Mii_matrix_key, Mii_matrix_input, true)) {
							std::vector<input_value> Mii_matrix_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value
							(InputValueType::IVType_REAL, Mii_matrix_key, Mii_matrix_input, true);
							if (Mii_matrix_value.size() != ConRegions::instance().region_con_number(region_index)) {
								WriteDebugFile("# ERROR: Size of Mii should be "
									+ std::to_string(ConRegions::instance().region_con_number(region_index)) + " \n");
								SYS_PROGRAM_STOP;
							}
							for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
								parameters::Mii(phi_property, ConRegions::instance().region_con(region_index, index)) = Mii_matrix_value[index].REAL_value;
						}
					}
				}
				if (main_field::is_temp_field_on) {
					REAL Qii_const_value = 0;
					WriteDebugFile("# BulkMobility(T) = BulkMobility_0 * exp( - BulkMobilityQ / R / T) \n");
					WriteDebugFile("# Model.DDCCPAI.Con.BulkMobilityQ.const  = Qii_value \n");
					WriteDebugFile("#                                .PhaseName.matrix = ( Q11_value, ... ) \n");
					infile_reader::read_real_value("Model.DDCCPAI.Con.BulkMobilityQ.const", Qii_const_value, true);
					parameters::Qii.resize(PhiProperties::instance().phi_property_number(), main_field::con_number, Qii_const_value);
					for (size_t phi_property = 0; phi_property < PhiProperties::instance().phi_property_number(); phi_property++) {
						std::string Qii_matrix_key = "Model.DDCCPAI.Con.BulkMobilityQ." + PhiProperties::instance().phi_property_name(phi_property)
							+ ".matrix", Qii_matrix_input = "(";
						size_t region_index = ConRegions::instance().phi_property_region(phi_property);
						for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
							if (index < ConRegions::instance().region_con_number(region_index) - 1)
								Qii_matrix_input += "Qbulk_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ",";
							else
								Qii_matrix_input += "Qbulk_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ")";
						if (infile_reader::read_string_value(Qii_matrix_key, Qii_matrix_input, true)) {
							mob_with_temp = true;
							std::vector<input_value> Qii_matrix_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value
							(InputValueType::IVType_REAL, Qii_matrix_key, Qii_matrix_input, true);
							if (Qii_matrix_value.size() != ConRegions::instance().region_con_number(region_index)) {
								WriteDebugFile("# ERROR: Size of Qii should be "
									+ std::to_string(ConRegions::instance().region_con_number(region_index)) + " \n");
								SYS_PROGRAM_STOP;
							}
							for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
								parameters::Qii(phi_property, ConRegions::instance().region_con(region_index, index)) = Qii_matrix_value[index].REAL_value;
						}
					}
				}
				if (mob_with_temp)
					parameters::mobility.push_back(concentration_field_functions::bulk_diffusion_mobility_with_temperature);
				else
					parameters::mobility.push_back(concentration_field_functions::bulk_diffusion_mobility);
				{
					REAL Mii_surf_value = 0;
					bool mob_on_surf = false;
					WriteDebugFile("# Model.DDCCPAI.Con.SurfaceMobility.const = Mii_surf\n");
					WriteDebugFile("#                                  .RegionName.matrix = ( Msurf_Con1, ... ) \n");
					infile_reader::read_real_value("Model.DDCCPAI.Con.SurfaceMobility.const", Mii_surf_value, true);
					parameters::Mii_surf.resize(ConRegions::instance().region_number(), main_field::con_number, Mii_surf_value);
					for (size_t region_index = 0; region_index < ConRegions::instance().region_number(); region_index++) {
						std::string Mii_matrix_key = "Model.DDCCPAI.Con.SurfaceMobility." + ConRegions::instance().region_name(region_index)
							+ ".matrix", Mii_matrix_input = "(";
						for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
							if (index < ConRegions::instance().region_con_number(region_index) - 1)
								Mii_matrix_input += "Msurf_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ",";
							else
								Mii_matrix_input += "Msurf_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ")";
						if (infile_reader::read_string_value(Mii_matrix_key, Mii_matrix_input, true)) {
							mob_on_surf = true;
							std::vector<input_value> Mii_matrix_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value
							(InputValueType::IVType_REAL, Mii_matrix_key, Mii_matrix_input, true);
							if (Mii_matrix_value.size() != ConRegions::instance().region_con_number(region_index)) {
								WriteDebugFile("# ERROR: Size of Mii_surf should be "
									+ std::to_string(ConRegions::instance().region_con_number(region_index)) + " \n");
								SYS_PROGRAM_STOP;
							}
							for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
								parameters::Mii_surf(region_index, ConRegions::instance().region_con(region_index, index)) = Mii_matrix_value[index].REAL_value;
						}
					}
					if (mob_on_surf)
						parameters::mobility.push_back(concentration_field_functions::surface_diffusion_mobility);
				}
				{
					REAL Mii_int_value = 0;
					bool mob_on_int = false;
					WriteDebugFile("# Model.DDCCPAI.Con.IntMobility.const = Mii_int\n");
					WriteDebugFile("#                              .PhaseName1|PhaseName2.matrix = ( Mint_Con1, ... ) \n");
					infile_reader::read_real_value("Model.DDCCPAI.Con.IntMobility.const", Mii_int_value, true);
					parameters::Mii_grain.resize(PhiProperties::instance().phi_property_number(), 
						PhiProperties::instance().phi_property_number(), main_field::con_number, Mii_int_value);
					for (size_t region_index = 0; region_index < ConRegions::instance().region_number(); region_index++)
						for (size_t alpha_index = 0; alpha_index < ConRegions::instance().region_phi_property_number(region_index); alpha_index++)
							for (size_t beta_index = alpha_index + 1; beta_index < ConRegions::instance().region_phi_property_number(region_index); beta_index++) {
								size_t alpha_property = ConRegions::instance().region_phi_property(region_index, alpha_index),
									beta_property = ConRegions::instance().region_phi_property(region_index, beta_index);
								std::string Mii_matrix_key = "Model.DDCCPAI.Con.IntMobility." + PhiProperties::instance().phi_property_name(alpha_property)
									+ "|" + PhiProperties::instance().phi_property_name(beta_property) + ".matrix", Mii_matrix_input = "(";
								for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
									if (index < ConRegions::instance().region_con_number(region_index) - 1)
										Mii_matrix_input += "Mint_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ",";
									else
										Mii_matrix_input += "Mint_" + ConRegions::instance().con_name(ConRegions::instance().region_con(region_index, index)) + ")";
								if (infile_reader::read_string_value(Mii_matrix_key, Mii_matrix_input, true)) {
									mob_on_int = true;
									std::vector<input_value> Mii_matrix_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value
									(InputValueType::IVType_REAL, Mii_matrix_key, Mii_matrix_input, true);
									if (Mii_matrix_value.size() != ConRegions::instance().region_con_number(region_index)) {
										WriteDebugFile("# ERROR: Size of Mii_int should be "
											+ std::to_string(ConRegions::instance().region_con_number(region_index)) + " \n");
										SYS_PROGRAM_STOP;
									}
									for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++) {
										parameters::Mii_grain(alpha_property, beta_property, ConRegions::instance().region_con(region_index, index)) = Mii_matrix_value[index].REAL_value;
										parameters::Mii_grain(beta_property, alpha_property, ConRegions::instance().region_con(region_index, index)) = Mii_matrix_value[index].REAL_value;
									}
								}
							};
					if (mob_on_int)
						parameters::mobility.push_back(concentration_field_functions::interphase_diffusion_mobility);
				}
				{
					bool is_dfdcon_con = false;
					WriteDebugFile("# Model.DDCCPAI.Con.DiffusionPotential.con : \\nabla Mii \\cdot \\nabla con_i \n");
					if (infile_reader::read_bool_value("Model.DDCCPAI.Con.DiffusionPotential.con", is_dfdcon_con, true))
						if (is_dfdcon_con)
							parameters::delt_Fbulk_delt_con.push_back(concentration_field_functions::dfdcon_con);
				}
				{
					WriteDebugFile("# Model.DDCCPAI.Con.DiffusionPotential.Ki : \\nabla Mii \\cdot \\nabla [-K_i * laplace(con_i)] \n");
					WriteDebugFile("# Model.DDCCPAI.Con.DiffusionPotential.Ki = ( K1, ... )\n");
					parameters::Ki.resize(main_field::con_number, 0);
					std::string K_key = "Model.DDCCPAI.Con.DiffusionPotential.Ki", K_input = "(";
					for (size_t con_index = 0; con_index < main_field::con_number; con_index++)
						if (con_index == main_field::con_number - 1)
							K_input += "K_" + ConRegions::instance().con_name(con_index) + ")";
						else
							K_input += "K_" + ConRegions::instance().con_name(con_index) + ",";
					if (infile_reader::read_string_value(K_key, K_input, true)) {
						std::vector<input_value> K_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value
						(InputValueType::IVType_REAL, K_key, K_input, true);
						if (K_value.size() != main_field::con_number) {
							WriteDebugFile("# ERROR: Size of Ki should be " + std::to_string(main_field::con_number) + " \n");
							SYS_PROGRAM_STOP;
						}
						for (size_t index = 0; index < main_field::con_number; index++)
							parameters::Ki[index] = K_value[index].REAL_value;
						if (parameters::diff_method == DifferenceMethod::SEVEN_POINT)
							parameters::delt_Fbulk_delt_con.push_back(concentration_field_functions::dfdcon_lap_con_7P);
						else if(parameters::diff_method == DifferenceMethod::NINETEEN_POINT)
							parameters::delt_Fbulk_delt_con.push_back(concentration_field_functions::dfdcon_lap_con_19P);
					}
				}
			}
			// ==========================================================================================================================
			// - init temp euqation
			if (main_field::is_temp_field_on) {
				if (parameters::diff_method == DifferenceMethod::SEVEN_POINT) {
					parameters::cal_mob_temp_grad_lap
						= temperature_functions::cal_mob_temp_grad_lap_7P;
				}
				else if (parameters::diff_method == DifferenceMethod::NINETEEN_POINT) {
					parameters::cal_mob_temp_grad_lap
						= temperature_functions::cal_mob_temp_grad_lap_19P;
				}
				{
					REAL M_const_value = 0;
					WriteDebugFile("# Model.DDCCPAI.Temp.Mobility.const  = M_value \n");
					WriteDebugFile("#                            .PhaseName.matrix = M_value \n");
					infile_reader::read_real_value("Model.DDCCPAI.Temp.Mobility.const", M_const_value, true);
					parameters::Mtemp.resize(PhiProperties::instance().phi_property_number(), M_const_value);
					for (size_t phi_property = 0; phi_property < PhiProperties::instance().phi_property_number(); phi_property++) {
						std::string M_matrix_key = "Model.DDCCPAI.Temp.Mobility." + PhiProperties::instance().phi_property_name(phi_property)
							+ ".matrix";
						infile_reader::read_real_value(M_matrix_key, parameters::Mtemp[phi_property], true);
					}
					parameters::mobility_temp.push_back(temperature_functions::temperature_mobility_0);
				}
				{
					parameters::Ktemp.resize(PhiProperties::instance().phi_property_number(), 0);
					WriteDebugFile("# Model.DDCCPAI.Temp.Source.dphidtemp  = [(phase name, value), ... ] \n");
					std::string source_key = "Model.DDCCPAI.Temp.Source.dphidtemp", source_value = "[()]";
					if (infile_reader::read_string_value(source_key, source_value, true)) {
						std::vector<std::vector<input_value>> K_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value
						({ InputValueType::IVType_STRING, InputValueType::IVType_REAL }, source_key, source_value, true);
						for (size_t index = 0; index < K_value.size(); index++) {
							if (K_value[index].size() != 2) {
								WriteDebugFile("# ERROR: element in .Source.dphidtemp should be (phase name, value). \n");
								SYS_PROGRAM_STOP;
							}
							if (!PhiProperties::instance().is_phi_property(K_value[index][0].string_value)) {
								WriteDebugFile("# ERROR: phase name in .Source.dphidtemp has not been defined. \n");
								SYS_PROGRAM_STOP;
							}
							parameters::Ktemp[PhiProperties::instance().phi_property(K_value[index][0].string_value)] = K_value[index][1].REAL_value;
						}
						parameters::source_temp.push_back(temperature_functions::temperature_source_phi);
					}
				}
			}
			// ==========================================================================================================================
			// - init correlation terms in physics-informed equations
			{
				// - CALPHAD coupling
				parameters::is_energy_minimization.resize(ConRegions::instance().region_number(), false);
				WriteDebugFile("# Model.DDCCPAI.ThermoCalc.regions = ( region_name, ... )\n");
				std::string rg_key = "Model.DDCCPAI.ThermoCalc.regions", rg_input = "()";
				infile_reader::read_string_value(rg_key, rg_input, true);
				bool is_thermocalc = false;
				std::vector<input_value> rg_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, rg_key, rg_input, true);
				for(size_t index = 0; index < rg_value.size(); index++)
					if (ConRegions::instance().is_region(rg_value[index].string_value)) {
						size_t rg_index = ConRegions::instance().region_index(rg_value[index].string_value);
						parameters::is_energy_minimization[rg_index] = true;
						is_thermocalc = true;
					}
				if (is_thermocalc) {
					parameters::con_orders.resize(PhiProperties::instance().phi_property_number());
					parameters::params.resize(PhiProperties::instance().phi_property_number());
					parameters::temp_orders.resize(PhiProperties::instance().phi_property_number());
					parameters::terms_number.resize(PhiProperties::instance().phi_property_number(), 0);
					parameters::temp_terms_number.resize(PhiProperties::instance().phi_property_number());
					parameters::local_concentration_redistribution = ddc_calphad_ai_model::chemical_energy_functions::local_concentration_redistribution;
					for (size_t rg_index = 0; rg_index < ConRegions::instance().region_number(); rg_index++)
						if (parameters::is_energy_minimization[rg_index]) {
							for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(rg_index); index++) {
								WriteDebugFile("# Model.DDCCPAI.ThermoCalc.ChemEnergy.PhaseName = {[( Con_0_Order, Con_1_Order, ... ), ( Temp_0_Order, ... ), ( Temp_0_Param, ... )], ... }\n");
								size_t phi_property = ConRegions::instance().region_phi_property(rg_index, index);
								std::string chem_key = "Model.DDCCPAI.ThermoCalc.ChemEnergy." + PhiProperties::instance().phi_property_name(phi_property), 
									chem_input = "{[(),(),()]}";
								infile_reader::read_string_value(chem_key, chem_input, true);
								std::vector<std::vector<std::vector<input_value>>> chem_value = InputFileReader::get_instance()->
									trans_matrix_3d_const_array_const_to_input_value({ InputValueType::IVType_INT, InputValueType::IVType_INT, InputValueType::IVType_REAL }, chem_key, chem_input, true);
								parameters::terms_number[phi_property] = chem_value.size();
								parameters::con_orders[phi_property].resize(parameters::terms_number[phi_property], main_field::con_number, 0);
								parameters::temp_terms_number[phi_property].resize(parameters::terms_number[phi_property], 0);
								size_t max_temp_term = 0;
								for (size_t tindex = 0; tindex < parameters::terms_number[phi_property]; tindex++) {
									if (chem_value[tindex].size() != 3) {
										WriteDebugFile("# ERROR: Please follow format of " + chem_key + " to write the input value ! (terms " + std::to_string(tindex) + " )\n");
										SYS_PROGRAM_STOP;
									}
									if (chem_value[tindex][0].size() != main_field::con_number) {
										WriteDebugFile("# ERROR:  " + chem_key + " : Number of component should be " + std::to_string(main_field::con_number) + ", which defined before ! (terms " + std::to_string(tindex) + " ) \n");
										SYS_PROGRAM_STOP;
									}
									for (size_t cindex = 0; cindex < main_field::con_number; cindex++)
										parameters::con_orders[phi_property](tindex, cindex) = chem_value[tindex][0][cindex].int_value;
									if (chem_value[tindex][1].size() == 0) {
										WriteDebugFile("# ERROR:  " + chem_key + " : Temperature orders and parameters should be defined ! (terms " + std::to_string(tindex) + " ) \n");
										SYS_PROGRAM_STOP;
									}
									if (chem_value[tindex][1].size() != chem_value[tindex][2].size()) {
										WriteDebugFile("# ERROR:  " + chem_key + " : Number of Temperature orders and parameters should be equal ! (terms " + std::to_string(tindex) + " ) \n");
										SYS_PROGRAM_STOP;
									}
									parameters::temp_terms_number[phi_property][tindex] = chem_value[tindex][1].size();
									if (max_temp_term < parameters::temp_terms_number[phi_property][tindex])
										max_temp_term = parameters::temp_terms_number[phi_property][tindex];
								}
								parameters::params[phi_property].resize(parameters::terms_number[phi_property], max_temp_term, 0);
								parameters::temp_orders[phi_property].resize(parameters::terms_number[phi_property], max_temp_term, 0);
								for (size_t tindex = 0; tindex < parameters::terms_number[phi_property]; tindex++)
									for (size_t temp_index = 0; temp_index < parameters::temp_terms_number[phi_property][tindex]; temp_index++) {
										parameters::temp_orders[phi_property](tindex, temp_index) = chem_value[tindex][1][temp_index].int_value;
										parameters::params[phi_property](tindex, temp_index) = chem_value[tindex][2][temp_index].REAL_value;
									}
							}
						}
					infile_reader::read_int_value("Model.DDCCPAI.ThermoCalc.max_variation_step", parameters::max_phasecon_variation_step, true);
					infile_reader::read_real_value("Model.DDCCPAI.ThermoCalc.L", parameters::L_phasecon, true);
					infile_reader::read_real_value("Model.DDCCPAI.ThermoCalc.epsilon", parameters::phasecon_epsilon, true);
					parameters::fchem = chemical_energy_functions::fchem_polynomial;
					parameters::miu = chemical_energy_functions::miu_polynomial;
					parameters::delt_Fbulk_delt_phi.push_back(chemical_energy_functions::delt_Fchem_delt_phi);
					parameters::delt_Fbulk_delt_con.push_back(chemical_energy_functions::delt_Fchem_delt_con);
				}
			}
			// ==========================================================================================================================
			// - output a standalone thermodynamic energy-minimization scan
			{
				const std::string csv_prefix = "Model.DDCCPAI.Output.CSV.EnergyMinimization.";
				WriteDebugFile("# " + csv_prefix + "active = true/false\n");
				infile_reader::read_bool_value(csv_prefix + "active", parameters::is_write_thermo_calc_csv, true);
				if (parameters::is_write_thermo_calc_csv) {
					auto input_error = [](const std::string& message) {
						WriteLog("> ERROR: " + message + "\n");
						SYS_PROGRAM_STOP;
					};
					if (!main_field::is_con_field_on)
						input_error("DDCCPAI energy-minimization CSV output requires the concentration field.");

					WriteDebugFile("# " + csv_prefix + "file_name = thermo_calc_data.csv\n");
					infile_reader::read_string_value(csv_prefix + "file_name", parameters::thermo_calc_equi_csv_name, true);
					if (parameters::thermo_calc_equi_csv_name.empty())
						input_error(csv_prefix + "file_name cannot be empty.");

					std::string region_name;
					WriteDebugFile("# " + csv_prefix + "region = RegionName\n");
					if (!infile_reader::read_string_value(csv_prefix + "region", region_name, true) ||
						!ConRegions::instance().is_region(region_name))
						input_error(csv_prefix + "region is missing or has not been defined.");
					parameters::thermo_calc_region = ConRegions::instance().region_index(region_name);
					if (parameters::thermo_calc_region >= parameters::is_energy_minimization.size() ||
						!parameters::is_energy_minimization[parameters::thermo_calc_region])
						input_error(csv_prefix + "region must be enabled in Model.DDCCPAI.ThermoCalc.regions");

					parameters::thermo_calc_phases.clear();
					parameters::thermo_calc_phi_property.assign(PhiProperties::instance().phi_property_number(), false);
					for (size_t index = 0; index < ConRegions::instance().region_phi_property_number(parameters::thermo_calc_region); index++) {
						size_t phi_property = ConRegions::instance().region_phi_property(parameters::thermo_calc_region, index);
						parameters::thermo_calc_phi_property[phi_property] = true;
						parameters::thermo_calc_phases.push_back(phi_property);
					}
					if (parameters::thermo_calc_phases.empty())
						input_error(csv_prefix + "region does not contain a defined phase.");

					const std::string con_key = csv_prefix + "concentration";
					std::string con_input = "()";
					WriteDebugFile("# " + con_key + " = (start, end, step)   # common scan range for all components in the selected region\n");
					if (!infile_reader::read_string_value(con_key, con_input, true))
						input_error(con_key + " is required.");
					std::vector<input_value> con_values = InputFileReader::get_instance()->
						trans_matrix_1d_array_to_input_value(
							{ InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL },
							con_key, con_input, true);
					ThermoCalcScanRange con_range{ con_values[0].REAL_value, con_values[1].REAL_value, con_values[2].REAL_value };
					if (con_range.begin < 0 || con_range.end >= 1 || con_range.begin > con_range.end || con_range.step <= 0)
						input_error(con_key + " must satisfy 0 <= start <= end < 1 and step > 0.");
					parameters::thermo_calc_con_ranges.assign(main_field::con_number, ThermoCalcScanRange());
					parameters::thermo_calc_fix_con.assign(main_field::con_number, { false, 0 });
					for (size_t index = 0; index < ConRegions::instance().region_con_number(parameters::thermo_calc_region); index++) {
						size_t con_index = ConRegions::instance().region_con(parameters::thermo_calc_region, index);
						parameters::thermo_calc_con_ranges[con_index] = con_range;
						parameters::thermo_calc_fix_con[con_index] = { con_range.begin == con_range.end, con_range.begin };
					}

					const std::string phi_key = csv_prefix + "phase_fraction";
					std::string phi_input = "()";
					WriteDebugFile("# " + phi_key + " = (start, end, step)    # common scan range for the first N-1 phases; the last phase completes the sum to 1\n");
					if (!infile_reader::read_string_value(phi_key, phi_input, true))
						input_error(phi_key + " is required.");
					std::vector<input_value> phi_values = InputFileReader::get_instance()->
						trans_matrix_1d_array_to_input_value(
							{ InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL },
							phi_key, phi_input, true);
					ThermoCalcScanRange phi_range{ phi_values[0].REAL_value, phi_values[1].REAL_value, phi_values[2].REAL_value };
					if (phi_range.begin <= 0 || phi_range.end >= 1 || phi_range.begin > phi_range.end || phi_range.step <= 0)
						input_error(phi_key + " must satisfy 0 < start <= end < 1 and step > 0.");
					parameters::thermo_calc_phi_ranges.assign(PhiProperties::instance().phi_property_number(), ThermoCalcScanRange());
					for (size_t index = 0; index < parameters::thermo_calc_phases.size(); index++)
						parameters::thermo_calc_phi_ranges[parameters::thermo_calc_phases[index]] = phi_range;

					const std::string temp_key = csv_prefix + "temperature";
					std::string temp_input = "()";
					WriteDebugFile("# " + temp_key + " = (start, end, step)\n");
					if (!infile_reader::read_string_value(temp_key, temp_input, true))
						input_error(temp_key + " is required.");
					std::vector<input_value> temp_values = InputFileReader::get_instance()->
						trans_matrix_1d_array_to_input_value(
							{ InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL },
							temp_key, temp_input, true);
					parameters::thermo_calc_temp_range =
						{ temp_values[0].REAL_value, temp_values[1].REAL_value, temp_values[2].REAL_value };
					if (parameters::thermo_calc_temp_range.begin > parameters::thermo_calc_temp_range.end ||
						parameters::thermo_calc_temp_range.step <= 0)
						input_error(temp_key + " must satisfy start <= end and step > 0.");
				}
			}
			// ==========================================================================================================================
			// - output chemical-energy scans for selected phases
			{
				const std::string energy_key = "Model.DDCCPAI.Output.CSV.ChemEnergy";
				std::string energy_input = "()";
				WriteDebugFile("# " + energy_key + " = (PhaseName, ...)\n");
				parameters::thermo_calc_energy_phases.clear();
				if (infile_reader::read_string_value(energy_key, energy_input, true)) {
					auto input_error = [](const std::string& message) {
						WriteLog("> ERROR: " + message + "\n");
						SYS_PROGRAM_STOP;
					};
					std::vector<input_value> energy_values = InputFileReader::get_instance()->
						trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, energy_key, energy_input, true);
					std::vector<bool> phase_defined(PhiProperties::instance().phi_property_number(), false);
					std::vector<bool> required_con(main_field::con_number, false);
					for (const input_value& energy_value : energy_values) {
						const std::string& phase_name = energy_value.string_value;
						if (!PhiProperties::instance().is_phi_property(phase_name))
							input_error(energy_key + " contains an undefined phase: " + phase_name + ".");
						size_t phi_property = PhiProperties::instance().phi_property(phase_name);
						size_t region_index = ConRegions::instance().phi_property_region(phi_property);
						if (region_index >= parameters::is_energy_minimization.size() ||
							!parameters::is_energy_minimization[region_index])
							input_error("Phase " + phase_name + " is not in an energy-minimization region.");
						if (phase_defined[phi_property])
							input_error(energy_key + " contains duplicate phase: " + phase_name + ".");
						phase_defined[phi_property] = true;
						parameters::thermo_calc_energy_phases.push_back(phi_property);
						for (size_t index = 0; index < ConRegions::instance().region_con_number(region_index); index++)
							required_con[ConRegions::instance().region_con(region_index, index)] = true;
					}

					if (!parameters::thermo_calc_energy_phases.empty()) {
						WriteDebugFile("# " + energy_key + ".file_name = thermo_calc_energy_data.csv\n");
						infile_reader::read_string_value(energy_key + ".file_name", parameters::thermo_calc_energy_csv_name, true);
						if (parameters::thermo_calc_energy_csv_name.empty())
							input_error(energy_key + ".file_name cannot be empty.");

						const std::string con_key = energy_key + ".concentration";
						std::string con_input = "[()]";
						WriteDebugFile("# " + con_key + " = [(ConName, start, end, step), ...]   # all components used by selected phases\n");
						if (!infile_reader::read_string_value(con_key, con_input, true))
							input_error(con_key + " is required.");
						std::vector<std::vector<input_value>> con_values = InputFileReader::get_instance()->
							trans_matrix_2d_const_array_to_input_value(
								{ InputValueType::IVType_STRING, InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL },
								con_key, con_input, true);
						parameters::thermo_calc_energy_con_ranges.assign(main_field::con_number, ThermoCalcScanRange());
						std::vector<bool> con_defined(main_field::con_number, false);
						for (const std::vector<input_value>& con_value : con_values) {
							const std::string& con_name = con_value[0].string_value;
							if (!ConRegions::instance().is_con(con_name))
								input_error(con_key + " contains an undefined concentration: " + con_name + ".");
							size_t con_index = ConRegions::instance().con_index(con_name);
							if (!required_con[con_index])
								input_error("Concentration " + con_name + " is not used by a selected phase.");
							if (con_defined[con_index])
								input_error(con_key + " contains duplicate concentration: " + con_name + ".");
							ThermoCalcScanRange range{ con_value[1].REAL_value, con_value[2].REAL_value, con_value[3].REAL_value };
							if (range.begin < 0 || range.end > 1 || range.begin > range.end || range.step <= 0)
								input_error("Concentration range for " + con_name + " must satisfy 0 <= start <= end <= 1 and step > 0.");
							parameters::thermo_calc_energy_con_ranges[con_index] = range;
							con_defined[con_index] = true;
						}
						for (size_t con_index = 0; con_index < main_field::con_number; con_index++)
							if (required_con[con_index] && !con_defined[con_index])
								input_error(con_key + " must define concentration " + ConRegions::instance().con_name(con_index) + ".");

						const std::string temp_key = energy_key + ".temperature";
						std::string temp_input = "()";
						WriteDebugFile("# " + temp_key + " = (start, end, step)\n");
						if (!infile_reader::read_string_value(temp_key, temp_input, true))
							input_error(temp_key + " is required.");
						std::vector<input_value> temp_values = InputFileReader::get_instance()->
							trans_matrix_1d_array_to_input_value(
								{ InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL },
								temp_key, temp_input, true);
						parameters::thermo_calc_energy_temp_range =
							{ temp_values[0].REAL_value, temp_values[1].REAL_value, temp_values[2].REAL_value };
						if (parameters::thermo_calc_energy_temp_range.begin > parameters::thermo_calc_energy_temp_range.end ||
							parameters::thermo_calc_energy_temp_range.step <= 0)
							input_error(temp_key + " must satisfy start <= end and step > 0.");
					}
				}
			}
			// ==========================================================================================================================
			bool buff = false;
			InputFileReader::get_instance()->read_bool_value("Model.DDCCPAI.Output.VTS.active_phis", buff, true);
			if (buff)
				write_vts::load_vts_func(phase_field_functions::write_scalar_active_phi_number); buff = false;
			InputFileReader::get_instance()->read_bool_value("Model.DDCCPAI.Output.VTS.con_all", buff, true);
			if (buff)
				write_vts::load_vts_func(concentration_field_functions::write_scalar_con_all);
			{
				WriteDebugFile("# Model.DDCCPAI.Output.VTS.phase_con = [( phase_name, con_name ), ... ]\n");
				std::string vts_key = "Model.DDCCPAI.Output.VTS.phase_con", vts_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(vts_key, vts_input, true)) {
					std::vector<std::vector<input_value>> vts_value = InputFileReader::get_instance()->
						trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_STRING, vts_key, vts_input, true);
					write_vts::load_vts_func(chemical_energy_functions::write_scalar_phasecon);
					for (size_t index = 0; index < vts_value.size(); index++) {
						if (!PhiProperties::instance().is_phi_property(vts_value[index][0].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_con : The phase name " + vts_value[index][0].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						}
						size_t region_index = ConRegions::instance().phi_property_region(PhiProperties::instance().phi_property(vts_value[index][0].string_value));
						if (!parameters::is_energy_minimization[region_index]) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_con : The phase name " + vts_value[index][0].string_value + " do not need calphad calculation. \n");
							SYS_PROGRAM_STOP;
						}
						if (!ConRegions::instance().is_con(vts_value[index][1].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_con : The con name " + vts_value[index][1].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						}
						if (!ConRegions::instance().is_con_in_region(region_index, ConRegions::instance().con_index(vts_value[index][1].string_value))) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_con : The con name " + vts_value[index][1].string_value + " has not been defined in this phase & region. \n");
							SYS_PROGRAM_STOP;
						}
						parameters::is_write_phase_con.push_back(
							{ PhiProperties::instance().phi_property(vts_value[index][0].string_value), 
							  ConRegions::instance().con_index(vts_value[index][1].string_value)});
					}
				}
			}
			{
				WriteDebugFile("# Model.DDCCPAI.Output.VTS.con = ( con_name, ... ) \n");
				std::string vts_key = "Model.DDCCPAI.Output.VTS.con", vts_input = "()";
				if (InputFileReader::get_instance()->read_string_value(vts_key, vts_input, true)) {
					std::vector<input_value> vts_value = InputFileReader::get_instance()->
						trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, vts_key, vts_input, true);
					write_vts::load_vts_func(chemical_energy_functions::write_scalar_con);
					for (size_t index = 0; index < vts_value.size(); index++) {
						if (!ConRegions::instance().is_con(vts_value[index].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.con : The con name " + vts_value[index].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						}
						parameters::is_write_con.push_back(ConRegions::instance().con_index(vts_value[index].string_value));
					}
				}
			}
			{
				WriteDebugFile("# Model.DDCCPAI.Output.VTS.phase_miu = [( phase_name, con_name ), ... ]\n");
				std::string vts_key = "Model.DDCCPAI.Output.VTS.phase_miu", vts_input = "[()]";
				if (InputFileReader::get_instance()->read_string_value(vts_key, vts_input, true)) {
					std::vector<std::vector<input_value>> vts_value = InputFileReader::get_instance()->
						trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_STRING, vts_key, vts_input, true);
					write_vts::load_vts_func(chemical_energy_functions::write_scalar_phasemiu);
					for (size_t index = 0; index < vts_value.size(); index++) {
						if (!PhiProperties::instance().is_phi_property(vts_value[index][0].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_miu : The phase name " + vts_value[index][0].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						} 
						size_t region_index = ConRegions::instance().phi_property_region(PhiProperties::instance().phi_property(vts_value[index][0].string_value));
						if (!parameters::is_energy_minimization[region_index]) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_miu : The phase name " + vts_value[index][0].string_value + " do not need calphad calculation. \n");
							SYS_PROGRAM_STOP;
						}
						if (!ConRegions::instance().is_con(vts_value[index][1].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_miu : The con name " + vts_value[index][1].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						}
						if (!ConRegions::instance().is_con_in_region(region_index, ConRegions::instance().con_index(vts_value[index][1].string_value))) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.phase_miu : The con name " + vts_value[index][1].string_value + " has not been defined in this phase & region. \n");
							SYS_PROGRAM_STOP;
						}
						parameters::is_write_phase_miu.push_back(
							{ PhiProperties::instance().phi_property(vts_value[index][0].string_value),
							  ConRegions::instance().con_index(vts_value[index][1].string_value) });
					}
				}
			}
			{
				WriteDebugFile("# Model.DDCCPAI.Output.VTS.miu = ( con_name, ... ) \n");
				std::string vts_key = "Model.DDCCPAI.Output.VTS.miu", vts_input = "()";
				if (InputFileReader::get_instance()->read_string_value(vts_key, vts_input, true)) {
					std::vector<input_value> vts_value = InputFileReader::get_instance()->
						trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, vts_key, vts_input, true);
					write_vts::load_vts_func(chemical_energy_functions::write_scalar_miu);
					for (size_t index = 0; index < vts_value.size(); index++) {
						if (!ConRegions::instance().is_con(vts_value[index].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.miu : The con name " + vts_value[index].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						}
						parameters::is_write_miu.push_back(ConRegions::instance().con_index(vts_value[index].string_value));
					}
				}
			}
			{
				WriteDebugFile("# Model.DDCCPAI.Output.VTS.f_chemical = ( phase_name, ... ) \n");
				std::string vts_key = "Model.DDCCPAI.Output.VTS.f_chemical", vts_input = "()";
				if (InputFileReader::get_instance()->read_string_value(vts_key, vts_input, true)) {
					std::vector<input_value> vts_value = InputFileReader::get_instance()->
						trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, vts_key, vts_input, true);
					write_vts::load_vts_func(chemical_energy_functions::write_scalar_fchem);
					for (size_t index = 0; index < vts_value.size(); index++) {
						if (!PhiProperties::instance().is_phi_property(vts_value[index].string_value)) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.f_chemical : The phase name " + vts_value[index].string_value + " has not been defined. \n");
							SYS_PROGRAM_STOP;
						}
						size_t region_index = ConRegions::instance().phi_property_region(PhiProperties::instance().phi_property(vts_value[index].string_value));
						if (!parameters::is_energy_minimization[region_index]) {
							WriteDebugFile("# ERROR: Model.DDCCPAI.Output.VTS.f_chemical : The phase name " + vts_value[index].string_value + " do not need calphad calculation. \n");
							SYS_PROGRAM_STOP;
						}
						parameters::is_write_fchem.push_back(PhiProperties::instance().phi_property(vts_value[index].string_value));
					}
				}
			}
			// ==========================================================================================================================
			// infile_reader::read_real_value("SimulationModels.DataDrivenComplex.L", parameters::L, true);
			load_a_new_module(nullptr, exec_pre_ii, exec_pre_iii,  // exec_pre_i   exec_pre_ii    exec_pre_iii
				exec_i, nullptr, nullptr,  // exec_i   exec_ii   exec_iii
				exec_pos_i, nullptr, exec_pos_iii,   // exec_pos_i   exec_pos_ii   exec_pos_iii
				deinit);  // deinit
		}
	}
}
