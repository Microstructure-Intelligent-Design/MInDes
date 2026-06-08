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
				string error_report = "> ERROR, the PCT command line settings do not meet the requirements : PCT = (N>0,K,false/true) !\n";
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
			WriteDebugFile("# Model.DDCCPAI.Phi.IntMobility.Anisotropic.property = (phi_name, AnisotropicModelIndex)\n");
			WriteDebugFile("#          AnisotropicModelIndex : 0 - IEA_ISO, 1 - IEA_CUBIC, 2 - IEA_HEX_BOETTGER, 3 - IEA_HEX_SUN, 4 - IEA_HEX_YANG \n");
			std::string matrix_aniso_key = "Model.DDCCPAI.Phi.IntMobility.Anisotropic.property", matrix_aniso_input = "()";
			if (infile_reader::read_string_value(matrix_aniso_key, matrix_aniso_input, true)) {
				std::vector<input_value> matrix_aniso_value =
					InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(
						{ InputValueType::IVType_STRING , InputValueType::IVType_INT }, 
						matrix_aniso_key, matrix_aniso_input, true);
				parameters::intMobAniso_phi_property = PhiProperties::instance().phi_property(matrix_aniso_value[0].string_value);
				parameters::intMobAniso_model = parameters::Int_Mobility_Anisotropic(matrix_aniso_value[1].int_value);
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

			WriteDebugFile("# Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.property = (phi_name, AnisotropicModelIndex)\n");
			WriteDebugFile("#         AnisotropicModelIndex : 0 - IEA_ISO, 1 - IEA_CUBIC, 2 - IEA_HEX_BOETTGER, 3 - IEA_HEX_SUN, 4 - IEA_HEX_YANG \n");
			std::string matrix_aniso_key2 = "Model.DDCCPAI.Phi.InterfaceEnergy.Anisotropic.property", matrix_aniso_input2 = "()";
			if (infile_reader::read_string_value(matrix_aniso_key2, matrix_aniso_input2, true)) {
				std::vector<input_value> matrix_aniso_value2 = InputFileReader::get_instance()->
					trans_matrix_1d_array_to_input_value({ InputValueType::IVType_STRING , InputValueType::IVType_INT }, matrix_aniso_key2, matrix_aniso_input2, true);
				parameters::intEnAniso_phi_property = PhiProperties::instance().phi_property(matrix_aniso_value2[0].string_value);
				parameters::intEnAniso_model = parameters::Int_Energy_Anisotropic(matrix_aniso_value2[1].int_value);
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
			}
			// ==========================================================================================================================
			// - init correlation terms in physics-informed equations
			
			// ==========================================================================================================================
			bool buff = false;
			InputFileReader::get_instance()->read_bool_value("Model.DDCCPAI.Output.VTS.active_phis", buff, true);
			if (buff)
				write_vts::load_vts_func(phase_field_functions::write_scalar_active_phi_number);
			// ==========================================================================================================================
			// infile_reader::read_real_value("SimulationModels.DataDrivenComplex.L", parameters::L, true);
			load_a_new_module(nullptr, exec_pre_ii, exec_pre_iii,  // exec_pre_i   exec_pre_ii    exec_pre_iii
				exec_i, nullptr, nullptr,  // exec_i   exec_ii   exec_iii
				exec_pos_i, nullptr, nullptr,   // exec_pos_i   exec_pos_ii   exec_pos_iii
				deinit);  // deinit
		}
	}
}
