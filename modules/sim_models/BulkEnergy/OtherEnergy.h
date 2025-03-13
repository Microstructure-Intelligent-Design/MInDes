#pragma once
#include "../../Base.h"
namespace pf {
	namespace other_energy {
		enum dfdpctType { // phi , phase con , grand potential , temperature
			df_Const,
			df_dphi_supercooling
		};
		struct supercool_para_t {
			double alpha{}, gamma{}, T_eq{};
			double energy_A{}, energy_B{};
		}supercool_paras;

		static int model_type = dfdpctType::df_Const;
		static double dfother_dphi_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return 0.0;
		}
		static void dfother_dcon_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return;
		}
		static double dfother_dconi_const(pf::PhaseNode& node, int con_i) {
			return 0.0;
		}
		static double dfother_dphi_supercooling(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double m = supercool_paras.alpha / (4 * atan(1)) * std::atan(supercool_paras.gamma * (supercool_paras.T_eq - node.temperature.T));
			return supercool_paras.energy_A * m * phase.phi * phase.phi - supercool_paras.energy_B * m * phase.phi;
		};

		static double (*fother_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*dfother_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of this model

		static void (*dfother_dcon)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*dfother_dconi)(pf::PhaseNode& node, int con_i);

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			fother_density = dfother_dphi_const;
			dfother_dphi = dfother_dphi_const;
			dfother_dcon = dfother_dcon_const;
			dfother_dconi = dfother_dconi_const;
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard) {
				if (infile_debug)
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.PhiConT.OtherEnergy.type = 1 - Supercool Driving Force \n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.PhiConT.OtherEnergy.type", model_type, infile_debug);
				string in_key = "", in_string = ""; vector<input_value> in_value;
				switch (model_type)
				{
				case pf::other_energy::df_dphi_supercooling:
					if (infile_debug) {
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("\n# Supercool Driving Force : df_dphi = A m phi^2 + B m phi \n# m(T) = (alpha/pi) arctan (gamma(Teq-T))  \n", InputFileReader::get_instance()->debug_file);
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.OtherEnergy.Supercooling.driving_force_para = (alpha, gamma, T_eq) \n", InputFileReader::get_instance()->debug_file);
						InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.Phi.OtherEnergy.Supercooling.Energy_Coefficient = (A, B) \n", InputFileReader::get_instance()->debug_file);
					}
					in_string = "()";
					in_key = "ModelsManager.Phi.OtherEnergy.Supercooling.driving_force_para";
					InputFileReader::get_instance()->read_string_value(in_key, in_string, infile_debug);
					in_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, in_key, in_string, infile_debug);
					supercool_paras.alpha = in_value[0].double_value;
					supercool_paras.gamma = in_value[1].double_value;
					supercool_paras.T_eq = in_value[2].double_value;
					in_key = "ModelsManager.Phi.OtherEnergy.Supercooling.Energy_Coefficient";
					InputFileReader::get_instance()->read_string_value(in_key, in_string, infile_debug);
					in_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, in_key, in_string, infile_debug);
					supercool_paras.energy_A = in_value[0].double_value;
					supercool_paras.energy_B = in_value[1].double_value;
					dfother_dphi = dfother_dphi_supercooling;

					break;
				default:
					break;
				}

			}

		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}