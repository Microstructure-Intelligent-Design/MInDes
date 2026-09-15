#pragma once
#include "../../Module.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "DS_Params.h"
#include "DS_Functions.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace pf {
	namespace dendrite_solidification_model {
		inline void init_model_modules() {
			auto input_error = [](const std::string& message) {
				WriteLog("> ERROR, Dendrite Solidification: " + message + "\n");
				WriteDebugFile("> ERROR, Dendrite Solidification: " + message + "\n");
				SYS_PROGRAM_STOP;
			};
			if (main_field::phi_number < 2 || main_field::con_number % 2 != 0 || !main_field::is_temp_field_on)
				input_error("PCT must be (N>=2,0,true) or (N>=2,2K,true), with K>=0.");
			parameters::component_number = main_field::con_number / 2;
			parameters::callbacks_frozen = false;
			parameters::component_names.clear();
			parameters::liquid_con_equilibrium.clear();
			parameters::driving_force_dcon.clear();
			parameters::con_mobility_field.clear();
			parameters::con_mobility_const = 0;

			WriteLog("> \n> Simulation Model - Dendrite Solidification - is Activated !\n> \n");
			WriteDebugFile(parameters::is_thermal_only()
				? "# DS mode: thermal-only synchronous phi-temperature evolution, PCT=(N,0,true).\n"
				: "# DS mode: multicomponent synchronous phi-concentration-temperature evolution, PCT=(N,2K,true).\n");
			WriteDebugFile("# DS uses a seven-point stencil.\n");
			WriteDebugFile("# DS phase convention: phi[0] is liquid; phi[index>0] is a solid grain.\n");
			WriteDebugFile("# dphi_a/dt = sum_b L_ab/eta * (M_b-M_a) + S_ab\n");
			WriteDebugFile("# dT/dt = div(M_T grad(T)) + sum_a K_a*dphi_a/dt + S_T\n");
			if (parameters::has_components()) {
				WriteDebugFile("# DS concentration convention: the first K fields are liquid components and the last K fields are the matching solid components.\n");
				WriteDebugFile("# dc_j/dt = div(M_j grad(mu_j)) + S_j, j=0,...,2K-1\n");
				WriteDebugFile("# DeltaG = sum_i dGdc_liquid_i*(c_liquid_i-c_liquid_i_eq) + dGdT*(T-Teq); positive DeltaG promotes solidification.\n");
			}
			else {
				WriteDebugFile("# DS has no concentration fields or concentration equation in thermal-only mode.\n");
				WriteDebugFile("# Model.DS.Con.* and Model.DS.Phi.DrivingForce.components are not applicable and are not read.\n");
				WriteDebugFile("# DeltaG = dGdT*(T-Teq); positive DeltaG promotes solidification.\n");
			}
			WriteDebugFile("# DS main fields are the immutable step-start state until all right-hand sides are complete.\n");
			WriteDebugFile("# DS workspace and per-thread normalization scratch are allocated once in exec_pre_iii and remain fixed during exec_i.\n");

			if (parameters::has_components()) {
				const std::string component_names_key = "Model.DS.Con.ComponentNames";
				std::string component_names_input = "()";
				if (!infile_reader::read_string_value(component_names_key, component_names_input, true))
					input_error(component_names_key + " must define every independent component.");
				auto component_name_values = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(
					InputValueType::IVType_STRING, component_names_key, component_names_input, true);
				if (component_name_values.size() != parameters::component_number)
					input_error("Con.ComponentNames must contain exactly K names.");
				std::unordered_set<std::string> unique_component_names;
				for (const auto& value : component_name_values) {
					const std::string& name = value.string_value;
					if (name.empty() || !unique_component_names.insert(name).second)
						input_error("Con.ComponentNames contains an empty or duplicate name.");
					parameters::component_names.push_back(name);
				}
				for (size_t i = 0; i < parameters::component_number; ++i) {
					WriteDebugFile("# DS concentration_field[" + std::to_string(parameters::liquid_con_index(i)) + "] = liquid/" + parameters::component_names[i] + "\n");
					WriteDebugFile("# DS concentration_field[" + std::to_string(parameters::solid_con_index(i)) + "] = solid/" + parameters::component_names[i] + "\n");
				}
			}

			GrainsOrientations::instance().init();
			Matrix3x3 identity; identity.set_to_unity();
			parameters::grain_rotation_matrix.resize(main_field::phi_number, identity);
			for (size_t i = 0; i < main_field::phi_number; ++i)
				parameters::grain_rotation_matrix[i] = GrainsOrientations::instance().RotationMatrix(i);

			parameters::PAIRWISE_ACC_STOP = main_field::phi_number;
			parameters::PHI_ACC_NUMBER = main_field::phi_number;
			infile_reader::read_int_value("Model.DS.Phi.PairWiseAcc.container_size", parameters::PHI_ACC_NUMBER, true);
			if (parameters::PHI_ACC_NUMBER == 0) input_error("Phi.PairWiseAcc.container_size must be positive.");
			parameters::PHI_ACC_NUMBER = std::min(parameters::PHI_ACC_NUMBER, main_field::phi_number);
			infile_reader::read_real_value("Model.DS.Phi.cutoff", parameters::Phi_Cut_Off, true);
			infile_reader::read_real_value("Model.DS.Phi.con_cutoff", parameters::PhiCon_Cut_Off, true);
			infile_reader::read_bool_value("Model.DS.Phi.is_normalize", parameters::is_phi_normalized, true);
			if (parameters::Phi_Cut_Off < 0 || parameters::Phi_Cut_Off >= REAL(0.5)) input_error("Phi.cutoff must satisfy 0 <= cutoff < 0.5.");
			if (parameters::PhiCon_Cut_Off < 0 || parameters::PhiCon_Cut_Off > 1) input_error("Phi.con_cutoff must lie in [0,1].");
			parameters::Phi_Cut_Off_R = REAL(1) - parameters::Phi_Cut_Off;

			REAL mobility_const = 0;
			infile_reader::read_real_value("Model.DS.Phi.IntMobility.const", mobility_const, true);
			parameters::Lij.resize(main_field::phi_number, main_field::phi_number, mobility_const);
			std::string mobility_input = "[()]";
			if (infile_reader::read_string_value("Model.DS.Phi.IntMobility.matrix", mobility_input, true)) {
				auto values = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(
					{ InputValueType::IVType_INT, InputValueType::IVType_INT, InputValueType::IVType_REAL },
					"Model.DS.Phi.IntMobility.matrix", mobility_input, true);
				for (const auto& row : values) {
					const size_t a = row[0].int_value, b = row[1].int_value;
					if (a >= main_field::phi_number || b >= main_field::phi_number || a == b) input_error("IntMobility.matrix contains an invalid phase pair.");
					parameters::Lij(a, b) = parameters::Lij(b, a) = row[2].REAL_value;
				}
			}
			for (size_t a = 0; a < main_field::phi_number; ++a)
				for (size_t b = a + 1; b < main_field::phi_number; ++b)
					if (parameters::Lij(a, b) < 0) input_error("Interface mobilities must be non-negative.");
			parameters::interface_mobility_anisotropy_delta = 0;
			infile_reader::read_real_value("Model.DS.Phi.IntMobility.Anisotropic.delta",
				parameters::interface_mobility_anisotropy_delta, true);
			if (!std::isfinite(parameters::interface_mobility_anisotropy_delta)
				|| parameters::interface_mobility_anisotropy_delta <= REAL(-1)
				|| parameters::interface_mobility_anisotropy_delta >= REAL(1.5))
				input_error("Phi.IntMobility.Anisotropic.delta must be finite and satisfy -1 < delta < 1.5.");
			WriteDebugFile("# DS default interface mobility: L_ab(n) = L_ab^0 * [1-delta_L*(1.5-2.5*sum_i(n_i^c)^4)].\n");
			WriteDebugFile("# DS interface-mobility orientation uses canonical a<b and n^c=R_b*n; an undefined normal uses L_ab^0.\n");

			REAL sigma_const = 0;
			infile_reader::read_real_value("Model.DS.Phi.InterfaceEnergy.const", sigma_const, true);
			parameters::sigma_ab.resize(main_field::phi_number, main_field::phi_number, sigma_const);
			std::string sigma_input = "[()]";
			if (infile_reader::read_string_value("Model.DS.Phi.InterfaceEnergy.matrix", sigma_input, true)) {
				auto values = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(
					{ InputValueType::IVType_INT, InputValueType::IVType_INT, InputValueType::IVType_REAL },
					"Model.DS.Phi.InterfaceEnergy.matrix", sigma_input, true);
				for (const auto& row : values) {
					const size_t a = row[0].int_value, b = row[1].int_value;
					if (a >= main_field::phi_number || b >= main_field::phi_number || a == b) input_error("InterfaceEnergy.matrix contains an invalid phase pair.");
					parameters::sigma_ab(a, b) = parameters::sigma_ab(b, a) = row[2].REAL_value;
				}
			}
			for (size_t a = 0; a < main_field::phi_number; ++a)
				for (size_t b = a + 1; b < main_field::phi_number; ++b)
					if (parameters::sigma_ab(a, b) < 0) input_error("Interface energy coefficients must be non-negative.");

			infile_reader::read_real_value("Model.DS.Phi.InterfaceEnergy.int_width", parameters::interface_width, true);
			infile_reader::read_real_value("Model.DS.Phi.InterfaceEnergy.Anisotropic.delta", parameters::anisotropy_delta, true);
			infile_reader::read_real_value("Model.DS.Phi.TripleJunctionEnergy.const", parameters::triple_junction_energy, true);
			if (parameters::interface_width <= 0) input_error("Phi.InterfaceEnergy.int_width must be positive.");
			if (parameters::anisotropy_delta <= REAL(-1.5) || parameters::anisotropy_delta >= REAL(1.0))
				input_error("Phi.InterfaceEnergy.Anisotropic.delta must satisfy -1.5 < delta < 1 so sigma remains positive.");

			infile_reader::read_real_value("Model.DS.Phi.DrivingForce.Teq", parameters::equilibrium_temperature, true);
			infile_reader::read_real_value("Model.DS.Phi.DrivingForce.dGdT", parameters::driving_force_dtemp, true);
			if (!std::isfinite(parameters::equilibrium_temperature) || !std::isfinite(parameters::driving_force_dtemp))
				input_error("Driving-force temperature parameters must be finite.");

			if (parameters::has_components()) {
				const std::string driving_components_key = "Model.DS.Phi.DrivingForce.components";
				std::string driving_components_input = "[()]";
				if (!infile_reader::read_string_value(driving_components_key, driving_components_input, true))
					input_error(driving_components_key + " must define every independent component.");
				auto driving_component_values = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(
					{ InputValueType::IVType_STRING, InputValueType::IVType_REAL, InputValueType::IVType_REAL },
					driving_components_key, driving_components_input, true);
				parameters::liquid_con_equilibrium.assign(parameters::component_number, 0);
				parameters::driving_force_dcon.assign(parameters::component_number, 0);
				std::vector<bool> driving_component_defined(parameters::component_number, false);
				for (const auto& row : driving_component_values) {
					const auto found = std::find(parameters::component_names.begin(), parameters::component_names.end(), row[0].string_value);
					if (found == parameters::component_names.end())
						input_error("Phi.DrivingForce.components contains an unknown component name.");
					const size_t i = size_t(found - parameters::component_names.begin());
					if (driving_component_defined[i])
						input_error("Phi.DrivingForce.components contains a duplicate component.");
					const REAL equilibrium = row[1].REAL_value, derivative = row[2].REAL_value;
					if (!std::isfinite(equilibrium) || equilibrium < 0 || equilibrium > 1 || !std::isfinite(derivative))
						input_error("Driving-force equilibrium concentrations must lie in [0,1] and derivatives must be finite.");
					parameters::liquid_con_equilibrium[i] = equilibrium;
					parameters::driving_force_dcon[i] = derivative;
					driving_component_defined[i] = true;
				}
				for (size_t i = 0; i < parameters::component_number; ++i)
					if (!driving_component_defined[i])
						input_error("Phi.DrivingForce.components must define component " + parameters::component_names[i] + ".");
			}

			if (parameters::has_components()) {
				infile_reader::read_real_value("Model.DS.Con.Mobility.const", parameters::con_mobility_const, true);
				if (!std::isfinite(parameters::con_mobility_const) || parameters::con_mobility_const < 0)
					input_error("Concentration mobility must be finite and non-negative.");
				parameters::con_mobility_field.assign(main_field::con_number, parameters::con_mobility_const);
				std::string con_mobility_input = "[()]";
				if (infile_reader::read_string_value("Model.DS.Con.Mobility.matrix", con_mobility_input, true)) {
					auto values = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(
						{ InputValueType::IVType_STRING, InputValueType::IVType_STRING, InputValueType::IVType_REAL },
						"Model.DS.Con.Mobility.matrix", con_mobility_input, true);
					std::vector<bool> mobility_defined(main_field::con_number, false);
					for (const auto& row : values) {
						const std::string& state = row[0].string_value;
						const auto found = std::find(parameters::component_names.begin(), parameters::component_names.end(), row[1].string_value);
						if (found == parameters::component_names.end())
							input_error("Con.Mobility.matrix contains an unknown component name.");
						const size_t component = size_t(found - parameters::component_names.begin());
						size_t field_index = 0;
						if (state == "liquid") field_index = parameters::liquid_con_index(component);
						else if (state == "solid") field_index = parameters::solid_con_index(component);
						else input_error("Con.Mobility.matrix state must be liquid or solid.");
						if (mobility_defined[field_index])
							input_error("Con.Mobility.matrix contains a duplicate state/component pair.");
						if (!std::isfinite(row[2].REAL_value) || row[2].REAL_value < 0)
							input_error("Concentration mobilities must be finite and non-negative.");
						parameters::con_mobility_field[field_index] = row[2].REAL_value;
						mobility_defined[field_index] = true;
					}
				}
			}

			infile_reader::read_real_value("Model.DS.Temp.Mobility.const", parameters::temp_mobility_const, true);
			if (!std::isfinite(parameters::temp_mobility_const) || parameters::temp_mobility_const < 0)
				input_error("Temperature mobility must be finite and non-negative.");
			parameters::temp_mobility_phase.clear();
			std::string temp_mobility_input = "[()]";
			if (infile_reader::read_string_value("Model.DS.Temp.Mobility.matrix", temp_mobility_input, true)) {
				parameters::temp_mobility_phase.resize(main_field::phi_number, parameters::temp_mobility_const);
				auto values = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(
					{ InputValueType::IVType_INT, InputValueType::IVType_REAL }, "Model.DS.Temp.Mobility.matrix", temp_mobility_input, true);
				for (const auto& row : values) {
					const size_t i = row[0].int_value; if (i >= main_field::phi_number) input_error("Temp.Mobility.matrix contains an invalid phase index.");
					parameters::temp_mobility_phase[i] = row[1].REAL_value;
				}
				for (REAL value : parameters::temp_mobility_phase)
					if (value < 0) input_error("Temperature phase mobilities must be non-negative.");
			}

			parameters::latent_heat_phase.resize(main_field::phi_number, 0);
			std::string latent_input = "[()]";
			if (infile_reader::read_string_value("Model.DS.Temp.Source.dphidtemp", latent_input, true)) {
				auto values = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(
					{ InputValueType::IVType_INT, InputValueType::IVType_REAL }, "Model.DS.Temp.Source.dphidtemp", latent_input, true);
				for (const auto& row : values) {
					const size_t i = row[0].int_value; if (i >= main_field::phi_number) input_error("Temp.Source.dphidtemp contains an invalid phase index.");
					parameters::latent_heat_phase[i] = row[1].REAL_value;
				}
			}

			if (!parameters::interface_mobility) parameters::interface_mobility = field_functions::interface_mobility_const;
			if (parameters::has_components() && !parameters::con_mobility)
				parameters::con_mobility = field_functions::concentration_mobility_const;
			if (!parameters::temp_mobility) parameters::temp_mobility = field_functions::temperature_mobility_mixture;

			load_a_new_module(nullptr, nullptr, dendrite_solidification_model::exec_pre_iii,
				dendrite_solidification_model::exec_i, nullptr, nullptr,
				nullptr, nullptr, nullptr, dendrite_solidification_model::deinit);
		}
	}
}
