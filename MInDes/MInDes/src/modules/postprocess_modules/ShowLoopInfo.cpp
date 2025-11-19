#include "ShowLoopInfo.h"
namespace pf {
	namespace show_loop_information {

		std::vector<REAL> statistical_phi() {
			std::vector<REAL> phi_info(main_field::phi_number, 0);
			size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
				for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
					for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
						std::vector<REAL>& phi = main_field::phase_field(x, y, z);
						for (size_t index = 0; index < main_field::phi_number; index++)
							phi_info[index] += phi[index];
					}
			for (size_t index = 0; index < main_field::phi_number; index++)
				phi_info[index] /= SIZE;
			return phi_info;
		}

		std::vector<REAL> statistical_con() {
			std::vector<REAL> con_info(main_field::con_number, 0);
			size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (long long x = main_field::concentration_field.COMP_X_BGN(); x <= main_field::concentration_field.COMP_X_END(); x++)
				for (long long y = main_field::concentration_field.COMP_Y_BGN(); y <= main_field::concentration_field.COMP_Y_END(); y++)
					for (long long z = main_field::concentration_field.COMP_Z_BGN(); z <= main_field::concentration_field.COMP_Z_END(); z++) {
						std::vector<REAL>& con = main_field::concentration_field(x, y, z);
						for (size_t index = 0; index < main_field::con_number; index++)
							con_info[index] += con[index];
					}
			for (size_t index = 0; index < main_field::con_number; index++)
				con_info[index] /= SIZE;
			return con_info;
		}

		REAL statistical_temp() {
			REAL temp = 0;
			size_t SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (long long x = main_field::temperature_field.COMP_X_BGN(); x <= main_field::temperature_field.COMP_X_END(); x++)
				for (long long y = main_field::temperature_field.COMP_Y_BGN(); y <= main_field::temperature_field.COMP_Y_END(); y++)
					for (long long z = main_field::temperature_field.COMP_Z_BGN(); z <= main_field::temperature_field.COMP_Z_END(); z++) {
						temp += main_field::temperature_field(x, y, z);
					}
			return temp / SIZE;
		}

		void exec_pre_iii() {
			std::vector<REAL> phi_info, con_info; REAL temp_info = 0;
			if (main_field::is_phi_field_on)
				phi_info = statistical_phi();
			if (main_field::is_con_field_on)
				con_info = statistical_con();
			if (main_field::is_temp_field_on)
				temp_info = statistical_temp();
			stringstream log;
			log << "> PHI CON TEMP information :" << endl;
			if (main_field::is_phi_field_on)
				for (int index = 0; index < main_field::phi_number; index++)
					log << ">  Phi " << index << " = "
					<< setprecision(5) << phi_info[index] << endl;
			if (main_field::is_con_field_on)
				for (int index = 0; index < main_field::con_number; index++)
					log << ">  Con " << index << " = "
					<< setprecision(5) << con_info[index] << endl;
			if (main_field::is_temp_field_on)
				log << ">  Temp " << " = "
				<< setprecision(5) << temp_info << endl;
			WriteLog(log.str());
		}

		void exec_pos_i() {
			time_parameters::Real_Time += time_parameters::delt_t;
			if (screen_output_step == 0)
				return;
			if (main_iterator::Current_ITE_step % screen_output_step == 0) {
				std::vector<REAL> phi_info, con_info; REAL temp_info = 0;
				if (main_field::is_phi_field_on)
					phi_info = statistical_phi();
				if (main_field::is_con_field_on)
					con_info = statistical_con();
				if (main_field::is_temp_field_on)
					temp_info = statistical_temp();
				stringstream log;
				log << "#====================================================================================================" << endl;
				log << timer::return_cunrrent_time_by_string();
				log.setf(ios::fixed);
				log << "# Simulation step " << main_iterator::Current_ITE_step << " has been finished!" << endl;
				log << "# Real time is " << setprecision(3) << time_parameters::Real_Time << endl;
				log << "# This " << screen_output_step << " steps used " << setprecision(3) << timer::interval_end(main_iterator::t_interval_begin) << "(secs.), " << endl;
				log << "# Total simulation used " << setprecision(3) << timer::total_duration_sec(main_iterator::t_total_begin) << "(secs.)." << endl;
				log << "#----------------------------------------------------------------------------------------------------" << endl;
				if (main_field::is_phi_field_on)
					for (int index = 0; index < main_field::phi_number; index++)
						log << "#  Phi " << index << " = "
						<< setprecision(5) << phi_info[index] << endl;
				if (main_field::is_con_field_on)
					for (int index = 0; index < main_field::con_number; index++)
						log << "#  Con " << index << " = "
						<< setprecision(5) << con_info[index] << endl;
				if (main_field::is_temp_field_on)
					log << "#  Temp " << " = "
					<< setprecision(5) << temp_info << endl;
				WriteLog(log.str());
				timer::interval_begin(main_iterator::t_interval_begin);
			}
		}

		void exec_pos_ii() {
			bool output = false;
			if (screen_loop_step != 0)
				if (main_iterator::Current_ITE_step % screen_loop_step == 0)
					output = true;
			if (screen_output_step != 0)
				if (main_iterator::Current_ITE_step % screen_output_step == 0)
					output = true;
			if (output) {
				stringstream log;
				log << "#------------------------------------------ PCT Field -----------------------------------------------" << endl;
				log << "# CURRENT STEP = " << main_iterator::Current_ITE_step << ", REAL TIME = " << time_parameters::Real_Time << endl;
				if (main_field::is_phi_field_on)
					log << "# MAX PHI  INCREMENT = " << setprecision(5) << main_field::PHI_MAX_VARIATION << endl;
				if (main_field::is_con_field_on)
					log << "# MAX CON  INCREMENT = " << setprecision(5) << main_field::CON_MAX_VARIATION << endl;
				if (main_field::is_temp_field_on)
					log << "# MAX TEMP INCREMENT = " << setprecision(5) << main_field::TEMP_MAX_VARIATION << endl;
				log << "#----------------------------------------------------------------------------------------------------" << endl;
				main_field::PHI_MAX_VARIATION = 0;
				main_field::CON_MAX_VARIATION = 0;
				main_field::TEMP_MAX_VARIATION = 0;
				if (screen_output_step != 0)
					if (main_iterator::Current_ITE_step % screen_output_step != 0)
						log << endl << endl;
				WriteLog(log.str());
			}
		}

		// info end
		void exec_pos_iii() {
			if (screen_output_step == 0)
				return;
			if (main_iterator::Current_ITE_step % screen_output_step == 0) {
				stringstream log;
				log << "#====================================================================================================" << endl;
				log << endl << endl;
				WriteLog(log.str());
			}
		}

		void init_show_loop_information() {
			infile_reader::read_int_value("Solver.Output.LOG.loop_info_step", screen_loop_step, true);
			infile_reader::read_int_value("Solver.Output.LOG.screen_output_step", screen_output_step, true);
			load_a_new_module(nullptr, nullptr, exec_pre_iii,
				nullptr, nullptr, nullptr,
				exec_pos_i, exec_pos_ii, exec_pos_iii, nullptr);
		}
	}
}