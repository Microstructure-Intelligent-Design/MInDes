#pragma once
#include "../Module.h"
#include "inputfiles/QuickStartUp.h"
#include "inputfiles/infile_selector.h"
#include "ioFiles_Params.h"
#include "inputfiles/InputFileReader.h"
#include "inputfiles/UserStartUp.h"
#include "../../MainIterator_Params.h"
#include "../Modules_Params.h"
namespace pf {
	static void init_input_modules() {
		stringstream out; 
		out << "============================== Parameters Definition Format =============================" << endl;
		out << "                             Parameter_Name = Parameter_Value" << endl;
		out << "======================================= M a c r o =======================================" << endl;
		out << "1.tube               " << ", $TUBE[1,10,2]$ = 1,3,5,7,9" << endl;
		out << "2.rand               " << ", $RAND_INT[1,10]$ = 1 - 10" << endl;
		out << "                     " << ", $RAND_REAL[1,10]$ = 1.000000 - 10.000000" << endl;
		out << "========================= Define Custom Variables and Functions =========================" << endl;
		out << "Define.Var  = VarName,0.1" << endl;
		out << "Define.Func = FuncName@{[VarName*pow(VarName,2)]}@" << endl;
		out << "default functions      : \"pow(val, ord)\", \"sqrt(val)\", \"abs(val)\", \"exp(val)\", \"ln(val)\"," << endl;
		out << "                         \"log(base_val, val)\", \"sin(val)\", \"cos(val)\", \"tan(val)\"," << endl;
		out << "                         \"asin(val)\", \"acos(val)\", \"atan(val)\", \"cos(val)\", \"tan(val)\"," << endl;
		out << "                         \"asin(val)\", \"acos(val)\", \"atan(val)\", \"rand_int(min_val, max_val)\"," << endl;
		out << "                         \"rand_real(min_val, max_val)\"" << endl;
		out << "=========================================================================================" << endl;
		add_string_to_file(out.str(), input_output_files_parameters::DebugFile_Path);
		InputFileReader::get_instance()->debug_infile_and_valid_words();
		InputFileReader::get_instance()->debug_custom_variavle_and_funcs();
		out.str("");
		out << "=========================================================================================" << endl;
		out << "================================= Parameters Definition =================================" << endl;
		out << "=========================================================================================" << endl;
		add_string_to_file(out.str(), input_output_files_parameters::DebugFile_Path);
		// parallel
		infile_reader::read_int_value("Solver.Loop.OpenMP_Thread", main_iterator::OpenMP_Thread_Counts, true);
		// - mesh and time parameters
		infile_reader::read_int_value("Solver.Loop.begin_step", main_iterator::ITE_Begin_Step, true);
		infile_reader::read_int_value("Solver.Loop.end_step", main_iterator::ITE_End_Step, true);
		infile_reader::read_real_value("Solver.Loop.RealTime.init", time_parameters::Real_Time, true);
		infile_reader::read_real_value("Solver.Loop.dt", time_parameters::delt_t, true);
		infile_reader::read_int_value("Solver.Mesh.Nx", mesh_parameters::MESH_NX, true);
		infile_reader::read_int_value("Solver.Mesh.Ny", mesh_parameters::MESH_NY, true);
		infile_reader::read_int_value("Solver.Mesh.Nz", mesh_parameters::MESH_NZ, true);
		// dimension rules
		if (mesh_parameters::MESH_NX < 1 || mesh_parameters::MESH_NY < 1 || mesh_parameters::MESH_NZ < 1) {
			string error_report = "> ERROR, one edge length (Nx or Ny or Nz) of domain defined error, domain does not exist ! Please set Nx > 0 and Ny > 0 and Nz > 0 !\n";
			WriteLog(error_report);
			std::exit(0);
		}
		else if (mesh_parameters::MESH_NX == 1 && (mesh_parameters::MESH_NY > 1 || mesh_parameters::MESH_NZ > 1)) {
			string error_report = "> ERROR, edge length (Nx, Ny, Nz) of domain should be set (Nx > 1, Ny = 1, Nz = 1) in 1 dimension or (Nx > 1, Ny > 1, Nz = 1) in 2 dimension !\n";
			WriteLog(error_report);
			std::exit(0);
		}
		else if (mesh_parameters::MESH_NX > 1 && mesh_parameters::MESH_NY == 1 && mesh_parameters::MESH_NZ > 1) {
			string error_report = "> ERROR, edge length (Nx, Ny, Nz) of domain should be set as (Nx > 1, Ny > 1, Nz = 1) in 2 dimension !\n";
			WriteLog(error_report);
			std::exit(0);
		}
		if (mesh_parameters::MESH_NY == 1 && mesh_parameters::MESH_NZ == 1)
			mesh_parameters::dimention = Dimension::One_Dimension;
		else if (mesh_parameters::MESH_NZ == 1)
			mesh_parameters::dimention = Dimension::Two_Dimension;
		else
			mesh_parameters::dimention = Dimension::Three_Dimension;
		InputFileReader::get_instance()->read_REAL_value("Solver.Mesh.dr", mesh_parameters::delt_r, true);
		WriteDebugFile("# Solver.Mesh.BoundaryCondition : 0 - FIXED , 1 - PERIODIC , 2 - ZEROFLUX\n");
		int boundary_condition = 1;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_up", boundary_condition, true);
		mesh_parameters::x_up = BoundaryCondition(boundary_condition);
		boundary_condition = 1;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_down", boundary_condition, true);
		mesh_parameters::x_down = BoundaryCondition(boundary_condition);
		boundary_condition = 1;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_up", boundary_condition, true);
		mesh_parameters::y_up = BoundaryCondition(boundary_condition);
		boundary_condition = 1;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_down", boundary_condition, true);
		mesh_parameters::y_down = BoundaryCondition(boundary_condition);
		boundary_condition = 1;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_up", boundary_condition, true);
		mesh_parameters::z_up = BoundaryCondition(boundary_condition);
		boundary_condition = 1;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_dowm", boundary_condition, true);
		mesh_parameters::z_down = BoundaryCondition(boundary_condition);
		// - 
		string pct_field_key = "Solver.Mesh.PCT", pct_field_string = "(0,0,false)";
		WriteDebugFile("# Solver.Mesh.PCT = ( phi_number, con_number, is_temperature_on ) \n");
		if (InputFileReader::get_instance()->read_string_value(pct_field_key, pct_field_string, true)) {
			std::vector<input_value> pct_field_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(
				{ InputValueType::IVType_INT ,InputValueType::IVType_INT ,InputValueType::IVType_BOOL },
				pct_field_key, pct_field_string, true);
			main_field::phi_number = pct_field_value[0].int_value;
			if (main_field::phi_number > 0) {
				main_field::is_phi_field_on = true;
				main_field::init_phase_field();
			}
			main_field::con_number = pct_field_value[1].int_value;
			if (main_field::con_number > 0) {
				main_field::is_con_field_on = true;
				main_field::init_concentration_field();
			}
			main_field::is_temp_field_on = pct_field_value[2].bool_value;
			if (main_field::is_temp_field_on)
				main_field::init_temperature_field();
		}
		WriteLog("> MODULE INIT : Input File (.mindes) Ready !\n");
	}
}