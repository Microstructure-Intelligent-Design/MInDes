#pragma once
#include "../Module.h"
#include "Model_Params.h"
#include "grain_grows_spinodal/Model_Manager.h"
#include "../../MainIterator_Params.h"
namespace pf {
	enum SimulationModels { SM_None, SM_GrainGrowsSpinodal };
	inline void init_basic_time_mesh_parameters();
	inline void init_model_modules() {
		// parallel
		infile_reader::read_int_value("Solver.Loop.OpenMP_Thread", main_iterator::OpenMP_Thread_Counts, true);
		// - mesh and time parameters
		init_basic_time_mesh_parameters();
		// - init models
		WriteDebugFile("# SimulationModels.model =  0 - None \n");
		WriteDebugFile("#                           1 - Grain Grows Spinodal , PCT = (N,1,false) \n");
		int sm_model = SimulationModels::SM_None;
		infile_reader::read_int_value("SimulationModels.model", sm_model, true);
		switch (SimulationModels(sm_model)) {
		case SimulationModels::SM_GrainGrowsSpinodal: {
			// - model settings
			if (main_field::phi_number == 0 || main_field::con_number != 1 || main_field::is_temp_field_on != false) {
				string error_report = "> ERROR, the PCT command line settings do not meet the requirements : PCT = (N,1,false) !\n";
				WriteLog(error_report);
				std::exit(0);
			}
			grain_grows_spinodal_model::init_model_modules();
			break;
		}
		}
		// - 
		WriteLog("> MODULE INIT : Models Ready !\n");
	}

	inline void init_basic_time_mesh_parameters() {
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
		int boundary_condition = 0;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_up", boundary_condition, true);
		mesh_parameters::x_up = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_down", boundary_condition, true);
		mesh_parameters::x_down = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_up", boundary_condition, true);
		mesh_parameters::y_up = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_down", boundary_condition, true);
		mesh_parameters::y_down = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_up", boundary_condition, true);
		mesh_parameters::z_up = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_dowm", boundary_condition, true);
		mesh_parameters::z_down = BoundaryCondition(boundary_condition);
		// - 
		string pct_field_key = "Solver.Mesh.PCT", pct_field_string = "(0,0,false)";
		WriteDebugFile("# Solver.Mesh.PCT = ( phi_number, con_number, is_temperature_on ) \n");
		if (InputFileReader::get_instance()->read_string_value(pct_field_key, pct_field_string, true)) {
			vector<input_value> pct_field_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(
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
	}
}
