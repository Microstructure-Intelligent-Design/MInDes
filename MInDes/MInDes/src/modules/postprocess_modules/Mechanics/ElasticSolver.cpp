#include "ElasticSolver.h"
namespace pf {
	namespace elastic_solver {
		// boundary condition
		void change_fix_boundaty_condition_implicity() {
			if (mechanical_field_solver_im.AppStrainMask[0]) {
				for (auto vec = fix_boundary_x_change_rate.begin(); vec < fix_boundary_x_change_rate.end(); vec++)
					if (time_parameters::Real_Time >= (*vec)[0] && time_parameters::Real_Time < (*vec)[1]) {
						mechanical_field_solver_im.applied_strain[0] += (*vec)[2] * time_parameters::delt_t;
						break;
					}
			}
			else if (mechanical_field_solver_im.LoadStressMask[0]) {
				for (auto vec = fix_boundary_x_change_rate.begin(); vec < fix_boundary_x_change_rate.end(); vec++)
					if (time_parameters::Real_Time >= (*vec)[0] && time_parameters::Real_Time < (*vec)[1]) {
						mechanical_field_solver_im.applied_stress[0] += (*vec)[2] * time_parameters::delt_t;
						break;
					}
			}
			if (mechanical_field_solver_im.AppStrainMask[1]) {
				for (auto vec = fix_boundary_y_change_rate.begin(); vec < fix_boundary_y_change_rate.end(); vec++)
					if (time_parameters::Real_Time >= (*vec)[0] && time_parameters::Real_Time < (*vec)[1]) {
						mechanical_field_solver_im.applied_strain[1] += (*vec)[2] * time_parameters::delt_t;
						break;
					}
			}
			else if (mechanical_field_solver_im.LoadStressMask[1]) {
				for (auto vec = fix_boundary_y_change_rate.begin(); vec < fix_boundary_y_change_rate.end(); vec++)
					if (time_parameters::Real_Time >= (*vec)[0] && time_parameters::Real_Time < (*vec)[1]) {
						mechanical_field_solver_im.applied_stress[1] += (*vec)[2] * time_parameters::delt_t;
						break;
					}
			}
			if (mechanical_field_solver_im.AppStrainMask[2]) {
				for (auto vec = fix_boundary_z_change_rate.begin(); vec < fix_boundary_z_change_rate.end(); vec++)
					if (time_parameters::Real_Time >= (*vec)[0] && time_parameters::Real_Time < (*vec)[1]) {
						mechanical_field_solver_im.applied_strain[2] += (*vec)[2] * time_parameters::delt_t;
						break;
					}
			}
			else if (mechanical_field_solver_im.LoadStressMask[2]) {
				for (auto vec = fix_boundary_z_change_rate.begin(); vec < fix_boundary_z_change_rate.end(); vec++)
					if (time_parameters::Real_Time >= (*vec)[0] && time_parameters::Real_Time < (*vec)[1]) {
						mechanical_field_solver_im.applied_stress[2] += (*vec)[2] * time_parameters::delt_t;
						break;
					}
			}
		}
		// get infomation
		vStrain get_applied_strain() {
			if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach 
				|| MFType == MechanicalFieldType::MFType_Implicit_Khachaturyan)
				return mechanical_field_solver_im.applied_strain;
			return vStrain();
		}
		vStress get_applied_stress() {
			if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach 
				|| MFType == MechanicalFieldType::MFType_Implicit_Khachaturyan)
				return mechanical_field_solver_im.applied_stress;
			return vStress();
		}
		// eigenstrain and stiffness
		Matrix6x6 cal_stiffness(int x, int y, int z) {
			Matrix6x6 stiffness;
			stiffness_eigenstrain::stiffness(x, y, z, stiffness);
			return stiffness;
		}
		vStrain cal_eigenstrain(int x, int y, int z) {
			vStrain eigenstrain;
			stiffness_eigenstrain::eigen_strain(x, y, z, eigenstrain);
			return eigenstrain;
		}
		// solver loop for MFType_Implicit
		void exec_pre_im_steinbach() {
			std::stringstream output;
			output << "> Mechanical field solver debug:" << std::endl;
			WriteLog(output.str());
			// =========================================================
			mechanical_field_solver_im.SetMAXElasticConstants(stiffness_eigenstrain::phi_index_stiffness);
			mechanical_field_solver_im.cal_eigenstrain_stiffness();
			mechanical_field_solver_im.initStrainIncrements();
			int map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			std::string elastic_solver_output = "";
			do
			{
				map_steps++;
				output.str("");
				output << "> Map: " << map_steps << std::endl;
				elastic_solver_output = mechanical_field_solver_im.Solve(solver_strain_accuracy, 
					solver_max_iterate_times, bc_incre_rate, solver_debug, output_displacement_field);
				output << elastic_solver_output;
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain();
					output << "  dplastic strain:     " << dplastic_strain << std::endl
						<< "  dave plastic strain: " << dave_plastic_strain << std::endl;
				}
				WriteLog(output.str());
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy);
		}
		std::string exec_loop_im_steinbach() {
			std::stringstream output;
			change_fix_boundaty_condition_implicity();
			mechanical_field_solver_im.cal_eigenstrain_stiffness();
			if (restart_iterator_in_loop && main_iterator::Current_ITE_step % restart_iterator_in_loop_steps == 0)
				mechanical_field_solver_im.initStrainIncrements();
			int map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				elastic_solver_output = mechanical_field_solver_im.Solve(solver_strain_accuracy, solver_max_iterate_times, bc_incre_rate, false, output_displacement_field);
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain();
				}
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy);

			output << "  Mapping steps:       " << map_steps << std::endl
				<< elastic_solver_output
				<< "  dplastic strain:     " << dplastic_strain << std::endl
				<< "  dave plastic strain: " << dave_plastic_strain << std::endl;
			return output.str();
		}
		void exec_pre_im_khachaturyan() {
			stringstream output;
			if (solver_debug) {
				output << "> Mechanical field solver debug:" << std::endl;
				WriteLog(output.str());
			}
			mechanical_field_solver_im.SetMAXElasticConstants(stiffness_eigenstrain::phi_index_stiffness);
			mechanical_field_solver_im.cal_eigenstrain_stiffness();
			mechanical_field_solver_im.initVirtualEigenstrain();
			int map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				output.str("");
				output << "> Map: " << map_steps << std::endl;
				elastic_solver_output = mechanical_field_solver_im.Solve2(solver_strain_accuracy, solver_max_iterate_times, virtual_strain_iterate_rate, solver_debug, output_displacement_field);
				output << elastic_solver_output;
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain();
					output << "  dplastic strain:     " << dplastic_strain << std::endl
						<< "  dave plastic strain: " << dave_plastic_strain << std::endl;
				}
				WriteLog(output.str());
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy);
		}
		std::string exec_loop_im_khachaturyan() {
			stringstream output;
			change_fix_boundaty_condition_implicity();
			mechanical_field_solver_im.cal_eigenstrain_stiffness();
			if (restart_iterator_in_loop && main_iterator::Current_ITE_step % restart_iterator_in_loop_steps == 0)
				mechanical_field_solver_im.initVirtualEigenstrain();
			int map_steps = 0;
			double dplastic_strain = 0.0, dave_plastic_strain = 0.0;
			string elastic_solver_output = "";
			do
			{
				map_steps++;
				elastic_solver_output = mechanical_field_solver_im.Solve2(solver_strain_accuracy, solver_max_iterate_times, 
					virtual_strain_iterate_rate, false, output_displacement_field);
				if (is_plastic_on) {
					plastic_solver::solve_plastic_flow(dplastic_strain, dave_plastic_strain);
					mechanical_field_solver_im.recal_eigenstrain();
				}
				if (map_steps >= mechanic_map_steps)
					break;
			} while (dplastic_strain > solver_strain_accuracy || dave_plastic_strain > solver_strain_accuracy);

			output << "  Mapping steps:       " << map_steps << std::endl
				<< elastic_solver_output
				<< "  dplastic strain:     " << dplastic_strain << std::endl
				<< "  dave plastic strain: " << dave_plastic_strain << std::endl;
			return output.str();
		}

		void init() {
			const size_t AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2;
			WriteDebugFile("# Postprocess.SolidMechanics.momentum_balance = 0 - None , 1 - Implicit (Ingo Steinbach) , 2 - Implicit (Armen G. Khachaturyan) \n");
			int _MFType = MechanicalFieldType::MFType_None;
			infile_reader::read_int_value("Postprocess.SolidMechanics.momentum_balance", _MFType, true);
			MFType = MechanicalFieldType(_MFType);
			if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach) {
				pf::BoundaryCondition bc_x = pf::BoundaryCondition::PERIODIC, bc_y = pf::BoundaryCondition::PERIODIC, bc_z = pf::BoundaryCondition::PERIODIC;
				if (mesh_parameters::x_down != pf::BoundaryCondition::PERIODIC
					|| mesh_parameters::x_up != pf::BoundaryCondition::PERIODIC)
					bc_x = pf::BoundaryCondition::ZEROFLUX;
				if (mesh_parameters::y_down != pf::BoundaryCondition::PERIODIC
					|| mesh_parameters::y_up != pf::BoundaryCondition::PERIODIC)
					bc_y = pf::BoundaryCondition::ZEROFLUX;
				if (mesh_parameters::z_down != pf::BoundaryCondition::PERIODIC
					|| mesh_parameters::z_up != pf::BoundaryCondition::PERIODIC)
					bc_z = pf::BoundaryCondition::ZEROFLUX;
				mechanical_field_solver_im.init(int(mesh_parameters::MESH_NX), int(mesh_parameters::MESH_NY), int(mesh_parameters::MESH_NZ), 
					bc_x, bc_y, bc_z, main_field::elastic_field);
				main_field::is_mech_field_elas_on = true;
				mechanical_field_solver_im.cal_stiffness = cal_stiffness;
				mechanical_field_solver_im.cal_eigenstrain = cal_eigenstrain;
				if (infile_reader::read_bool_value("Postprocess.SolidMechanics.restart_iterator_in_loop", restart_iterator_in_loop, true))
					infile_reader::read_int_value("Postprocess.SolidMechanics.RestartIteratorInLoop.steps", restart_iterator_in_loop_steps, true);
				infile_reader::read_real_value("Postprocess.SolidMechanics.Implicit.bc_ite_rate", bc_incre_rate, true);
				fix_domain_boundary.resize(3);
				WriteDebugFile("# Postprocess.SolidMechanics.fix_boundary.type = (BC_X, BC_Y, BC_Z) , 0 - Average , 1 - Strain , 2 - Stress \n");
				string fix_boundary_key = "Postprocess.SolidMechanics.fix_boundary.type", fix_boundary_input = "(0,0,0)";
				if (infile_reader::read_string_value(fix_boundary_key, fix_boundary_input, true)) {
					vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, fix_boundary_key, fix_boundary_input, true);
					fix_domain_boundary[AXIS_X] = FixBoundaryCondition(fix_boundary_value[0].int_value);
					fix_domain_boundary[AXIS_Y] = FixBoundaryCondition(fix_boundary_value[1].int_value);
					fix_domain_boundary[AXIS_Z] = FixBoundaryCondition(fix_boundary_value[2].int_value);
					for (int direction = 0; direction < 3; direction++) {
						if (fix_domain_boundary[direction] == 1) {
							mechanical_field_solver_im.AppStrainMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (direction)
							{
							case AXIS_X:
								_d = "x";
								break;
							case AXIS_Y:
								_d = "y";
								break;
							case AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d;
							double strain = 0.0;
							infile_reader::read_real_value(fix_boundary_val_key, strain, true);
							mechanical_field_solver_im.applied_strain[direction] = strain;
							// [(time_begin,time_end,change rate), ... ]
							WriteDebugFile("# Postprocess.SolidMechanics.fix_boundary.strain_? = [(real_time_begin, real_time_end, dstrain_dt), ... ] \n");
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (infile_reader::read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, true)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, fix_boundary_rate_key, fix_boundary_rate_input, true);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].REAL_value;
									vec[1] = fix_boundary_rate_value[index][1].REAL_value;
									vec[2] = fix_boundary_rate_value[index][2].REAL_value;
									switch (direction)
									{
									case AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
						else if (fix_domain_boundary[direction] == 2) {
							mechanical_field_solver_im.LoadStressMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (direction)
							{
							case AXIS_X:
								_d = "x";
								break;
							case AXIS_Y:
								_d = "y";
								break;
							case AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d;
							double stress = 0.0;
							infile_reader::read_real_value(fix_boundary_val_key, stress, true);
							mechanical_field_solver_im.applied_stress[direction] = stress;
							// [(time_begin,time_end,change rate), ... ]
							WriteDebugFile("# Postprocess.SolidMechanics.fix_boundary.stress_? = [(real_time_begin, real_time_end, dstress_dt), ... ] \n");
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (infile_reader::read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, true)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, fix_boundary_rate_key, fix_boundary_rate_input, true);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].REAL_value;
									vec[1] = fix_boundary_rate_value[index][1].REAL_value;
									vec[2] = fix_boundary_rate_value[index][2].REAL_value;
									switch (direction)
									{
									case AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
					}
				}
			}
			else if (MechanicalFieldType(MFType) == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				pf::BoundaryCondition bc_x = pf::BoundaryCondition::PERIODIC, bc_y = pf::BoundaryCondition::PERIODIC, bc_z = pf::BoundaryCondition::PERIODIC;
				if (mesh_parameters::x_down != pf::BoundaryCondition::PERIODIC
					|| mesh_parameters::x_up != pf::BoundaryCondition::PERIODIC)
					bc_x = pf::BoundaryCondition::ZEROFLUX;
				if (mesh_parameters::y_down != pf::BoundaryCondition::PERIODIC
					|| mesh_parameters::y_up != pf::BoundaryCondition::PERIODIC)
					bc_y = pf::BoundaryCondition::ZEROFLUX;
				if (mesh_parameters::z_down != pf::BoundaryCondition::PERIODIC
					|| mesh_parameters::z_up != pf::BoundaryCondition::PERIODIC)
					bc_z = pf::BoundaryCondition::ZEROFLUX;
				mechanical_field_solver_im.init(int(mesh_parameters::MESH_NX), int(mesh_parameters::MESH_NY), int(mesh_parameters::MESH_NZ)
					, bc_x, bc_y, bc_z, main_field::elastic_field);
				main_field::is_mech_field_elas_on = true;
				mechanical_field_solver_im.cal_stiffness = cal_stiffness;
				mechanical_field_solver_im.cal_eigenstrain = cal_eigenstrain;
				infile_reader::read_real_value("Postprocess.SolidMechanics.VisualStrain.L_ijkl", virtual_strain_iterate_rate, true);
				if (infile_reader::read_bool_value("Postprocess.SolidMechanics.restart_iterator_in_loop", restart_iterator_in_loop, true))
					infile_reader::read_int_value("Postprocess.SolidMechanics.RestartIteratorInLoop.steps", restart_iterator_in_loop_steps, true);
				fix_domain_boundary.resize(3);
				WriteDebugFile("# Postprocess.SolidMechanics.fix_boundary.type = (BC_X, BC_Y, BC_Z) , 0 - Average , 1 - Strain , 2 - Stress \n");
				string fix_boundary_key = "Postprocess.SolidMechanics.fix_boundary.type", fix_boundary_input = "(0,0,0)";
				if (infile_reader::read_string_value(fix_boundary_key, fix_boundary_input, true)) {
					vector<input_value> fix_boundary_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, fix_boundary_key, fix_boundary_input, true);
					fix_domain_boundary[AXIS_X] = FixBoundaryCondition(fix_boundary_value[0].int_value);
					fix_domain_boundary[AXIS_Y] = FixBoundaryCondition(fix_boundary_value[1].int_value);
					fix_domain_boundary[AXIS_Z] = FixBoundaryCondition(fix_boundary_value[2].int_value);
					for (int direction = 0; direction < 3; direction++) {
						if (fix_domain_boundary[direction] == 1) {
							mechanical_field_solver_im.AppStrainMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (direction) 
							{
							case AXIS_X:
								_d = "x";
								break;
							case AXIS_Y:
								_d = "y";
								break;
							case AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d;
							double strain = 0.0;
							infile_reader::read_real_value(fix_boundary_val_key, strain, true);
							mechanical_field_solver_im.applied_strain[direction] = strain;
							// [(time_begin,time_end,change rate), ... ]
							WriteDebugFile("# Postprocess.SolidMechanics.fix_boundary.strain_? = [(real_time_begin, real_time_end, dstrain_dt), ... ] \n");
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.strain_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (infile_reader::read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, true)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, fix_boundary_rate_key, fix_boundary_rate_input, true);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].REAL_value;
									vec[1] = fix_boundary_rate_value[index][1].REAL_value;
									vec[2] = fix_boundary_rate_value[index][2].REAL_value;
									switch (direction)
									{
									case AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
						else if (fix_domain_boundary[direction] == 2) {
							mechanical_field_solver_im.LoadStressMask[direction] = true;
							mechanical_field_solver_im.AvgStrainMask[direction] = false;
							string _d = "x";
							switch (direction)
							{
							case AXIS_X:
								_d = "x";
								break;
							case AXIS_Y:
								_d = "y";
								break;
							case AXIS_Z:
								_d = "z";
								break;
							}
							string fix_boundary_val_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d;
							double stress = 0.0;
							infile_reader::read_real_value(fix_boundary_val_key, stress, true);
							mechanical_field_solver_im.applied_stress[direction] = stress;
							// [(time_begin,time_end,change rate), ... ]
							WriteDebugFile("# Postprocess.SolidMechanics.fix_boundary.stress_? = [(real_time_begin, real_time_end, dstress_dt), ... ] \n");
							string fix_boundary_rate_key = "Postprocess.SolidMechanics.fix_boundary.stress_" + _d + ".rate", fix_boundary_rate_input = "[()]";
							if (infile_reader::read_string_value(fix_boundary_rate_key, fix_boundary_rate_input, true)) {
								vector<vector<input_value>> fix_boundary_rate_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, fix_boundary_rate_key, fix_boundary_rate_input, true);
								for (int index = 0; index < fix_boundary_rate_value.size(); index++) {
									Vector3 vec;
									vec[0] = fix_boundary_rate_value[index][0].REAL_value;
									vec[1] = fix_boundary_rate_value[index][1].REAL_value;
									vec[2] = fix_boundary_rate_value[index][2].REAL_value;
									switch (direction)
									{
									case AXIS_X:
										fix_boundary_x_change_rate.push_back(vec);
										break;
									case AXIS_Y:
										fix_boundary_y_change_rate.push_back(vec);
										break;
									case AXIS_Z:
										fix_boundary_z_change_rate.push_back(vec);
										break;
									}
								}
							}
						}
					}
				}
			}
			infile_reader::read_bool_value("Postprocess.SolidMechanics.write_displacement_field", output_displacement_field, true);
			infile_reader::read_int_value("Postprocess.SolidMechanics.max_iteration_steps", solver_max_iterate_times, true);
			infile_reader::read_bool_value("Postprocess.SolidMechanics.debug", solver_debug, true);
			infile_reader::read_real_value("Postprocess.SolidMechanics.strain_accuracy", solver_strain_accuracy, true);
			infile_reader::read_bool_value("Postprocess.SolidMechanics.plasticity", is_plastic_on, true);
			if (is_plastic_on) {
				infile_reader::read_int_value("Postprocess.SolidMechanics.Elastoplasticity.mapping_steps", mechanic_map_steps, true);
				if (mechanic_map_steps < 1)
					mechanic_map_steps = 1;
			}
			WriteLog("> MODULE INIT : Elastic Solver \n");
		}

		void exec_pre() {
			if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach) {
				exec_pre_im_steinbach();
			}
			else if (MFType == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				exec_pre_im_khachaturyan();
			}
		}

		std::string exec_loop() {
			if (MFType == MechanicalFieldType::MFType_Implicit_Steinbach) {
				return exec_loop_im_steinbach();
			}
			else if (MFType == MechanicalFieldType::MFType_Implicit_Khachaturyan) {
				return exec_loop_im_khachaturyan();
			}
			return "";
		}

		void deinit() {
			mechanical_field_solver_im.free();
		}

		void write_vts(std::ofstream& fout) {
			vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
			for (int ele_index = 0; ele_index < 6; ele_index++)
			{
				string compname = "\"stress_" + compNameV[ele_index] + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							if (i > 0 && i <= mesh_parameters::MESH_NX && j > 0 && j <= mesh_parameters::MESH_NY
								&& k > 0 && k <= mesh_parameters::MESH_NZ) {
								ElasticPoint& point = main_field::elastic_field(i - 1, j - 1, k - 1);
								fout << point.Stress[ele_index] << std::endl;
							}
							else {
								fout << 0 << std::endl;
							}
						}
				fout << "</DataArray>" << std::endl;
			}
			for (int ele_index = 0; ele_index < 6; ele_index++)
			{
				string compname = "\"strain_" + compNameV[ele_index] + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							if (i > 0 && i <= mesh_parameters::MESH_NX && j > 0 && j <= mesh_parameters::MESH_NY
								&& k > 0 && k <= mesh_parameters::MESH_NZ) {
								ElasticPoint& point = main_field::elastic_field(i - 1, j - 1, k - 1);
								fout << point.Strain[ele_index] << std::endl;
							}
							else {
								fout << 0 << std::endl;
							}
						}
				fout << "</DataArray>" << std::endl;
			}
			for (int ele_index = 0; ele_index < 6; ele_index++)
			{
				string compname = "\"nonElasticStrain_" + compNameV[ele_index] + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << compname <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							if (i > 0 && i <= mesh_parameters::MESH_NX && j > 0 && j <= mesh_parameters::MESH_NY
								&& k > 0 && k <= mesh_parameters::MESH_NZ) {
								ElasticPoint& point = main_field::elastic_field(i - 1, j - 1, k - 1);
								fout << point.EffectiveEigenStrain[ele_index] << std::endl;
							}
							else {
								fout << 0 << std::endl;
							}
						}
				fout << "</DataArray>" << std::endl;
			}
			fout << "<DataArray type = \"Float64\" Name = \"" << "stress_J1" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
			for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						if (i > 0 && i <= mesh_parameters::MESH_NX && j > 0 && j <= mesh_parameters::MESH_NY
							&& k > 0 && k <= mesh_parameters::MESH_NZ) {
							ElasticPoint& point = main_field::elastic_field(i - 1, j - 1, k - 1);
							fout << point.Stress.J1() << std::endl;
						}
						else {
							fout << 0 << std::endl;
						}
					}
			fout << "</DataArray>" << std::endl;
			fout << "<DataArray type = \"Float64\" Name = \"" << "stress_vMises" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
			for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						if (i > 0 && i <= mesh_parameters::MESH_NX && j > 0 && j <= mesh_parameters::MESH_NY
							&& k > 0 && k <= mesh_parameters::MESH_NZ) {
							ElasticPoint& point = main_field::elastic_field(i - 1, j - 1, k - 1);
							fout << point.Stress.Mises() << std::endl;
						}
						else {
							fout << 0 << std::endl;
						}
					}
			fout << "</DataArray>" << std::endl;
			if (output_displacement_field) {
				string name = "mech_displacement";
				fout << "<DataArray type = \"Float64\" Name = \"" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							if (i > 0 && i <= mesh_parameters::MESH_NX && j > 0 && j <= mesh_parameters::MESH_NY
								&& k > 0 && k <= mesh_parameters::MESH_NZ) {
								fout << mechanical_field_solver_im.get_u_main_node(i - 1, j - 1, k - 1) << " "
									<< mechanical_field_solver_im.get_v_main_node(i - 1, j - 1, k - 1) << " "
									<< mechanical_field_solver_im.get_w_main_node(i - 1, j - 1, k - 1) << std::endl;
							}
							else {
								fout << 0 << " "
									<< 0 << " "
									<< 0 << std::endl;
							}
						}
				fout << "</DataArray>" << std::endl;
			}
		}
	}
}