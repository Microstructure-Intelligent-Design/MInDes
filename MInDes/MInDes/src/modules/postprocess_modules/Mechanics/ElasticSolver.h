#pragma once
#include "GoverningEquationsM.h"
#include "StiffnessEigenStrain.h"
#include "PlasticSolver.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../Modules_Params.h"
#include "../../../MainIterator_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
namespace pf {
	enum MechanicalFieldType { MFType_None, MFType_Implicit_Steinbach, MFType_Implicit_Khachaturyan };
	enum FixBoundaryCondition { FixBC_Average, FixBC_Strain, FixBC_Stress };
	namespace elastic_solver {
		inline MechanicalField_Implicit mechanical_field_solver_im;
		inline std::vector<size_t> calculation_step;
		// General Parameters
		inline MechanicalFieldType MFType = MechanicalFieldType::MFType_None;
		inline double bc_incre_rate = 1.0;
		//-
		inline std::vector<int> fix_domain_boundary;
		inline int solver_max_iterate_times = 0;
		inline double solver_strain_accuracy = 1e-3;
		inline bool solver_debug = false;
		inline std::vector<Vector3> fix_boundary_x_change_rate;
		inline std::vector<Vector3> fix_boundary_y_change_rate;
		inline std::vector<Vector3> fix_boundary_z_change_rate;
		// MFType_Implicit Parameters
		inline double virtual_strain_iterate_rate = 1.0;
		inline bool restart_iterator_in_loop = false;
		inline int restart_iterator_in_loop_steps = 1;
		inline bool is_displacement_field_output = false;
		// plastic solver
		inline int mechanic_map_steps = 1;
		// - statistic
		const std::pair<std::string, std::string> statistic_app_strain = { "app_strain", "applied strain on each direction" };
		const std::pair<std::string, std::string> statistic_app_stress = { "app_stress", "applied stress on each direction" };
		const std::pair<std::string, std::string> statistic_ave_strain = { "ave_strain", "average strain on each direction" };
		const std::pair<std::string, std::string> statistic_ave_stress = { "ave_stress", "average stress on each direction" };
		const std::pair<std::string, std::string> statistic_max_vMises = { "max_vMises_stress", "max von Mises stress" };
		const std::pair<std::string, std::string> statistic_ave_plas_strain = { "ave_plas_strain", "average cumulative plastic strain" };
		inline bool is_app_strain_statistic = false;
		inline bool is_app_stress_statistic = false;
		inline bool is_ave_strain_statistic = false;
		inline bool is_ave_stress_statistic = false;
		inline bool is_ave_plas_strain_statistic = false;
		inline bool is_max_vMises_stress_statistic = false;
		// boundary condition
		void change_fix_boundaty_condition_implicity();
		// get infomation
		vStrain get_applied_strain();
		vStress get_applied_stress();
		// eigenstrain and stiffness
		Matrix6x6 cal_stiffness(long long x, long long y, long long z);
		vStrain cal_eigenstrain(long long x, long long y, long long z);
		// solver loop for MFType_Implicit
		void exec_pre_im_steinbach();
		void exec_loop_im_steinbach();
		void exec_pre_im_khachaturyan();
		void exec_loop_im_khachaturyan();

		void init();

		void exec_pre_i();

		void exec_pre();

		void exec_loop();

		void deinit();

		void write_vts_displacement(std::ofstream& fout);

	}
}