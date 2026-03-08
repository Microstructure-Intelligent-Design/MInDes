#pragma once
#include "GoverningEquations.h"
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
		inline bool output_displacement_field = false;
		inline double virtual_strain_iterate_rate = 1.0;
		inline bool restart_iterator_in_loop = false;
		inline int restart_iterator_in_loop_steps = 1;
		// plastic solver
		inline bool is_plastic_on = false;
		inline int mechanic_map_steps = 1;
		// boundary condition
		void change_fix_boundaty_condition_implicity();
		// get infomation
		vStrain get_applied_strain();
		vStress get_applied_stress();
		// eigenstrain and stiffness
		Matrix6x6 cal_stiffness(int x, int y, int z);
		vStrain cal_eigenstrain(int x, int y, int z);
		// solver loop for MFType_Implicit
		void exec_pre_im_steinbach();
		std::string exec_loop_im_steinbach();
		void exec_pre_im_khachaturyan();
		std::string exec_loop_im_khachaturyan();

		void init();

		void exec_pre();

		std::string exec_loop();

		void deinit();

		void write_field(std::ofstream& fout);

	}
}