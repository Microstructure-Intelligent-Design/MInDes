#pragma once
#include "../../base/VectorMatrix.h"
#include "../../base/Mesh_0.h"
#include "../../Modules_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../model_modules/GrainsOrientations.h"
#include "../../model_modules/PhiProperties.h"
namespace pf {
	namespace stiffness_eigenstrain {
		enum StiffnessEigenStrainType { ESType_None, ESType_PhaseDependent, ESType_GrainDependent };
		inline size_t phi_number = 0;
		inline size_t phi_property_number = 0;
		inline std::vector<Matrix6x6> phi_property_stiffness;
		inline std::vector<Matrix6x6> phi_index_stiffness;
		inline std::vector<vStrain> phi_property_eigen_strain;
		inline std::vector<vStrain> phi_index_eigen_strain;
		// - function list, which is called in elastic solver, and more eigenstain model can be defined and push in this list
		inline std::vector<void(*)(long long, long long, long long, vStrain&)> eigenstrain_list;
		// - 
		inline Matrix6x6 get_phase_stiffness(size_t phi_property) {
			return phi_property_stiffness[phi_property];
		}
		inline Matrix6x6 get_phi_stiffness(size_t phi_index) {
			return phi_index_stiffness[phi_index];
		}
		inline vStrain get_phase_eigen_strain(size_t phi_property) {
			return phi_property_eigen_strain[phi_property];
		}
		inline vStrain get_phi_eigen_strain(size_t phi_index) {
			return phi_index_eigen_strain[phi_index];
		}
		void do_phi_index_orientation(size_t phi_index, size_t phi_property);
		// stiffness and eigenstain
		inline void(*stiffness)(long long, long long, long long, Matrix6x6&);
		inline void(*eigen_strain)(long long, long long, long long, vStrain&);
		// 
		void init(size_t _phi_number, size_t _phi_property_number);
	}
}