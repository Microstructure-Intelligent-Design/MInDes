#pragma once
#include "../../base/Mesh_0.h"
#include "../../base/RotationMatrix.h"
#include "../../Modules_Params.h"
namespace pf {
	namespace data_driven_complex_model {
		struct FIELD_VAR
		{
			std::vector<size_t> active_index;
			std::vector<bool> interphase;
			std::vector<REAL> old_phi;
			std::vector<REAL> new_phi;
			std::vector<REAL> lap_phi;
			std::vector<Vector3> grad_phi;
			std::vector<REAL> old_con;
			std::vector<REAL> new_con;
			REAL old_temp;
			REAL new_temp;
			void init(size_t phi_number, size_t con_number, size_t acc_number) {
				active_index.resize(acc_number, 0);
				interphase.resize(phi_number, false);
				old_phi.resize(phi_number, 0);
				new_phi.resize(phi_number, 0);
				lap_phi.resize(phi_number, 0);
				grad_phi.resize(phi_number, Vector3(0, 0, 0));
				old_con.resize(con_number, 0);
				new_con.resize(con_number, 0);
			}
			FIELD_VAR() {
				old_temp = 0;
				new_temp = 0;
			}
			FIELD_VAR& operator=(const FIELD_VAR& n) {
				active_index = n.active_index;
				interphase = n.interphase;
				old_phi = n.old_phi;
				new_phi = n.new_phi;
				lap_phi = n.lap_phi;
				grad_phi = n.grad_phi;
				old_con = n.old_con;
				new_con = n.new_con;
				old_temp = n.old_temp;
				new_temp = n.new_temp;
			}
		};
		struct FIELD_VAR_BC
		{
			std::vector<REAL> miu;
			std::vector<REAL> mob;
			void init(size_t con_numberr) {
				miu.resize(con_numberr, 0);
				mob.resize(con_numberr, 0);
			}
			FIELD_VAR_BC() { }
			FIELD_VAR_BC& operator=(const FIELD_VAR_BC& n) {
				miu = n.miu;
				mob = n.mob;
			}
		};
		class GrainsOrientations
		{
		public:
			GrainsOrientations() = default;
			~GrainsOrientations() {
				orientations.clear();
			}
			GrainsOrientations& operator=(const GrainsOrientations& n) {
				rotation_gauge = n.rotation_gauge;
				orientations = n.orientations;
				return *this;
			}
			void init(RotationGauge _rotation_gauge, size_t phi_number) {
				rotation_gauge = _rotation_gauge;
				orientations.resize(phi_number);
				for (size_t index = 0; index < phi_number; index++) {
					orientations[index][0] = REAL(0.0);
					orientations[index][1] = REAL(0.0);
					orientations[index][2] = REAL(0.0);
				}
			}
			void set_phi_orientation(size_t phi_index, Vector3 radian) {
				orientations[phi_index] = radian;
			}
			Vector3 get_phi_orientation(size_t phi_index) {
				return orientations[phi_index];
			}
			Matrix3x3 get_phi_rotationMatrix(size_t phi_index) {
				return RotationMatrix::rotationMatrix(orientations[phi_index], rotation_gauge);
			}
			RotationGauge rotation_gauge;
			std::vector<Vector3> orientations;
		};
		namespace parameters {
			// phases name/property in this simulation
			inline std::vector<std::string> PHASES;
			inline size_t phase_property(std::string phi_name) {
				for (size_t index = 0; index < PHASES.size(); index++)
					if (phi_name.compare(PHASES[index]) == 0)
						return index;
				std::cout << "ERROR, phase name: " << phi_name << " has not been defined !";
				exit(0);
			}
			// components name/index in this simulation
			inline std::vector<std::string> COMPONENTS;
			inline size_t comp_index(std::string con_name) {
				for (size_t index = 0; index < COMPONENTS.size(); index++)
					if (con_name.compare(COMPONENTS[index]) == 0)
						return index;
				std::cout << "ERROR, component name: " << con_name << " has not been defined !";
				exit(0);
			}
			// - property/index/region parameters
			inline size_t phi_property_number = 0;
			inline std::vector<size_t> phi_property; // [index] -> property
			inline size_t con_region_number = 0;
			inline std::vector<size_t> con_region; // [index] -> region
			inline std::vector<std::vector<size_t>> rg_con_index; // [region] -> { con_index1 , con_index2 , ... }
			inline std::vector<std::vector<size_t>> rg_phi_property; // [region] -> { phi_property1 , phi_property2 , ... }
			inline std::vector<std::vector<size_t>> pp_con_index; // [property] -> { con_index1 , con_index2 , ... }
			inline std::vector<std::vector<size_t>> rg_phi_index; // [region] -> { phi_index1 , phi_index2 , ... }
			// - pairwise accelerate
			inline size_t PAIRWISE_ACC_STOP = 0;
			inline size_t PHI_ACC_NUMBER = 0;
			// - field
			inline Mesh_Boundry<FIELD_VAR> field_variable;
			inline Mesh_Boundry<FIELD_VAR_BC> field_variable_bc;
			// ===================================================================================================
			// - parameters
			// Statement
			enum class Int_Gradient { Steinbach_1996, Steinbach_1999, Steinbach_G2009 };
			enum class Int_Potential { Nestler_Well, Nestler_Obstacle, Steinbach_P2009 };
			enum class DifferenceMethod { FIVE_POINT, NINE_POINT };
			inline bool is_phi_normalized = false;
			inline bool is_con_normalized = false;
			inline bool is_temp_normalized = false;
			inline DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT;
			// - driving force functions
			inline std::vector<REAL(*)(long long x, long long y, long long z, size_t phi_index)> delt_Fbulk_delt_phi;
			// - source functions
			inline std::vector<REAL(*)(long long x, long long y, long long z, size_t alpha_index, size_t beta_index)> source_alpha_beta;
			// - anisotropy
			// - mobility
			inline std::vector<std::vector<REAL>> Lij; // <- phi index
			inline std::vector<std::vector<REAL>> Qij; // <- phi index
			const REAL R = REAL(8.314);
			// - interface mobility anisotropy
			inline GrainsOrientations grains_orientation;
			enum Int_Mobility_Anisotropic { IMA_ISO, IMA_CUBIC, IMA_HEX_BOETTGER, IMA_HEX_SUN, IMA_HEX_YANG };
			inline std::vector<std::vector<int>>  intMob_anisotropic_model_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intMobAniso1_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intMobAniso2_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intMobAniso3_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intMobAniso4_matrix; // <- phi property
			// - interface energy models
			inline Int_Gradient interface_gradient = Int_Gradient::Steinbach_G2009;
			inline Int_Potential interface_potential = Int_Potential::Steinbach_P2009;
			// interface energy
			inline REAL interface_width = REAL(4.0);
			inline std::vector<std::vector<REAL>> xi_ab; // <- phi property
			inline std::vector<std::vector<std::vector<REAL>>> xi_abc; // <- phi property
			// interface energy anisotropy
			enum Int_Energy_Anisotropic { IEA_ISO, IEA_CUBIC, IEA_HEX_BOETTGER, IEA_HEX_SUN, IEA_HEX_YANG };
			inline std::vector<std::vector<int>>  intEn_anisotropic_model_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intEnAniso1_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intEnAniso2_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intEnAniso3_matrix; // <- phi property
			inline std::vector<std::vector<REAL>> intEnAniso4_matrix; // <- phi property
		}
	}
}