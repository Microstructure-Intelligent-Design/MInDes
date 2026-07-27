#pragma once
#include "../../base/Mesh_0.h"
#include "../../base/RotationMatrix.h"
#include "../../Modules_Params.h"
#include "../PhiProperties.h"
#include "../ConRegions.h"
#include "../GrainsOrientations.h"
namespace pf {
	namespace ddc_calphad_ai_model {
		// Statement
		enum Int_Gradient { Steinbach_1996, Steinbach_1999, Steinbach_G2009 };
		enum Int_Potential { Nestler_Well, Nestler_Obstacle, Steinbach_P2009 };
		enum DifferenceMethod { SEVEN_POINT, NINETEEN_POINT };
		enum InterfaceFlag { IF_BULK, IF_NEAR_INTERFACE, IF_INTERFACE };
		struct ThermoCalcScanRange {
			REAL begin = 0;
			REAL end = 0;
			REAL step = 1;
		};
		struct FIELD_PhiTemp
		{
			// - phi
			std::vector<size_t> active_index;
			size_t active_number = 0;
			std::vector<int> intflag;
			std::vector<REAL> old_phi;
			std::vector<REAL> new_phi;
			std::vector<REAL> lap_phi;
			std::vector<Vector3> grad_phi;
			// - temp
			REAL old_temp = 0;
			REAL new_temp = 0;
			REAL mob_temp = 0;
			// - 
			void init_phi(size_t phi_number, size_t acc_number) {
				active_index.resize(acc_number, 0);
				intflag.resize(phi_number, InterfaceFlag::IF_BULK);
				old_phi.resize(phi_number, 0);
				new_phi.resize(phi_number, 0);
				lap_phi.resize(phi_number, 0);
				grad_phi.resize(phi_number, Vector3(0, 0, 0));
			}
			// 矩阵加法 (+)
			FIELD_PhiTemp operator+(const FIELD_PhiTemp& other) const {
				FIELD_PhiTemp result = *this;
				result += other;
				return result;
			}
			// 矩阵减法 (-)
			FIELD_PhiTemp operator-(const FIELD_PhiTemp& other) const {
				FIELD_PhiTemp result = *this;
				result -= other;
				return result;
			}
			// 矩阵除法 (*) - 逐元素相乘
			FIELD_PhiTemp operator*(const REAL& other) const {
				FIELD_PhiTemp result = *this;
				result *= other;
				return result;
			}
			// 矩阵除法 (/) - 逐元素相除
			FIELD_PhiTemp operator/(const REAL& other) const {
				FIELD_PhiTemp result = *this;
				result /= other;
				return result;
			}
			// 复合赋值运算符 (+=)，可以提升连续运算时的性能
			FIELD_PhiTemp& operator+=(const FIELD_PhiTemp& other) {
				size_t phi_number = old_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					old_phi[i] += other.old_phi[i];
					new_phi[i] += other.new_phi[i];
				}
				old_temp += other.old_temp;
				new_temp += other.new_temp;
				mob_temp += other.mob_temp;
				return *this;
			}
			// 复合赋值运算符 (-=)，可以提升连续运算时的性能
			FIELD_PhiTemp& operator-=(const FIELD_PhiTemp& other) {
				size_t phi_number = old_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					old_phi[i] -= other.old_phi[i];
					new_phi[i] -= other.new_phi[i];
				}
				old_temp -= other.old_temp;
				new_temp -= other.new_temp;
				mob_temp -= other.mob_temp;
				return *this;
			}
			// 复合赋值运算符 (*=)，可以提升连续运算时的性能
			FIELD_PhiTemp& operator*=(const REAL& other) {
				size_t phi_number = old_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					old_phi[i] *= other;
					new_phi[i] *= other;
				}
				old_temp *= other;
				new_temp *= other;
				mob_temp *= other;
				return *this;
			}
			// 复合赋值运算符 (/=)，可以提升连续运算时的性能
			FIELD_PhiTemp& operator/=(const REAL& other) {
				size_t phi_number = old_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					old_phi[i] /= other;
					new_phi[i] /= other;
				}
				old_temp /= other;
				new_temp /= other;
				mob_temp /= other;
				return *this;
			}
		};
		struct FIELD_Con
		{
			// - con
			std::vector<REAL> old_region;
			std::vector<REAL> new_region;
			std::vector<REAL> old_con;
			std::vector<REAL> new_con;
			std::vector<REAL> miu_con;
			std::vector<REAL> mob_con;
			// - 
			void init_con(size_t con_number, size_t region_number) {
				old_con.resize(con_number, 0);
				new_con.resize(con_number, 0);
				miu_con.resize(con_number, 0);
				mob_con.resize(con_number, 0);
				old_region.resize(region_number, 0);
				new_region.resize(region_number, 0);
			}
			// 矩阵加法 (+)
			FIELD_Con operator+(const FIELD_Con& other) const {
				FIELD_Con result = *this;
				result += other;
				return result;
			}
			// 矩阵减法 (-)
			FIELD_Con operator-(const FIELD_Con& other) const {
				FIELD_Con result = *this;
				result -= other;
				return result;
			}
			// 矩阵除法 (*) - 逐元素相乘
			FIELD_Con operator*(const REAL& other) const {
				FIELD_Con result = *this;
				result *= other;
				return result;
			}
			// 矩阵除法 (/) - 逐元素相除
			FIELD_Con operator/(const REAL& other) const {
				FIELD_Con result = *this;
				result /= other;
				return result;
			}
			// 复合赋值运算符 (+=)，可以提升连续运算时的性能
			FIELD_Con& operator+=(const FIELD_Con& other) {
				size_t con_number = old_con.size(), region_number = old_region.size();
				for (size_t i = 0; i < con_number; ++i) {
					old_con[i] += other.old_con[i];
					new_con[i] += other.new_con[i];
					miu_con[i] += other.miu_con[i];
					mob_con[i] += other.mob_con[i];
				}
				for (size_t i = 0; i < region_number; ++i) {
					old_region[i] += other.old_region[i];
					new_region[i] += other.new_region[i];
				}
				return *this;
			}
			// 复合赋值运算符 (-=)，可以提升连续运算时的性能
			FIELD_Con& operator-=(const FIELD_Con& other) {
				size_t con_number = old_con.size(), region_number = old_region.size();
				for (size_t i = 0; i < con_number; ++i) {
					old_con[i] -= other.old_con[i];
					new_con[i] -= other.new_con[i];
					miu_con[i] -= other.miu_con[i];
					mob_con[i] -= other.mob_con[i];
				}
				for (size_t i = 0; i < region_number; ++i) {
					old_region[i] -= other.old_region[i];
					new_region[i] -= other.new_region[i];
				}
				return *this;
			}
			// 复合赋值运算符 (*=)，可以提升连续运算时的性能
			FIELD_Con& operator*=(const REAL& other) {
				size_t con_number = old_con.size(), region_number = old_region.size();
				for (size_t i = 0; i < con_number; ++i) {
					old_con[i] *= other;
					new_con[i] *= other;
					miu_con[i] *= other;
					mob_con[i] *= other;
				}
				for (size_t i = 0; i < region_number; ++i) {
					old_region[i] *= other;
					new_region[i] *= other;
				}
				return *this;
			}
			// 复合赋值运算符 (/=)，可以提升连续运算时的性能
			FIELD_Con& operator/=(const REAL& other) {
				size_t con_number = old_con.size(), region_number = old_region.size();
				for (size_t i = 0; i < con_number; ++i) {
					old_con[i] /= other;
					new_con[i] /= other;
					miu_con[i] /= other;
					mob_con[i] /= other;
				}
				for (size_t i = 0; i < region_number; ++i) {
					old_region[i] /= other;
					new_region[i] /= other;
				}
				return *this;
			}
		};
		struct FIELD_PhaseCon
		{
			// - phase phi
			std::vector<REAL> phase_old_phi;
			std::vector<REAL> phase_phi;
			// - phase con
			std::vector<std::vector<REAL>> phase_con;
			// - phase miu
			std::vector<std::vector<REAL>> phase_miu;
			// - 
			void init_con(size_t phase_number, size_t con_number) {
				phase_old_phi.resize(phase_number, 0);
				phase_phi.resize(phase_number, 0);
				phase_con.resize(phase_number);
				phase_miu.resize(phase_number);
				for (size_t index = 0; index < phase_number; index++) {
					phase_con[index].resize(con_number, 0);
					phase_miu[index].resize(con_number, 0);
				}
			}
			// 矩阵加法 (+)
			FIELD_PhaseCon operator+(const FIELD_PhaseCon& other) const {
				FIELD_PhaseCon result = *this;
				result += other;
				return result;
			}
			// 矩阵减法 (-)
			FIELD_PhaseCon operator-(const FIELD_PhaseCon& other) const {
				FIELD_PhaseCon result = *this;
				result -= other;
				return result;
			}
			// 矩阵除法 (*) - 逐元素相乘
			FIELD_PhaseCon operator*(const REAL& other) const {
				FIELD_PhaseCon result = *this;
				result *= other;
				return result;
			}
			// 矩阵除法 (/) - 逐元素相除
			FIELD_PhaseCon operator/(const REAL& other) const {
				FIELD_PhaseCon result = *this;
				result /= other;
				return result;
			}
			// 复合赋值运算符 (+=)，可以提升连续运算时的性能
			FIELD_PhaseCon& operator+=(const FIELD_PhaseCon& other) {
				size_t con_number = phase_con[0].size(), phi_number = phase_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					phase_phi[i] += other.phase_phi[i];
					phase_old_phi[i] += other.phase_old_phi[i];
					for (size_t j = 0; j < con_number; ++j) {
						phase_con[i][j] += other.phase_con[i][j];
						phase_miu[i][j] += other.phase_miu[i][j];
					}
				}
				return *this;
			}
			// 复合赋值运算符 (-=)，可以提升连续运算时的性能
			FIELD_PhaseCon& operator-=(const FIELD_PhaseCon& other) {
				size_t con_number = phase_con[0].size(), phi_number = phase_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					phase_phi[i] -= other.phase_phi[i];
					phase_old_phi[i] -= other.phase_old_phi[i];
					for (size_t j = 0; j < con_number; ++j) {
						phase_con[i][j] -= other.phase_con[i][j];
						phase_miu[i][j] -= other.phase_miu[i][j];
					}
				}
				return *this;
			}
			// 复合赋值运算符 (*=)，可以提升连续运算时的性能
			FIELD_PhaseCon& operator*=(const REAL& other) {
				size_t con_number = phase_con[0].size(), phi_number = phase_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					phase_phi[i] *= other;
					phase_old_phi[i] *= other;
					for (size_t j = 0; j < con_number; ++j) {
						phase_con[i][j] *= other;
						phase_miu[i][j] *= other;
					}
				}
				return *this;
			}
			// 复合赋值运算符 (/=)，可以提升连续运算时的性能
			FIELD_PhaseCon& operator/=(const REAL& other) {
				size_t con_number = phase_con[0].size(), phi_number = phase_phi.size();
				for (size_t i = 0; i < phi_number; ++i) {
					phase_phi[i] /= other;
					phase_old_phi[i] /= other;
					for (size_t j = 0; j < con_number; ++j) {
						phase_con[i][j] /= other;
						phase_miu[i][j] /= other;
					}
				}
				return *this;
			}
		};
		namespace parameters {
			inline REAL Phi_Cut_Off = 0.001;
			inline REAL Phi_Cut_Off_R = 0.999;
			inline REAL PhiCon_Cut_Off = 0.1;
			inline REAL PhiCon_Cut_Off_R = 0.9;
			// ===================================================================================================
			// - buff field for Phi Con Temp
			inline Mesh_Boundry<FIELD_PhiTemp>  PhiTemp_field;
			inline Mesh_Boundry<FIELD_Con>  Con_field;
			inline Mesh_Boundry<FIELD_PhaseCon>  PhaseCon_field;
			inline bool is_phi_normalized = false;
			inline bool is_con_normalized = false;
			inline bool is_temp_normalized = false;
			inline bool is_phasecon_on = false;
			inline DifferenceMethod diff_method = DifferenceMethod::SEVEN_POINT;
			// - statistic
			const std::pair<std::string, std::string> statistic_phase_volume = { "vol_phase", "volume fraction of each phase" };
			const std::pair<std::string, std::string> statistic_region_volume = { "vol_region", "volume fraction of each region" };
			const std::pair<std::string, std::string> statistic_average_con = { "ave_con", "average concentration of each component" };
			const std::pair<std::string, std::string> statistic_total_con = { "tot_con", "total concentration of each component" };
			inline bool is_phase_volume_statistic = false;
			inline bool is_region_volume_statistic = false;
			inline bool is_average_con_statistic = false;
			inline bool is_total_con_statistic = false;
			// ===================================================================================================
			// - parameters phi
			// - pairwise accelerate
			inline size_t PAIRWISE_ACC_STOP = 10;
			inline size_t PHI_ACC_NUMBER = 10;
			inline size_t MAX_ACTIVE_PHI_NUMBER = 0;
			// - driving force functions
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z, size_t phi_index)> delt_Fbulk_delt_phi;
			// - source functions
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z, size_t alpha_index, size_t beta_index)> source_alpha_beta;
			// - anisotropy
			inline std::vector<Matrix3x3> grain_rotation_matrix;
			// - mobility
			inline Matrix2D<REAL> Lij; // <- phi index i j
			inline Matrix2D<REAL> Qij; // <- phi index i j
			const REAL R = REAL(8.314);
			// - interface mobility anisotropy
			enum Int_Mobility_Anisotropic { IMA_ISO, IMA_CUBIC, IMA_HEX_BOETTGER, IMA_HEX_SUN, IMA_HEX_YANG, IMA_DENDRITE_YANG };
			inline Int_Mobility_Anisotropic intMobAniso_model = Int_Mobility_Anisotropic::IMA_ISO;
			inline REAL intMobAniso_param1;
			inline REAL intMobAniso_param2;
			inline REAL intMobAniso_param3;
			inline REAL intMobAniso_param4;
			// - interface energy models
			inline Int_Gradient interface_gradient = Int_Gradient::Steinbach_G2009;
			inline Int_Potential interface_potential = Int_Potential::Steinbach_P2009;
			// interface energy
			inline REAL interface_width = REAL(4.0);
			inline Matrix2D<REAL> xi_ab; // <- phi property 
			inline Matrix3D<REAL> xi_abc; // <- phi property 
			// interface energy anisotropy 
			enum Int_Energy_Anisotropic { IEA_ISO, IEA_CUBIC, IEA_HEX_BOETTGER, IEA_HEX_SUN, IEA_HEX_YANG, IEA_DENDRITE_YANG	};
			inline Int_Energy_Anisotropic intEnAniso_model = Int_Energy_Anisotropic::IEA_ISO;
			inline REAL intEnAniso_param1;
			inline REAL intEnAniso_param2;
			inline REAL intEnAniso_param3;
			inline REAL intEnAniso_param4;
			// const bulk energy density
			inline std::vector<REAL> f_bulk_0; // <- phi property 
			// ===================================================================================================
			// - parameters con 
			inline std::pair<REAL, size_t>(*local_concentration_redistribution)(size_t x, size_t y, size_t z);  // { MAX_VARIATION , MAX_ITERATION_STEP }
			// - moving region method 
			inline void (*init_con_in_moving_region)(size_t x, size_t y, size_t z, size_t region_index);
			inline void (*deinit_con_in_moving_region)(size_t x, size_t y, size_t z, size_t region_index);
			// - driving force functions 
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z, size_t region_index, size_t con_index)> delt_Fbulk_delt_con;
			// - mobility functions 
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z, size_t region_index, size_t con_index)> mobility;
			// - 
			inline void (*cal_mob_miu_grad_lap)(size_t x, size_t y, size_t z, size_t region_index,
				Vector3& grad_region, std::vector<Vector3>& grad_miu, std::vector<Vector3>& grad_mob, std::vector<REAL>& lap_miu);
			// - mobility parameters 
			inline Matrix2D<REAL> Mii;             // (phi property, con_index) -> Mii
			inline Matrix2D<REAL> Qii;             // (phi property, con_index) -> Qii
			// - anisotropic diffusion 
			inline Matrix2D<REAL> Mii_surf;             // (region index, con_index) -> Mii_surf 
			inline Matrix3D<REAL> Mii_grain;             // (phi property_1, phi property_2, con_index) -> Mii_grain 
			// - reaction functions 
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z, size_t region_index, size_t con_index)> bulk_reaction;
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z, size_t region_index, size_t con_index)> surface_reaction;
			// - con laplace term 
			inline std::vector<REAL> Ki;                // [con_index] -> Ki 
			// ===================================================================================================
			// - parameters temperature
			// - mobility functions 
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z)> mobility_temp;
			// - source functions 
			inline std::vector<REAL(*)(size_t x, size_t y, size_t z)> source_temp;
			// - 
			inline void (*cal_mob_temp_grad_lap)(size_t x, size_t y, size_t z, Vector3& grad_temp, Vector3& grad_mob, REAL& lap_temp);
			// - mobility parameters 
			inline std::vector<REAL> Mtemp;             // [phi property] -> Mtemp
			// - source term
			inline std::vector<REAL> Ktemp;             // [phi property] -> Ktemp     Ktemp * dphi / dt
			// ===================================================================================================
			// - parameters bulk energy
			inline std::vector<bool> is_energy_minimization;      // [region index] -> true/false
			inline size_t max_phasecon_variation_step = 100;
			inline REAL L_phasecon = 1.0;
			inline REAL phasecon_epsilon = 1e-4;
			inline std::pair<REAL, size_t> MAX_MINIMIZATION_RESULTS = {0, 0};
			// - chemical energy
			inline REAL (*fchem)(std::vector<REAL> con, REAL temperature, size_t phi_property);
			inline REAL (*miu)(std::vector<REAL> con, REAL temperature, size_t phi_property, size_t con_index);
			// - polynomial
			inline std::vector<Matrix2D<int>> con_orders;      // [phi_property](term index, con index) -> con order
			inline std::vector<Matrix2D<REAL>> params;         // [phi_property](term index, temp term index) -> param
			inline std::vector<Matrix2D<int>> temp_orders;     // [phi_property](term index, temp term index) -> temp order
			inline std::vector<size_t> terms_number;           // [phi_property] -> term number
			inline std::vector<std::vector<size_t>> temp_terms_number;       // [phi_property][term index] -> temp term number
			// - AI

			// ===================================================================================================
			// - output
			inline std::vector<std::pair<size_t, size_t>> is_write_phase_con;  // ( phase_property , con_index )
			inline std::vector<size_t> is_write_con;  // ( con_index )
			inline std::vector<std::pair<size_t, size_t>> is_write_phase_miu;  // ( phase_property , con_index )
			inline std::vector<size_t> is_write_miu;                           // ( con_index )
			inline std::vector<size_t> is_write_fchem;                         // ( phase_property )
			// - write thermodynamic energy
			inline std::vector<size_t> thermo_calc_energy_phases;
			inline std::vector<ThermoCalcScanRange> thermo_calc_energy_con_ranges;
			inline ThermoCalcScanRange thermo_calc_energy_temp_range;
			inline std::string thermo_calc_energy_csv_name = "thermo_calc_energy_data.csv";
			// - write thermodynamic calculation data
			inline bool is_write_thermo_calc_csv = false;
			inline size_t thermo_calc_region = 0;
			inline std::vector<bool> thermo_calc_phi_property;
			inline std::vector<std::pair<bool, REAL>> thermo_calc_fix_con;
			inline std::vector<size_t> thermo_calc_phases;
			inline std::vector<ThermoCalcScanRange> thermo_calc_con_ranges;
			inline std::vector<ThermoCalcScanRange> thermo_calc_phi_ranges;
			inline ThermoCalcScanRange thermo_calc_temp_range;
			inline std::string thermo_calc_equi_csv_name = "thermo_calc_equi_data.csv";
		}
	}
}
