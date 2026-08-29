#include "DDCCPAI_Phi_Functions.h"
#include <cmath>
namespace pf {
	namespace ddc_calphad_ai_model {
		namespace phase_field_functions {
			void interphase_gradient_lapace_calculation_7P(size_t x, size_t y, size_t z, size_t phi_index) {
				if (!main_field::phase_field.IS_COMP_POINT(x, y, z))
					return;
				FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
				field_var.lap_phi[phi_index] =
					(main_field::phase_field(x + 1, y, z)[phi_index] + main_field::phase_field(x - 1, y, z)[phi_index]
						+ main_field::phase_field(x, y + 1, z)[phi_index] + main_field::phase_field(x, y - 1, z)[phi_index]
						+ main_field::phase_field(x, y, z + 1)[phi_index] + main_field::phase_field(x, y, z - 1)[phi_index]
						- 6 * main_field::phase_field(x, y, z)[phi_index]) / mesh_parameters::delt_r / mesh_parameters::delt_r;
				field_var.grad_phi[phi_index][0] = (main_field::phase_field(x + 1, y, z)[phi_index]
					- main_field::phase_field(x - 1, y, z)[phi_index]) / 2 / mesh_parameters::delt_r;
				field_var.grad_phi[phi_index][1] = (main_field::phase_field(x, y + 1, z)[phi_index]
					- main_field::phase_field(x, y - 1, z)[phi_index]) / 2 / mesh_parameters::delt_r;
				field_var.grad_phi[phi_index][2] = (main_field::phase_field(x, y, z + 1)[phi_index]
					- main_field::phase_field(x, y, z - 1)[phi_index]) / 2 / mesh_parameters::delt_r;
			}
			void interphase_gradient_lapace_calculation_19P(size_t x, size_t y, size_t z, size_t phi_index) {
				if (!main_field::phase_field.IS_COMP_POINT(x, y, z))
					return;
				FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
				field_var.lap_phi[phi_index] =
					((main_field::phase_field(x + 1, y, z)[phi_index] + main_field::phase_field(x - 1, y, z)[phi_index]
						+ main_field::phase_field(x, y + 1, z)[phi_index] + main_field::phase_field(x, y - 1, z)[phi_index]
						+ main_field::phase_field(x, y, z + 1)[phi_index] + main_field::phase_field(x, y, z - 1)[phi_index]) * 4
						+ main_field::phase_field(x - 1, y - 1, z)[phi_index] + main_field::phase_field(x - 1, y + 1, z)[phi_index]
						+ main_field::phase_field(x + 1, y - 1, z)[phi_index] + main_field::phase_field(x + 1, y + 1, z)[phi_index]
						+ main_field::phase_field(x - 1, y, z - 1)[phi_index] + main_field::phase_field(x - 1, y, z + 1)[phi_index]
						+ main_field::phase_field(x + 1, y, z - 1)[phi_index] + main_field::phase_field(x + 1, y, z + 1)[phi_index]
						+ main_field::phase_field(x, y - 1, z - 1)[phi_index] + main_field::phase_field(x, y - 1, z + 1)[phi_index]
						+ main_field::phase_field(x, y + 1, z - 1)[phi_index] + main_field::phase_field(x, y + 1, z + 1)[phi_index]
						- 36 * main_field::phase_field(x, y, z)[phi_index]) / 6 / mesh_parameters::delt_r / mesh_parameters::delt_r;
				field_var.grad_phi[phi_index][0] = (main_field::phase_field(x + 1, y, z)[phi_index]
					- main_field::phase_field(x - 1, y, z)[phi_index]) / 2 / mesh_parameters::delt_r;
				field_var.grad_phi[phi_index][1] = (main_field::phase_field(x, y + 1, z)[phi_index]
					- main_field::phase_field(x, y - 1, z)[phi_index]) / 2 / mesh_parameters::delt_r;
				field_var.grad_phi[phi_index][2] = (main_field::phase_field(x, y, z + 1)[phi_index]
					- main_field::phase_field(x, y, z - 1)[phi_index]) / 2 / mesh_parameters::delt_r;
			}
			InterfaceFlag currentFlag_7P(size_t x, size_t y, size_t z, size_t phi_index) {
				if (!main_field::phase_field.IS_COMP_POINT(x, y, z))
					return InterfaceFlag::IF_BULK;
				Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
				if (phi[phi_index] >= parameters::Phi_Cut_Off && phi[phi_index] <= parameters::Phi_Cut_Off_R)
					return InterfaceFlag::IF_INTERFACE;
				if (phi[phi_index] < parameters::Phi_Cut_Off) {
					if (main_field::phase_field(x + 1, y, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z + 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z - 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
				}
				else if (phi[phi_index] > parameters::Phi_Cut_Off_R) {
					if (main_field::phase_field(x + 1, y, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z + 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z - 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
				}
				return InterfaceFlag::IF_BULK;
			}
			InterfaceFlag currentFlag_19P(size_t x, size_t y, size_t z, size_t phi_index) {
				if (!main_field::phase_field.IS_COMP_POINT(x, y, z))
					return InterfaceFlag::IF_BULK;
				Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
				if (phi[phi_index] >= parameters::Phi_Cut_Off && phi[phi_index] <= parameters::Phi_Cut_Off_R)
					return InterfaceFlag::IF_INTERFACE;
				if (phi[phi_index] < parameters::Phi_Cut_Off) {
					if (main_field::phase_field(x + 1, y, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z + 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z - 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y + 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y - 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y, z + 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y, z - 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y + 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y - 1, z)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z + 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z - 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z + 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z - 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z + 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z - 1)[phi_index] >= parameters::Phi_Cut_Off)
						return InterfaceFlag::IF_NEAR_INTERFACE;
				}
				else if (phi[phi_index] > parameters::Phi_Cut_Off_R) {
					if (main_field::phase_field(x + 1, y, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z + 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y, z - 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y + 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y - 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y, z + 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x + 1, y, z - 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y + 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y - 1, z)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z + 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x - 1, y, z - 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z + 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y + 1, z - 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z + 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
					else if (main_field::phase_field(x, y - 1, z - 1)[phi_index] <= parameters::Phi_Cut_Off_R)
						return InterfaceFlag::IF_NEAR_INTERFACE;
				}
				return InterfaceFlag::IF_BULK;
			}
			void upgradeFlag_7P(size_t x, size_t y, size_t z, size_t phi_index) {
				REAL phi_fraction = 0;
				if (main_field::phase_field.IS_COMP_POINT(x - 1, y, z)) {
					phi_fraction = main_field::phase_field(x - 1, y, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x - 1, y, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x + 1, y, z)) {
					phi_fraction = main_field::phase_field(x + 1, y, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x + 1, y, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y - 1, z)) {
					phi_fraction = main_field::phase_field(x, y - 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y - 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y + 1, z)) {
					phi_fraction = main_field::phase_field(x, y + 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y + 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y, z - 1)) {
					phi_fraction = main_field::phase_field(x, y, z - 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y, z - 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y, z + 1)) {
					phi_fraction = main_field::phase_field(x, y, z + 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y, z + 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
			}
			void upgradeFlag_19P(size_t x, size_t y, size_t z, size_t phi_index) {
				REAL phi_fraction = 0;
				if (main_field::phase_field.IS_COMP_POINT(x - 1, y, z)) {
					phi_fraction = main_field::phase_field(x - 1, y, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x - 1, y, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x + 1, y, z)) {
					phi_fraction = main_field::phase_field(x + 1, y, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x + 1, y, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y - 1, z)) {
					phi_fraction = main_field::phase_field(x, y - 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y - 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y + 1, z)) {
					phi_fraction = main_field::phase_field(x, y + 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y + 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y, z - 1)) {
					phi_fraction = main_field::phase_field(x, y, z - 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y, z - 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y, z + 1)) {
					phi_fraction = main_field::phase_field(x, y, z + 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y, z + 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x - 1, y - 1, z)) {
					phi_fraction = main_field::phase_field(x - 1, y - 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x - 1, y - 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x - 1, y + 1, z)) {
					phi_fraction = main_field::phase_field(x - 1, y + 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x - 1, y + 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x - 1, y, z - 1)) {
					phi_fraction = main_field::phase_field(x - 1, y, z - 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x - 1, y, z - 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x - 1, y, z + 1)) {
					phi_fraction = main_field::phase_field(x - 1, y, z + 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x - 1, y, z + 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x + 1, y - 1, z)) {
					phi_fraction = main_field::phase_field(x + 1, y - 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x + 1, y - 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x + 1, y + 1, z)) {
					phi_fraction = main_field::phase_field(x + 1, y + 1, z)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x + 1, y + 1, z).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x + 1, y, z - 1)) {
					phi_fraction = main_field::phase_field(x + 1, y, z - 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x + 1, y, z - 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x + 1, y, z + 1)) {
					phi_fraction = main_field::phase_field(x + 1, y, z + 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x + 1, y, z + 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y - 1, z - 1)) {
					phi_fraction = main_field::phase_field(x, y - 1, z - 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y - 1, z - 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y - 1, z + 1)) {
					phi_fraction = main_field::phase_field(x, y - 1, z + 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y - 1, z + 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y + 1, z - 1)) {
					phi_fraction = main_field::phase_field(x, y + 1, z - 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y + 1, z - 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
				if (main_field::phase_field.IS_COMP_POINT(x, y + 1, z + 1)) {
					phi_fraction = main_field::phase_field(x, y + 1, z + 1)[phi_index];
					if (phi_fraction < parameters::Phi_Cut_Off || phi_fraction > parameters::Phi_Cut_Off_R)
						parameters::PhiTemp_field(x, y + 1, z + 1).intflag[phi_index] = InterfaceFlag::IF_NEAR_INTERFACE;
				}
			}
			static void normalize_new_phi(FIELD_PhiTemp& var) {
				REAL all_phi = 0;
				for (size_t index = 0; index < main_field::phi_number; index++)
					all_phi += var.new_phi[index];
				if (all_phi < SYS_EPSILON) {
					return;
				}
				else {
					for (size_t index = 0; index < main_field::phi_number; index++)
						var.new_phi[index] /= all_phi;
				}
			}
			// - 
			Vector3 normals(REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				Vector3 normal;
				normal[0] = alpha_grad[0] * beta_phi - alpha_phi * beta_grad[0];
				normal[1] = alpha_grad[1] * beta_phi - alpha_phi * beta_grad[1];
				normal[2] = alpha_grad[2] * beta_phi - alpha_phi * beta_grad[2];
				REAL length = sqrt(normal * normal);
				if (length < SYS_EPSILON) {
					normal[0] = 0.0;
					normal[1] = 0.0;
					normal[2] = 0.0;
				}
				else {
					normal[0] = normal[0] / length;
					normal[1] = normal[1] / length;
					normal[2] = normal[2] / length;
				};
				return normal;
			};
			// - interface mobility
			REAL Lij(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				return parameters::Lij(alpha_index, beta_index);
			};
			REAL Lij_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1 *
					(1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return Lij;
			};
			REAL Lij_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 + parameters::intMobAniso_param1 * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
				return Lij;
			};
			REAL Lij_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1 * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- parameters::intMobAniso_param2 * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- parameters::intMobAniso_param3 * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- parameters::intMobAniso_param4 * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
				return Lij;
			};
			REAL Lij_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1 * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- parameters::intMobAniso_param2 * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + parameters::intMobAniso_param3), 2.0));
				return Lij;
			};
			REAL Lij_dendrite_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad)/*beta_grad*/;
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1 *
					(3.0 + 4.0 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return Lij;
			};
			REAL Lij_temp(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				return parameters::Lij(alpha_index, beta_index) * std::exp(-parameters::Qij(alpha_index, beta_index) / parameters::R / temperature);
			};
			REAL Lij_temp_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index) * std::exp(-parameters::Qij(alpha_index, beta_index) / parameters::R / temperature);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return Lij;
			};
			REAL Lij_temp_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index) * std::exp(-parameters::Qij(alpha_index, beta_index) / parameters::R / temperature);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 + parameters::intMobAniso_param1 * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
				return Lij;
			};
			REAL Lij_temp_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index) * std::exp(-parameters::Qij(alpha_index, beta_index) / parameters::R / temperature);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1 * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- parameters::intMobAniso_param2 * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- parameters::intMobAniso_param3 * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- parameters::intMobAniso_param4 * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
				return Lij;
			};
			REAL Lij_temp_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index) * std::exp(-parameters::Qij(alpha_index, beta_index) / parameters::R / temperature);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				Lij = Lij * REAL(1.0 - parameters::intMobAniso_param1 * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- parameters::intMobAniso_param2 * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + parameters::intMobAniso_param3), 2.0));
				return Lij;
			};
			REAL Lij_temp_dendrite_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad, REAL temperature) {
				// to be defined
				REAL Lij = parameters::Lij(alpha_index, beta_index) * std::exp(-parameters::Qij(alpha_index, beta_index) / parameters::R / temperature);
				// - anisotropic
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad)/*beta_grad*/;
				Lij = Lij * REAL(1.0 + parameters::intMobAniso_param1 *
					(3.0 + 4.0 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return Lij;
			}
			// - interface energy
			REAL xi_ab(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				return parameters::xi_ab(PhiProperties::instance()[alpha_index], PhiProperties::instance()[beta_index]);
			};
			REAL xi_ab_cubic(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = PhiProperties::instance()[alpha_index], phi_beta_property = PhiProperties::instance()[beta_index];
				REAL loc_xi_ab = parameters::xi_ab(phi_alpha_property, phi_beta_property);
				// - 
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				loc_xi_ab = loc_xi_ab * REAL(1.0 + parameters::intEnAniso_param1 *
					(1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return loc_xi_ab;
			};
			REAL xi_ab_hex_boettger(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = PhiProperties::instance()[alpha_index], phi_beta_property = PhiProperties::instance()[beta_index];
				REAL loc_xi_ab = parameters::xi_ab(phi_alpha_property, phi_beta_property);
				// - 
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				loc_xi_ab = loc_xi_ab * REAL(1.0 - parameters::intEnAniso_param1 * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] +
					15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
				return loc_xi_ab;
			};
			REAL xi_ab_hex_sun(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = PhiProperties::instance()[alpha_index], phi_beta_property = PhiProperties::instance()[beta_index];
				REAL loc_xi_ab = parameters::xi_ab(phi_alpha_property, phi_beta_property);
				// - 
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				loc_xi_ab = loc_xi_ab * REAL(1.0 + parameters::intEnAniso_param1 * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					+ parameters::intEnAniso_param2 * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0) - 30.0 * norm[2] * norm[2] + 3.0)
					+ parameters::intEnAniso_param3 * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0) - 315.0 * std::pow(norm[2], 4.0)
						+ 105.0 * norm[2] * norm[2] - 5.0)
					+ parameters::intEnAniso_param4 * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0) - std::pow(norm[1], 6.0)));
				return loc_xi_ab;
			};
			REAL xi_ab_hex_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = PhiProperties::instance()[alpha_index], phi_beta_property = PhiProperties::instance()[beta_index];
				REAL loc_xi_ab = parameters::xi_ab(phi_alpha_property, phi_beta_property);
				// - 
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad);
				loc_xi_ab = loc_xi_ab * REAL(1.0 + parameters::intEnAniso_param1 * std::pow((3.0 * norm[2] * norm[2] - 1.0), 2.0)
					+ parameters::intEnAniso_param2 * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + parameters::intEnAniso_param3), 2.0));
				return loc_xi_ab;
			};
			REAL xi_ab_dendrite_yang(size_t alpha_index, size_t beta_index, REAL alpha_phi, REAL beta_phi, Vector3& alpha_grad, Vector3& beta_grad) {
				size_t phi_alpha_property = PhiProperties::instance()[alpha_index], phi_beta_property = PhiProperties::instance()[beta_index];
				REAL loc_xi_ab = parameters::xi_ab(phi_alpha_property, phi_beta_property);
				// - 
				Vector3 norm = parameters::grain_rotation_matrix[beta_index] * normals(alpha_phi, beta_phi, alpha_grad, beta_grad)/*beta_grad*/;
				loc_xi_ab = loc_xi_ab * REAL(1.0 - parameters::intEnAniso_param1 *
					(3.0 + 4.0 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
				return loc_xi_ab;
			};
			REAL xi_abc(size_t alpha_index, size_t beta_index, size_t gamma_index) {
				return parameters::xi_abc(PhiProperties::instance()[alpha_index], PhiProperties::instance()[beta_index], PhiProperties::instance()[gamma_index]);
			};
			// - source
			REAL noise_pairwise_acc(size_t x, size_t y, size_t z, size_t alpha_index, size_t beta_index) {
				REAL noise = 0.0;
				if (main_iterator::Current_ITE_step < parameters::phi_noise_begin ||
					main_iterator::Current_ITE_step > parameters::phi_noise_end)
					return noise;
				if (main_iterator::Current_ITE_step % parameters::phi_noise_frequency != 0)
					return noise;
				if (parameters::is_phi_noise[alpha_index] && parameters::is_phi_noise[beta_index]) {
#ifdef _OPENMP
#pragma omp critical(phi_noise_random_number)
#endif
					{
						noise = parameters::phi_noise_amplitude * parameters::phi_noise_real_dist(parameters::phi_noise_gen);
					}
				}
				return noise;
			}
			// - interface energy model
			REAL dfint_dphi_grad_S1996_acc(FIELD_PhiTemp& point, size_t phi_index) {
				REAL grad = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && point.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (point.intflag[phi_bIndex] && phi_bIndex != phi_index) {
						grad += _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index], point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex])
							* (point.grad_phi[phi_bIndex] * point.grad_phi[phi_bIndex] * point.old_phi[phi_index] * 2
								+ point.old_phi[phi_index] * point.old_phi[phi_bIndex] * point.lap_phi[phi_bIndex]
								- point.grad_phi[phi_index] * point.grad_phi[phi_bIndex] * point.old_phi[phi_bIndex] * 2
								- point.old_phi[phi_bIndex] * point.old_phi[phi_bIndex] * point.lap_phi[phi_index]);
					}
				}
				return 2 * parameters::interface_width * grad;
			};
			REAL dfint_dphi_grad_S1999_acc(FIELD_PhiTemp& point, size_t phi_index) {
				REAL grad = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && point.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (point.intflag[phi_bIndex] && phi_bIndex != phi_index) {
						grad += _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index], point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex]) * point.lap_phi[phi_bIndex];
					}
				}
				return parameters::interface_width * grad;
			};
			REAL dfint_dphi_pot_Nobstacle_acc(FIELD_PhiTemp& point, size_t phi_index) {
				REAL pot = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && point.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (point.intflag[phi_bIndex] && phi_bIndex != phi_index) {
						pot += 16 / PI2 * _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index],point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex]) * point.old_phi[phi_bIndex];
						for (size_t index2 = index + 1; index2 < parameters::PHI_ACC_NUMBER && point.active_index[index2] != parameters::PAIRWISE_ACC_STOP; index2++) {
							size_t phi_gIndex = point.active_index[index2];
							if (point.intflag[phi_gIndex] && phi_gIndex != phi_index)
								pot += _xi_abc(phi_index, phi_bIndex, phi_gIndex) * point.old_phi[phi_bIndex] * point.old_phi[phi_gIndex];
						}
					}
				}
				return pot / parameters::interface_width;
			};
			REAL dfint_dphi_pot_Nwell_acc(FIELD_PhiTemp& point, size_t phi_index) {
				REAL pot = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && point.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t phi_bIndex = point.active_index[index];
					if (point.intflag[phi_bIndex] && phi_bIndex != phi_index) {
						pot += 18 * _xi_ab(phi_index, phi_bIndex, point.old_phi[phi_index],point.old_phi[phi_bIndex], point.grad_phi[phi_index], point.grad_phi[phi_bIndex]) 
							* point.old_phi[phi_bIndex] * point.old_phi[phi_bIndex];
						for (size_t index2 = index + 1; index2 < parameters::PHI_ACC_NUMBER && point.active_index[index2] != parameters::PAIRWISE_ACC_STOP; index2++) {
							size_t phi_gIndex = point.active_index[index2];
							if (point.intflag[phi_gIndex] && phi_gIndex != phi_index)
								pot += 6 * _xi_abc(phi_index, phi_bIndex, phi_gIndex)
								* point.old_phi[phi_bIndex] * point.old_phi[phi_gIndex] * point.old_phi[phi_bIndex] * point.old_phi[phi_gIndex];
						}
					}
				}
				return point.old_phi[phi_index] / parameters::interface_width * pot;
			};
			// - 
			REAL dfint_dphi_pairwise_S2009(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index) {
				return _xi_ab(alpha_index, beta_index, point.old_phi[alpha_index], point.old_phi[beta_index], point.grad_phi[alpha_index], point.grad_phi[beta_index])
					* (PI * PI / 2 / parameters::interface_width * (point.old_phi[alpha_index] - point.old_phi[beta_index])
						+ parameters::interface_width * (point.old_phi[beta_index] * point.lap_phi[alpha_index] - point.old_phi[alpha_index] * point.lap_phi[beta_index]));
			};
			REAL dfint_dphi_pairwise_S1996_NO(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1996_acc(point, beta_index) + dfint_dphi_pot_Nobstacle_acc(point, beta_index) -
					dfint_dphi_grad_S1996_acc(point, alpha_index) - dfint_dphi_pot_Nobstacle_acc(point, alpha_index);
			};
			REAL dfint_dphi_pairwise_S1996_NW(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1996_acc(point, beta_index) + dfint_dphi_pot_Nwell_acc(point, beta_index) -
					dfint_dphi_grad_S1996_acc(point, alpha_index) - dfint_dphi_pot_Nwell_acc(point, alpha_index);
			};
			REAL dfint_dphi_pairwise_S1999_NO(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1999_acc(point, beta_index) + dfint_dphi_pot_Nobstacle_acc(point, beta_index) -
					dfint_dphi_grad_S1999_acc(point, alpha_index) - dfint_dphi_pot_Nobstacle_acc(point, alpha_index);
			};
			REAL dfint_dphi_pairwise_S1999_NW(FIELD_PhiTemp& point, size_t alpha_index, size_t beta_index) {
				return dfint_dphi_grad_S1999_acc(point, beta_index) + dfint_dphi_pot_Nwell_acc(point, beta_index) -
					dfint_dphi_grad_S1999_acc(point, alpha_index) - dfint_dphi_pot_Nwell_acc(point, alpha_index);
			};
			REAL dfbulk_dphi_pairwise_acc(long long x, long long y, long long z, size_t phi_index) {
				REAL dfbulk_sum = 0.0;
				for (auto dfbulk : parameters::delt_Fbulk_delt_phi)
					dfbulk_sum += dfbulk(x, y, z, PhiProperties::instance()[phi_index]);
				return dfbulk_sum;
			}
			REAL source_pairwise_acc(long long x, long long y, long long z, size_t alpha_index, size_t beta_index) {
				REAL source_sum = 0.0;
				for (auto source : parameters::source_alpha_beta)
					source_sum += source(x, y, z, PhiProperties::instance()[alpha_index], PhiProperties::instance()[beta_index]);
				return source_sum;
			}
			// - 
			REAL dfbulk_dphi_0(size_t x, size_t y, size_t z, size_t phi_index) {
				return parameters::f_bulk_0[PhiProperties::instance()[phi_index]];
			}
			// - 
			void pairwise_normalize_acc(std::vector<size_t>& active_index, std::vector<int>& intflag,
				const Matrix1D<REAL>& old_phi, std::vector<REAL>& phi_increment) {
				REAL scale = 1.0, increment = 0.0;
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t acc_index = active_index[index];
					if (intflag[acc_index]) {
						increment = time_parameters::delt_t * phi_increment[acc_index];
						if (isTwoREALEquality(increment, 0))
							continue;
						else if ((old_phi[acc_index] + increment) > 1) {
							REAL p_scale = (1 - old_phi[acc_index]) / increment;
							if (p_scale < scale)
								scale = p_scale;
						}
						else if ((old_phi[acc_index] + increment) < 0) {
							REAL p_scale = abs(old_phi[acc_index] / increment);
							if (p_scale < scale)
								scale = p_scale;
						}
					}
				}
				for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
					size_t acc_index = active_index[index];
					if (intflag[acc_index])
						phi_increment[acc_index] *= scale;
				}
			}
			//=========================================================================================================================================
			void init_phi_pair_wise() {
				if (parameters::is_phi_noise_rand)
					parameters::phi_noise_gen.seed(parameters::phi_noise_rd());
				else
					parameters::phi_noise_gen.seed(static_cast<std::mt19937::result_type>(parameters::phi_noise_seed));

				for (long long x = 0; x < main_field::phase_field.Nx(); x++)
					for (long long y = 0; y < main_field::phase_field.Ny(); y++)
						for (long long z = 0; z < main_field::phase_field.Nz(); z++) {
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
							field_var.init_phi(main_field::phi_number, parameters::PHI_ACC_NUMBER);
							if (parameters::is_phi_normalized) {
								REAL sum_phi = 0.0;
								for (size_t index = 0; index < main_field::phi_number; index++)
									sum_phi += phi[index];
								if (sum_phi > SYS_EPSILON) {
									for (size_t index = 0; index < main_field::phi_number; index++)
										phi[index] /= sum_phi;
								}
							}
							for (size_t index = 0; index < main_field::phi_number; index++) {
								field_var.old_phi[index] = phi[index];
								field_var.new_phi[index] = phi[index];
								field_var.intflag[index] = currentFlag(x, y, z, index);
								if (field_var.intflag[index])
									interphase_gradient_lapace_calculation(x, y, z, index);
							}
						}
			}
			void pre_calculation_phi_pair_wise() {
				parameters::MAX_ACTIVE_PHI_NUMBER = 0;
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
							field_var.active_number = 0;
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER; index++)
								field_var.active_index[index] = parameters::PAIRWISE_ACC_STOP;
							for (size_t index = 0; index < main_field::phi_number; index++) {
								field_var.old_phi[index] = phi[index]; // save old phi for other physical field
								field_var.new_phi[index] = 0; // increment during calculation
								if ((field_var.intflag[index] || phi[index] > parameters::PhiCon_Cut_Off) && field_var.active_number < parameters::PHI_ACC_NUMBER) {
									field_var.active_index[field_var.active_number] = index; // save the active phi index 
									field_var.active_number++;
								}
								if (field_var.intflag[index]) {
									interphase_gradient_lapace_calculation(x, y, z, index);
								}
								else {
									field_var.lap_phi[index] = 0;
									field_var.grad_phi[index][0] = 0;
									field_var.grad_phi[index][1] = 0;
									field_var.grad_phi[index][2] = 0;
								}
							}
							if (field_var.active_number > parameters::MAX_ACTIVE_PHI_NUMBER)
								parameters::MAX_ACTIVE_PHI_NUMBER = field_var.active_number;
						}
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
							if (field_var.active_number <= 1)
								continue;
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							std::vector<REAL>& phi_increment = field_var.new_phi;
							// interface energy
							for (size_t aIndex = 0; aIndex < parameters::PHI_ACC_NUMBER && field_var.active_index[aIndex] != parameters::PAIRWISE_ACC_STOP; aIndex++) {
								size_t alpha_index = field_var.active_index[aIndex];
								if (!field_var.intflag[alpha_index])
									continue;
								for (size_t bIndex = aIndex + 1; bIndex < parameters::PHI_ACC_NUMBER && field_var.active_index[bIndex] != parameters::PAIRWISE_ACC_STOP; bIndex++) {
									size_t beta_index = field_var.active_index[bIndex];
									if (!field_var.intflag[beta_index])
										continue;
									REAL int_incre_b_a = _Lij(alpha_index, beta_index, phi[alpha_index], phi[beta_index],
										field_var.grad_phi[alpha_index], field_var.grad_phi[beta_index], field_var.new_temp) / parameters::interface_width
										* dfint_dphi_pairwise_acc(field_var, alpha_index, beta_index);
									if (parameters::is_phi_normalized) {
										if ((int_incre_b_a > SYS_EPSILON && (phi[alpha_index] > SYS_EPSILON_R || phi[beta_index] < SYS_EPSILON))
											|| (int_incre_b_a < -SYS_EPSILON && (phi[alpha_index] < SYS_EPSILON || phi[beta_index] > SYS_EPSILON_R)))
											int_incre_b_a = 0.0;
									}
									phi_increment[alpha_index] += int_incre_b_a;
									phi_increment[beta_index] -= int_incre_b_a;
#ifdef _DEBUG
									if (std::isnan(int_incre_b_a)) {
										std::cout << "DEBUG: interface energy (pair-wise functions) error !" << std::endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
							// driving force and source term
							for (size_t aIndex = 0; aIndex < parameters::PHI_ACC_NUMBER && field_var.active_index[aIndex] != parameters::PAIRWISE_ACC_STOP; aIndex++) {
								size_t alpha_index = field_var.active_index[aIndex];
								if (!field_var.intflag[alpha_index])
									continue;
								for (size_t bIndex = aIndex + 1; bIndex < parameters::PHI_ACC_NUMBER && field_var.active_index[bIndex] != parameters::PAIRWISE_ACC_STOP; bIndex++) {
									size_t beta_index = field_var.active_index[bIndex];
									if (!field_var.intflag[beta_index])
										continue;
									if (phi[alpha_index] > parameters::PhiCon_Cut_Off && phi[beta_index] > parameters::PhiCon_Cut_Off) {
										REAL bulk_increment_b_a = _Lij(alpha_index, beta_index, phi[alpha_index], phi[beta_index],
											field_var.grad_phi[alpha_index], field_var.grad_phi[beta_index], field_var.new_temp) / parameters::interface_width
											* (dfbulk_dphi_pairwise_acc(x, y, z, beta_index) - dfbulk_dphi_pairwise_acc(x, y, z, alpha_index))
											+ source_pairwise_acc(x, y, z, alpha_index, beta_index);
										if (parameters::is_phi_normalized) {
											if ((bulk_increment_b_a > SYS_EPSILON && (phi[alpha_index] > SYS_EPSILON_R || phi[beta_index] < SYS_EPSILON))
												|| (bulk_increment_b_a < -SYS_EPSILON && (phi[alpha_index] < SYS_EPSILON || phi[beta_index] > SYS_EPSILON_R)))
												bulk_increment_b_a = 0.0;
										}
										phi_increment[alpha_index] += bulk_increment_b_a;
										phi_increment[beta_index] -= bulk_increment_b_a;
									}
								}
							}
							// numerical treatment to nomalize the phi
							if (parameters::is_phi_normalized)
								pairwise_normalize_acc(field_var.active_index, field_var.intflag, phi, phi_increment);
						}
			}
			REAL solve_phi_pair_wise() {
				REAL MAX_PHI_INCREMENT = 0.0;
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
							Matrix1D<REAL>& phi = main_field::phase_field(x, y, z);
							if (field_var.active_number <= 1) {
								// save all phase field variables to new phi
								for (size_t index = 0; index < main_field::phi_number; index++)
									field_var.new_phi[index] = phi[index];
								continue;
							}
							std::vector<REAL>& phi_increment = field_var.new_phi;
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_var.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
								size_t phi_index = field_var.active_index[index];
								if (field_var.intflag[phi_index]) {
									field_var.new_phi[phi_index] = field_var.old_phi[phi_index] + (time_parameters::delt_t)*phi_increment[phi_index];
#ifdef _OPENMP
#pragma omp critical
#endif
									{
										if (abs(field_var.old_phi[phi_index] - field_var.new_phi[phi_index]) > MAX_PHI_INCREMENT)
											MAX_PHI_INCREMENT = abs(field_var.old_phi[phi_index] - field_var.new_phi[phi_index]);
									}
#ifdef _DEBUG
									if (std::isnan(field_var.new_phi[phi_index])) {
										std::string error_report = "ERROR : Phase fraction is NaN in x = " + std::to_string(x) + " , y = " + std::to_string(y)
											+ " , z = " + std::to_string(z) + " position, the phase index is " + std::to_string(phi_index) + "\n";
										std::cout << error_report << std::endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
							// normalize new phi
							if (parameters::is_phi_normalized) {
								normalize_new_phi(field_var);
							}
							// save changed phis to phase field
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_var.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
								size_t phi_index = field_var.active_index[index];
								if (field_var.intflag[phi_index])
									phi[phi_index] = field_var.new_phi[phi_index];
							}
							// save all phase field variables to new phi
							// for (size_t index = 0; index < main_field::phi_number; index++)
							// 	field_var.new_phi[index] = phi[index];
						}
				// update phi flag
#pragma omp parallel for
				for (long long x = main_field::phase_field.COMP_X_BGN(); x <= main_field::phase_field.COMP_X_END(); x++)
					for (long long y = main_field::phase_field.COMP_Y_BGN(); y <= main_field::phase_field.COMP_Y_END(); y++)
						for (long long z = main_field::phase_field.COMP_Z_BGN(); z <= main_field::phase_field.COMP_Z_END(); z++) {
							FIELD_PhiTemp& field_var = parameters::PhiTemp_field(x, y, z);
							if (field_var.active_number <= 1)
								continue;
							for (size_t index = 0; index < parameters::PHI_ACC_NUMBER && field_var.active_index[index] != parameters::PAIRWISE_ACC_STOP; index++) {
								size_t phi_index = field_var.active_index[index];
								if (field_var.intflag[phi_index]) {
									int new_flag = currentFlag(x, y, z, phi_index);
									if (new_flag > field_var.intflag[phi_index]) {
										// this phi change from near_interface to interface
										upgradeFlag(x, y, z, phi_index);
									}
									field_var.intflag[phi_index] = new_flag;
								}
							}
						}
				main_field::phase_field.do_boundary_condition();
				return MAX_PHI_INCREMENT;
			}
			//=========================================================================================================================================
			// - output
			void write_scalar_active_phi_number(std::ofstream& fout) {
				fout << "<DataArray type = \"Float64\" Name = \"" << "active_phis" <<
					"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
				for (size_t k = write_vts::z_begin; k <= write_vts::z_end; ++k)
					for (size_t j = write_vts::y_begin; j <= write_vts::y_end; ++j)
						for (size_t i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
							fout << parameters::PhiTemp_field(i, j, k).active_number << std::endl;
						}
				fout << "</DataArray>" << std::endl;
			}
		}
	}
}
