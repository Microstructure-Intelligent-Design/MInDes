/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "Solver_Mechanics.h"
using namespace std;

namespace pf {

	double MechanicalField_Explicit::get_u_main_node(int _x, int _y, int _z) {
		int now_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW];
		VectorNode& u = U(_x + Minus_0_5, _y, _z);
		return (u.vals[now_domain] + u.get_neighbor_node(Direction::x_up).vals[now_domain]) / 2.0;
	}
	double MechanicalField_Explicit::get_v_main_node(int _x, int _y, int _z) {
		int now_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW];
		VectorNode& v = V(_x, _y + Minus_0_5, _z);
		return (v.vals[now_domain] + v.get_neighbor_node(Direction::y_up).vals[now_domain]) / 2.0;
	}
	double MechanicalField_Explicit::get_w_main_node(int _x, int _y, int _z) {
		int now_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW];
		VectorNode& w = W(_x, _y, _z + Minus_0_5);
		return (w.vals[now_domain] + w.get_neighbor_node(Direction::z_up).vals[now_domain]) / 2.0;
	}

	void MechanicalField_Explicit::cal_parameters_before_calculation() {
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					cal_parameters_for_main_domain(node, phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z);
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        average_stiffness += node.customMatrix6x6s[ExternalFields::MECH_stiffness];
                    }
				}
        average_stiffness /= double(phaseMesh->limit_x * phaseMesh->limit_y * phaseMesh->limit_z);
        average_compliences = average_stiffness.get_inverted_matrix();
	}

	void MechanicalField_Explicit::init_velocity_field() {
		int now_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW], old_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_OLD];
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
					u.vals[now_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][0];
					u.vals[old_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][0];
					v.vals[now_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][1];
					v.vals[old_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][1];
					w.vals[now_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][2];
					w.vals[old_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][2];
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					u.vals[now_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][0];
					u.vals[old_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][0];
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					v.vals[now_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][1];
					v.vals[old_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][1];
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					w.vals[now_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][2];
					w.vals[old_domain] = node.customVec6s[ExternalFields::MECH_eigen_strain][2];
				}
		}
	}

    double MechanicalField_Explicit::boundary_condition(vStress& AverageStress, vStrain& AverageStrain, double incre_rate) {
        vStress TargetStress, oldTargetStrain;
        TargetStress.set_to_zero();
        oldTargetStrain.set_to_zero();

        TargetStress[0] = (-AverageStress[0] + applied_stress[0]) * incre_rate;
        TargetStress[1] = (-AverageStress[1] + applied_stress[1]) * incre_rate;
        TargetStress[2] = (-AverageStress[2] + applied_stress[2]) * incre_rate;

        oldTargetStrain = AverageStrain;

        AverageStrain.set_to_zero();

        for (int n = 0; n < 3; n++)
            for (int m = 0; m < 6; m++)
            {
                AverageStrain[n] += TargetStress[m] * average_compliences(n, m);
            }

        double MAXTargetStrainDifference = 0.0;
        for (int n = 0; n < 3; n++)
            if (fabs(AverageStrain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
            {
                MAXTargetStrainDifference = fabs(AverageStrain[n] - oldTargetStrain[n]);
            }

        if (average_stiffness(0, 0) == 0 or average_stiffness(1, 1) == 0 or average_stiffness(2, 2) == 0
            or average_stiffness(1, 2) == 0 or average_stiffness(0, 2) == 0 or average_stiffness(0, 1) == 0)
        {
            cout << "> ERROR ! In SetElasticBoundaryConditions(), Zero component in AverageElasticConstants." << endl;
            SYS_PROGRAM_STOP;
        }

        if (AppStrainMask[0] and !AppStrainMask[1] and !AppStrainMask[2])   /// XX
        {
            AverageStrain[0] = applied_strain[0];
            AverageStrain[1] = (-average_stiffness(0, 2) * average_stiffness(1, 2) * applied_strain[0] +
                average_stiffness(0, 1) * average_stiffness(2, 2) * applied_strain[0] -
                average_stiffness(2, 2) * TargetStress[1] + average_stiffness(1, 2) * TargetStress[2]) /
                (average_stiffness(1, 2) * average_stiffness(1, 2) - average_stiffness(1, 1) * average_stiffness(2, 2));
            AverageStrain[2] = (-average_stiffness(1, 1) * AverageStrain[1] - average_stiffness(0, 1) * applied_strain[0] + TargetStress[1]) /
                average_stiffness(1, 2);
        }

        if (!AppStrainMask[0] and AppStrainMask[1] and !AppStrainMask[2])   /// YY
        {
            AverageStrain[0] = (-average_stiffness(0, 2) * average_stiffness(1, 2) * applied_strain[1] +
                average_stiffness(0, 1) * average_stiffness(2, 2) * applied_strain[1] -
                average_stiffness(2, 2) * TargetStress[0] + average_stiffness(0, 2) * TargetStress[2]) /
                (average_stiffness(0, 2) * average_stiffness(0, 2) - average_stiffness(0, 0) * average_stiffness(2, 2));
            AverageStrain[1] = applied_strain[1];
            AverageStrain[2] = (-average_stiffness(0, 0) * AverageStrain[0] - average_stiffness(0, 1) * applied_strain[1] + TargetStress[0]) /
                average_stiffness(0, 2);
        }

        if (!AppStrainMask[0] and !AppStrainMask[1] and AppStrainMask[2])   /// ZZ
        {
            AverageStrain[0] = (average_stiffness(0, 2) * average_stiffness(1, 1) * applied_strain[2] -
                average_stiffness(0, 1) * average_stiffness(1, 2) * applied_strain[2] -
                average_stiffness(1, 1) * TargetStress[0] + average_stiffness(0, 1) * TargetStress[1]) /
                (average_stiffness(0, 1) * average_stiffness(0, 1) - average_stiffness(0, 0) * average_stiffness(1, 1));

            AverageStrain[1] = (-average_stiffness(0, 0) * AverageStrain[0] - average_stiffness(0, 2) * applied_strain[2] + TargetStress[0]) /
                average_stiffness(0, 1);
            AverageStrain[2] = applied_strain[2];
        }

        if (AppStrainMask[0] and AppStrainMask[1] and !AppStrainMask[2])    /// XX & YY
        {
            AverageStrain[0] = applied_strain[0];
            AverageStrain[1] = applied_strain[1];
            AverageStrain[2] = (-average_stiffness(0, 2) * applied_strain[0] - average_stiffness(1, 2) * applied_strain[1] + TargetStress[2]) /
                average_stiffness(2, 2);
        }

        if (AppStrainMask[0] and !AppStrainMask[1] and AppStrainMask[2])    /// XX & ZZ
        {
            AverageStrain[0] = applied_strain[0];
            AverageStrain[1] = (-average_stiffness(0, 1) * applied_strain[0] - average_stiffness(1, 2) * applied_strain[2] + TargetStress[1]) /
                average_stiffness(1, 1);
            AverageStrain[2] = applied_strain[2];
        }

        if (!AppStrainMask[0] and AppStrainMask[1] and AppStrainMask[2])    /// YY & ZZ
        {
            AverageStrain[0] = (-average_stiffness(0, 1) * applied_strain[1] - average_stiffness(0, 2) * applied_strain[2] + TargetStress[0]) /
                average_stiffness(0, 0);
            AverageStrain[1] = applied_strain[1];
            AverageStrain[2] = applied_strain[2];
        }

        if (AppStrainMask[0] and AppStrainMask[1] and AppStrainMask[2])    /// XX & YY & ZZ
        {
            AverageStrain[0] = applied_strain[0];
            AverageStrain[1] = applied_strain[1];
            AverageStrain[2] = applied_strain[2];
        }
        return MAXTargetStrainDifference;
    }

	void MechanicalField_Explicit::cal_stress(vStress& AverageStress, bool is_plastic) {
        AverageStress.set_to_zero();
        if (is_plastic) {
#pragma omp parallel for
            for (int x = 0; x < phaseMesh->limit_x; x++)
                for (int y = 0; y < phaseMesh->limit_y; y++)
                    for (int z = 0; z < phaseMesh->limit_z; z++) {
                        PhaseNode& node = (*phaseMesh)(x, y, z);
                        node.customVec6s[ExternalFields::MECH_stress] = node.customMatrix6x6s[ExternalFields::MECH_stiffness] * (node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain] - node.customVec6s[ExternalFields::MECH_plastic_strain]);
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            AverageStress += node.customVec6s[ExternalFields::MECH_stress];
                        }
                    }
        }
        else {
#pragma omp parallel for
            for (int x = 0; x < phaseMesh->limit_x; x++)
                for (int y = 0; y < phaseMesh->limit_y; y++)
                    for (int z = 0; z < phaseMesh->limit_z; z++) {
                        PhaseNode& node = (*phaseMesh)(x, y, z);
                        node.customVec6s[ExternalFields::MECH_stress] = node.customMatrix6x6s[ExternalFields::MECH_stiffness] * (node.customVec6s[ExternalFields::MECH_strain] - node.customVec6s[ExternalFields::MECH_eigen_strain]);
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            AverageStress += node.customVec6s[ExternalFields::MECH_stress];
                        }
                    }
        }
        AverageStress /= phaseMesh->limit_x * phaseMesh->limit_y * phaseMesh->limit_z;
	}

	double MechanicalField_Explicit::evolve_momentum_equation(double mass_density, double mechanic_dt, vStrain& AverageStrain) {
		double MAX_VARIATION = 0.0;
		double dr = phaseMesh->dr;
		int now_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW], old_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_OLD], future_domain = velocity_doamin_index[VelocityDomainIndex::VDIndex_FUTURE];
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5); 
					Vector6& stress_xdown = node.get_neighbor_node(Direction::x_down).customVec6s[ExternalFields::MECH_stress],
						stress_ydown = node.get_neighbor_node(Direction::y_down).customVec6s[ExternalFields::MECH_stress],
						stress_zdown = node.get_neighbor_node(Direction::z_down).customVec6s[ExternalFields::MECH_stress],
						stress = node.customVec6s[ExternalFields::MECH_stress];
					// momentum equation for u
					{
						double A = (stress[0] - stress_xdown[0]) / dr + (stress[5] - stress_ydown[5]) / dr + (stress[4] - stress_zdown[4]) / dr;

						u.vals[future_domain] = A * mechanic_dt * mechanic_dt / mass_density - u.vals[old_domain] + 2.0 * u.vals[now_domain];
					}
					// momentum equation for v
					{
						double B = (stress[5] - stress_xdown[5]) / dr + (stress[1] - stress_ydown[1]) / dr + (stress[3] - stress_zdown[3]) / dr;

						v.vals[future_domain] = B * mechanic_dt * mechanic_dt / mass_density - v.vals[old_domain] + 2.0 * v.vals[now_domain];
					}
					// momentum equation for w
					{
						double C = (stress[4] - stress_xdown[4]) / dr + (stress[3] - stress_ydown[3]) / dr + (stress[2] - stress_zdown[2]) / dr;

						w.vals[future_domain] = C * mechanic_dt * mechanic_dt / mass_density - w.vals[old_domain] + 2.0 * w.vals[now_domain];
					}
				}
		if (U.x_up_bc != PERIODIC) {
			int x = phaseMesh->limit_x - 1;
#pragma omp parallel for
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Plus_0_5, y, z);
					Vector6& stress_xup = node.get_neighbor_node(Direction::x_up).customVec6s[ExternalFields::MECH_stress],
						stress_yup = node.get_neighbor_node(Direction::y_up).customVec6s[ExternalFields::MECH_stress],
						stress_zup = node.get_neighbor_node(Direction::z_up).customVec6s[ExternalFields::MECH_stress],
						stress = node.customVec6s[ExternalFields::MECH_stress];
					// momentum equation for u
					{
						double A = (stress_xup[0] - stress[0]) / dr + (stress_yup[5] - stress[5]) / dr + (stress_zup[4] - stress[4]) / dr;

						u.vals[future_domain] = A * mechanic_dt * mechanic_dt / mass_density - u.vals[old_domain] + 2.0 * u.vals[now_domain];
					}
				}
		}
		if (V.y_up_bc != PERIODIC) {
			int y = phaseMesh->limit_y - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& v = V(x, y + Plus_0_5, z);
					Vector6& stress_xup = node.get_neighbor_node(Direction::x_up).customVec6s[ExternalFields::MECH_stress],
						stress_yup = node.get_neighbor_node(Direction::y_up).customVec6s[ExternalFields::MECH_stress],
						stress_zup = node.get_neighbor_node(Direction::z_up).customVec6s[ExternalFields::MECH_stress],
						stress = node.customVec6s[ExternalFields::MECH_stress];
					// momentum equation for v
					{
						double B = (stress_xup[5] - stress[5]) / dr + (stress_yup[1] - stress[1]) / dr + (stress_zup[3] - stress[3]) / dr;

						v.vals[future_domain] = B * mechanic_dt * mechanic_dt / mass_density - v.vals[old_domain] + 2.0 * v.vals[now_domain];
					}
				}
		}
		if (W.z_up_bc != PERIODIC) {
			int z = phaseMesh->limit_z - 1;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& w = W(x, y, z + Plus_0_5);
					Vector6& stress_xup = node.get_neighbor_node(Direction::x_up).customVec6s[ExternalFields::MECH_stress],
						stress_yup = node.get_neighbor_node(Direction::y_up).customVec6s[ExternalFields::MECH_stress],
						stress_zup = node.get_neighbor_node(Direction::z_up).customVec6s[ExternalFields::MECH_stress],
						stress = node.customVec6s[ExternalFields::MECH_stress];
					// momentum equation for w
					{
						double C = (stress_xup[4] - stress[4]) / dr + (stress_yup[3] - stress[3]) / dr + (stress_zup[2] - stress[2]) / dr;

						w.vals[future_domain] = C * mechanic_dt * mechanic_dt / mass_density - w.vals[old_domain] + 2.0 * w.vals[now_domain];
					}
				}
		}

		// calculate total strain
#pragma omp parallel for
		for (int x = 0; x < phaseMesh->limit_x; x++)
			for (int y = 0; y < phaseMesh->limit_y; y++)
				for (int z = 0; z < phaseMesh->limit_z; z++) {
					PhaseNode& node = (*phaseMesh)(x, y, z);
					VectorNode& u = U(x + Minus_0_5, y, z);
					VectorNode& v = V(x, y + Minus_0_5, z);
					VectorNode& w = W(x, y, z + Minus_0_5);
                    Vector6 vec_old = node.customVec6s[ExternalFields::MECH_strain];
					node.customVec6s[ExternalFields::MECH_strain][0] = (u.get_neighbor_node(Direction::x_up).vals[future_domain] - u.vals[future_domain]) / dr + AverageStrain[0];
					node.customVec6s[ExternalFields::MECH_strain][1] = (v.get_neighbor_node(Direction::y_up).vals[future_domain] - v.vals[future_domain]) / dr + AverageStrain[1];
					node.customVec6s[ExternalFields::MECH_strain][2] = (w.get_neighbor_node(Direction::z_up).vals[future_domain] - w.vals[future_domain]) / dr + AverageStrain[2];
					node.customVec6s[ExternalFields::MECH_strain][3] = 0.5 * ((v.get_neighbor_node(Direction::z_up).vals[future_domain] - v.vals[future_domain]) 
																			+ (w.get_neighbor_node(Direction::y_up).vals[future_domain] - w.vals[future_domain])) / dr + AverageStrain[3];
					node.customVec6s[ExternalFields::MECH_strain][4] = 0.5 * ((u.get_neighbor_node(Direction::z_up).vals[future_domain] - u.vals[future_domain])
																			+ (w.get_neighbor_node(Direction::x_up).vals[future_domain] - w.vals[future_domain])) / dr + AverageStrain[4];
					node.customVec6s[ExternalFields::MECH_strain][5] = 0.5 * ((u.get_neighbor_node(Direction::y_up).vals[future_domain] - u.vals[future_domain])
																			+ (v.get_neighbor_node(Direction::x_up).vals[future_domain] - v.vals[future_domain])) / dr + AverageStrain[5];
                    for (int n = 0; n < 6; n++)
                    {
                        double locStrainDifference = fabs(node.customVec6s[ExternalFields::MECH_strain][n] - vec_old[n]);  //difference for locStrain
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (locStrainDifference > MAX_VARIATION)
                            {
                                MAX_VARIATION = locStrainDifference;
                            }
                        }
                    }
				}
		
		// exchange domain
		int buff = velocity_doamin_index[VelocityDomainIndex::VDIndex_OLD];
		velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW] = velocity_doamin_index[VelocityDomainIndex::VDIndex_FUTURE];
		velocity_doamin_index[VelocityDomainIndex::VDIndex_OLD] = velocity_doamin_index[VelocityDomainIndex::VDIndex_NOW];
		velocity_doamin_index[VelocityDomainIndex::VDIndex_FUTURE] = buff;
		return MAX_VARIATION;
	}

    PhaseNode& MechanicalField_Implicit::get_phaseNode_in_elasticField(int elas_x, int elas_y, int elas_z) {
        int mirror_x = elas_x, mirror_y = elas_y, mirror_z = elas_z;
        if (elas_x >= phaseMesh->limit_x) {
            mirror_x = mechanicalField.Nx - 1 - elas_x;
        }
        if (elas_y >= phaseMesh->limit_y) {
            mirror_y = mechanicalField.Ny - 1 - elas_y;
        }
        if (elas_z >= phaseMesh->limit_z) {
            mirror_z = mechanicalField.Nz - 1 - elas_z;
        }
        return (*phaseMesh)(mirror_x, mirror_y, mirror_z);
    }

    void MechanicalField_Implicit::SetMAXElasticConstants(vector<Matrix6x6> Cijs) {
        MAX_ElasticConstants.set_to_zero();
        for (auto Cij = Cijs.begin(); Cij < Cijs.end(); Cij++) {
            for (int n = 0; n < 6; n++)
                for (int m = 0; m < 6; m++)
                    if ((*Cij)(n, m) > MAX_ElasticConstants(n, m))
                        MAX_ElasticConstants(n, m) = (*Cij)(n, m);
        }
    }
    void MechanicalField_Implicit::SetMAXElasticConstants(tensor1_matrix6 Cijs) {
        MAX_ElasticConstants.set_to_zero();
        for (auto Cij = Cijs.begin(); Cij < Cijs.end(); Cij++) {
            for (int n = 0; n < 6; n++)
                for (int m = 0; m < 6; m++)
                    if (Cij->val(n, m) > MAX_ElasticConstants(n, m))
                        MAX_ElasticConstants(n, m) = Cij->val(n, m);
        }
    }

    void MechanicalField_Implicit::cal_parameters_before_calculation(bool is_plastic) {
        average_stiffness.set_to_zero();
        if (is_plastic) {
#pragma omp parallel for
            for (int i = 0; i < mechanicalField.Nx; i++)
                for (int j = 0; j < mechanicalField.Ny; j++)
                    for (int k = 0; k < mechanicalField.Nz; k++) {
                        MechanicNode& mech = mechanicalField(i, j, k);
                        PhaseNode& node = get_phaseNode_in_elasticField(i, j, k);
                        //mech.StrainIncrements.set_to_zero();
                        cal_parameters_for_main_domain(node, phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z);
                        mech.EffectiveEigenStrains = node.customVec6s[ExternalFields::MECH_eigen_strain] + node.customVec6s[ExternalFields::MECH_plastic_strain];
                        mech.EffectiveElasticConstants = node.customMatrix6x6s[ExternalFields::MECH_stiffness];
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            average_stiffness += mech.EffectiveElasticConstants;
                        }
#ifdef _DEBUG
                        if (mech.EffectiveEigenStrains.is_nan_val_exist()) {
                            cout << "DEBUG: mech.EffectiveEigenStrains error !" << endl;
                            SYS_PROGRAM_STOP;
                        }
#endif
                    }
        }
        else {
#pragma omp parallel for
            for (int i = 0; i < mechanicalField.Nx; i++)
                for (int j = 0; j < mechanicalField.Ny; j++)
                    for (int k = 0; k < mechanicalField.Nz; k++) {
                        MechanicNode& mech = mechanicalField(i, j, k);
                        PhaseNode& node = get_phaseNode_in_elasticField(i, j, k);
                        //mech.StrainIncrements.set_to_zero();
                        cal_parameters_for_main_domain(node, phaseMesh->limit_x, phaseMesh->limit_y, phaseMesh->limit_z);
                        mech.EffectiveEigenStrains = node.customVec6s[ExternalFields::MECH_eigen_strain];
                        mech.EffectiveElasticConstants = node.customMatrix6x6s[ExternalFields::MECH_stiffness];
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            average_stiffness += mech.EffectiveElasticConstants;
                        }
#ifdef _DEBUG
                        if (mech.EffectiveEigenStrains.is_nan_val_exist()) {
                            cout << "DEBUG: mech.EffectiveEigenStrains error !" << endl;
                            SYS_PROGRAM_STOP;
                        }
#endif
                    }
        }
        average_stiffness /= double(mechanicalField.Nx * mechanicalField.Ny * mechanicalField.Nz);
        average_compliences = average_stiffness.get_inverted_matrix();
    }

    void MechanicalField_Implicit::recal_eigenstrain_with_plasticity() {
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    PhaseNode& node = get_phaseNode_in_elasticField(i, j, k);
                    mech.EffectiveEigenStrains = node.customVec6s[ExternalFields::MECH_eigen_strain] + node.customVec6s[ExternalFields::MECH_plastic_strain];
#ifdef _DEBUG
                    if (mech.EffectiveEigenStrains.is_nan_val_exist()) {
                        cout << "DEBUG: mech.EffectiveEigenStrains error !" << endl;
                        SYS_PROGRAM_STOP;
                    }
#endif
                }
    }

    void MechanicalField_Implicit::initStrainIncrements() {
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    for (int n = 0; n < 6; n++) {
                        mech.StrainIncrements[n] = 0.0;
                    }
                }
    }

    string MechanicalField_Implicit::Solve(double StrainAccuracy, int MAXIterations, WriteToFile& writer, double incre_rate, bool is_dvStraindt_output, bool getU)
    {
        stringstream output;
        Matrix6x6 Cij;
        Cij = MAX_ElasticConstants * 1.1;

        vStress  AverageStress;
        vStress  TargetStress;

        AverageStress.set_to_zero();
        TargetStress.set_to_zero();

        vStrain   oldTargetStrain;
        oldTargetStrain.set_to_zero();
        average_strain.set_to_zero();

        int    IterationCount = 0;
        double MAXStrainDifference = 0.0;
        double MAXTargetStrainDifference = 0.0;

        do // Iteration loop begin
        {
#pragma omp parallel for
            for (int i = 0; i < rlSIZE; i++)
                for (int n = 0; n < 9; n++)
                    rlDefGrad[n][i] = 0.0;

            IterationCount++;

            MAXStrainDifference = 0.0;
            MAXTargetStrainDifference = 0.0;
            AverageStress.set_to_zero();

            CalculateRHS(Cij);
            ExecuteForwardFFT();
            CalculateFourierSolution(Cij, getU);
            ExecuteBackwardFFT(getU);
            SetElasticProperties1(MAXStrainDifference, AverageStress);

            TargetStress[0] = (-AverageStress[0] + applied_stress[0]) * incre_rate;
            TargetStress[1] = (-AverageStress[1] + applied_stress[1]) * incre_rate;
            TargetStress[2] = (-AverageStress[2] + applied_stress[2]) * incre_rate;

            oldTargetStrain = average_strain;

            average_strain.set_to_zero();

            for (int n = 0; n < 3; n++)
                for (int m = 0; m < 6; m++)
                {
                    average_strain[n] += TargetStress[m] * average_compliences(n, m);
                }

            SetElasticBoundaryConditions(TargetStress);

            SetElasticProperties2(MAXStrainDifference);

            for (int n = 0; n < 3; n++)
                if (fabs(average_strain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
                {
                    MAXTargetStrainDifference = fabs(average_strain[n] - oldTargetStrain[n]);
                }

            if (is_dvStraindt_output) {
                output.str("");
                output << "(Elastic solver) iterate step:                     " << IterationCount << endl
                    << "                 MAX dvstrain/dt:                  " << MAXStrainDifference << endl
                    << "                 MAX dTargetStrainDifference/dt:   " << MAXTargetStrainDifference << endl
                    << "                 Average Strain:                   " << "( " << average_strain[0] << ", "
                    << average_strain[1] << ", "
                    << average_strain[2] << ", "
                    << average_strain[3] << ", "
                    << average_strain[4] << ", "
                    << average_strain[5] << " )" << endl
                    << "                 Average Stress:                   " << "( " << AverageStress[0] << ", "
                    << AverageStress[1] << ", "
                    << AverageStress[2] << ", "
                    << AverageStress[3] << ", "
                    << AverageStress[4] << ", "
                    << AverageStress[5] << " )" << endl;
                writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
            }

            for (int n = 0; n < 3; n++)
                if (fabs(average_strain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
                {
                    MAXTargetStrainDifference = fabs(average_strain[n] - oldTargetStrain[n]);
                }
        } // Iteration loop end
        while ((MAXStrainDifference > StrainAccuracy ||
            MAXTargetStrainDifference > StrainAccuracy) && IterationCount < MAXIterations);

        assignment_to_phaseField();

        output.str("");
        output << "> Elastic Solver" << endl
            << "  iterate step:                     " << IterationCount << endl
            << "  MAX dvstrain/dt:                  " << MAXStrainDifference << endl
            << "  MAX dTargetStrainDifference/dt:   " << MAXTargetStrainDifference << endl
            << "  Average Strain:                   " << "( " << average_strain[0] << ", "
            << average_strain[1] << ", "
            << average_strain[2] << ", "
            << average_strain[3] << ", "
            << average_strain[4] << ", "
            << average_strain[5] << " )" << endl
            << "  Average Stress:                   " << "( " << AverageStress[0] << ", "
            << AverageStress[1] << ", "
            << AverageStress[2] << ", "
            << AverageStress[3] << ", "
            << AverageStress[4] << ", "
            << AverageStress[5] << " )" << endl;

        return output.str();
    }

    void MechanicalField_Implicit::CalculateRHS(Matrix6x6 Cij) {
        //Matrix3x3 unity;
        //unity.set_to_unity();
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    for (int n = 0; n < 6; n++) {
                        rlRHSide[n][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] = 0.0;
                        for (int m = 0; m < 6; m++) {
                            rlRHSide[n][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] +=
                                mech.EffectiveElasticConstants(n, m) * (mech.EffectiveEigenStrains[m]) + (Cij(n, m) - mech.EffectiveElasticConstants(n, m)) * (mech.StrainIncrements[m]);
                        }
                    }
                }
    }

    void MechanicalField_Implicit::CalculateFourierSolution(Matrix6x6 Cij, bool getU) {
        // XYZ from 0 to rcSIZE
#pragma omp parallel for
        for (int ite = 0; ite < rcSIZE; ite++)
        {
            double Qx = Q[0][ite], Qy = Q[1][ite], Qz = Q[2][ite];

            complex<double> rhsX = -I * (Qx * rcRHSide[0][ite] +
                Qy * rcRHSide[5][ite] +
                Qz * rcRHSide[4][ite]);
            complex<double> rhsY = -I * (Qx * rcRHSide[5][ite] +
                Qy * rcRHSide[1][ite] +
                Qz * rcRHSide[3][ite]);
            complex<double> rhsZ = -I * (Qx * rcRHSide[4][ite] +
                Qy * rcRHSide[3][ite] +
                Qz * rcRHSide[2][ite]);

            //< green function tensor inverse to the tensor (a)
            double a11 = (Cij(0, 0) * Qx * Qx + 2.0 * Cij(0, 5) * Qx * Qy + Cij(5, 5) * Qy * Qy +
                2.0 * Cij(0, 4) * Qx * Qz + 2.0 * Cij(4, 5) * Qy * Qz + Cij(4, 4) * Qz * Qz);

            double a21 = (Cij(0, 5) * Qx * Qx + Cij(0, 1) * Qx * Qy + Cij(5, 5) * Qx * Qy +
                Cij(1, 5) * Qy * Qy + Cij(0, 3) * Qx * Qz + Cij(4, 5) * Qx * Qz +
                Cij(1, 4) * Qy * Qz + Cij(3, 5) * Qy * Qz + Cij(3, 4) * Qz * Qz);

            double a31 = (Cij(0, 4) * Qx * Qx + Cij(0, 3) * Qx * Qy + Cij(4, 5) * Qx * Qy +
                Cij(3, 5) * Qy * Qy + Cij(0, 2) * Qx * Qz + Cij(4, 4) * Qx * Qz +
                Cij(2, 5) * Qy * Qz + Cij(3, 4) * Qy * Qz + Cij(2, 4) * Qz * Qz);

            double a12 = a21;

            double a22 = (Cij(5, 5) * Qx * Qx + 2.0 * Cij(1, 5) * Qx * Qy + Cij(1, 1) * Qy * Qy +
                2.0 * Cij(3, 5) * Qx * Qz + 2.0 * Cij(1, 3) * Qy * Qz + Cij(3, 3) * Qz * Qz);

            double a32 = (Cij(4, 5) * Qx * Qx + Cij(1, 4) * Qx * Qy + Cij(3, 5) * Qx * Qy +
                Cij(1, 3) * Qy * Qy + Cij(2, 5) * Qx * Qz + Cij(3, 4) * Qx * Qz +
                Cij(1, 2) * Qy * Qz + Cij(3, 3) * Qy * Qz + Cij(2, 3) * Qz * Qz);

            double a13 = a31;

            double a23 = a32;

            double a33 = (Cij(4, 4) * Qx * Qx + 2.0 * Cij(3, 4) * Qx * Qy + Cij(3, 3) * Qy * Qy +
                2.0 * Cij(2, 4) * Qx * Qz + 2.0 * Cij(2, 3) * Qy * Qz + Cij(2, 2) * Qz * Qz);

            // |a| = 
            double denominator = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32 -
                a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33);

            if (fabs(denominator) > DBL_EPSILON)
            {
                denominator = 1.0 / denominator;
            }
            else
            {
                denominator = 0.0;
            }
            // locUrc = a^-1 * rhs
            complex<double> locUrcX = (-a23 * a32 * rhsX + a22 * a33 * rhsX + a13 * a32 * rhsY -
                a12 * a33 * rhsY - a13 * a22 * rhsZ + a12 * a23 * rhsZ) * denominator * Norm;

            complex<double> locUrcY = (a23 * a31 * rhsX - a21 * a33 * rhsX - a13 * a31 * rhsY +
                a11 * a33 * rhsY + a13 * a21 * rhsZ - a11 * a23 * rhsZ) * denominator * Norm;

            complex<double> locUrcZ = (-a22 * a31 * rhsX + a21 * a32 * rhsX + a12 * a31 * rhsY -
                a11 * a32 * rhsY - a12 * a21 * rhsZ + a11 * a22 * rhsZ) * denominator * Norm;

            if (getU)  // Displacements in reciprocal space
            {
                rcU[0][ite] = locUrcX;
                rcU[1][ite] = locUrcY;
                rcU[2][ite] = locUrcZ;
            }

            //  Deformation gradient entries in Fourier space
            rcDefGrad[0][ite] = I * (Qx * locUrcX);
            rcDefGrad[1][ite] = I * (Qy * locUrcX);
            rcDefGrad[2][ite] = I * (Qz * locUrcX);
            rcDefGrad[3][ite] = I * (Qx * locUrcY);
            rcDefGrad[4][ite] = I * (Qy * locUrcY);
            rcDefGrad[5][ite] = I * (Qz * locUrcY);
            rcDefGrad[6][ite] = I * (Qx * locUrcZ);
            rcDefGrad[7][ite] = I * (Qy * locUrcZ);
            rcDefGrad[8][ite] = I * (Qz * locUrcZ);
        }
    }

    void MechanicalField_Implicit::SetElasticProperties1(double& MAXStrainDifference, vStress& AverageStress) {
        Matrix3x3 unity;
        unity.set_to_unity();
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);

                    Matrix3x3 locDefGrad;
                    locDefGrad(0, 0) = rlDefGrad[0][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    locDefGrad(1, 1) = rlDefGrad[4][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    locDefGrad(2, 2) = rlDefGrad[8][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;

                    locDefGrad(0, 1) = rlDefGrad[1][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(0, 2) = rlDefGrad[2][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(1, 2) = rlDefGrad[5][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    locDefGrad(1, 0) = rlDefGrad[3][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(2, 0) = rlDefGrad[6][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(2, 1) = rlDefGrad[7][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    Matrix3x3 locStrainTensor = (locDefGrad.get_transposed() * locDefGrad - unity) * 0.5;
                    vStrain locStrain = locStrainTensor.VoigtStrain();

                    for (int n = 0; n < 6; n++)
                    {
                        double locStrainDifference = fabs(mech.StrainIncrements[n] - locStrain[n] - average_strain[n]);  //difference for locStrain
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (locStrainDifference > MAXStrainDifference)
                            {
                                MAXStrainDifference = locStrainDifference;
                            }
                        }
                    }
                    for (int n = 0; n < 6; n++)
                    {
                        double locStress = 0.0;
                        for (int m = 0; m < 6; m++)
                        {
                            locStress += mech.EffectiveElasticConstants(n, m) *
                                (locStrain[m] - mech.EffectiveEigenStrains[m]/* - mech.PlasticStrains[m]*/);// + EP.RemeshedStrain[m]);
                        }
#ifdef _OPENMP
#pragma omp atomic
#endif
                        AverageStress[n] += locStress;
                    }
                }
        for (int n = 0; n < 6; n++)
        {
            AverageStress[n] *= Norm;
        }
    }

    void MechanicalField_Implicit::SetElasticProperties2(double& MAXStrainDifference) {
        Matrix3x3 unity;
        unity.set_to_unity();
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);

                    Matrix3x3 locDefGrad;
                    locDefGrad(0, 0) = rlDefGrad[0][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    locDefGrad(1, 1) = rlDefGrad[4][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    locDefGrad(2, 2) = rlDefGrad[8][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;

                    locDefGrad(0, 1) = rlDefGrad[1][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(0, 2) = rlDefGrad[2][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(1, 2) = rlDefGrad[5][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    locDefGrad(1, 0) = rlDefGrad[3][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(2, 0) = rlDefGrad[6][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    locDefGrad(2, 1) = rlDefGrad[7][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    Matrix3x3 locStrainTensor = (locDefGrad.get_transposed() * locDefGrad - unity) * 0.5;
                    //dMatrix3x3 locStrainTensor = (locDefGrad.transposed() + locDefGrad)*0.5 - unity;
                    vStrain locStrain = locStrainTensor.VoigtStrain();

                    for (int n = 0; n < 6; n++)
                    {
                        locStrain[n] += average_strain[n];

                        double locStrainDifference = fabs(mech.StrainIncrements[n] - locStrain[n]);   // difference for locStrain and average_strain

#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (locStrainDifference > MAXStrainDifference)
                            {
                                MAXStrainDifference = locStrainDifference;
                            }
                        }
                        mech.StrainIncrements[n] = locStrain[n];
                        // if (!EP.LargeDeformations)
                        mech.Strains[n] = locStrain[n];
                    }
                    for (int n = 0; n < 6; n++)
                    {
                        double locStress = 0.0;
                        for (int m = 0; m < 6; m++)
                        {
                            locStress += mech.EffectiveElasticConstants(n, m) *
                                (locStrain[m] - mech.EffectiveEigenStrains[m]/* - mech.PlasticStrains[m]*/);// + EP.RemeshedStrain[m]);
                        }
                        //double locStressDifference = fabs(mech.Stresses[n] - locStress);
                        mech.Stresses[n] = locStress;

                    }
                }
    }

    void MechanicalField_Implicit::SetElasticBoundaryConditions(vStress TargetStress) {
        if (average_stiffness(0, 0) == 0 or average_stiffness(1, 1) == 0 or average_stiffness(2, 2) == 0
            or average_stiffness(1, 2) == 0 or average_stiffness(0, 2) == 0 or average_stiffness(0, 1) == 0)
        {
            cout << "> ERROR ! In SetElasticBoundaryConditions(), Zero component in AverageElasticConstants." << endl;
            SYS_PROGRAM_STOP;
        }

        if (AppStrainMask[0] and !AppStrainMask[1] and !AppStrainMask[2])   /// XX
        {
            average_strain[0] = applied_strain[0];
            average_strain[1] = (-average_stiffness(0, 2) * average_stiffness(1, 2) * applied_strain[0] +
                average_stiffness(0, 1) * average_stiffness(2, 2) * applied_strain[0] -
                average_stiffness(2, 2) * TargetStress[1] + average_stiffness(1, 2) * TargetStress[2]) /
                (average_stiffness(1, 2) * average_stiffness(1, 2) - average_stiffness(1, 1) * average_stiffness(2, 2));
            average_strain[2] = (-average_stiffness(1, 1) * average_strain[1] - average_stiffness(0, 1) * applied_strain[0] + TargetStress[1]) /
                average_stiffness(1, 2);
        }

        if (!AppStrainMask[0] and AppStrainMask[1] and !AppStrainMask[2])   /// YY
        {
            average_strain[0] = (-average_stiffness(0, 2) * average_stiffness(1, 2) * applied_strain[1] +
                average_stiffness(0, 1) * average_stiffness(2, 2) * applied_strain[1] -
                average_stiffness(2, 2) * TargetStress[0] + average_stiffness(0, 2) * TargetStress[2]) /
                (average_stiffness(0, 2) * average_stiffness(0, 2) - average_stiffness(0, 0) * average_stiffness(2, 2));
            average_strain[1] = applied_strain[1];
            average_strain[2] = (-average_stiffness(0, 0) * average_strain[0] - average_stiffness(0, 1) * applied_strain[1] + TargetStress[0]) /
                average_stiffness(0, 2);
        }

        if (!AppStrainMask[0] and !AppStrainMask[1] and AppStrainMask[2])   /// ZZ
        {
            average_strain[0] = (average_stiffness(0, 2) * average_stiffness(1, 1) * applied_strain[2] -
                average_stiffness(0, 1) * average_stiffness(1, 2) * applied_strain[2] -
                average_stiffness(1, 1) * TargetStress[0] + average_stiffness(0, 1) * TargetStress[1]) /
                (average_stiffness(0, 1) * average_stiffness(0, 1) - average_stiffness(0, 0) * average_stiffness(1, 1));

            average_strain[1] = (-average_stiffness(0, 0) * average_strain[0] - average_stiffness(0, 2) * applied_strain[2] + TargetStress[0]) /
                average_stiffness(0, 1);
            average_strain[2] = applied_strain[2];
        }

        if (AppStrainMask[0] and AppStrainMask[1] and !AppStrainMask[2])    /// XX & YY
        {
            average_strain[0] = applied_strain[0];
            average_strain[1] = applied_strain[1];
            average_strain[2] = (-average_stiffness(0, 2) * applied_strain[0] - average_stiffness(1, 2) * applied_strain[1] + TargetStress[2]) /
                average_stiffness(2, 2);
        }

        if (AppStrainMask[0] and !AppStrainMask[1] and AppStrainMask[2])    /// XX & ZZ
        {
            average_strain[0] = applied_strain[0];
            average_strain[1] = (-average_stiffness(0, 1) * applied_strain[0] - average_stiffness(1, 2) * applied_strain[2] + TargetStress[1]) /
                average_stiffness(1, 1);
            average_strain[2] = applied_strain[2];
        }

        if (!AppStrainMask[0] and AppStrainMask[1] and AppStrainMask[2])    /// YY & ZZ
        {
            average_strain[0] = (-average_stiffness(0, 1) * applied_strain[1] - average_stiffness(0, 2) * applied_strain[2] + TargetStress[0]) /
                average_stiffness(0, 0);
            average_strain[1] = applied_strain[1];
            average_strain[2] = applied_strain[2];
        }

        if (AppStrainMask[0] and AppStrainMask[1] and AppStrainMask[2])    /// XX & YY & ZZ
        {
            average_strain[0] = applied_strain[0];
            average_strain[1] = applied_strain[1];
            average_strain[2] = applied_strain[2];
        }
    }

    void MechanicalField_Implicit::ExecuteForwardFFT() {
#pragma omp parallel sections// OMP BEGIN
        {
#pragma omp section
            {
                fftw_execute(ForwardPlanRHS[0]);
            }
#pragma omp section
            {
                fftw_execute(ForwardPlanRHS[1]);
            }
#pragma omp section
            {
                fftw_execute(ForwardPlanRHS[2]);
            }
#pragma omp section
            {
                fftw_execute(ForwardPlanRHS[3]);
            }
#pragma omp section
            {
                fftw_execute(ForwardPlanRHS[4]);
            }
#pragma omp section
            {
                fftw_execute(ForwardPlanRHS[5]);
            }
        }//OMP END
    }

    void MechanicalField_Implicit::ExecuteBackwardFFT(bool getU) {
#pragma omp parallel sections // OMP BEGIN
        {
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[0]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[1]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[2]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[3]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[4]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[5]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[6]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[7]);
            }
#pragma omp section
            {
                fftw_execute(BackwardPlanDefGrad[8]);
            }
#pragma omp section
            {
                if (getU) fftw_execute(BackwardPlanU[0]);
            }
#pragma omp section
            {
                if (getU) fftw_execute(BackwardPlanU[1]);
            }
#pragma omp section
            {
                if (getU) fftw_execute(BackwardPlanU[2]);
            }
        }
    }

    void MechanicalField_Implicit::initVirtualEigenstrain() {
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    mech.VirtualEigenStrains.set_to_zero();
                    //mech.VirtualEigenStrains = mech.EffectiveEigenStrains;
                }
    }

    string MechanicalField_Implicit::Solve2(double StrainAccuracy, int MAXIterations, double iterate_rate, WriteToFile& writer, bool is_dvStraindt_output, bool getU) {
        stringstream output;

        Matrix6x6 Cij;                                                             ///< assuming homogeneous elastic constants values
        Cij = MAX_ElasticConstants * 1.1;
        Matrix6x6 Sij;
        Sij = Cij.get_inverted_matrix();

        double MAXStrainDifference = 0.0;
        int IterationCount = 0;

        for (IterationCount = 1; IterationCount <= MAXIterations; IterationCount++) {
            MAXStrainDifference = 0.0;

            CalculateRHS2(Cij);
            ExecuteForwardFFT();
            CalculateFourierSolution(Cij, getU);
            ExecuteBackwardFFT(getU);

            SetElasticBoundaryConditions2(Cij, Sij);

            evaluate_virtualEigenstrain(Cij, Sij, MAXStrainDifference, iterate_rate);

            if (is_dvStraindt_output) {
                output.str("");
                output << "(Elastic solver) iterate step:                     " << IterationCount << endl
                    << "                 MAX dVirtualStrain/dt:                  " << MAXStrainDifference << endl
                    << "                 Average Strain:                   " << "( " << average_strain[0] << ", "
                    << average_strain[1] << ", "
                    << average_strain[2] << ", "
                    << average_strain[3] << ", "
                    << average_strain[4] << ", "
                    << average_strain[5] << " )" << endl;
                writer.add_string_to_txt_and_screen(output.str(), LOG_FILE_NAME);
            }

            if (MAXStrainDifference < StrainAccuracy)
                break;
        }

        assignment_to_phaseField();

        output.str("");
        output << "> Elastic Solver" << endl
            << "  iterate step:                     " << IterationCount << endl
            << "  MAX dVirtualStrain/dt:                  " << MAXStrainDifference << endl
            << "  Average Strain:                   " << "( " << average_strain[0] << ", "
            << average_strain[1] << ", "
            << average_strain[2] << ", "
            << average_strain[3] << ", "
            << average_strain[4] << ", "
            << average_strain[5] << " )" << endl;

        return output.str();
    }

    void MechanicalField_Implicit::CalculateRHS2(Matrix6x6& Cij) {
        average_virtual_strain.set_to_zero();
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    for (int n = 0; n < 6; n++) {
                        rlRHSide[n][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] = 0.0;
                        for (int m = 0; m < 6; m++) {
                            rlRHSide[n][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] +=
                                Cij(n, m) * mech.VirtualEigenStrains[m];
                        }
                    }
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        average_virtual_strain += mech.VirtualEigenStrains;
                    }
                }
        for (int n = 0; n < 6; n++)
        {
            average_virtual_strain[n] *= Norm;
        }
    }

    void MechanicalField_Implicit::evaluate_virtualEigenstrain(Matrix6x6 Cij, Matrix6x6 Sij, double& MAXvStrainDifference, double iterate_rate) {
        int VoigtIndex[6][2] = { {0,0},{4,4},{8,8},{5,7},{2,6},{1,3} };
        double VoigtFactor[6] = { 0.5, 0.5, 0.5, 1.0, 1.0, 1.0 };
        average_strain.set_to_zero();
#pragma omp parallel for
        for (int i = 0; i < mechanicalField.Nx; i++)
            for (int j = 0; j < mechanicalField.Ny; j++)
                for (int k = 0; k < mechanicalField.Nz; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    vStrain strain1, strain2, strain3, vStrain_increment;
                    Matrix6x6 deltS;
                    deltS = Cij - mech.EffectiveElasticConstants;
                    deltS.do_invert();
                    for (int n = 0; n < 6; n++)
                    {
                        strain1[n] = VoigtFactor[n] * (rlDefGrad[VoigtIndex[n][0]][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] +
                            rlDefGrad[VoigtIndex[n][1]][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)]);
                        mech.Strains[n] = applied_strain[n] + strain1[n];
                    }
                    strain2 = deltS * (Cij * (mech.VirtualEigenStrains - mech.EffectiveEigenStrains));
                    strain3 = Sij * applied_stress;
                    for (int n = 0; n < 6; n++)
                    {
                        vStrain_increment[n] = 0.0;
                        mech.Stresses[n] = 0.0;
                        for (int m = 0; m < 6; m++) {
                            vStrain_increment[n] += iterate_rate * Cij(n, m) * (strain1[m] - strain2[m] - mech.EffectiveEigenStrains[m] + average_virtual_strain[m] + strain3[m]);
                            mech.Stresses[n] += mech.EffectiveElasticConstants(n, m) * (mech.Strains[m] - mech.EffectiveEigenStrains[m]); // important
                        }
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (fabs(vStrain_increment[n]) > MAXvStrainDifference)
                            {
                                MAXvStrainDifference = fabs(vStrain_increment[n]);
                            }
                            average_strain[n] += mech.Strains[n];
                        }
                        mech.VirtualEigenStrains[n] += vStrain_increment[n];
                    }
                }
        average_strain *= Norm;
    }

    void MechanicalField_Implicit::SetElasticBoundaryConditions2(Matrix6x6& Cij, Matrix6x6& Sij) {
        pf::vStress exStress;
        exStress = Cij * (applied_strain - average_virtual_strain);
        if (AppStrainMask[0]) {
            applied_stress[0] = exStress[0];
        }
        if (AppStrainMask[1]) {
            applied_stress[1] = exStress[1];
        }
        if (AppStrainMask[2]) {
            applied_stress[2] = exStress[2];
        }
        pf::vStrain exStrain;
        exStrain = Sij * applied_stress + average_virtual_strain;
        if (LoadStressMask[0]) {
            applied_strain[0] = exStrain[0];
        }
        if (LoadStressMask[1]) {
            applied_strain[1] = exStrain[1];
        }
        if (LoadStressMask[2]) {
            applied_strain[2] = exStrain[2];
        }
    }

    double MechanicalField_Implicit::get_u_main_node(int _x, int _y, int _z) {
        return rlU[0][_z + mechanicalField.Nz * (_y + mechanicalField.Ny * _x)];
    }
    double MechanicalField_Implicit::get_v_main_node(int _x, int _y, int _z) {
        return rlU[1][_z + mechanicalField.Nz * (_y + mechanicalField.Ny * _x)];
    }
    double MechanicalField_Implicit::get_w_main_node(int _x, int _y, int _z) {
        return rlU[2][_z + mechanicalField.Nz * (_y + mechanicalField.Ny * _x)];
    }

    void MechanicalField_Implicit::assignment_to_phaseField() {
#pragma omp parallel for
        for (int i = 0; i < phaseMesh->limit_x; i++)
            for (int j = 0; j < phaseMesh->limit_y; j++)
                for (int k = 0; k < phaseMesh->limit_z; k++) {
                    MechanicNode& mech = mechanicalField(i, j, k);
                    PhaseNode& node = (*phaseMesh)(i, j, k);
                    node.customVec6s[ExternalFields::MECH_strain] = mech.Strains;
                    node.customVec6s[ExternalFields::MECH_stress] = mech.Stresses;
                }
    }

}