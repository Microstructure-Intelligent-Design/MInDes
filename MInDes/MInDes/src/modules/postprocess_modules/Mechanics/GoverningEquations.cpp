#include"GoverningEquations.h"
namespace pf {

    void MechanicalField_Implicit::SetMAXElasticConstants(std::vector<Matrix6x6> Cijs) {
        MAX_ElasticConstants.set_to_zero();
        for (auto Cij = Cijs.begin(); Cij < Cijs.end(); Cij++) {
            for (int n = 0; n < 6; n++)
                for (int m = 0; m < 6; m++)
                    if ((*Cij)(n, m) > MAX_ElasticConstants(n, m))
                        MAX_ElasticConstants(n, m) = (*Cij)(n, m);
        }
    }

    void MechanicalField_Implicit::cal_eigenstrain_stiffness() {
        average_stiffness.set_to_zero();
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    int mirror_x = i, mirror_y = j, mirror_z = k;
                    if (i >= Nx)
                        mirror_x = mech_Nx - 1 - i;
                    if (j >= Ny)
                        mirror_y = mech_Ny - 1 - j;
                    if (k >= Nz)
                        mirror_z = mech_Nz - 1 - k;
                    mech.EffectiveEigenStrain = cal_eigenstrain(mirror_x, mirror_y, mirror_z);
                    mech.EffectiveElasticConstant = cal_stiffness(mirror_x, mirror_y, mirror_z);
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        average_stiffness += mech.EffectiveElasticConstant;
                    }
#ifdef _DEBUG
                    if (mech.EffectiveEigenStrain.is_nan_val_exist()) {
                        std::cout << "DEBUG: mech.EffectiveEigenStrain error !" << std::endl;
                        SYS_PROGRAM_STOP;
                    }
#endif
                }
        average_stiffness /= REAL(mech_Nx * mech_Ny * mech_Nz);
        average_compliences = average_stiffness.get_inverted_matrix();
    }

    void MechanicalField_Implicit::recal_eigenstrain() {
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    int mirror_x = i, mirror_y = j, mirror_z = k;
                    if (i >= Nx)
                        mirror_x = mech_Nx - 1 - i;
                    if (j >= Ny)
                        mirror_y = mech_Ny - 1 - j;
                    if (k >= Nz)
                        mirror_z = mech_Nz - 1 - k;
                    mech.EffectiveEigenStrain = cal_eigenstrain(mirror_x, mirror_y, mirror_z);
#ifdef _DEBUG
                    if (mech.EffectiveEigenStrain.is_nan_val_exist()) {
                        std::cout << "DEBUG: mech.EffectiveEigenStrain error !" << std::endl;
                        SYS_PROGRAM_STOP;
                    }
#endif
                }
    }

    void MechanicalField_Implicit::initStrainIncrements() {
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    for (int n = 0; n < 6; n++) {
                        mech.StrainIncrement[n] = 0.0;
                    }
                }
    }

    std::string MechanicalField_Implicit::Solve(double StrainAccuracy, int MAXIterations, double incre_rate, bool is_dvStraindt_output, bool getU)
    {
        std::stringstream output;
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

#ifdef _DEBUG
            if (average_strain.is_nan_val_exist()) {
                std::cout << "DEBUG: average_strain error, nan value exists !" << std::endl;
                SYS_PROGRAM_STOP;
            }
#endif

            SetElasticBoundaryConditions(TargetStress);

#ifdef _DEBUG
            if (average_strain.is_nan_val_exist()) {
                std::cout << "DEBUG: average_strain error, nan value exists !" << std::endl;
                SYS_PROGRAM_STOP;
            }
            if (TargetStress.is_nan_val_exist()) {
                std::cout << "DEBUG: TargetStress error, nan value exists !" << std::endl;
                SYS_PROGRAM_STOP;
            }
#endif

            SetElasticProperties2(MAXStrainDifference);

            for (int n = 0; n < 3; n++)
                if (fabs(average_strain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
                {
                    MAXTargetStrainDifference = fabs(average_strain[n] - oldTargetStrain[n]);
                }

            if (is_dvStraindt_output) {
                output.str("");
                output << "(Elastic solver) iterate step:                     " << IterationCount << std::endl
                    << "                 MAX dvstrain/dt:                  " << MAXStrainDifference << std::endl
                    << "                 MAX dTargetStrainDifference/dt:   " << MAXTargetStrainDifference << std::endl
                    << "                 Average Strain:                   " << "( " << average_strain[0] << ", "
                    << average_strain[1] << ", "
                    << average_strain[2] << ", "
                    << average_strain[3] << ", "
                    << average_strain[4] << ", "
                    << average_strain[5] << " )" << std::endl
                    << "                 Average Stress:                   " << "( " << AverageStress[0] << ", "
                    << AverageStress[1] << ", "
                    << AverageStress[2] << ", "
                    << AverageStress[3] << ", "
                    << AverageStress[4] << ", "
                    << AverageStress[5] << " )" << std::endl;
                WriteLog(output.str());
            }

            for (int n = 0; n < 3; n++)
                if (fabs(average_strain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
                {
                    MAXTargetStrainDifference = fabs(average_strain[n] - oldTargetStrain[n]);
                }
        } // Iteration loop end
        while ((MAXStrainDifference > StrainAccuracy ||
            MAXTargetStrainDifference > StrainAccuracy) && IterationCount < MAXIterations);

        output.str("");
        output << "> Elastic Solver" << std::endl
            << "  iterate step:                     " << IterationCount << std::endl
            << "  MAX dvstrain/dt:                  " << MAXStrainDifference << std::endl
            << "  MAX dTargetStrainDifference/dt:   " << MAXTargetStrainDifference << std::endl
            << "  Average Strain:                   " << "( " << average_strain[0] << ", "
            << average_strain[1] << ", "
            << average_strain[2] << ", "
            << average_strain[3] << ", "
            << average_strain[4] << ", "
            << average_strain[5] << " )" << std::endl
            << "  Average Stress:                   " << "( " << AverageStress[0] << ", "
            << AverageStress[1] << ", "
            << AverageStress[2] << ", "
            << AverageStress[3] << ", "
            << AverageStress[4] << ", "
            << AverageStress[5] << " )" << std::endl;

        return output.str();
    }

    void MechanicalField_Implicit::CalculateRHS(Matrix6x6 Cij) {
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    for (int n = 0; n < 6; n++) {
                        rlRHSide[n][k + mech_Nz * (j + mech_Ny * i)] = 0.0;
                        for (int m = 0; m < 6; m++) {
                            rlRHSide[n][k + mech_Nz * (j + mech_Ny * i)] +=
                                mech.EffectiveElasticConstant(n, m) * (mech.EffectiveEigenStrain[m]) + (Cij(n, m) - mech.EffectiveElasticConstant(n, m)) * (mech.StrainIncrement[m]);
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

            std::complex<double> rhsX = -I * (Qx * rcRHSide[0][ite] +
                Qy * rcRHSide[5][ite] +
                Qz * rcRHSide[4][ite]);
            std::complex<double> rhsY = -I * (Qx * rcRHSide[5][ite] +
                Qy * rcRHSide[1][ite] +
                Qz * rcRHSide[3][ite]);
            std::complex<double> rhsZ = -I * (Qx * rcRHSide[4][ite] +
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
            std::complex<double> locUrcX = (-a23 * a32 * rhsX + a22 * a33 * rhsX + a13 * a32 * rhsY -
                a12 * a33 * rhsY - a13 * a22 * rhsZ + a12 * a23 * rhsZ) * denominator * Norm;

            std::complex<double> locUrcY = (a23 * a31 * rhsX - a21 * a33 * rhsX - a13 * a31 * rhsY +
                a11 * a33 * rhsY + a13 * a21 * rhsZ - a11 * a23 * rhsZ) * denominator * Norm;

            std::complex<double> locUrcZ = (-a22 * a31 * rhsX + a21 * a32 * rhsX + a12 * a31 * rhsY -
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
        int VoigtIndex[6][2] = { {0,0},{4,4},{8,8},{5,7},{2,6},{1,3} };
        double VoigtFactor[6] = { 0.5, 0.5, 0.5, 1.0, 1.0, 1.0 };
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    // - OpenPhase.3999.Mar2018
                    vStrain locStrain;
                    for (int n = 0; n < 6; n++)
                    {
                        locStrain[n] = VoigtFactor[n] * (rlDefGrad[VoigtIndex[n][0]][k + mech_Nz * (j + mech_Ny * i)] +
                            rlDefGrad[VoigtIndex[n][1]][k + mech_Nz * (j + mech_Ny * i)]);
                    }
                    // - OpenPhase.V0.9.2
                    //Matrix3x3 locDefGrad;
                    //locDefGrad(0, 0) = rlDefGrad[0][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    //locDefGrad(1, 1) = rlDefGrad[4][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    //locDefGrad(2, 2) = rlDefGrad[8][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;

                    //locDefGrad(0, 1) = rlDefGrad[1][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(0, 2) = rlDefGrad[2][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(1, 2) = rlDefGrad[5][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    //locDefGrad(1, 0) = rlDefGrad[3][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(2, 0) = rlDefGrad[6][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(2, 1) = rlDefGrad[7][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    //Matrix3x3 locStrainTensor = (locDefGrad.get_transposed() * locDefGrad - unity) * 0.5;
                    //// Matrix3x3 locStrainTensor = (locDefGrad.do_transpose() + locDefGrad) * 0.5 - unity;
                    //vStrain locStrain = locStrainTensor.VoigtStrain();

                    for (int n = 0; n < 6; n++)
                    {
                        double locStrainDifference = fabs(mech.StrainIncrement[n] - locStrain[n] - average_strain[n]);  //difference for locStrain
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
                            locStress += mech.EffectiveElasticConstant(n, m) *
                                (locStrain[m] - mech.EffectiveEigenStrain[m]/* - mech.PlasticStrains[m]*/);// + EP.RemeshedStrain[m]);
                        }
#ifdef _OPENMP
#pragma omp atomic
#endif
                        AverageStress[n] += locStress;
                    }
#ifdef _DEBUG
                    if (AverageStress.is_nan_val_exist()) {
                        std::cout << "DEBUG: AverageStress error, nan value exists !" << std::endl;
                        SYS_PROGRAM_STOP;
                    }
#endif
                }
        for (int n = 0; n < 6; n++)
        {
            AverageStress[n] *= Norm;
        }
    }

    void MechanicalField_Implicit::SetElasticProperties2(double& MAXStrainDifference) {
        Matrix3x3 unity;
        unity.set_to_unity();
        int VoigtIndex[6][2] = { {0,0},{4,4},{8,8},{5,7},{2,6},{1,3} };
        double VoigtFactor[6] = { 0.5, 0.5, 0.5, 1.0, 1.0, 1.0 };
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    // - OpenPhase.3999.Mar2018
                    vStrain locStrain;
                    for (int n = 0; n < 6; n++)
                    {
                        locStrain[n] = VoigtFactor[n] * (rlDefGrad[VoigtIndex[n][0]][k + mech_Nz * (j + mech_Ny * i)] +
                            rlDefGrad[VoigtIndex[n][1]][k + mech_Nz * (j + mech_Ny * i)]);
                    }
                    // - OpenPhase.V0.9.2
                    //Matrix3x3 locDefGrad;
                    //locDefGrad(0, 0) = rlDefGrad[0][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    //locDefGrad(1, 1) = rlDefGrad[4][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;
                    //locDefGrad(2, 2) = rlDefGrad[8][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)] + 1.0;

                    //locDefGrad(0, 1) = rlDefGrad[1][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(0, 2) = rlDefGrad[2][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(1, 2) = rlDefGrad[5][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    //locDefGrad(1, 0) = rlDefGrad[3][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(2, 0) = rlDefGrad[6][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];
                    //locDefGrad(2, 1) = rlDefGrad[7][k + mechanicalField.Nz * (j + mechanicalField.Ny * i)];

                    //Matrix3x3 locStrainTensor = (locDefGrad.get_transposed() * locDefGrad - unity) * 0.5;
                    //// Matrix3x3 locStrainTensor = (locDefGrad.do_transpose() + locDefGrad)*0.5 - unity;
                    //vStrain locStrain = locStrainTensor.VoigtStrain();

                    for (int n = 0; n < 6; n++)
                    {
                        locStrain[n] += average_strain[n];

                        double locStrainDifference = fabs(mech.StrainIncrement[n] - locStrain[n]);   // difference for locStrain and average_strain

#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (locStrainDifference > MAXStrainDifference)
                            {
                                MAXStrainDifference = locStrainDifference;
                            }
                        }
                        mech.StrainIncrement[n] = locStrain[n];
                        // if (!EP.LargeDeformations)
                        mech.Strain[n] = locStrain[n];
                    }
                    for (int n = 0; n < 6; n++)
                    {
                        double locStress = 0.0;
                        for (int m = 0; m < 6; m++)
                        {
                            locStress += mech.EffectiveElasticConstant(n, m) *
                                (locStrain[m] - mech.EffectiveEigenStrain[m]/* - mech.PlasticStrains[m]*/);// + EP.RemeshedStrain[m]);
                        }
                        //double locStressDifference = fabs(mech.Stresses[n] - locStress);
                        mech.Stress[n] = locStress;

                    }
#ifdef _DEBUG
                    if (locStrain.is_nan_val_exist()) {
                        std::cout << "DEBUG: locStrain error, nan value exists !" << std::endl;
                        SYS_PROGRAM_STOP;
                    }
#endif
                }
    }

    void MechanicalField_Implicit::SetElasticBoundaryConditions(vStress TargetStress) {
        if (average_stiffness(0, 0) == 0 or average_stiffness(1, 1) == 0 or average_stiffness(2, 2) == 0
            or average_stiffness(1, 2) == 0 or average_stiffness(0, 2) == 0 or average_stiffness(0, 1) == 0)
        {
            std::cout << "> ERROR ! In SetElasticBoundaryConditions(), Zero component in AverageElasticConstants." << std::endl;
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
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    mech.VirtualEigenStrain.set_to_zero();
                    //mech.VirtualEigenStrain = mech.EffectiveEigenStrain;
                }
    }

    std::string MechanicalField_Implicit::Solve2(double StrainAccuracy, int MAXIterations, double iterate_rate, bool is_dvStraindt_output, bool getU) {
        std::stringstream output;

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
                output << "(Elastic solver) iterate step:                     " << IterationCount << std::endl
                    << "                 MAX dVirtualStrain/dt:                  " << MAXStrainDifference << std::endl
                    << "                 Average Strain:                   " << "( " << average_strain[0] << ", "
                    << average_strain[1] << ", "
                    << average_strain[2] << ", "
                    << average_strain[3] << ", "
                    << average_strain[4] << ", "
                    << average_strain[5] << " )" << std::endl;
                WriteLog(output.str());
            }

            if (MAXStrainDifference < StrainAccuracy)
                break;
        }

        output.str("");
        output << "> Elastic Solver" << std::endl
            << "  iterate step:                     " << IterationCount << std::endl
            << "  MAX dVirtualStrain/dt:                  " << MAXStrainDifference << std::endl
            << "  Average Strain:                   " << "( " << average_strain[0] << ", "
            << average_strain[1] << ", "
            << average_strain[2] << ", "
            << average_strain[3] << ", "
            << average_strain[4] << ", "
            << average_strain[5] << " )" << std::endl;

        return output.str();
    }

    void MechanicalField_Implicit::CalculateRHS2(Matrix6x6& Cij) {
        average_virtual_strain.set_to_zero();
#pragma omp parallel for
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    for (int n = 0; n < 6; n++) {
                        rlRHSide[n][k + mech_Nz * (j + mech_Ny * i)] = 0.0;
                        for (int m = 0; m < 6; m++) {
                            rlRHSide[n][k + mech_Nz * (j + mech_Ny * i)] +=
                                Cij(n, m) * mech.VirtualEigenStrain[m];
                        }
                    }
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        average_virtual_strain += mech.VirtualEigenStrain;
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
        for (int i = 0; i < mech_Nx; i++)
            for (int j = 0; j < mech_Ny; j++)
                for (int k = 0; k < mech_Nz; k++) {
                    ElasticPoint& mech = elastic_field->at(i, j, k);
                    vStrain strain1, strain2, strain3, vStrain_increment;
                    Matrix6x6 deltS;
                    deltS = Cij - mech.EffectiveElasticConstant;
                    deltS.do_invert();
                    for (int n = 0; n < 6; n++)
                    {
                        strain1[n] = VoigtFactor[n] * (rlDefGrad[VoigtIndex[n][0]][k + mech_Nz * (j + mech_Ny * i)] +
                            rlDefGrad[VoigtIndex[n][1]][k + mech_Nz * (j + mech_Ny * i)]);
                        mech.Strain[n] = applied_strain[n] + strain1[n];
                    }
                    strain2 = deltS * (Cij * (mech.VirtualEigenStrain - mech.EffectiveEigenStrain));
                    strain3 = Sij * applied_stress;
                    for (int n = 0; n < 6; n++)
                    {
                        vStrain_increment[n] = 0.0;
                        mech.Stress[n] = 0.0;
                        for (int m = 0; m < 6; m++) {
                            vStrain_increment[n] += iterate_rate * Cij(n, m) * (strain1[m] - strain2[m] - mech.EffectiveEigenStrain[m] + average_virtual_strain[m] + strain3[m]);
                            mech.Stress[n] += mech.EffectiveElasticConstant(n, m) * (mech.Strain[m] - mech.EffectiveEigenStrain[m]); // important
                        }
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (fabs(vStrain_increment[n]) > MAXvStrainDifference)
                            {
                                MAXvStrainDifference = fabs(vStrain_increment[n]);
                            }
                            average_strain[n] += mech.Strain[n];
                        }
                        mech.VirtualEigenStrain[n] += vStrain_increment[n];
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

    double MechanicalField_Implicit::get_u_main_node(size_t _x, size_t _y, size_t _z) {
        return rlU[0][_z + mech_Nz * (_y + mech_Ny * _x)];
    }
    double MechanicalField_Implicit::get_v_main_node(size_t _x, size_t _y, size_t _z) {
        return rlU[1][_z + mech_Nz * (_y + mech_Ny * _x)];
    }
    double MechanicalField_Implicit::get_w_main_node(size_t _x, size_t _y, size_t _z) {
        return rlU[2][_z + mech_Nz * (_y + mech_Ny * _x)];
    }
}