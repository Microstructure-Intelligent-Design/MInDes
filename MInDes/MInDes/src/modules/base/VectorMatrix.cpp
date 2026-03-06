#include "VectorMatrix.h"
namespace pf {

    Matrix3x3::Matrix3x3() {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                storage[i][j] = 0.0;
    };

    Matrix3x3::Matrix3x3(const Matrix3x3& rhs)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                storage[i][j] = rhs(i, j);
    };

    REAL& Matrix3x3::operator()(const int i, const int j)
    {
#ifdef DEBUG
        if (i > 2 or j > 2 or i < 0 or j < 0)
        {
            std::cout << "Error in dMatrix3x3::operator()\n"
                << "Access beyond storage range. (i, j) = "
                << i << ", " << j << " > (2, 2)"
                << "\nTerminating!!!" << std::endl;
            SYS_PROGRAM_STOP;
        }
#endif
        return storage[i][j];
    };
    REAL const& Matrix3x3::operator()(const int i, const int j) const
    {
#ifdef DEBUG
        if (i > 2 or j > 2 or i < 0 or j < 0)
        {
            std::cout << "Error in dMatrix3x3::operator()\n"
                << "Access beyond storage range. (i, j) = "
                << i << ", " << j << " > (2, 2)"
                << "\nTerminating!!!" << std::endl;
            SYS_PROGRAM_STOP;
        }
#endif
        return storage[i][j];
    };
    void Matrix3x3::set_to_zero(void)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                storage[i][j] = 0.0;
    };
    void Matrix3x3::set_to_unity(void)
    {
        set_to_zero();
        storage[0][0] = 1.0;
        storage[1][1] = 1.0;
        storage[2][2] = 1.0;
    };
    Matrix3x3& Matrix3x3::operator=(const Matrix3x3& rhs)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                storage[i][j] = rhs(i, j);
        return *this;
    };
    Matrix3x3& Matrix3x3::operator=(const REAL rhs[3][3])
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                storage[i][j] = rhs[i][j];
        return *this;
    };
    bool Matrix3x3::operator==(Matrix3x3& rhs)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                if (storage[i][j] != rhs(i, j))
                    return false;
            }
        return true;
    };
    bool Matrix3x3::operator!=(Matrix3x3& rhs)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                if (storage[i][j] != rhs(i, j))
                    return true;
            }
        return false;
    };

    Matrix3x3 Matrix3x3::operator*(REAL m)
    {
        Matrix3x3 tmp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp(i, j) = storage[i][j] * m;
            }
        return tmp;
    };
    Matrix3x3 Matrix3x3::operator/(REAL m)
    {
        Matrix3x3 tmp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp(i, j) = storage[i][j] / m;
            }
        return tmp;
    };
    Matrix3x3 Matrix3x3::operator*(const Matrix3x3& rhs)
    {
        Matrix3x3 tmp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                {
                    tmp(i, j) += storage[i][k] * rhs(k, j);
                }
        return tmp;
    };
    Matrix3x3 Matrix3x3::operator+(Matrix3x3& rhs)
    {
        Matrix3x3 tmp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp(i, j) = storage[i][j] + rhs(i, j);
            }
        return tmp;
    };
    Matrix3x3 Matrix3x3::operator-(Matrix3x3& rhs)
    {
        Matrix3x3 tmp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp(i, j) = storage[i][j] - rhs(i, j);
            }
        return tmp;
    };

    Matrix3x3& Matrix3x3::operator+=(Matrix3x3& rhs)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                storage[i][j] += rhs(i, j);
            }
        return *this;
    };
    Matrix3x3& Matrix3x3::operator-=(Matrix3x3& rhs)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                storage[i][j] -= rhs(i, j);
            }
        return *this;
    };
    Matrix3x3& Matrix3x3::operator*=(REAL m)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                storage[i][j] *= m;
            }
        return *this;
    };
    Matrix3x3& Matrix3x3::operator/=(REAL m)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                storage[i][j] /= m;
            }
        return *this;
    };

    inline REAL Matrix3x3::cal_determinant(void) const
    {
        return (storage[0][0] * storage[1][1] * storage[2][2] +
            storage[0][1] * storage[1][2] * storage[2][0] +
            storage[0][2] * storage[1][0] * storage[2][1] -
            storage[0][2] * storage[1][1] * storage[2][0] -
            storage[0][1] * storage[1][0] * storage[2][2] -
            storage[0][0] * storage[1][2] * storage[2][1]);
    };
    Matrix3x3& Matrix3x3::do_invert(void)
    {
        Matrix3x3 tmp;
        REAL detInv = cal_determinant();

        if (detInv != 0.0) detInv = 1 / detInv;
        else
        {
            std::cout << "dMatrix3x3:\n" << this->print_matrix() << "Is not Invertible!" << std::endl;
            SYS_PROGRAM_STOP;
        }

        tmp(0, 0) = (storage[1][1] * storage[2][2] -
            storage[1][2] * storage[2][1]) * detInv;
        tmp(1, 0) = -(storage[1][0] * storage[2][2] -
            storage[1][2] * storage[2][0]) * detInv;
        tmp(2, 0) = (storage[1][0] * storage[2][1] -
            storage[1][1] * storage[2][0]) * detInv;
        tmp(0, 1) = -(storage[0][1] * storage[2][2] -
            storage[0][2] * storage[2][1]) * detInv;
        tmp(1, 1) = (storage[0][0] * storage[2][2] -
            storage[0][2] * storage[2][0]) * detInv;
        tmp(2, 1) = -(storage[0][0] * storage[2][1] -
            storage[0][1] * storage[2][0]) * detInv;
        tmp(0, 2) = (storage[0][1] * storage[1][2] -
            storage[1][1] * storage[0][2]) * detInv;
        tmp(1, 2) = -(storage[0][0] * storage[1][2] -
            storage[0][2] * storage[1][0]) * detInv;
        tmp(2, 2) = (storage[0][0] * storage[1][1] -
            storage[0][1] * storage[1][0]) * detInv;

        (*this) = tmp;
        return *this;
    };
    Matrix3x3 Matrix3x3::get_inverted_Matrix(void) const
    {
        Matrix3x3 tmp;

        REAL detInv = cal_determinant();

        if (detInv != 0.0) detInv = 1 / detInv;
        else
        {
            std::cout << "dMatrix3x3:\n" << this->print_matrix() << "Is not Invertible!" << std::endl;
            SYS_PROGRAM_STOP;
        }

        tmp(0, 0) = (storage[1][1] * storage[2][2] -
            storage[1][2] * storage[2][1]) * detInv;
        tmp(1, 0) = -(storage[1][0] * storage[2][2] -
            storage[1][2] * storage[2][0]) * detInv;
        tmp(2, 0) = (storage[1][0] * storage[2][1] -
            storage[1][1] * storage[2][0]) * detInv;
        tmp(0, 1) = -(storage[0][1] * storage[2][2] -
            storage[0][2] * storage[2][1]) * detInv;
        tmp(1, 1) = (storage[0][0] * storage[2][2] -
            storage[0][2] * storage[2][0]) * detInv;
        tmp(2, 1) = -(storage[0][0] * storage[2][1] -
            storage[0][1] * storage[2][0]) * detInv;
        tmp(0, 2) = (storage[0][1] * storage[1][2] -
            storage[1][1] * storage[0][2]) * detInv;
        tmp(1, 2) = -(storage[0][0] * storage[1][2] -
            storage[0][2] * storage[1][0]) * detInv;
        tmp(2, 2) = (storage[0][0] * storage[1][1] -
            storage[0][1] * storage[1][0]) * detInv;

        return tmp;
    };
    Matrix3x3& Matrix3x3::do_transpose(void)
    {
        REAL tmp[3][3] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp[i][j] = storage[j][i];
            }
        (*this) = tmp;
        return *this;
    };
    Matrix3x3 Matrix3x3::get_transposed(void) const
    {
        Matrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                TempMat(i, j) = storage[j][i];
            }
        return TempMat;
    };
    Matrix3x3& Matrix3x3::do_rotate(Matrix3x3& RotationMatrix)
    {
        Matrix3x3 TempMat;

        TempMat = (*this);

        (*this) = RotationMatrix * (TempMat * RotationMatrix.get_transposed());

        return *this;
    };
    Matrix3x3 Matrix3x3::get_rotated_Matrix(Matrix3x3& RotationMatrix) const
    {
        Matrix3x3 Out;
        Out = (*this);

        Out = RotationMatrix * (Out * RotationMatrix.get_transposed());

        return Out;
    };
    REAL Matrix3x3::REAL_contract(Matrix3x3& rHS) const
    {
        REAL tmp = 0.0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp += storage[i][j] * rHS(i, j);
            }
        return tmp;
    };
    REAL Matrix3x3::trace(void) const
    {
        return storage[0][0] + storage[1][1] + storage[2][2];
    };
    Matrix3x3 Matrix3x3::get_sym(void) const
    {
        Matrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                TempMat(i, j) = storage[i][j] + storage[j][i];
            }
        TempMat *= 0.5;

        return TempMat;
    };
    Matrix3x3 Matrix3x3::get_skew(void) const
    {
        Matrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                TempMat(i, j) = storage[i][j] - storage[j][i];
            }
        TempMat *= 0.5;

        return TempMat;
    };
    REAL Matrix3x3::norm(void) const  /// Frobenius norm
    {
        REAL tmp = 0.0;
        for (int i = 0; i < 3; i++)
        {
            tmp += storage[i][0] * storage[i][0]
                + storage[i][1] * storage[i][1]
                + storage[i][2] * storage[i][2];
        }

        return std::sqrt(tmp);
    };
    std::string Matrix3x3::print_matrix(void) const
    {
        std::stringstream out;
        for (int i = 0; i < 3; i++)
        {
            out << "||" << std::setprecision(6) << std::right
                << std::setw(8) << storage[i][0] << " "
                << std::setw(8) << storage[i][1] << " "
                << std::setw(8) << storage[i][2] << "||\n";
        }
        return out.str();
    };
    Vector3::Vector3() {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
    };
    Vector3::Vector3(REAL i, REAL j, REAL k) {
        storage[0] = i;
        storage[1] = j;
        storage[2] = k;
    }

    Vector3::Vector3(const Vector3& rhs)
    {
        for (int i = 0; i < 3; i++)
            storage[i] = rhs[i];
    };

    Vector3::Vector3(const REAL vecinit[3])
    {
#ifdef DEBUG
        if (vecinit.size() != 3)
        {
            std::cout << "Error in dVector3::constructor()\n"
                << "Initialization list size beyond storage range."
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        storage[0] = vecinit[0];
        storage[1] = vecinit[1];
        storage[2] = vecinit[2];
    }

    REAL& Vector3::operator[](const int i)
    {
#ifdef DEBUG
        if (i > 2)
        {
            std::cout << "Error in dVector3::operator[]\n"
                << "Access beyond storage range. i = "
                << i << " > 2"
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    REAL const& Vector3::operator[](const int i) const
    {
#ifdef DEBUG
        if (i > 2)
        {
            std::cout << "Error in dVector3::operator[]\n"
                << "Access beyond storage range. i = "
                << i << " > 2"
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };

    REAL Vector3::getX(void) const
    {
        return storage[0];
    };

    void Vector3::setX(const REAL newX)
    {
        storage[0] = newX;
    };

    REAL Vector3::getY(void) const
    {
        return storage[1];
    };

    void Vector3::setY(const REAL newY)
    {
        storage[1] = newY;
    };

    REAL Vector3::getZ(void) const
    {
        return storage[2];
    };

    void Vector3::setZ(const REAL newX)
    {
        storage[2] = newX;
    };

    void Vector3::set_to_zero(void)
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
    };
    void Vector3::set_to_unitX(void)
    {
        storage[0] = 1.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
    };
    void Vector3::set_to_unitY(void)
    {
        storage[0] = 0.0;
        storage[1] = 1.0;
        storage[2] = 0.0;
    };
    void Vector3::set_to_unitZ(void)
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 1.0;
    };
    bool Vector3::operator==(const Vector3& rhs) const
    {
        for (int i = 0; i < 3; i++) {
            if (storage[i] != rhs[i])
                return false;
        }
        return true;
    };
    bool Vector3::operator!=(const Vector3& rhs) const
    {
        for (int i = 0; i < 3; i++) {
            if (storage[i] != rhs[i])
                return true;
        }
        return false;
    };
    Vector3 Vector3::operator*(const REAL m) const
    {
        Vector3 tmp;
        tmp[0] = storage[0] * m;
        tmp[1] = storage[1] * m;
        tmp[2] = storage[2] * m;
        return tmp;
    };
    Vector3 Vector3::operator/(const REAL m) const
    {
        Vector3 tmp;
        tmp[0] = storage[0] / m;
        tmp[1] = storage[1] / m;
        tmp[2] = storage[2] / m;
        return tmp;
    };
    REAL Vector3::operator*(const Vector3& rhs) const
    {
        return storage[0] * rhs[0] + storage[1] * rhs[1] + storage[2] * rhs[2];
    };
    REAL Vector3::abs() const
    {
        return std::sqrt(storage[0] * storage[0] +
            storage[1] * storage[1] +
            storage[2] * storage[2]);
    };
    Vector3 Vector3::cross(const Vector3& rhs) const
    {
        Vector3 tmp;
        tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
        tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
        tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
        return tmp;
    };
    Vector3 Vector3::operator+(const Vector3& rhs) const
    {
        Vector3 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        return tmp;
    };
    Vector3 Vector3::operator-(const Vector3& rhs) const
    {
        Vector3 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        return tmp;
    };
    Vector3& Vector3::operator*=(const REAL m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        return *this;
    };
    Vector3& Vector3::operator/=(const REAL m)
    {
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        return *this;
    };
    Vector3& Vector3::operator-=(const Vector3& rhs)
    {
        storage[0] = storage[0] - rhs[0];
        storage[1] = storage[1] - rhs[1];
        storage[2] = storage[2] - rhs[2];
        return *this;
    };
    Vector3& Vector3::operator+=(const Vector3& rhs)
    {
        storage[0] = storage[0] + rhs[0];
        storage[1] = storage[1] + rhs[1];
        storage[2] = storage[2] + rhs[2];
        return *this;
    };
    Vector3& Vector3::operator=(const Vector3& rhs)
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        return *this;
    };
    Vector3& Vector3::operator=(const REAL rhs[3])
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        return *this;
    };
    REAL Vector3::length(void) const
    {
        return std::sqrt(storage[0] * storage[0] +
            storage[1] * storage[1] +
            storage[2] * storage[2]);
    };
    Vector3& Vector3::normalize(void)
    {
        REAL norm = std::sqrt(storage[0] * storage[0] +
            storage[1] * storage[1] +
            storage[2] * storage[2]);
        if (norm != 0.0)
        {
            storage[0] /= norm;
            storage[1] /= norm;
            storage[2] /= norm;
        }
        return *this;
    };
    Vector3 Vector3::normalized(void) const
    {
        REAL norm = std::sqrt(storage[0] * storage[0] +
            storage[1] * storage[1] +
            storage[2] * storage[2]);
        Vector3 tmp;
        if (norm != 0.0)
        {
            tmp[0] = storage[0] / norm;
            tmp[1] = storage[1] / norm;
            tmp[2] = storage[2] / norm;
        }
        else
        {
            tmp.set_to_zero();
        }
        return tmp;
    };
    Vector3& Vector3::do_rotate(const Matrix3x3& RotationMatrix)
    {
        Vector3 Out;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                Out[i] += RotationMatrix(i, j) * storage[j];
            }
        (*this) = Out;
        return *this;
    };
    Vector3 Vector3::get_rotated_vec3(const Matrix3x3& RotationMatrix) const
    {
        Vector3 Out;
        Out.set_to_zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                Out[i] += RotationMatrix(i, j) * storage[j];
            }
        return Out;
    };
    Vector3 Vector3::Xreflected(void) const
    {
        Vector3 Out;
        for (int i = 0; i < 3; ++i)
        {
            Out[i] = storage[i];
        }
        Out[0] *= -1.0;
        return Out;
    };
    Vector3 Vector3::Yreflected(void) const
    {
        Vector3 Out;
        for (int i = 0; i < 3; ++i)
        {
            Out[i] = storage[i];
        }
        Out[1] *= -1.0;
        return Out;
    };
    Vector3 Vector3::Zreflected(void) const
    {
        Vector3 Out;
        for (int i = 0; i < 3; ++i)
        {
            Out[i] = storage[i];
        }
        Out[2] *= -1.0;
        return Out;
    };
    std::string Vector3::print(void) const
    {
        std::stringstream out;

        out << "(" << storage[0] << ", "
            << storage[1] << ", "
            << storage[2] << ")";
        return out.str();
    };
    REAL* Vector3::data(void)
    {
        return storage;
    };
    const REAL* Vector3::const_data(void) const
    {
        return storage;
    };
    Vector3 Matrix3x3::operator*(const Vector3& rhs) const
    {
        Vector3 tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    }
    Vector6::Vector6() {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
    }

    Vector6::Vector6(REAL v1, REAL v2, REAL v3, REAL v4, REAL v5, REAL v6) {
        storage[0] = v1;
        storage[1] = v2;
        storage[2] = v3;
        storage[3] = v4;
        storage[4] = v5;
        storage[5] = v6;
    }

    Vector6::Vector6(std::initializer_list<REAL> vecinit)
    {
#ifdef DEBUG
        if (vecinit.size() != 6)
        {
            std::cout << "Error in Vector6::constructor()\n"
                << "Initialization list size beyond storage range."
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }

    REAL& Vector6::operator[](const int i)
    {
#ifdef DEBUG
        if (i > 5)
        {
            std::cout << "Error in Vector6::operator[]\n"
                << "Access beyond storage range. i = "
                << i << " > 5"
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    }
    REAL const& Vector6::operator[](const int i) const
    {
#ifdef DEBUG
        if (i > 5)
        {
            std::cout << "Error in dVector6::operator[]\n"
                << "Access beyond storage range. i = "
                << i << " > 5"
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    bool Vector6::is_nan_val_exist() {
        for (int i = 0; i < 6; i++)
            if (std::isnan(storage[i]))
                return true;
        return false;
    }
    void Vector6::set_to_zero(void)
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
    };
    void Vector6::set_to_unity(void)
    {
        storage[0] = 1.0;
        storage[1] = 1.0;
        storage[2] = 1.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
    };
    void Vector6::do_abs(void)
    {
        storage[0] = abs(storage[0]);
        storage[1] = abs(storage[1]);
        storage[2] = abs(storage[2]);
        storage[3] = abs(storage[3]);
        storage[4] = abs(storage[4]);
        storage[5] = abs(storage[5]);
    };

    REAL Vector6::norm(void) const  /// Frobenius norm
    {
        REAL tmp = storage[0] * storage[0]
            + storage[1] * storage[1]
            + storage[2] * storage[2]
            + storage[3] * storage[3] * 2
            + storage[4] * storage[4] * 2
            + storage[5] * storage[5] * 2;
        return std::sqrt(tmp);
    };

    REAL Vector6::trace(void) const
    {
        return storage[0] + storage[1] + storage[2];
    };
    REAL Vector6::max_abs(void) const
    {
        // Returns maximum absolute value
        REAL tempmax = 0.0;
        for (int i = 0; i < 6; i++)
        {
            if (abs(storage[i]) > tempmax)
            {
                tempmax = abs(storage[i]);
            }
        }
        return tempmax;
    };
    Vector6 Vector6::operator*(const REAL m) const
    {
        Vector6 tmp;
        tmp[0] = storage[0] * m;
        tmp[1] = storage[1] * m;
        tmp[2] = storage[2] * m;
        tmp[3] = storage[3] * m;
        tmp[4] = storage[4] * m;
        tmp[5] = storage[5] * m;
        return tmp;
    };
    bool Vector6::operator==(const Vector6& rhs) const
    {
        for (int i = 0; i < 6; i++) {
            if (storage[i] != rhs[i])
                return false;
        }
        return true;
    };
    bool Vector6::operator!=(const Vector6& rhs) const
    {
        for (int i = 0; i < 6; i++) {
            if (storage[i] != rhs[i])
                return true;
        }
        return false;
    };
    Vector6 Vector6::operator/(const REAL m) const
    {
        Vector6 tmp;
        tmp[0] = storage[0] / m;
        tmp[1] = storage[1] / m;
        tmp[2] = storage[2] / m;
        tmp[3] = storage[3] / m;
        tmp[4] = storage[4] / m;
        tmp[5] = storage[5] / m;
        return tmp;
    };
    Vector6& Vector6::operator*=(const REAL m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        storage[3] *= m;
        storage[4] *= m;
        storage[5] *= m;
        return *this;
    };

    REAL Vector6::operator*(const Vector6& rhs)
    {
        REAL tmp = 0.0;
        for (int i = 0; i < 6; i++)
        {
            tmp += storage[i] * rhs[i];
        }
        return tmp;
    };

    Vector6& Vector6::operator/=(const REAL m)
    {
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        storage[3] /= m;
        storage[4] /= m;
        storage[5] /= m;

        return *this;
    };
    Vector6 Vector6::operator+(const Vector6& rhs) const
    {
        Vector6 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        tmp[3] = storage[3] + rhs[3];
        tmp[4] = storage[4] + rhs[4];
        tmp[5] = storage[5] + rhs[5];
        return tmp;
    };
    Vector6& Vector6::operator+=(const Vector6& rhs)
    {
        storage[0] = storage[0] + rhs[0];
        storage[1] = storage[1] + rhs[1];
        storage[2] = storage[2] + rhs[2];
        storage[3] = storage[3] + rhs[3];
        storage[4] = storage[4] + rhs[4];
        storage[5] = storage[5] + rhs[5];
        return *this;
    };
    Vector6 Vector6::operator-(const Vector6& rhs) const
    {
        Vector6 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    Vector6& Vector6::operator-=(const Vector6& rhs)
    {
        storage[0] = storage[0] - rhs[0];
        storage[1] = storage[1] - rhs[1];
        storage[2] = storage[2] - rhs[2];
        storage[3] = storage[3] - rhs[3];
        storage[4] = storage[4] - rhs[4];
        storage[5] = storage[5] - rhs[5];
        return *this;
    };
    Vector6& Vector6::operator=(const Vector6& rhs)
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        storage[3] = rhs[3];
        storage[4] = rhs[4];
        storage[5] = rhs[5];
        return *this;
    };
    Vector6& Vector6::operator=(const REAL rhs[6])
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        storage[3] = rhs[3];
        storage[4] = rhs[4];
        storage[5] = rhs[5];
        return *this;
    };
    Vector6 Vector6::get_rotated_vec6(const Matrix3x3& RotationMatrix) const
    {
        REAL In[3][3] = {};
        REAL Out[3][3] = {};

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                Out[p][q] = 0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                    }
            }
        Vector6 out;

        out[0] = Out[0][0];
        out[5] = Out[0][1];
        out[4] = Out[0][2];
        out[1] = Out[1][1];
        out[3] = Out[1][2];
        out[2] = Out[2][2];

        return out;
    };
    Vector6& Vector6::do_rotate(Matrix3x3 RotationMatrix)
    {
        REAL In[3][3] = {};
        REAL Out[3][3] = {};

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                Out[p][q] = 0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                    }
            }

        storage[0] = Out[0][0];
        storage[5] = Out[0][1];
        storage[4] = Out[0][2];
        storage[1] = Out[1][1];
        storage[3] = Out[1][2];
        storage[2] = Out[2][2];

        return *this;
    };
    REAL* Vector6::data(void)
    {
        return storage;
    };
    const REAL* Vector6::const_data(void) const
    {
        return storage;
    };
    Vector6& Vector6::H_product(Vector6 vector)
    {
        storage[0] *= vector[0];
        storage[1] *= vector[1];
        storage[2] *= vector[2];
        storage[3] *= vector[3];
        storage[4] *= vector[4];
        storage[5] *= vector[5];

        return *this;
    };

    Matrix3x3 Vector6::sym_tensor(void) const
    {
        Matrix3x3 tmp;
        tmp(0, 0) = storage[0];
        tmp(0, 1) = storage[5];
        tmp(0, 2) = storage[4];
        tmp(1, 0) = storage[5];
        tmp(1, 1) = storage[1];
        tmp(1, 2) = storage[3];
        tmp(2, 0) = storage[4];
        tmp(2, 1) = storage[3];
        tmp(2, 2) = storage[2];
        return tmp;
    };

    Matrix3x3 Vector6::skew_tensor(void) const
    {
        Matrix3x3 tmp;
        tmp(0, 0) = storage[0];
        tmp(0, 1) = storage[5];
        tmp(0, 2) = storage[4];
        tmp(1, 0) = -storage[5];
        tmp(1, 1) = storage[1];
        tmp(1, 2) = storage[3];
        tmp(2, 0) = -storage[4];
        tmp(2, 1) = -storage[3];
        tmp(2, 2) = storage[2];
        return tmp;
    };

    std::string Vector6::print(void) const
    {
        std::stringstream out;
        out << "< | ";
        for (int i = 0; i < 6; i++)
        {
            out << storage[i] << " " << " | ";
        }
        out << " >";
        return out.str();
    };

    Matrix6x6::Matrix6x6() {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                storage[i][j] = 0.0;
    }
    REAL& Matrix6x6::operator()(const int i, const int j)
    {
#ifdef DEBUG
        if (i > 5 or j > 5)
        {
            std::cout << "Error in Matrix6x6::operator()\n"
                << "Access beyond storage range. (i, j) = "
                << i << ", " << j << " > (5, 5)"
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i][j];
    };
    bool Matrix6x6::is_nan_val_exist() {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                if (std::isnan(storage[i][j]))
                    return true;
        return false;
    }
    REAL const& Matrix6x6::operator()(const int i, const int j) const
    {
#ifdef DEBUG
        if (i > 5 or j > 5)
        {
            std::cout << "Error in Matrix6x6::operator()\n"
                << "Access beyond storage range. (i, j) = "
                << i << ", " << j << " > (5, 5)"
                << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i][j];
    };
    void Matrix6x6::set_to_zero(void)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                storage[i][j] = 0.0;
    };
    void Matrix6x6::set_to_unity(void)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                if (i == j)
                    storage[i][j] = 1.0;
                else
                    storage[i][j] = 0.0;
    };

    REAL Matrix6x6::norm(void) const
    {
        REAL tmp = 0.0;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp += storage[i][j] * storage[i][j];
            }
        return std::sqrt(tmp);
    };

    Matrix6x6 Matrix6x6::operator*(const REAL m) const
    {
        Matrix6x6 tmp;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp(i, j) = storage[i][j] * m;
            }
        return tmp;
    };
    Matrix6x6 Matrix6x6::operator/(const REAL m) const
    {
        Matrix6x6 tmp;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp(i, j) = storage[i][j] / m;
            }
        return tmp;
    };
    bool Matrix6x6::operator==(const Matrix6x6& rhs) const
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
                if (storage[i][j] != rhs(i, j))
                    return false;
            }
        return true;
    };
    bool Matrix6x6::operator!=(const Matrix6x6& rhs) const
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
                if (storage[i][j] != rhs(i, j))
                    return true;
            }
        return false;
    };
    Vector6 Matrix6x6::operator*(Vector6 rhs) const
    {
        Vector6 tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    };
    Matrix6x6 Matrix6x6::operator*(const Matrix6x6& rhs) const
    {
        Matrix6x6 tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                for (int k = 0; k < 6; k++)
                {
                    tmp(i, j) += storage[i][k] * rhs(k, j);
                }
        return tmp;
    };
    Matrix6x6 Matrix6x6::operator+(const Matrix6x6& rhs) const
    {
        Matrix6x6 tmp;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp(i, j) = storage[i][j] + rhs(i, j);
            }
        return tmp;
    };
    Matrix6x6 Matrix6x6::operator-(const Matrix6x6& rhs) const
    {
        Matrix6x6 tmp;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp(i, j) = storage[i][j] - rhs(i, j);
            }
        return tmp;
    };
    Matrix6x6& Matrix6x6::operator+=(const Matrix6x6& rhs)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                storage[i][j] += rhs(i, j);
            }
        return *this;
    };
    Matrix6x6& Matrix6x6::operator/=(const REAL m)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                storage[i][j] /= m;
            }
        return *this;
    };
    Matrix6x6& Matrix6x6::operator*=(const REAL m)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                storage[i][j] *= m;
            }
        return *this;
    };
    Matrix6x6& Matrix6x6::operator-=(const Matrix6x6& rhs)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                storage[i][j] -= rhs(i, j);
            }
        return *this;
    };
    Matrix6x6& Matrix6x6::operator=(const Matrix6x6& rhs)
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                storage[i][j] = rhs(i, j);
        return *this;
    };
    Matrix6x6& Matrix6x6::operator=(const REAL rhs[6][6])
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                storage[i][j] = rhs[i][j];
        return *this;
    };
    Matrix6x6& Matrix6x6::do_invert(void)
    {
        REAL Out[6][6] = {};

        int indxc[6] = {};
        int indxr[6] = {};
        int ipiv[6] = { 0, 0, 0, 0, 0, 0 };
        int icol = 0;
        int irow = 0;
        REAL pivinv;
        REAL dum;
        REAL Uni[6][6] = { {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} };


        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                Out[i][j] = storage[i][j];
            }

        for (int i = 0; i < 6; i++)
        {
            REAL big = 0.0;
            for (int j = 0; j < 6; j++)
                if (ipiv[j] != 1)
                    for (int k = 0; k < 6; k++)
                    {
                        if (ipiv[k] == 0)
                        {
                            if (fabs(Out[j][k]) >= big)
                            {
                                big = fabs(Out[j][k]);
                                irow = j;
                                icol = k;
                            };
                        }
                        else if (ipiv[k] > 1)
                        {
                            std::cout << "Matrix6x6: Can Not Compute Inverse Matrix."
                                << "Matrix:\n" << this->print()
                                << "is Singular 1!!!" << std::endl;
                            SYS_PROGRAM_STOP;
                        }
                    };
            ++(ipiv[icol]);
            if (irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    REAL temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < 6; l++)
                {
                    REAL temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cout << "dMatrix6x6: Can Not Compute Inverse Matrix. Matrix:\n" << this->print() << "is Singular 2!!!" << std::endl;
                SYS_PROGRAM_STOP;
            }
            pivinv = 1 / Out[icol][icol];
            Out[icol][icol] = 1.0;
            for (int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
            for (int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
            for (int ll = 0; ll < 6; ll++)
                if (ll != icol)
                {
                    dum = Out[ll][icol];
                    Out[ll][icol] = 0.0;
                    for (int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l] * dum;
                    for (int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l] * dum;
                }
        }
        for (int l = 5; l >= 0; l--)
        {
            if (indxr[l] != indxc[l])
                for (int k = 0; k < 6; k++)
                {
                    REAL temp = Out[k][indxr[l]];
                    Out[k][indxr[l]] = Out[k][indxc[l]];
                    Out[k][indxc[l]] = temp;
                };
        }

        (*this) = Out;
        return *this;
    };

    Matrix6x6 Matrix6x6::get_inverted_matrix(void) const
    {
        REAL Out[6][6] = {};

        int indxc[6] = {};
        int indxr[6] = {};
        int ipiv[6] = { 0, 0, 0, 0, 0, 0 };
        int icol = 0;
        int irow = 0;
        REAL pivinv;
        REAL dum;
        REAL Uni[6][6] = { {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} };

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                Out[i][j] = storage[i][j];
            }

        for (int i = 0; i < 6; i++)
        {
            REAL big = 0.0;
            for (int j = 0; j < 6; j++)
                if (ipiv[j] != 1)
                    for (int k = 0; k < 6; k++)
                    {
                        if (ipiv[k] == 0)
                        {
                            if (fabs(Out[j][k]) >= big)
                            {
                                big = fabs(Out[j][k]);
                                irow = j;
                                icol = k;
                            };
                        }
                        else if (ipiv[k] > 1)
                        {
                            std::cout << "Matrix6x6: Can Not Compute Inverse Matrix.\n"
                                << this->print() << "Matrix is Singular 2!!!"
                                << std::endl;
                            SYS_PROGRAM_STOP;
                        }
                    };
            ++(ipiv[icol]);
            if (irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    REAL temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < 6; l++)
                {
                    REAL temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cout << "Matrix6x6:\n" << this->print()
                    << " Is not Invertible!" << std::endl;
                SYS_PROGRAM_STOP;
            }
            pivinv = 1 / Out[icol][icol];
            Out[icol][icol] = 1.0;
            for (int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
            for (int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
            for (int ll = 0; ll < 6; ll++)
                if (ll != icol)
                {
                    dum = Out[ll][icol];
                    Out[ll][icol] = 0.0;
                    for (int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l] * dum;
                    for (int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l] * dum;
                }
        }
        for (int l = 5; l >= 0; l--)
        {
            if (indxr[l] != indxc[l])
                for (int k = 0; k < 6; k++)
                {
                    REAL temp = Out[k][indxr[l]];
                    Out[k][indxr[l]] = Out[k][indxc[l]];
                    Out[k][indxc[l]] = temp;
                };
        }

        Matrix6x6 tmp;
        tmp = Out;
        return tmp;
    };

    Matrix6x6& Matrix6x6::do_transpose(void)
    {
        REAL tmp[6][6] = {};
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp[i][j] = storage[j][i];
            }
        (*this) = tmp;
        return *this;
    };

    Matrix6x6 Matrix6x6::get_transposed_matrix(void) const
    {
        Matrix6x6 tmp;

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp(i, j) = storage[j][i];
            }
        return tmp;
    };

    Matrix6x6& Matrix6x6::do_rotate(const Matrix3x3& RotationMatrix)
    {
        REAL In[3][3][3][3] = {};
        REAL Out[3][3][3][3] = {};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        In[i][j][k][l] = tensor(i, j, k, l);
                    }

        for (int m = 0; m < 3; m++)
            for (int n = 0; n < 3; n++)
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                    {
                        Out[m][n][p][q] = 0.0;

                        for (int i = 0; i < 3; i++)
                            for (int j = 0; j < 3; j++)
                                for (int k = 0; k < 3; k++)
                                    for (int l = 0; l < 3; l++)
                                    {
                                        //                Former version
                                        //                Out[m][n][p][q] += RotationMatrix(i,m)*
                                        //                                   RotationMatrix(j,n)*
                                        //                                   In[i][j][k][l]*
                                        //                                   RotationMatrix(k,p)*
                                        //                                   RotationMatrix(l,q);
                                        //              New version (active rotation)
                                        Out[m][n][p][q] += RotationMatrix(m, i) *
                                            RotationMatrix(n, j) *
                                            In[i][j][k][l] *
                                            RotationMatrix(p, k) *
                                            RotationMatrix(q, l);
                                    }
                    }
        int VoigtIndex[6][2] = { {0,0},{1,1},{2,2},{1,2},{0,2},{0,1} };
        for (int m = 0; m < 6; m++)
            for (int n = 0; n < 6; n++)
            {
                int i = VoigtIndex[m][0];
                int j = VoigtIndex[m][1];
                int k = VoigtIndex[n][0];
                int l = VoigtIndex[n][1];

                storage[m][n] = Out[i][j][k][l];
            }
        return *this;
    };
    Matrix6x6 Matrix6x6::get_rotated_matrix(const Matrix3x3& RotationMatrix) const
    {
        REAL In[3][3][3][3] = {};
        REAL Out[3][3][3][3] = {};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        In[i][j][k][l] = tensor(i, j, k, l);
                    }

        for (int m = 0; m < 3; m++)
            for (int n = 0; n < 3; n++)
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                    {
                        Out[m][n][p][q] = 0.0;

                        for (int i = 0; i < 3; i++)
                            for (int j = 0; j < 3; j++)
                                for (int k = 0; k < 3; k++)
                                    for (int l = 0; l < 3; l++)
                                    {
                                        //                Former version
                                        //                Out[m][n][p][q] += RotationMatrix(i,m)*
                                        //                                   RotationMatrix(j,n)*
                                        //                                   In[i][j][k][l]*
                                        //                                   RotationMatrix(k,p)*
                                        //                                   RotationMatrix(l,q);
                                        //              New version (active rotation)
                                        Out[m][n][p][q] += RotationMatrix(m, i) *
                                            RotationMatrix(n, j) *
                                            In[i][j][k][l] *
                                            RotationMatrix(p, k) *
                                            RotationMatrix(q, l);
                                    }
                    }
        Matrix6x6 _OUT;
        int VoigtIndex[6][2] = { {0,0},{1,1},{2,2},{1,2},{0,2},{0,1} };

        for (int m = 0; m < 6; m++)
            for (int n = 0; n < 6; n++)
            {
                int i = VoigtIndex[m][0];
                int j = VoigtIndex[m][1];
                int k = VoigtIndex[n][0];
                int l = VoigtIndex[n][1];

                _OUT(m, n) = Out[i][j][k][l];
            }
        return _OUT;
    };
    Matrix6x6& Matrix6x6::do_rotate_Voigt(const Matrix3x3& RotationMatrix) {
        // Build 6x6 Bond transformation matrix T
        Matrix6x6 T;
        auto R = [&](int i, int j) { return RotationMatrix(i, j); };

        T(0, 0) = R(0, 0) * R(0, 0); T(0, 1) = R(0, 1) * R(0, 1); T(0, 2) = R(0, 2) * R(0, 2);
        T(0, 3) = 2 * R(0, 1) * R(0, 2); T(0, 4) = 2 * R(0, 0) * R(0, 2); T(0, 5) = 2 * R(0, 0) * R(0, 1);

        T(1, 0) = R(1, 0) * R(1, 0); T(1, 1) = R(1, 1) * R(1, 1); T(1, 2) = R(1, 2) * R(1, 2);
        T(1, 3) = 2 * R(1, 1) * R(1, 2); T(1, 4) = 2 * R(1, 0) * R(1, 2); T(1, 5) = 2 * R(1, 0) * R(1, 1);

        T(2, 0) = R(2, 0) * R(2, 0); T(2, 1) = R(2, 1) * R(2, 1); T(2, 2) = R(2, 2) * R(2, 2);
        T(2, 3) = 2 * R(2, 1) * R(2, 2); T(2, 4) = 2 * R(2, 0) * R(2, 2); T(2, 5) = 2 * R(2, 0) * R(2, 1);

        T(3, 0) = R(1, 0) * R(2, 0); T(3, 1) = R(1, 1) * R(2, 1); T(3, 2) = R(1, 2) * R(2, 2);
        T(3, 3) = R(1, 1) * R(2, 2) + R(1, 2) * R(2, 1);
        T(3, 4) = R(1, 0) * R(2, 2) + R(1, 2) * R(2, 0);
        T(3, 5) = R(1, 0) * R(2, 1) + R(1, 1) * R(2, 0);

        T(4, 0) = R(0, 0) * R(2, 0); T(4, 1) = R(0, 1) * R(2, 1); T(4, 2) = R(0, 2) * R(2, 2);
        T(4, 3) = R(0, 1) * R(2, 2) + R(0, 2) * R(2, 1);
        T(4, 4) = R(0, 0) * R(2, 2) + R(0, 2) * R(2, 0);
        T(4, 5) = R(0, 0) * R(2, 1) + R(0, 1) * R(2, 0);

        T(5, 0) = R(0, 0) * R(1, 0); T(5, 1) = R(0, 1) * R(1, 1); T(5, 2) = R(0, 2) * R(1, 2);
        T(5, 3) = R(0, 1) * R(1, 2) + R(0, 2) * R(1, 1);
        T(5, 4) = R(0, 0) * R(1, 2) + R(0, 2) * R(1, 0);
        T(5, 5) = R(0, 0) * R(1, 1) + R(0, 1) * R(1, 0);
        (*this) = T * (*this) * T.get_transposed_matrix();

        return (*this);
    }
    Matrix6x6 Matrix6x6::get_rotated_matrix_Voigt(const Matrix3x3& RotationMatrix) const {
        // Build 6x6 Bond transformation matrix T
        Matrix6x6 T;
        auto R = [&](int i, int j) { return RotationMatrix(i, j); };

        T(0, 0) = R(0, 0) * R(0, 0); T(0, 1) = R(0, 1) * R(0, 1); T(0, 2) = R(0, 2) * R(0, 2);
        T(0, 3) = 2 * R(0, 1) * R(0, 2); T(0, 4) = 2 * R(0, 0) * R(0, 2); T(0, 5) = 2 * R(0, 0) * R(0, 1);

        T(1, 0) = R(1, 0) * R(1, 0); T(1, 1) = R(1, 1) * R(1, 1); T(1, 2) = R(1, 2) * R(1, 2);
        T(1, 3) = 2 * R(1, 1) * R(1, 2); T(1, 4) = 2 * R(1, 0) * R(1, 2); T(1, 5) = 2 * R(1, 0) * R(1, 1);

        T(2, 0) = R(2, 0) * R(2, 0); T(2, 1) = R(2, 1) * R(2, 1); T(2, 2) = R(2, 2) * R(2, 2);
        T(2, 3) = 2 * R(2, 1) * R(2, 2); T(2, 4) = 2 * R(2, 0) * R(2, 2); T(2, 5) = 2 * R(2, 0) * R(2, 1);

        T(3, 0) = R(1, 0) * R(2, 0); T(3, 1) = R(1, 1) * R(2, 1); T(3, 2) = R(1, 2) * R(2, 2);
        T(3, 3) = R(1, 1) * R(2, 2) + R(1, 2) * R(2, 1);
        T(3, 4) = R(1, 0) * R(2, 2) + R(1, 2) * R(2, 0);
        T(3, 5) = R(1, 0) * R(2, 1) + R(1, 1) * R(2, 0);

        T(4, 0) = R(0, 0) * R(2, 0); T(4, 1) = R(0, 1) * R(2, 1); T(4, 2) = R(0, 2) * R(2, 2);
        T(4, 3) = R(0, 1) * R(2, 2) + R(0, 2) * R(2, 1);
        T(4, 4) = R(0, 0) * R(2, 2) + R(0, 2) * R(2, 0);
        T(4, 5) = R(0, 0) * R(2, 1) + R(0, 1) * R(2, 0);

        T(5, 0) = R(0, 0) * R(1, 0); T(5, 1) = R(0, 1) * R(1, 1); T(5, 2) = R(0, 2) * R(1, 2);
        T(5, 3) = R(0, 1) * R(1, 2) + R(0, 2) * R(1, 1);
        T(5, 4) = R(0, 0) * R(1, 2) + R(0, 2) * R(1, 0);
        T(5, 5) = R(0, 0) * R(1, 1) + R(0, 1) * R(1, 0);

        return T * (*this) * T.get_transposed_matrix();
    }

    std::string Matrix6x6::print(void) const
    {
        std::stringstream out;
        for (int i = 0; i < 6; i++)
        {
            out << "||" << std::setprecision(6)
                << std::setw(8) << storage[i][0] << " "
                << std::setw(8) << storage[i][1] << " "
                << std::setw(8) << storage[i][2] << " "
                << std::setw(8) << storage[i][3] << " "
                << std::setw(8) << storage[i][4] << " "
                << std::setw(8) << storage[i][5] << "||\n";
        }
        return out.str();
    };
    const REAL& Matrix6x6::tensor(const int i, const int j, const int k, const int l) const
    {
        return storage[(i == j) ? (i) : (6 - (i + j))][(k == l) ? (k) : (6 - (k + l))];
    };
    REAL* Matrix6x6::data(void)
    {
        return &storage[0][0];
    };
    const REAL* Matrix6x6::const_data(void) const
    {
        return &storage[0][0];
    };
    Matrix6x6 Matrix3x3::outer(const Matrix3x3& rhs) const
    {
        Matrix6x6 tmp;

        tmp(0, 0) += storage[0][0] * rhs(0, 0);
        tmp(0, 5) += storage[0][0] * rhs(0, 1);
        tmp(0, 4) += storage[0][0] * rhs(0, 2);
        tmp(0, 5) += storage[0][0] * rhs(1, 0);
        tmp(0, 1) += storage[0][0] * rhs(1, 1);
        tmp(0, 3) += storage[0][0] * rhs(1, 2);
        tmp(0, 4) += storage[0][0] * rhs(2, 0);
        tmp(0, 3) += storage[0][0] * rhs(2, 1);
        tmp(0, 2) += storage[0][0] * rhs(2, 2);
        tmp(5, 0) += storage[0][1] * rhs(0, 0);
        tmp(5, 5) += storage[0][1] * rhs(0, 1);
        tmp(5, 4) += storage[0][1] * rhs(0, 2);
        tmp(5, 5) += storage[0][1] * rhs(1, 0);
        tmp(5, 1) += storage[0][1] * rhs(1, 1);
        tmp(5, 3) += storage[0][1] * rhs(1, 2);
        tmp(5, 4) += storage[0][1] * rhs(2, 0);
        tmp(5, 3) += storage[0][1] * rhs(2, 1);
        tmp(5, 2) += storage[0][1] * rhs(2, 2);
        tmp(4, 0) += storage[0][2] * rhs(0, 0);
        tmp(4, 5) += storage[0][2] * rhs(0, 1);
        tmp(4, 4) += storage[0][2] * rhs(0, 2);
        tmp(4, 5) += storage[0][2] * rhs(1, 0);
        tmp(4, 1) += storage[0][2] * rhs(1, 1);
        tmp(4, 3) += storage[0][2] * rhs(1, 2);
        tmp(4, 4) += storage[0][2] * rhs(2, 0);
        tmp(4, 3) += storage[0][2] * rhs(2, 1);
        tmp(4, 2) += storage[0][2] * rhs(2, 2);
        tmp(5, 0) += storage[1][0] * rhs(0, 0);
        tmp(5, 5) += storage[1][0] * rhs(0, 1);
        tmp(5, 4) += storage[1][0] * rhs(0, 2);
        tmp(5, 5) += storage[1][0] * rhs(1, 0);
        tmp(5, 1) += storage[1][0] * rhs(1, 1);
        tmp(5, 3) += storage[1][0] * rhs(1, 2);
        tmp(5, 4) += storage[1][0] * rhs(2, 0);
        tmp(5, 3) += storage[1][0] * rhs(2, 1);
        tmp(5, 2) += storage[1][0] * rhs(2, 2);
        tmp(1, 0) += storage[1][1] * rhs(0, 0);
        tmp(1, 5) += storage[1][1] * rhs(0, 1);
        tmp(1, 4) += storage[1][1] * rhs(0, 2);
        tmp(1, 5) += storage[1][1] * rhs(1, 0);
        tmp(1, 1) += storage[1][1] * rhs(1, 1);
        tmp(1, 3) += storage[1][1] * rhs(1, 2);
        tmp(1, 4) += storage[1][1] * rhs(2, 0);
        tmp(1, 3) += storage[1][1] * rhs(2, 1);
        tmp(1, 2) += storage[1][1] * rhs(2, 2);
        tmp(3, 0) += storage[1][2] * rhs(0, 0);
        tmp(3, 5) += storage[1][2] * rhs(0, 1);
        tmp(3, 4) += storage[1][2] * rhs(0, 2);
        tmp(3, 5) += storage[1][2] * rhs(1, 0);
        tmp(3, 1) += storage[1][2] * rhs(1, 1);
        tmp(3, 3) += storage[1][2] * rhs(1, 2);
        tmp(3, 4) += storage[1][2] * rhs(2, 0);
        tmp(3, 3) += storage[1][2] * rhs(2, 1);
        tmp(3, 2) += storage[1][2] * rhs(2, 2);
        tmp(4, 0) += storage[2][0] * rhs(0, 0);
        tmp(4, 5) += storage[2][0] * rhs(0, 1);
        tmp(4, 4) += storage[2][0] * rhs(0, 2);
        tmp(4, 5) += storage[2][0] * rhs(1, 0);
        tmp(4, 1) += storage[2][0] * rhs(1, 1);
        tmp(4, 3) += storage[2][0] * rhs(1, 2);
        tmp(4, 4) += storage[2][0] * rhs(2, 0);
        tmp(4, 3) += storage[2][0] * rhs(2, 1);
        tmp(4, 2) += storage[2][0] * rhs(2, 2);
        tmp(3, 0) += storage[2][1] * rhs(0, 0);
        tmp(3, 5) += storage[2][1] * rhs(0, 1);
        tmp(3, 4) += storage[2][1] * rhs(0, 2);
        tmp(3, 5) += storage[2][1] * rhs(1, 0);
        tmp(3, 1) += storage[2][1] * rhs(1, 1);
        tmp(3, 3) += storage[2][1] * rhs(1, 2);
        tmp(3, 4) += storage[2][1] * rhs(2, 0);
        tmp(3, 3) += storage[2][1] * rhs(2, 1);
        tmp(3, 2) += storage[2][1] * rhs(2, 2);
        tmp(2, 0) += storage[2][2] * rhs(0, 0);
        tmp(2, 5) += storage[2][2] * rhs(0, 1);
        tmp(2, 4) += storage[2][2] * rhs(0, 2);
        tmp(2, 5) += storage[2][2] * rhs(1, 0);
        tmp(2, 1) += storage[2][2] * rhs(1, 1);
        tmp(2, 3) += storage[2][2] * rhs(1, 2);
        tmp(2, 4) += storage[2][2] * rhs(2, 0);
        tmp(2, 3) += storage[2][2] * rhs(2, 1);
        tmp(2, 2) += storage[2][2] * rhs(2, 2);

        return tmp;
    }

    Vector6 Matrix3x3::VoigtVector() const
    {
        Vector6 tmp;
        tmp[0] = storage[0][0];
        tmp[1] = storage[1][1];
        tmp[2] = storage[2][2];
        tmp[3] = storage[1][2];
        tmp[4] = storage[0][2];
        tmp[5] = storage[0][1];
        return tmp;
    }
    vStress::vStress()
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
    };
    vStress::vStress(const Vector6& rhs)
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        storage[3] = rhs[3];
        storage[4] = rhs[4];
        storage[5] = rhs[5];
    };
    vStress& vStress::operator=(const Vector6& rhs)
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        storage[3] = rhs[3];
        storage[4] = rhs[4];
        storage[5] = rhs[5];
        return *this;
    };
    vStress vStress::operator-(const vStress& rhs) const
    {
        vStress tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    REAL vStress::norm(void) const  /// Frobenius norm
    {
        REAL tmp = storage[0] * storage[0]
            + storage[1] * storage[1]
            + storage[2] * storage[2]
            + storage[3] * storage[3] * 2
            + storage[4] * storage[4] * 2
            + storage[5] * storage[5] * 2;
        return std::sqrt(tmp);
    };

    REAL vStress::REALcontract(const vStress& Bstress) const  // "REAL-dot product"
    {
        REAL tmp =
            storage[0] * Bstress[0] +
            storage[1] * Bstress[1] +
            storage[2] * Bstress[2] +
            2 * storage[3] * Bstress[3] +
            2 * storage[4] * Bstress[4] +
            2 * storage[5] * Bstress[5];
        return tmp;
    };

    REAL vStress::REALcontract(const Vector6& symTensorV) const  // "REAL-dot product"
    {
        REAL tmp =
            storage[0] * symTensorV[0] +
            storage[1] * symTensorV[1] +
            storage[2] * symTensorV[2] +
            2 * storage[3] * symTensorV[3] +
            2 * storage[4] * symTensorV[4] +
            2 * storage[5] * symTensorV[5];
        return tmp;
    };

    inline REAL vStress::Pressure() const
    {
        return -(storage[0] + storage[1] + storage[2]) / 3;
    };

    inline REAL vStress::J1() const
    {
        return storage[0] + storage[1] + storage[2];
    };

    inline REAL vStress::Trace() const
    {
        return storage[0] + storage[1] + storage[2];
    };

    inline REAL vStress::Determinant() const
    {
        return storage[0] * storage[1] * storage[2] -
            storage[0] * storage[3] * storage[3] -
            storage[1] * storage[4] * storage[4] -
            storage[2] * storage[5] * storage[5] +
            storage[3] * storage[4] * storage[5] * 2;
    };

    Vector3 vStress::Invariants(void) const
    {
        Vector3 vec;
        vec[0] = Trace();
        vec[1] = REAL(0.5) * (Trace() * Trace() -
            (storage[0] * storage[0] +
                storage[1] * storage[1] +
                storage[2] * storage[2] +
                storage[3] * storage[3] * 2 +
                storage[4] * storage[4] * 2 +
                storage[5] * storage[5] * 2));
        vec[2] = Determinant();
        return vec;
    }

    REAL vStress::Mises(void) const
    {
        REAL vMises = 0.0;
        vMises = (storage[0] - storage[1]) * (storage[0] - storage[1]) +
            (storage[1] - storage[2]) * (storage[1] - storage[2]) +
            (storage[2] - storage[0]) * (storage[2] - storage[0]) +
            6 * (storage[3] * storage[3] + storage[4] * storage[4] + storage[5] * storage[5]);

        return std::sqrt(REAL(0.5) * vMises);
    }

    vStress vStress::get_rotated_matrix(const Matrix3x3& RotationMatrix) const
    {
        REAL In[3][3] = {};
        REAL Out[3][3] = {};

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                Out[p][q] = 0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                    }
            }
        vStress _OUT;

        _OUT[0] = Out[0][0];
        _OUT[5] = Out[0][1];
        _OUT[4] = Out[0][2];
        _OUT[1] = Out[1][1];
        _OUT[3] = Out[1][2];
        _OUT[2] = Out[2][2];

        return _OUT;
    };
    vStress& vStress::do_rotate(const Matrix3x3& RotationMatrix)
    {
        REAL In[3][3] = {};
        REAL Out[3][3] = {};

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                Out[p][q] = 0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                    }
            }

        storage[0] = Out[0][0];
        storage[5] = Out[0][1];
        storage[4] = Out[0][2];
        storage[1] = Out[1][1];
        storage[3] = Out[1][2];
        storage[2] = Out[2][2];

        return *this;
    };

    REAL vStress::get_tensor(const int i, const int j) const
    {
        return storage[(i == j) ? (i) : (6 - (i + j))];
    };

    Matrix3x3 vStress::tensor(void) const
    {
        Matrix3x3 tmp;
        tmp(0, 0) = storage[0];
        tmp(0, 1) = storage[5];
        tmp(0, 2) = storage[4];
        tmp(1, 0) = storage[5];
        tmp(1, 1) = storage[1];
        tmp(1, 2) = storage[3];
        tmp(2, 0) = storage[4];
        tmp(2, 1) = storage[3];
        tmp(2, 2) = storage[2];
        return tmp;
    };
    vStrain::vStrain()
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
    };
    vStrain::vStrain(const Vector6& rhs)
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        storage[3] = rhs[3];
        storage[4] = rhs[4];
        storage[5] = rhs[5];
    };
    vStrain& vStrain::operator=(const Vector6& rhs)
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        storage[3] = rhs[3];
        storage[4] = rhs[4];
        storage[5] = rhs[5];
        return *this;
    };
    bool vStrain::operator==(const Vector6& rhs)
    {
        for (int i = 0; i < 6; i++)
        {
            if (storage[i] != rhs[i]) { return false; };
        }
        return true;
    };
    vStrain vStrain::operator-(const vStrain& rhs) const
    {
        vStrain tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    REAL vStrain::norm(void) const  /// Frobenius norm
    {
        /*
         * multiplication by 0.5 of the off diagonal elements before taking
         * square and doubling the result due to REAL appearance of the off
         * diagonal elements => 0.5^2 / 2 = 0.5
         */
        REAL tmp = storage[0] * storage[0]
            + storage[1] * storage[1]
            + storage[2] * storage[2]
            + storage[3] * storage[3] / 2
            + storage[4] * storage[4] / 2
            + storage[5] * storage[5] / 2;
        return std::sqrt(tmp);
    };

    REAL vStrain::REALcontract(const vStrain& Bstrain) const // "REAL-dot product"
    {
        REAL tmp =
            storage[0] * Bstrain[0] +
            storage[1] * Bstrain[1] +
            storage[2] * Bstrain[2] +
            2 * (storage[3] / 2) * (Bstrain[3] / 2) +
            2 * (storage[4] / 2) * (Bstrain[4] / 2) +
            2 * (storage[5] / 2) * (Bstrain[5] / 2);
        return tmp;
    };

    vStrain vStrain::get_rotated_matrix(const Matrix3x3& RotationMatrix) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix and multiplying by two during translating back to
         * a Voigt strain vector
         */
        REAL In[3][3] = {};
        REAL Out[3][3] = {};

        In[0][0] = storage[0];
        In[0][1] = storage[5] / 2;
        In[0][2] = storage[4] / 2;
        In[1][0] = storage[5] / 2;
        In[1][1] = storage[1];
        In[1][2] = storage[3] / 2;
        In[2][0] = storage[4] / 2;
        In[2][1] = storage[3] / 2;
        In[2][2] = storage[2];

        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                Out[p][q] = 0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                    }
            }
        vStrain _OUT;

        _OUT[0] = Out[0][0];
        _OUT[5] = Out[0][1] * 2;
        _OUT[4] = Out[0][2] * 2;
        _OUT[1] = Out[1][1];
        _OUT[3] = Out[1][2] * 2;
        _OUT[2] = Out[2][2];

        return _OUT;
    };

    vStrain& vStrain::do_rotate(const Matrix3x3& RotationMatrix)
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix and multiplying by two during translating back to
         * a Voigt strain vector
         */
        REAL In[3][3] = {};
        REAL Out[3][3] = {};

        In[0][0] = storage[0];
        In[0][1] = storage[5] / 2;
        In[0][2] = storage[4] / 2;
        In[1][0] = storage[5] / 2;
        In[1][1] = storage[1];
        In[1][2] = storage[3] / 2;
        In[2][0] = storage[4] / 2;
        In[2][1] = storage[3] / 2;
        In[2][2] = storage[2];

        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                Out[p][q] = 0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                    }
            }

        storage[0] = Out[0][0];
        storage[5] = Out[0][1] * 2;
        storage[4] = Out[0][2] * 2;
        storage[1] = Out[1][1];
        storage[3] = Out[1][2] * 2;
        storage[2] = Out[2][2];

        return *this;
    };

    REAL vStrain::get_tensor(const int i, const int j) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix
         */
        return (storage[(i == j) ? (i) : (6 - (i + j))] * ((i == j) ? REAL(1.0) : REAL(0.5)));
    };

    Matrix3x3 vStrain::tensor(void) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix
         */
        Matrix3x3 tmp;
        tmp(0, 0) = storage[0];
        tmp(0, 1) = storage[5] / 2;
        tmp(0, 2) = storage[4] / 2;
        tmp(1, 0) = storage[5] / 2;
        tmp(1, 1) = storage[1];
        tmp(1, 2) = storage[3] / 2;
        tmp(2, 0) = storage[4] / 2;
        tmp(2, 1) = storage[3] / 2;
        tmp(2, 2) = storage[2];
        return tmp;
    };
    vStrain Matrix3x3::VoigtStrain() const
    {
        vStrain tmp;
        tmp[0] = storage[0][0];
        tmp[1] = storage[1][1];
        tmp[2] = storage[2][2];
        tmp[3] = storage[1][2] * 2;
        tmp[4] = storage[0][2] * 2;
        tmp[5] = storage[0][1] * 2;
        return tmp;
    }

    vStress Matrix3x3::VoigtStress() const
    {
        vStress tmp;
        tmp[0] = storage[0][0];
        tmp[1] = storage[1][1];
        tmp[2] = storage[2][2];
        tmp[3] = storage[1][2];
        tmp[4] = storage[0][2];
        tmp[5] = storage[0][1];
        return tmp;
    }

    vStrain Matrix6x6::operator*(const vStress& rhs) const
    {
        vStrain tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    }

    vStress Matrix6x6::operator*(const vStrain& rhs) const
    {
        vStress tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    }
}