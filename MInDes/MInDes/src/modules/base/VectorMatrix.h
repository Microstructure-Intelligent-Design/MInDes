#pragma once
#include "MACRO_DEF.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <vector>
namespace pf {
    class Vector3;
    class Matrix3x3;
    class Vector6;
    class Matrix6x6;
    class vStrain;
    class vStress;
    class Matrix3x3
    {
    public:

        Matrix3x3();

        Matrix3x3(const Matrix3x3& rhs);

        REAL& operator()(const int i, const int j);
        REAL const& operator()(const int i, const int j) const;
        REAL& operator()(const size_t i, const size_t j);
        REAL const& operator()(const size_t i, const size_t j) const;
        void set_to_zero(void);
        void set_to_unity(void);
        Matrix3x3& operator=(const Matrix3x3& rhs);
        Matrix3x3& operator=(const REAL rhs[3][3]);
        bool operator==(Matrix3x3& rhs);
        bool operator!=(Matrix3x3& rhs);

        Matrix3x3 operator*(REAL m);
        Matrix3x3 operator/(REAL m);
        Vector3 operator*(const Vector3& rhs) const;
        Matrix3x3 operator*(const Matrix3x3& rhs);
        Matrix3x3 operator+(Matrix3x3& rhs);
        Matrix3x3 operator-(Matrix3x3& rhs);

        Matrix6x6 outer(const Matrix3x3& rhs) const;

        Matrix3x3& operator+=(Matrix3x3& rhs);
        Matrix3x3& operator-=(Matrix3x3& rhs);
        Matrix3x3& operator*=(REAL m);
        Matrix3x3& operator/=(REAL m);

        inline REAL cal_determinant(void) const;
        Matrix3x3& do_invert(void);
        Matrix3x3 get_inverted_Matrix(void) const;
        Matrix3x3& do_transpose(void);
        Matrix3x3 get_transposed(void) const;
        Matrix3x3& do_rotate(Matrix3x3& RotationMatrix);
        Matrix3x3 get_rotated_Matrix(Matrix3x3& RotationMatrix) const;
        REAL REAL_contract(Matrix3x3& rHS) const;
        REAL trace(void) const;
        Matrix3x3 get_sym(void) const;
        Matrix3x3 get_skew(void) const;
        REAL norm(void) const;
        std::string print_matrix(void) const;

        inline Vector6 VoigtVector() const;
        inline vStrain  VoigtStrain() const;
        inline vStress  VoigtStress() const;

        REAL storage[3][3];
    };

    class Vector3
    {
    public:
        Vector3();
        Vector3(REAL i, REAL j, REAL k);

        Vector3(const Vector3& rhs);

        Vector3(const REAL vecinit[3]);

        REAL& operator[](const int i);
        REAL const& operator[](const int i) const;

        REAL& operator[](const size_t i);
        REAL const& operator[](const size_t i) const;

        REAL getX(void) const;

        void setX(const REAL newX);

        REAL getY(void) const;

        void setY(const REAL newY);

        REAL getZ(void) const;

        void setZ(const REAL newX);

        void set_to_zero(void);
        void set_to_unitX(void);
        void set_to_unitY(void);
        void set_to_unitZ(void);
        bool operator==(const Vector3& rhs) const;
        bool operator!=(const Vector3& rhs) const;
        Vector3 operator*(const REAL m) const;
        Vector3 operator/(const REAL m) const;
        REAL operator*(const Vector3& rhs) const;
        REAL abs() const;
        Vector3 cross(const Vector3& rhs) const;
        Vector3 operator+(const Vector3& rhs) const;
        Vector3 operator-(const Vector3& rhs) const;
        Vector3& operator*=(const REAL m);
        Vector3& operator/=(const REAL m);
        Vector3& operator-=(const Vector3& rhs);
        Vector3& operator+=(const Vector3& rhs);
        Vector3& operator=(const Vector3& rhs);
        Vector3& operator=(const REAL rhs[3]);
        REAL length(void) const;
        Vector3& normalize(void);
        Vector3 normalized(void) const;
        Vector3& do_rotate(const Matrix3x3& RotationMatrix);
        Vector3 get_rotated_vec3(const Matrix3x3& RotationMatrix) const;
        Vector3 Xreflected(void) const;
        Vector3 Yreflected(void) const;
        Vector3 Zreflected(void) const;
        std::string print(void) const;
        REAL* data(void);
        const REAL* const_data(void) const;
    protected:
    private:
        REAL storage[3];
    };

    class Vector6
    {
    public:

        Vector6();

        Vector6(REAL v1, REAL v2, REAL v3, REAL v4, REAL v5, REAL v6);

        Vector6(std::initializer_list<REAL> vecinit);

        REAL& operator[](const int i);
        REAL const& operator[](const int i) const;
        REAL& operator[](const size_t i);
        REAL const& operator[](const size_t i) const;
        bool is_nan_val_exist();
        void set_to_zero(void);
        void set_to_unity(void);
        void do_abs(void);

        REAL norm(void) const;
        REAL trace(void) const;
        REAL max_abs(void) const;
        Vector6 operator*(const REAL m) const;
        bool operator==(const Vector6& rhs) const;
        bool operator!=(const Vector6& rhs) const;
        Vector6 operator/(const REAL m) const;
        Vector6& operator*=(const REAL m);

        REAL operator*(const Vector6& rhs);

        Vector6& operator/=(const REAL m);
        Vector6 operator+(const Vector6& rhs) const;
        Vector6& operator+=(const Vector6& rhs);
        Vector6 operator-(const Vector6& rhs) const;
        Vector6& operator-=(const Vector6& rhs);
        Vector6& operator=(const Vector6& rhs);
        Vector6& operator=(const REAL rhs[6]);
        Vector6 get_rotated_vec6(const Matrix3x3& RotationMatrix) const;
        Vector6& do_rotate(Matrix3x3 RotationMatrix);
        REAL* data(void);
        const REAL* const_data(void) const;
        Vector6& H_product(Vector6 vector);

        Matrix3x3 sym_tensor(void) const;

        Matrix3x3 skew_tensor(void) const;

        std::string print(void) const;
    protected:
        REAL storage[6] = {};
    private:
    };

    class Matrix6x6
    {
    public:
        Matrix6x6();
        REAL& operator()(const int i, const int j);
        REAL const& operator()(const int i, const int j) const;
        REAL& operator()(const size_t i, const size_t j);
        REAL const& operator()(const size_t i, const size_t j) const;
        bool is_nan_val_exist();
        void set_to_zero(void);
        void set_to_unity(void);

        REAL norm(void) const;

        Matrix6x6 operator*(const REAL m) const;
        Matrix6x6 operator/(const REAL m) const;
        bool operator==(const Matrix6x6& rhs) const;
        bool operator!=(const Matrix6x6& rhs) const;
        Vector6 operator*(Vector6 rhs) const;
        Matrix6x6 operator*(const Matrix6x6& rhs) const;

        vStrain operator*(const vStress& rhs) const;

        vStress operator*(const vStrain& rhs) const;

        Matrix6x6 operator+(const Matrix6x6& rhs) const;
        Matrix6x6 operator-(const Matrix6x6& rhs) const;
        Matrix6x6& operator+=(const Matrix6x6& rhs);
        Matrix6x6& operator/=(const REAL m);
        Matrix6x6& operator*=(const REAL m);
        Matrix6x6& operator-=(const Matrix6x6& rhs);
        Matrix6x6& operator=(const Matrix6x6& rhs);
        Matrix6x6& operator=(const REAL rhs[6][6]);
        Matrix6x6& do_invert(void);

        Matrix6x6 get_inverted_matrix(void) const;

        Matrix6x6& do_transpose(void);

        Matrix6x6 get_transposed_matrix(void) const;

        Matrix6x6& do_rotate(const Matrix3x3& RotationMatrix);
        Matrix6x6 get_rotated_matrix(const Matrix3x3& RotationMatrix) const;

        Matrix6x6& do_rotate_Voigt(const Matrix3x3& RotationMatrix);
        Matrix6x6 get_rotated_matrix_Voigt(const Matrix3x3& RotationMatrix) const;

        std::string print(void) const;
        const REAL& tensor(const int i, const int j, const int k, const int l) const;
        const REAL& tensor(const size_t i, const size_t j, const size_t k, const size_t l) const;
        REAL* data(void);
        const REAL* const_data(void) const;
    protected:
    private:

        REAL storage[6][6];
    };

    class vStress : public Vector6
    {
        /*
         * Stress Voigt vector
         */
    public:
        vStress();
        vStress(const Vector6& rhs);
        vStress& operator=(const Vector6& rhs);
        vStress operator-(const vStress& rhs) const;
        REAL norm(void) const;

        REAL REALcontract(const vStress& Bstress) const;

        REAL REALcontract(const Vector6& symTensorV) const;

        REAL Pressure() const;

        REAL J1() const;

        REAL Trace() const;

        REAL Determinant() const;

        Vector3 Invariants(void) const;

        REAL Mises(void) const;

        vStress get_rotated_matrix(const Matrix3x3& RotationMatrix) const;
        vStress& do_rotate(const Matrix3x3& RotationMatrix);

        REAL get_tensor(const int i, const int j) const;
        REAL get_tensor(const size_t i, const size_t j) const;

        Matrix3x3 tensor(void) const;
    protected:
    private:
    };

    class vStrain : public Vector6
    {
        /*
         *  This is a special version of the six component Voigt vector:
         *  the off diagonal elements of the strain matrix are multiplied by 2(!)
         *  while transforming the 3x3 strain matrix to a 6-component vector,
         *  and have to be divided by two if they are transformed to a 3x3 matrix
         *  again. Therefore the transformation functions and some other functions
         *  have to be specified in different way then for the stress-type
         *  Voigt vector. The differences are specified in the functions.
         */
    public:
        vStrain();
        vStrain(const Vector6& rhs);
        vStrain& operator=(const Vector6& rhs);
        bool operator==(const Vector6& rhs);
        vStrain operator-(const vStrain& rhs) const;
        REAL norm(void) const;

        REAL REALcontract(const vStrain& Bstrain) const;

        vStrain get_rotated_matrix(const Matrix3x3& RotationMatrix) const;

        vStrain& do_rotate(const Matrix3x3& RotationMatrix);

        REAL get_tensor(const int i, const int j) const;
        REAL get_tensor(const size_t i, const size_t j) const;

        Matrix3x3 tensor(void) const;
    };

}