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
    // ---------------------------------------------------------------------------
    // 1D 矩阵封装 (行优先存储)
    // ---------------------------------------------------------------------------
    template <typename T>
    class Matrix1D {
    private:
        std::vector<T> data_;
        size_t length_;
    public:
        Matrix1D() : length_(0) {}
        Matrix1D(size_t length, T value = T{})
            : data_(length, value), length_(length) {
        }
        void resize(size_t length, T value = T{}) {
            length_ = length;
            data_.assign(length, value);
        }
        void resize(int length, T value = T{}) {
            length_ = length;
            data_.assign(size_t(length), value);
        }
        inline Matrix1D& operator=(const Matrix1D& rhs) {
            length_ = rhs.length_;
            data_ = rhs.data_;
            return *this;
        }
        inline T& operator()(size_t i) {
            return data_[i];
        }
        inline T& operator()(int i) {
            return data_[i];
        }
        inline const T& operator()(size_t i) const {
            return data_[i];
        }
        inline const T& operator()(int i) const {
            return data_[i];
        }
        inline T& operator[](size_t i) {
            return data_[i];
        }
        inline T& operator[](int i) {
            return data_[i];
        }
        inline const T& operator[](size_t i) const {
            return data_[i];
        }
        inline const T& operator[](int i) const {
            return data_[i];
        }
        // 矩阵加法 (+)
        Matrix1D operator+(const Matrix1D& other) const {
            if (length_ != other.length_) {
                throw std::invalid_argument("Matrix1D + error: lengths do not match!");
            }
            Matrix1D result(length_);
            for (size_t i = 0; i < length_; ++i)
                result[i] = data_[i] + other[i];
            return result;
        }
        // 矩阵减法 (-)
        Matrix1D operator-(const Matrix1D& other) const {
            if (length_ != other.length_) {
                throw std::invalid_argument("Matrix1D - error: lengths do not match!");
            }
            Matrix1D result(length_);
            for (size_t i = 0; i < length_; ++i)
                result[i] = data_[i] - other[i];
            return result;
        }

        // 矩阵除法 (*) - 逐元素相乘
        Matrix1D operator*(const REAL& other) const {
            Matrix1D result(length_);
            for (size_t i = 0; i < length_; ++i)
                result[i] = data_[i] * other;
            return result;
        }

        // 矩阵除法 (/) - 逐元素相除
        Matrix1D operator/(const REAL& other) const {
            Matrix1D result(length_);
            for (size_t i = 0; i < length_; ++i)
                result[i] = data_[i] / other;
            return result;
        }

        // 复合赋值运算符 (+=)，可以提升连续运算时的性能
        Matrix1D& operator+=(const Matrix1D& other) {
            if (length_ != other.length_) {
                throw std::invalid_argument("Matrix1D += error: lengths do not match!");
            }
            for (size_t i = 0; i < length_; ++i) {
                data_[i] += other.data_[i];
            }
            return *this;
        }
        // 复合赋值运算符 (-=)，可以提升连续运算时的性能
        Matrix1D& operator-=(const Matrix1D& other) {
            if (length_ != other.length_) {
                throw std::invalid_argument("Matrix1D -= error: lengths do not match!");
            }
            for (size_t i = 0; i < length_; ++i) {
                data_[i] -= other.data_[i];
            }
            return *this;
        }

        // 复合赋值运算符 (*=)，可以提升连续运算时的性能
        Matrix1D& operator*=(const REAL& other) {
            for (size_t i = 0; i < length_; ++i) {
                data_[i] *= other;
            }
            return *this;
        }

        // 复合赋值运算符 (/=)，可以提升连续运算时的性能
        Matrix1D& operator/=(const REAL& other) {
            for (size_t i = 0; i < length_; ++i) {
                data_[i] /= other;
            }
            return *this;
        }
        size_t size() const { return length_; }
        T* data() { return data_.data(); }
        const T* data() const { return data_.data(); }
    };
    // ---------------------------------------------------------------------------
    // 2D 矩阵封装 (行优先存储)
    // ---------------------------------------------------------------------------
    template <typename T>
    class Matrix2D {
    private:
        std::vector<T> data_;
        size_t rows_;
        size_t cols_;
    public:
        Matrix2D() : rows_(0), cols_(0) {}
        Matrix2D(size_t rows, size_t cols, T value = T{})
            : data_(rows* cols, value), rows_(rows), cols_(cols) {
        }
        void resize(size_t rows, size_t cols, T value = T{}) {
            rows_ = rows;
            cols_ = cols;
            data_.assign(rows * cols, value);
        }
        void resize(int rows, int cols, T value = T{}) {
            rows_ = rows;
            cols_ = cols;
            data_.assign(size_t(rows * cols), value);
        }
        inline Matrix2D& operator=(const Matrix2D& rhs) {
            rows_ = rhs.rows_;
            cols_ = rhs.cols_;
            data_ = rhs.data_;
            return *this;
        }
        inline T& operator()(size_t i, size_t j) {
            return data_[i * cols_ + j];
        }
        inline T& operator()(int i, int j) {
            return data_[i * cols_ + j];
        }
        inline const T& operator()(size_t i, size_t j) const {
            return data_[i * cols_ + j];
        }
        inline const T& operator()(int i, int j) const {
            return data_[i * cols_ + j];
        }
        size_t rows() const { return rows_; }
        size_t cols() const { return cols_; }
        T* data() { return data_.data(); }
        const T* data() const { return data_.data(); }
    };
    // ---------------------------------------------------------------------------
    // 3D 矩阵封装 (行优先存储)
    // 内存布局：[z][y][x] -> 线性索引 = z * (Y*X) + y * X + x
    // ---------------------------------------------------------------------------
    template <typename T>
    class Matrix3D {
    private:
        std::vector<T> data_;
        size_t dim_x_;
        size_t dim_y_;
        size_t dim_z_;
    public:
        Matrix3D() : dim_x_(0), dim_y_(0), dim_z_(0) {}
        Matrix3D(size_t x, size_t y, size_t z, T value = T{})
            : data_(x* y* z, value), dim_x_(x), dim_y_(y), dim_z_(z) {
        }
        void resize(size_t x, size_t y, size_t z, T value = T{}) {
            dim_x_ = x;
            dim_y_ = y;
            dim_z_ = z;
            data_.assign(x * y * z, value);
        }
        void resize(int x, int y, int z, T value = T{}) {
            dim_x_ = x;
            dim_y_ = y;
            dim_z_ = z;
            data_.assign(x * y * z, value);
        }
        inline T& operator()(size_t i, size_t j, size_t k) {
            return data_[k * (dim_y_ * dim_x_) + j * dim_x_ + i];
        }
        inline T& operator()(int i, int j, int k) {
            return data_[k * (dim_y_ * dim_x_) + j * dim_x_ + i];
        }
        inline const T& operator()(size_t i, size_t j, size_t k) const {
            return data_[k * (dim_y_ * dim_x_) + j * dim_x_ + i];
        }
        inline const T& operator()(int i, int j, int k) const {
            return data_[k * (dim_y_ * dim_x_) + j * dim_x_ + i];
        }
        inline Matrix3D& operator=(const Matrix3D& rhs) {
            dim_x_ = rhs.dim_x_;
            dim_y_ = rhs.dim_y_;
            dim_z_ = rhs.dim_z_;
            data_ = rhs.data_;
            return *this;
        }
        size_t x() const { return dim_x_; }
        size_t y() const { return dim_y_; }
        size_t z() const { return dim_z_; }
    };
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