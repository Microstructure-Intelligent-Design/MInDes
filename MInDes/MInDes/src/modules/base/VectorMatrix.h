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


    struct vec3_elem {
        int index;
        Vector3 vec;
        vec3_elem() {
            index = 0;
            vec[0] = 0.0;
            vec[1] = 0.0;
            vec[2] = 0.0;
        }
        vec3_elem& operator=(const vec3_elem& n) {
            index = n.index;
            vec[0] = n.vec[0];
            vec[1] = n.vec[1];
            vec[2] = n.vec[2];
            return *this;
        }
    };

    class vec3_box {
    public:
        vec3_box() {
            _vec_box.reserve(0);
        }
        ~vec3_box() {
            clear();
        }
        std::vector<vec3_elem> _vec_box;
        typedef std::vector<vec3_elem>::iterator iterator;
        typedef std::vector<vec3_elem>::const_iterator citerator;
        iterator  begin() { return _vec_box.begin(); };
        iterator  end() { return _vec_box.end(); };
        Vector3& operator[](const int index) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                if (i->index == index) return i->vec;
            }
            std::cout << "vec3_box error, can't find the vec index : " << index << std::endl;
            SYS_PROGRAM_STOP;
        }
        REAL& operator()(const int index, const int index2) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                if (i->index == index) return i->vec[index2];
            }
            std::cout << "vec3_box error, can't find the vec index : " << index << std::endl;
            SYS_PROGRAM_STOP;
        }
        vec3_box& operator=(const vec3_box& n) {
            _vec_box = n._vec_box;
            return *this;
        }
        void add_vec(int _index, REAL* vec) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = vec[0];
                    i->vec[1] = vec[1];
                    i->vec[2] = vec[2];
                    return;
                }
            _vec_box.reserve(_vec_box.size() + 1);
            vec3_elem elem;
            elem.index = _index;
            elem.vec[0] = vec[0];
            elem.vec[1] = vec[1];
            elem.vec[2] = vec[2];
            _vec_box.push_back(elem);
        }
        void add_vec(int _index, Vector3 vec) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = vec[0];
                    i->vec[1] = vec[1];
                    i->vec[2] = vec[2];
                    return;
                }
            _vec_box.reserve(_vec_box.size() + 1);
            vec3_elem elem;
            elem.index = _index;
            elem.vec[0] = vec[0];
            elem.vec[1] = vec[1];
            elem.vec[2] = vec[2];
            _vec_box.push_back(elem);
        }
        void erase(int index) {
            for (auto i = _vec_box.begin(); i < _vec_box.end();) {
                if (i->index == index) {
                    i = _vec_box.erase(i);
                }
                else
                    ++i;
            }
        }
        void clear() {
            _vec_box.clear();
        }
        int size() const {
            return int(_vec_box.size());
        }
    };

    class Vector6
    {
    public:

        Vector6();

        Vector6(REAL v1, REAL v2, REAL v3, REAL v4, REAL v5, REAL v6);

        Vector6(std::initializer_list<REAL> vecinit);

        REAL& operator[](const int i);
        REAL const& operator[](const int i) const;
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

    struct vec6_elem {
        int index;
        Vector6 vec;
        vec6_elem() {
            index = 0;
            vec[0] = 0.0;
            vec[1] = 0.0;
            vec[2] = 0.0;
            vec[3] = 0.0;
            vec[4] = 0.0;
            vec[5] = 0.0;
        }
        vec6_elem& operator=(const vec6_elem& n) {
            index = n.index;
            vec[0] = n.vec[0];
            vec[1] = n.vec[1];
            vec[2] = n.vec[2];
            vec[3] = n.vec[3];
            vec[4] = n.vec[4];
            vec[5] = n.vec[5];
            return *this;
        }
    };

    class vec6_box {
    public:
        vec6_box() {
            _vec_box.reserve(0);
        }
        ~vec6_box() {
            clear();
        }
        std::vector<vec6_elem> _vec_box;
        typedef std::vector<vec6_elem>::iterator iterator;
        typedef std::vector<vec6_elem>::const_iterator citerator;
        iterator  begin() { return _vec_box.begin(); };
        iterator  end() { return _vec_box.end(); };
        Vector6& operator[](const int index) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                if (i->index == index) return i->vec;
            }
            std::cout << "vec6_box error, can't find the vec index : " << index << std::endl;
            SYS_PROGRAM_STOP;
        }
        REAL& operator()(const int index, const int index2) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                if (i->index == index) return i->vec[index2];
            }
            std::cout << "vec6_box error, can't find the vec index : " << index << std::endl;
            SYS_PROGRAM_STOP;
        }
        vec6_box& operator=(const vec6_box& n) {
            _vec_box = n._vec_box;
            return *this;
        }
        void add_vec(int _index, REAL* vec) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = vec[0];
                    i->vec[1] = vec[1];
                    i->vec[2] = vec[2];
                    i->vec[3] = vec[3];
                    i->vec[4] = vec[4];
                    i->vec[5] = vec[5];
                    return;
                }
            _vec_box.reserve(_vec_box.size() + 1);
            vec6_elem elem;
            elem.index = _index;
            elem.vec[0] = vec[0];
            elem.vec[1] = vec[1];
            elem.vec[2] = vec[2];
            elem.vec[3] = vec[3];
            elem.vec[4] = vec[4];
            elem.vec[5] = vec[5];
            _vec_box.push_back(elem);
        }
        void add_vec(int _index, Vector6 vec) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = vec[0];
                    i->vec[1] = vec[1];
                    i->vec[2] = vec[2];
                    i->vec[3] = vec[3];
                    i->vec[4] = vec[4];
                    i->vec[5] = vec[5];
                    return;
                }
            _vec_box.reserve(_vec_box.size() + 1);
            vec6_elem elem;
            elem.index = _index;
            elem.vec[0] = vec[0];
            elem.vec[1] = vec[1];
            elem.vec[2] = vec[2];
            elem.vec[3] = vec[3];
            elem.vec[4] = vec[4];
            elem.vec[5] = vec[5];
            _vec_box.push_back(elem);
        }
        void erase(int index) {
            for (auto i = _vec_box.begin(); i < _vec_box.end();) {
                if (i->index == index) {
                    i = _vec_box.erase(i);
                }
                else
                    ++i;
            }
        }
        void clear() {
            _vec_box.clear();
        }
        int size() const {
            return int(_vec_box.size());
        }
    };

    class Matrix6x6
    {
    public:
        Matrix6x6();
        REAL& operator()(const int i, const int j);
        bool is_nan_val_exist();
        REAL const& operator()(const int i, const int j) const;
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
        REAL* data(void);
        const REAL* const_data(void) const;
    protected:
    private:

        REAL storage[6][6];
    };

    struct matrix6x6_elem {
        int index;
        Matrix6x6 matrix;
        matrix6x6_elem() {
            index = 0;
            matrix(0, 0) = 0.0; matrix(0, 1) = 0.0; matrix(0, 2) = 0.0; matrix(0, 3) = 0.0; matrix(0, 4) = 0.0; matrix(0, 5) = 0.0;
            matrix(1, 0) = 0.0; matrix(1, 1) = 0.0; matrix(1, 2) = 0.0; matrix(1, 3) = 0.0; matrix(1, 4) = 0.0; matrix(1, 5) = 0.0;
            matrix(2, 0) = 0.0; matrix(2, 1) = 0.0; matrix(2, 2) = 0.0; matrix(2, 3) = 0.0; matrix(2, 4) = 0.0; matrix(2, 5) = 0.0;
            matrix(3, 0) = 0.0; matrix(3, 1) = 0.0; matrix(3, 2) = 0.0; matrix(3, 3) = 0.0; matrix(3, 4) = 0.0; matrix(3, 5) = 0.0;
            matrix(4, 0) = 0.0; matrix(4, 1) = 0.0; matrix(4, 2) = 0.0; matrix(4, 3) = 0.0; matrix(4, 4) = 0.0; matrix(4, 5) = 0.0;
            matrix(5, 0) = 0.0; matrix(5, 1) = 0.0; matrix(5, 2) = 0.0; matrix(5, 3) = 0.0; matrix(5, 4) = 0.0; matrix(5, 5) = 0.0;
        }
        matrix6x6_elem& operator=(const matrix6x6_elem& n) {
            index = n.index;
            matrix(0, 0) = n.matrix(0, 0); matrix(0, 1) = n.matrix(0, 1); matrix(0, 2) = n.matrix(0, 2); matrix(0, 3) = n.matrix(0, 3); matrix(0, 4) = n.matrix(0, 4); matrix(0, 5) = n.matrix(0, 5);
            matrix(1, 0) = n.matrix(1, 0); matrix(1, 1) = n.matrix(1, 1); matrix(1, 2) = n.matrix(1, 2); matrix(1, 3) = n.matrix(1, 3); matrix(1, 4) = n.matrix(1, 4); matrix(1, 5) = n.matrix(1, 5);
            matrix(2, 0) = n.matrix(2, 0); matrix(2, 1) = n.matrix(2, 1); matrix(2, 2) = n.matrix(2, 2); matrix(2, 3) = n.matrix(2, 3); matrix(2, 4) = n.matrix(2, 4); matrix(2, 5) = n.matrix(2, 5);
            matrix(3, 0) = n.matrix(3, 0); matrix(3, 1) = n.matrix(3, 1); matrix(3, 2) = n.matrix(3, 2); matrix(3, 3) = n.matrix(3, 3); matrix(3, 4) = n.matrix(3, 4); matrix(3, 5) = n.matrix(3, 5);
            matrix(4, 0) = n.matrix(4, 0); matrix(4, 1) = n.matrix(4, 1); matrix(4, 2) = n.matrix(4, 2); matrix(4, 3) = n.matrix(4, 3); matrix(4, 4) = n.matrix(4, 4); matrix(4, 5) = n.matrix(4, 5);
            matrix(5, 0) = n.matrix(5, 0); matrix(5, 1) = n.matrix(5, 1); matrix(5, 2) = n.matrix(5, 2); matrix(5, 3) = n.matrix(5, 3); matrix(5, 4) = n.matrix(5, 4); matrix(5, 5) = n.matrix(5, 5);
            return *this;
        }
    };

    class matrix6x6_box {
    public:
        matrix6x6_box() {
            _matrix_box.reserve(0);
        }
        ~matrix6x6_box() {
            clear();
        }
        std::vector<matrix6x6_elem> _matrix_box;
        typedef std::vector<matrix6x6_elem>::iterator iterator;
        typedef std::vector<matrix6x6_elem>::const_iterator citerator;
        iterator  begin() { return _matrix_box.begin(); };
        iterator  end() { return _matrix_box.end(); };
        Matrix6x6& operator[](const int index) {
            for (auto i = _matrix_box.begin(); i < _matrix_box.end(); ++i) {
                if (i->index == index) return i->matrix;
            }
            std::cout << "Matrix6x6 error, can't find the vec index : " << index << std::endl;
            SYS_PROGRAM_STOP;
        }
        REAL& operator()(const int index, const int index_i, const int index_j) {
            for (auto i = _matrix_box.begin(); i < _matrix_box.end(); ++i) {
                if (i->index == index) return i->matrix(index_i, index_j);
            }
            std::cout << "Matrix6x6 error, can't find the vec index : " << index << std::endl;
            SYS_PROGRAM_STOP;
        }
        matrix6x6_box& operator=(const matrix6x6_box& n) {
            _matrix_box = n._matrix_box;
            return *this;
        }
        void add_matrix(int _index, Matrix6x6 _matrix) {
            for (auto i = _matrix_box.begin(); i < _matrix_box.end(); ++i)
                if (i->index == _index) {
                    for (int ii = 0; ii < 6; ii++)
                        for (int jj = 0; jj < 6; jj++)
                            i->matrix(ii, jj) = _matrix(ii, jj);
                    return;
                }
            _matrix_box.reserve(_matrix_box.size() + 1);
            matrix6x6_elem elem;
            elem.index = _index;
            for (int ii = 0; ii < 6; ii++)
                for (int jj = 0; jj < 6; jj++)
                    elem.matrix(ii, jj) = _matrix(ii, jj);
            _matrix_box.push_back(elem);
        }
        void erase(int index) {
            for (auto i = _matrix_box.begin(); i < _matrix_box.end();) {
                if (i->index == index) {
                    i = _matrix_box.erase(i);
                }
                else
                    ++i;
            }
        }
        void clear() {
            _matrix_box.clear();
        }
        int size() const {
            return int(_matrix_box.size());
        }
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

        inline REAL Pressure() const;

        inline REAL J1() const;

        inline REAL Trace() const;

        inline REAL Determinant() const;

        Vector3 Invariants(void) const;

        REAL Mises(void) const;

        vStress get_rotated_matrix(const Matrix3x3& RotationMatrix) const;
        vStress& do_rotate(const Matrix3x3& RotationMatrix);

        REAL get_tensor(const int i, const int j) const;

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

        Matrix3x3 tensor(void) const;
    };

}