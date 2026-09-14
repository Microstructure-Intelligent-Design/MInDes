#pragma once
#include "../base/RotationMatrix.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Modules_Params.h"
namespace pf {

    class GrainsOrientations
    {
    public:
        GrainsOrientations(const GrainsOrientations&) = delete;
        GrainsOrientations& operator=(const GrainsOrientations&) = delete;

        static GrainsOrientations& instance();

        void init();

        Vector3 get_phi_orientation(size_t phi_index) const
        {
            if (phi_index >= _grain_number) {
                return Vector3(0, 0, 0);
            }
            return orientations[phi_index];
        }
        Vector3 get_phi_orientation(int phi_index) const
        {
            if (phi_index >= _grain_number) {
                return Vector3(0, 0, 0);
            }
            return orientations[phi_index];
        }
        Matrix3x3 get_phi_rotationMatrix(size_t phi_index) const
        {
            if (phi_index >= _grain_number) {
                Matrix3x3 unity;
                unity.set_to_unity();
                return unity;
            }
            return RotationMatrix::rotationMatrix(orientations[phi_index], _rotation_gauge);
        }
        Matrix3x3 get_phi_rotationMatrix(int phi_index) const
        {
            if (phi_index >= _grain_number) {
                Matrix3x3 unity;
                unity.set_to_unity();
                return unity;
            }
            return RotationMatrix::rotationMatrix(orientations[phi_index], _rotation_gauge);
        }
        // return RotationMatrix directly
        Matrix3x3 RotationMatrix(size_t phi_index) {
            return RotationMatrix::rotationMatrix(orientations[phi_index], _rotation_gauge);
        }
        // return RotationMatrix directly
        Matrix3x3 RotationMatrix(int phi_index) {
            return RotationMatrix::rotationMatrix(orientations[phi_index], _rotation_gauge);
        }
    private:
        GrainsOrientations();
        ~GrainsOrientations() = default;
        bool _is_init = false;
        size_t _grain_number = 0;
        RotationGauge _rotation_gauge = RotationGauge::RG_XZX;
        std::vector<Vector3> orientations;
    };

}