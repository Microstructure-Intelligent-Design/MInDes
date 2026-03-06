#include "GrainsOrientations.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Modules_Params.h"
namespace pf {

    GrainsOrientations::GrainsOrientations()
        : _is_init(false)
        , _grain_number(0)
        , _rotation_gauge(RotationGauge::RG_XZX) {
    };

    void GrainsOrientations::init() {
		if (!_is_init && main_field::is_phi_field_on && main_field::phi_number > 0) {
			// init grain orientation
			WriteDebugFile("# Model.GrainsOrientations.rotation_gauge = 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
			WriteDebugFile("#                                            6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
			int rotation_gauge = RotationGauge::RG_ZXZ;
			infile_reader::read_int_value("Model.GrainsOrientations.rotation_gauge", rotation_gauge, true);
			_rotation_gauge = RotationGauge(rotation_gauge);
			_grain_number = main_field::phi_number;
			orientations.clear();
			orientations.resize(main_field::phi_number, Vector3(0, 0, 0));
			WriteDebugFile("# Model.GrainsOrientations = {[(phi_index_0, phi_index_2, ... ),(rotation_angle_1, rotation_angle_2, rotation_angle_3)],  ... } \n");
			string grains_key = "Model.GrainsOrientations", grains_input = "{}";
			infile_reader::read_string_value(grains_key, grains_input, true);
			std::vector<std::vector<std::vector<input_value>>> grains_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value
			({ InputValueType::IVType_INT , InputValueType::IVType_REAL }, grains_key, grains_input, true);
			for (size_t i = 0; i < grains_value.size(); i++)
				for (size_t j = 0; j < grains_value[i][0].size(); j++) {
					int phi_index = grains_value[i][0][j].int_value;
					if (phi_index < _grain_number)
						orientations[phi_index] = Vector3(AngleToRadians(grains_value[i][1][0].REAL_value),
							AngleToRadians(grains_value[i][1][1].REAL_value), AngleToRadians(grains_value[i][1][2].REAL_value));
				}
            _is_init = true;
        }
    }

    GrainsOrientations& GrainsOrientations::instance()
    {
        static GrainsOrientations inst;
        return inst;
    };

}