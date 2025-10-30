#include "Bmp24Structure.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
namespace pf {
	namespace bmp24_structure {

		void generate_structure_from_BMP_pic(size_t MESH_NX, size_t MESH_NY, size_t MESH_NZ,
			bool is_phi_field_on, bool is_con_field_on, bool is_temp_field_on) {
			if (MESH_NZ > 1) {
				std::cout << "> Warnning : generate structure error : Nz > 1, cant init structure by a BMP picture !" << std::endl;
				return;
			}
			// using stb_image read BMP
			int img_width, img_height, channels;
			unsigned char* image_data = stbi_load(bmp24file_path.c_str(), &img_width, &img_height, &channels, 3); //  3 channel (RGB)

			if (!image_data) {
				std::cout << "> Error : Failed to load BMP image: " << bmp24file_path << std::endl;
				return;
			}
			// make sure 24 bit (3 channel)
			if (channels != 3) {
				std::cout << "> Error : BMP must be 24-bit (3 channels), but got " << std::to_string(channels) << " channels." << std::endl;
				stbi_image_free(image_data);
				return;
			}

			for (size_t layer_index = 0; layer_index < bmp24_layer; layer_index++) {
				geometry_structure::PointSet set;
				for (int y = 0; y < MESH_NY; y++)
					for (int x = 0; x < MESH_NX; x++) {
						int xx = x * img_width / int(MESH_NX),
							yy = y * img_height / int(MESH_NY);
						xx = std::max(0, std::min(xx, img_width - 1));
						yy = std::max(0, std::min(yy, img_height - 1));

						int pixel_index = (yy * img_width + xx) * 3;
						unsigned char r = image_data[pixel_index + 0];
						unsigned char g = image_data[pixel_index + 1];
						unsigned char b = image_data[pixel_index + 2];
						double graypercent = (r * 0.299 + g * 0.587 + b * 0.114) / 255.0;

						if (graypercent > bmp24_threshold[layer_index][0] &&
							graypercent < bmp24_threshold[layer_index][1]) {
							set.add_point(x + 1, y + 1, 1, bmp24_phi_value[layer_index]);
						}
					}
				set.generate_step = 0;
				if (is_phi_field_on) {
					set.phaseIndex = bmp24_phi_index[layer_index];
					set.is_normalized = bmp24_phi_normalized[layer_index];
				}
				if (is_con_field_on)
					set.con = bmp24_con[layer_index];
				if (is_temp_field_on)
					set.temperature = bmp24_temperature[layer_index];
				geometry_structure::nucleation_box.point_set_box.push_back(set);
			}
			// free
			stbi_image_free(image_data);
		}
	}
}