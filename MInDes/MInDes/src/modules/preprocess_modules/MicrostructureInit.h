#pragma once
#include "../input_modules/ioFiles_Params.h"
#include "../input_modules/inputfiles/selectFile.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../base/RotationMatrix.h"
#include "../../MainIterator_Params.h"
#include "../Modules_Params.h"
#include "../Module.h"
#include "MicrostructureInit/WriteMeshData.h"
#include "MicrostructureInit/GeometryStructure.h"
#include "MicrostructureInit/Bmp24Structure.h"
#include "MicrostructureInit/PorousStructure.h"
#include "MicrostructureInit/VoronoiStructure.h"
namespace pf {
	namespace microstructure_init {
		inline std::vector<size_t> porous_phis_indexs;
		inline std::vector<size_t> voronoi_phis_indexs;
		namespace functions {
			void init_mesh_with_datafile(pf::Data_MeshInfo& report, std::string dat_path);
			void definiteNucleation(size_t Current_ITE_step);
		}
		void check_phi_index(size_t phi_index);
		void check_con_index(size_t con_index);
		void check_con_size(size_t con_size);
		// -
		void init_microstructure_pre_i();
		// - 
		void write_data_pre_iii();
		void write_data_pos_iii();
		void dinit();
		// -
		void init_microstructure();
	}
}