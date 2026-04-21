#include "MachineLearning.h"
#include "MachineLearning/field_predictor.h"
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
namespace pf {
	namespace machine_learning {
		void init_machine_learning() {
			WriteDebugFile("# Postprocess.PCT.MachineLearning.field = (is_phi, is_con, is_temp) \n");
			std::string field_key = "Postprocess.PCT.MachineLearning.field", field_in = "(false,false,false)";
			if (infile_reader::read_string_value(field_key, field_in, true)) {
				std::vector<input_value> field_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(
					InputValueType::IVType_BOOL, field_key, field_in, true);
				bool is_phi_acc = field_value[0].bool_value, is_con_acc = field_value[1].bool_value,
					is_temp_acc = field_value[2].bool_value;
				if (is_phi_acc || is_con_acc || is_temp_acc) {
					load_a_new_module(field_predictor::exec_i, nullptr, nullptr, nullptr, nullptr, nullptr,
						nullptr, nullptr, nullptr, nullptr);
				}
			}
		}
	}
}