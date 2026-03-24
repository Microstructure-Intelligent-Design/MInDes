#include "MechanicalField.h"
namespace pf {
	namespace mechanical_field {
		void init() {
			PhiProperties::instance().init();
			GrainsOrientations::instance().init();
			stiffness_eigenstrain::init(main_field::phi_number, PhiProperties::instance().phi_property_number());
			elastic_solver::init();
			plastic_solver::init();
			WriteLog("> MODULE INIT : MechanicalField !\n");
		}
		void exec_pre_iii() {

		}
		std::string exec_loop_iii() {
			std::string mech_report = "";

			return mech_report;
		}
		void deinit() {
			elastic_solver::deinit();
			plastic_solver::deinit();
		}
	}
}