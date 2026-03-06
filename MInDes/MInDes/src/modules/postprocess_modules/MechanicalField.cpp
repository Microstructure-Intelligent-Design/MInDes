#include "../input_modules/ioFiles_Params.h"
#include "../Modules_Params.h"
#include "WriteVTS.h"
#include "MechanicalField.h"
#include "Mechanics/StiffnessEigenStrain.h"
#include "Mechanics/ElasticSolver.h"
#include "Mechanics/PlasticSolver.h"
#include "../model_modules/GrainsOrientations.h"
#include "../model_modules/PhiProperties.h"
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
			string mech_report = "";

			return mech_report;
		}
		void deinit() {
			elastic_solver::deinit();
			plastic_solver::deinit();
		}
	}
}