#include "FluidField.h"
namespace pf {
	namespace fluid_field {
		// statement functons
		void init() {
			lattice_boltzmann::init();
			WriteLog("> MODULE INIT : FluidField !\n");
		}
		void exec_pre() {
			lattice_boltzmann::exec_pre();
		}
		std::string exec_loop() {
			std::string report = "";
			report = lattice_boltzmann::exec_loop();
			return report;
		}
		void deinit() {
			lattice_boltzmann::deinit();
		}
	}
}