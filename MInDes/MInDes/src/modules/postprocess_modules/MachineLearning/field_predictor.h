#include "FieldPredictor/shm_manager_win.h"
namespace pf {
    namespace field_predictor {
        bool wait_for_state(SharedHeader* header, ShmState target, uint32_t timeout_ms);

        void exec_i();
    }
}