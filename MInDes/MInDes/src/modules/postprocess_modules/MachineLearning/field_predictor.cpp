#include "field_predictor.h"
#include <windows.h>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>
#include <iomanip>
namespace pf {

    namespace field_predictor {
        bool wait_for_state(SharedHeader* header, ShmState target, uint32_t timeout_ms) {
            auto start = std::chrono::steady_clock::now();

            while (true) {
                auto current = static_cast<ShmState>(header->state);
                if (current == target) {
                    return true;
                }
                if (current == ShmState::Error) {
                    return false;
                }

                auto now = std::chrono::steady_clock::now();
                auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
                if (elapsed_ms > timeout_ms) {
                    return false;
                }

                std::this_thread::sleep_for(std::chrono::milliseconds(20));
            }
        }

        void exec_i() {
            const std::wstring shm_name = L"Local\\SimML_Shared_01";

            ShmManagerWin shm;
            if (!shm.create_or_open(shm_name)) {
                return;
            }

            SharedHeader* header = shm.header();
            float* input = shm.input_data();
            float* output = shm.output_data();

            // C++ internal data uses double
            std::vector<double> solver_data = { 1.0, 2.0, 3.5, 4.25, 10.0 };
            if (solver_data.size() > kMaxElems) {
                std::cerr << "Too many elements.\n";
                return;
            }

            // Wait until shared memory is idle
            if (static_cast<ShmState>(header->state) != ShmState::Idle) {
                std::cout << "Shared memory not idle. Resetting to Idle for demo.\n";
                header->state = static_cast<uint32_t>(ShmState::Idle);
            }

            // Convert double -> float and write to shared memory
            header->input_count = static_cast<uint64_t>(solver_data.size());
            header->output_count = 0;
            header->request_id += 1;
            header->result_code = 0;

            for (size_t i = 0; i < solver_data.size(); ++i) {
                input[i] = static_cast<float>(solver_data[i]);
            }

            std::cout << "C++ wrote input: ";
            for (double v : solver_data) {
                std::cout << v << " ";
            }
            std::cout << "\n";

            // Publish request
            header->state = static_cast<uint32_t>(ShmState::InputReady);
            std::cout << "C++ set state = InputReady, request_id = " << header->request_id << "\n";
            std::cout << "Start python service now if not already running.\n";

            // Wait for Python to finish
            if (!wait_for_state(header, ShmState::OutputReady, 30000)) {
                auto current = static_cast<ShmState>(header->state);
                std::cerr << "Timeout or error while waiting for output. current state = " << static_cast<uint32_t>(current) << "\n";
                return;
            }

            if (header->result_code != 100) {
                std::cerr << "Python returned non-success result_code = " << header->result_code << "\n";
                header->state = static_cast<uint32_t>(ShmState::Idle);
                return;
            }

            std::vector<double> result;
            result.resize(static_cast<size_t>(header->output_count));

            for (size_t i = 0; i < result.size(); ++i) {
                result[i] = static_cast<double>(output[i]); // float -> double
            }

            std::cout << "C++ received output: ";
            for (double v : result) {
                std::cout << std::fixed << std::setprecision(3) << v << " ";
            }
            std::cout << "\n";

            // Reset for next round
            header->state = static_cast<uint32_t>(ShmState::Idle);
            std::cout << "C++ reset state = Idle\n";

            return;
        }
    }
}