#pragma once
#include <windows.h>
#include <cstdint>
#include <string>
#include "shm_layout.h"
namespace pf {
    class ShmManagerWin {
    public:
        ShmManagerWin() = default;
        ~ShmManagerWin();

        bool create_or_open(const std::wstring& name);
        void close();

        SharedHeader* header() const { return header_; }
        float* input_data() const { return input_; }
        float* output_data() const { return output_; }
        void* base() const { return base_; }

    private:
        HANDLE mapping_ = nullptr;
        void* base_ = nullptr;
        SharedHeader* header_ = nullptr;
        float* input_ = nullptr;
        float* output_ = nullptr;
    };
}