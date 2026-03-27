#include "shm_manager_win.h"
#include <iostream>

namespace pf {
    ShmManagerWin::~ShmManagerWin() {
        close();
    }

    bool ShmManagerWin::create_or_open(const std::wstring& name) {
        mapping_ = CreateFileMappingW(
            INVALID_HANDLE_VALUE,
            nullptr,
            PAGE_READWRITE,
            static_cast<DWORD>((kShmTotalBytes >> 32) & 0xFFFFFFFF),
            static_cast<DWORD>(kShmTotalBytes & 0xFFFFFFFF),
            name.c_str()
        );

        if (!mapping_) {
            std::cerr << "CreateFileMappingW failed: " << GetLastError() << "\n";
            return false;
        }

        base_ = MapViewOfFile(
            mapping_,
            FILE_MAP_ALL_ACCESS,
            0, 0,
            kShmTotalBytes
        );

        if (!base_) {
            std::cerr << "MapViewOfFile failed: " << GetLastError() << "\n";
            CloseHandle(mapping_);
            mapping_ = nullptr;
            return false;
        }

        header_ = reinterpret_cast<SharedHeader*>(base_);
        input_ = reinterpret_cast<float*>(reinterpret_cast<std::uint8_t*>(base_) + input_offset_bytes());
        output_ = reinterpret_cast<float*>(reinterpret_cast<std::uint8_t*>(base_) + output_offset_bytes());

        // Always initialize for this simple prototype.
        header_->magic = kMagic;
        header_->version = kVersion;
        header_->state = static_cast<uint32_t>(ShmState::Idle);
        header_->input_dtype = static_cast<uint32_t>(ShmDType::Float32);
        header_->output_dtype = static_cast<uint32_t>(ShmDType::Float32);
        header_->input_count = 0;
        header_->output_count = 0;
        header_->input_offset = input_offset_bytes();
        header_->output_offset = output_offset_bytes();
        header_->request_id = 0;
        header_->result_code = 0;

        ZeroMemory(header_->reserved, sizeof(header_->reserved));
        ZeroMemory(input_, kInputBytes);
        ZeroMemory(output_, kOutputBytes);

        return true;
    }

    void ShmManagerWin::close() {
        if (base_) {
            UnmapViewOfFile(base_);
            base_ = nullptr;
        }
        if (mapping_) {
            CloseHandle(mapping_);
            mapping_ = nullptr;
        }
        header_ = nullptr;
        input_ = nullptr;
        output_ = nullptr;
    }
}