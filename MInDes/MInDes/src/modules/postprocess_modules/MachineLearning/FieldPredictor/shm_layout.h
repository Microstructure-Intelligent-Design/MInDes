#pragma once
#include <cstdint>
namespace pf {
#pragma pack(push, 1)

    enum class ShmState : uint32_t {
        Idle = 0,
        InputReady = 1,
        Processing = 2,
        OutputReady = 3,
        Error = 4
    };

    enum class ShmDType : uint32_t {
        Float32 = 1
    };

    struct SharedHeader {
        uint32_t magic;              // "SIML"
        uint32_t version;            // protocol version
        volatile uint32_t state;     // ShmState

        uint32_t input_dtype;        // ShmDType
        uint32_t output_dtype;       // ShmDType

        uint64_t input_count;        // number of float elements
        uint64_t output_count;       // number of float elements

        uint64_t input_offset;       // byte offset from base
        uint64_t output_offset;      // byte offset from base

        uint64_t request_id;         // request sequence id
        uint64_t result_code;        // 0=unset, 100=ok, 500=error

        char reserved[64];           // future use
    };

#pragma pack(pop)

    constexpr uint32_t kMagic = 0x4C4D4953; // 'SIML' little-endian
    constexpr uint32_t kVersion = 1;
    constexpr uint64_t kMaxElems = 1024;

    constexpr uint64_t kInputBytes = kMaxElems * sizeof(float);
    constexpr uint64_t kOutputBytes = kMaxElems * sizeof(float);
    constexpr uint64_t kShmTotalBytes = sizeof(SharedHeader) + kInputBytes + kOutputBytes;

    inline uint64_t input_offset_bytes() {
        return sizeof(SharedHeader);
    }

    inline uint64_t output_offset_bytes() {
        return sizeof(SharedHeader) + kInputBytes;
    }
}