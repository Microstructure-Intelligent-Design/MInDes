import ctypes

MAGIC = 0x4C4D4953
VERSION = 1
MAX_ELEMS = 1024

STATE_IDLE = 0
STATE_INPUT_READY = 1
STATE_PROCESSING = 2
STATE_OUTPUT_READY = 3
STATE_ERROR = 4

DTYPE_FLOAT32 = 1

INPUT_BYTES = MAX_ELEMS * 4
OUTPUT_BYTES = MAX_ELEMS * 4


class SharedHeader(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("magic", ctypes.c_uint32),
        ("version", ctypes.c_uint32),
        ("state", ctypes.c_uint32),

        ("input_dtype", ctypes.c_uint32),
        ("output_dtype", ctypes.c_uint32),

        ("input_count", ctypes.c_uint64),
        ("output_count", ctypes.c_uint64),

        ("input_offset", ctypes.c_uint64),
        ("output_offset", ctypes.c_uint64),

        ("request_id", ctypes.c_uint64),
        ("result_code", ctypes.c_uint64),

        ("reserved", ctypes.c_char * 64),
    ]


HEADER_SIZE = ctypes.sizeof(SharedHeader)
SHM_TOTAL_BYTES = HEADER_SIZE + INPUT_BYTES + OUTPUT_BYTES


def input_offset_bytes() -> int:
    return HEADER_SIZE


def output_offset_bytes() -> int:
    return HEADER_SIZE + INPUT_BYTES