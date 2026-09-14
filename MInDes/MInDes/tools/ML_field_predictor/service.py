import ctypes
import mmap
import time
import numpy as np
from multiprocessing import shared_memory

from shm_layout import (
    MAGIC,
    VERSION,
    MAX_ELEMS,
    STATE_IDLE,
    STATE_INPUT_READY,
    STATE_PROCESSING,
    STATE_OUTPUT_READY,
    STATE_ERROR,
    DTYPE_FLOAT32,
    SharedHeader,
    SHM_TOTAL_BYTES,
)

SHM_NAME = "SimML_Shared_01"


def get_header(buf) -> SharedHeader:
    return SharedHeader.from_buffer(buf, 0)


def get_input_view(buf, header: SharedHeader) -> np.ndarray:
    return np.ndarray(
        shape=(header.input_count,),
        dtype=np.float32,
        buffer=buf,
        offset=header.input_offset,
    )


def get_output_view(buf, header: SharedHeader) -> np.ndarray:
    return np.ndarray(
        shape=(MAX_ELEMS,),
        dtype=np.float32,
        buffer=buf,
        offset=header.output_offset,
    )


def main():
    print("Python service starting...")
    print(f"Connecting to shared memory: {SHM_NAME}")

    shm = shared_memory.SharedMemory(name=SHM_NAME, create=False, size=SHM_TOTAL_BYTES)
    buf = shm.buf
    header = get_header(buf)

    if header.magic != MAGIC:
        raise RuntimeError(f"Bad magic: {header.magic:#x}")
    if header.version != VERSION:
        raise RuntimeError(f"Bad version: {header.version}")

    print("Python attached shared memory successfully.")

    try:
        while True:
            state = header.state

            if state == STATE_INPUT_READY:
                print(f"[Python] request_id={header.request_id}, input_count={header.input_count}")
                header.state = STATE_PROCESSING

                x = get_input_view(buf, header)
                x_local = np.array(x, copy=True)  # first prototype, safe copy for processing

                # Demo processing: y = 2*x + 1
                y_local = x_local * np.float32(2.0) + np.float32(1.0)

                out = get_output_view(buf, header)
                out[:len(y_local)] = y_local

                header.output_count = len(y_local)
                header.result_code = 100
                header.output_dtype = DTYPE_FLOAT32
                header.state = STATE_OUTPUT_READY

                print(f"[Python] output ready, first values = {y_local[:min(5, len(y_local))]}")

            elif state == STATE_ERROR:
                print("[Python] detected error state")
                time.sleep(0.1)

            else:
                time.sleep(0.05)

    except KeyboardInterrupt:
        print("Python service stopped by user.")
    finally:
        del header
        shm.close()


if __name__ == "__main__":
    main()