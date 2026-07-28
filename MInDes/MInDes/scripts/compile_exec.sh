#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)

CC=${CC:-gcc}
CXX=${CXX:-g++}
BUILD_TYPE=${BUILD_TYPE:-Release}
TPM_MODE=${TPM_MODE:-auto}
NPROC=${NPROC:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}

case "$(uname -s)" in
    Linux) ;;
    *)
        echo "> This build script only supports Linux." >&2
        exit 2
        ;;
esac

case "$(uname -m)" in
    x86_64|amd64) ;;
    *)
        echo "> Unsupported Linux architecture: $(uname -m). MInDes supports Linux x86_64 only." >&2
        exit 2
        ;;
esac

case "$BUILD_TYPE" in
    Release|Debug) ;;
    *)
        echo "> Invalid BUILD_TYPE: $BUILD_TYPE (expected Release or Debug)." >&2
        exit 2
        ;;
esac

case "$TPM_MODE" in
    auto|required|off) ;;
    *)
        echo "> Invalid TPM_MODE: $TPM_MODE (expected auto, required, or off)." >&2
        exit 2
        ;;
esac

if ! [[ $NPROC =~ ^[1-9][0-9]*$ ]]; then
    echo "> Invalid NPROC value: $NPROC" >&2
    exit 2
fi

for tool in "$CC" "$CXX" pkg-config find sort; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "> Required build tool was not found: $tool" >&2
        exit 2
    fi
done

if ! pkg-config --exists libcrypto; then
    echo "> OpenSSL development files were not found (pkg-config module: libcrypto)." >&2
    echo "> On Ubuntu/Debian, install the libssl-dev package." >&2
    exit 2
fi

FFTW_INCLUDE_DIR="$PROJECT_DIR/lib/linux"
FFTW_LIBRARY_DIR="$PROJECT_DIR/lib/linux/lib"
FFTW_HEADER="$FFTW_INCLUDE_DIR/fftw3.h"
FFTW_LIBRARY="$FFTW_LIBRARY_DIR/libfftw3.a"
if [[ ! -f "$FFTW_HEADER" || ! -f "$FFTW_LIBRARY" ]]; then
    echo "> Bundled Linux FFTW files are missing:" >&2
    echo ">   $FFTW_HEADER" >&2
    echo ">   $FFTW_LIBRARY" >&2
    exit 2
fi

if ! "$CXX" -std=c++17 -fopenmp -x c++ -fsyntax-only /dev/null >/dev/null 2>&1; then
    echo "> The selected C++ compiler does not support the required OpenMP option (-fopenmp): $CXX" >&2
    exit 2
fi

declare -a CRYPTO_CFLAGS CRYPTO_LIBS TPM_CFLAGS TPM_LIBS
read -r -a CRYPTO_CFLAGS <<< "$(pkg-config --cflags libcrypto)"
read -r -a CRYPTO_LIBS <<< "$(pkg-config --libs libcrypto)"
TPM_CFLAGS=()
TPM_LIBS=()
TPM_ENABLED=0

if [[ "$TPM_MODE" != off ]]; then
    if pkg-config --exists tss2-fapi; then
        read -r -a TPM_CFLAGS <<< "$(pkg-config --cflags tss2-fapi)"
        read -r -a TPM_LIBS <<< "$(pkg-config --libs tss2-fapi)"
        TPM_ENABLED=1
    elif [[ "$TPM_MODE" == required ]]; then
        echo "> TPM_MODE=required, but TPM2-TSS FAPI development files were not found." >&2
        echo "> On Ubuntu/Debian, install the libtss2-dev package." >&2
        exit 2
    else
        echo "> Warning: tss2-fapi was not found; building an SMBIOS-only License V2 client." >&2
    fi
fi

BUILD_DIR="$PROJECT_DIR/.build/linux/$BUILD_TYPE"
OBJECT_DIR="$BUILD_DIR/objects"
OUTPUT_FILE="$BUILD_DIR/MInDes"

COMMON_C_FLAGS=(-std=c11 -march=x86-64 -mtune=generic -pthread -I"$FFTW_INCLUDE_DIR")
COMMON_CXX_FLAGS=(-std=c++17 -march=x86-64 -mtune=generic -fopenmp -pthread -I"$FFTW_INCLUDE_DIR")
COMMON_LINK_FLAGS=(-fopenmp -pthread)

if [[ "$BUILD_TYPE" == Release ]]; then
    COMMON_C_FLAGS+=(-O3 -DNDEBUG -flto)
    COMMON_CXX_FLAGS+=(-O3 -DNDEBUG -flto)
    COMMON_LINK_FLAGS+=(-flto)
else
    COMMON_C_FLAGS+=(-O0 -g -D_DEBUG)
    COMMON_CXX_FLAGS+=(-O0 -g -D_DEBUG)
fi

COMMON_CXX_FLAGS+=("${CRYPTO_CFLAGS[@]}")
if (( TPM_ENABLED != 0 )); then
    COMMON_CXX_FLAGS+=(-DMINDES_HAS_TPM2_TSS=1 "${TPM_CFLAGS[@]}")
fi

if [[ -n ${CFLAGS:-} ]]; then
    read -r -a EXTRA_C_FLAGS <<< "$CFLAGS"
    COMMON_C_FLAGS+=("${EXTRA_C_FLAGS[@]}")
fi
if [[ -n ${CXXFLAGS:-} ]]; then
    read -r -a EXTRA_CXX_FLAGS <<< "$CXXFLAGS"
    COMMON_CXX_FLAGS+=("${EXTRA_CXX_FLAGS[@]}")
fi
if [[ -n ${LDFLAGS:-} ]]; then
    read -r -a EXTRA_LINK_FLAGS <<< "$LDFLAGS"
    COMMON_LINK_FLAGS+=("${EXTRA_LINK_FLAGS[@]}")
fi

rm -rf -- "$BUILD_DIR"
mkdir -p -- "$OBJECT_DIR"

mapfile -d '' SOURCES < <(find "$PROJECT_DIR/src" -type f \( -name '*.cpp' -o -name '*.c' \) -print0 | sort -z)
if (( ${#SOURCES[@]} == 0 )); then
    echo "> No source files found under $PROJECT_DIR/src." >&2
    exit 2
fi

OBJECTS=()
BATCH_PIDS=()
BATCH_SOURCES=()

wait_for_batch() {
    local failed=0
    local index
    for index in "${!BATCH_PIDS[@]}"; do
        if ! wait "${BATCH_PIDS[$index]}"; then
            echo "> Compilation failed: ${BATCH_SOURCES[$index]}" >&2
            failed=1
        fi
    done
    BATCH_PIDS=()
    BATCH_SOURCES=()
    if (( failed != 0 )); then
        exit 1
    fi
}

echo "> Configuration: $BUILD_TYPE"
echo "> TPM2-TSS FAPI: $([[ $TPM_ENABLED == 1 ]] && echo enabled || echo disabled)"
echo "> Compiling ${#SOURCES[@]} source files with $NPROC parallel jobs..."
for source in "${SOURCES[@]}"; do
    relative=${source#"$PROJECT_DIR/src/"}
    object="$OBJECT_DIR/${relative%.*}.o"
    mkdir -p -- "$(dirname -- "$object")"
    OBJECTS+=("$object")

    if [[ $source == *.c ]]; then
        ("$CC" "${COMMON_C_FLAGS[@]}" -c "$source" -o "$object") &
    else
        ("$CXX" "${COMMON_CXX_FLAGS[@]}" -c "$source" -o "$object") &
    fi
    BATCH_PIDS+=("$!")
    BATCH_SOURCES+=("$source")

    if (( ${#BATCH_PIDS[@]} >= NPROC )); then
        wait_for_batch
    fi
done
if (( ${#BATCH_PIDS[@]} != 0 )); then
    wait_for_batch
fi

echo "> Linking $OUTPUT_FILE..."
"$CXX" "${OBJECTS[@]}" -o "$OUTPUT_FILE" "${COMMON_LINK_FLAGS[@]}" \
    -L"$FFTW_LIBRARY_DIR" -lfftw3 "${CRYPTO_LIBS[@]}" "${TPM_LIBS[@]}" -ldl -lm

echo "> MInDes has been built successfully: $OUTPUT_FILE"
