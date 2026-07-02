#!/usr/bin/env bash
set -euo pipefail

CC=${CC:-gcc}
CXX=${CXX:-g++}
NPROC=${NPROC:-$(nproc 2>/dev/null || echo 4)}
BUILD_DIR=.build/linux

C_FLAGS=(-std=c11 -O3 -march=x86-64 -mtune=generic -flto -pthread -Ilib/linux)
CXX_FLAGS=(-std=c++17 -O3 -march=x86-64 -mtune=generic -fopenmp -flto -pthread -Ilib/linux)
LINK_FLAGS=(-fopenmp -flto -pthread)

if [[ -n ${CFLAGS:-} ]]; then
    read -r -a EXTRA_C_FLAGS <<< "$CFLAGS"
    C_FLAGS+=("${EXTRA_C_FLAGS[@]}")
fi
if [[ -n ${CXXFLAGS:-} ]]; then
    read -r -a EXTRA_CXX_FLAGS <<< "$CXXFLAGS"
    CXX_FLAGS+=("${EXTRA_CXX_FLAGS[@]}")
fi
if [[ -n ${LDFLAGS:-} ]]; then
    read -r -a EXTRA_LINK_FLAGS <<< "$LDFLAGS"
    LINK_FLAGS+=("${EXTRA_LINK_FLAGS[@]}")
fi

if ! [[ $NPROC =~ ^[1-9][0-9]*$ ]]; then
    echo "> Invalid NPROC value: $NPROC" >&2
    exit 2
fi

rm -f -- ./*.o
rm -rf -- "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

mapfile -d '' SOURCES < <(find src -type f \( -name '*.cpp' -o -name '*.c' \) -print0 | sort -z)
if (( ${#SOURCES[@]} == 0 )); then
    echo "> No source files found under src." >&2
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

echo "> Compiling ${#SOURCES[@]} source files with $NPROC parallel jobs..."
for source in "${SOURCES[@]}"; do
    relative=${source#src/}
    object="$BUILD_DIR/src/${relative%.*}.o"
    mkdir -p "$(dirname "$object")"
    OBJECTS+=("$object")

    if [[ $source == *.c ]]; then
        ("$CC" "${C_FLAGS[@]}" -c "$source" -o "$object") &
    else
        ("$CXX" "${CXX_FLAGS[@]}" -c "$source" -o "$object") &
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

echo "> Linking MInDes..."
"$CXX" "${OBJECTS[@]}" -o MInDes "${LINK_FLAGS[@]}" \
    -Llib/linux/lib -lfftw3 -ldl -lm

echo "> MInDes has been built from source."
