#!/usr/bin/env bash
# Minimal giac 2.0.0 icas build (tommath, no FLTK/GMP/QuickJS).
# Verified flags — feed giac/CMakeLists.txt from this recipe.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${ROOT}/giac/giac-2.0.0"
TOMMATH_INC="${ROOT}/libtommath/libtommath-0.39"
TOMMATH_LIB="${ROOT}/build/lib"

[[ -f "${TOMMATH_LIB}/libtommath.a" ]] || {
  echo "Build libtommath first: cmake -S ${ROOT} -B ${ROOT}/build && cmake --build ${ROOT}/build --target tommath"
  exit 1
}

cd "${SRC}"
autoreconf -iv

export CPPFLAGS="-I${TOMMATH_INC}"
export LDFLAGS="-L${TOMMATH_LIB}"
export LIBS="-ltommath"

./configure \
  --enable-tommath \
  --disable-shared \
  --enable-static \
  --disable-fltk \
  --disable-pari \
  --disable-cocoa \
  --disable-ntl \
  --disable-ecm \
  --disable-bernmm \
  --disable-gsl \
  --disable-lapack \
  --disable-ao \
  --disable-glpk \
  --disable-samplerate \
  --disable-curl \
  --disable-micropy \
  --disable-quickjs \
  --disable-png

make -C src -j"$(nproc)" icas

# Install config.h for CMake (must match src/ and top-level copies).
install -m 644 config.h "${ROOT}/giac/config.h"
install -m 644 config.h "${SRC}/src/config.h"
# CMake CAS build: no readline/gettext in CI
sed -i 's/^#define ENABLE_NLS 1$/\/* #undef ENABLE_NLS *\//' "${ROOT}/giac/config.h"
sed -i 's/^#define HAVE_LIBREADLINE 1$/\/* #undef HAVE_LIBREADLINE *\//' "${ROOT}/giac/config.h"
sed -i 's/^#define HAVE_READLINE_HISTORY_H 1$/\/* #undef HAVE_READLINE_HISTORY_H *\//' "${ROOT}/giac/config.h"
sed -i 's/^#define HAVE_READLINE_READLINE_H 1$/\/* #undef HAVE_READLINE_READLINE_H *\//' "${ROOT}/giac/config.h"
sed -i 's/^#define DEBUG_SUPPORT \/\*\*\/$/\/* #undef DEBUG_SUPPORT *\//' "${ROOT}/giac/config.h"
cp "${ROOT}/giac/config.h" "${SRC}/src/config.h"
cp "${ROOT}/giac/config.h" "${SRC}/config.h"

echo "Built: ${SRC}/src/icas"
"${SRC}/src/icas" --version
