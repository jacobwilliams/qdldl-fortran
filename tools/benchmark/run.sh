#!/usr/bin/env bash
#
# Times upstream C QDLDL and the matching low-level Fortran kernels on identical
# 2-D Laplacian matrices. Run from any directory with `pixi run benchmark`.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

cc -O3 -I"$root/tools/crosscheck" -I"$root/reference/qdldl/include" \
    -c "$here/benchmark_c.c" -o "$work/benchmark_c.o"
cc -O3 -I"$root/tools/crosscheck" -I"$root/reference/qdldl/include" \
    -c "$root/reference/qdldl/src/qdldl.c" -o "$work/qdldl.o"
gfortran -funroll-loops -O3 -ffree-line-length-none -J"$work" -I"$work" \
    "$root/src/qdldl_kinds.F90" "$root/src/qdldl_core.f90" "$here/benchmark.f90" \
    "$work/benchmark_c.o" "$work/qdldl.o" -o "$work/benchmark"

"$work/benchmark"
