#!/usr/bin/env bash
#
# Cross-checks the Fortran port against upstream's C QDLDL (in reference/qdldl):
# both factor and solve the same random quasi-definite matrices, and their
# elimination trees, L, D, and solutions are compared (the structure must be
# identical, the values equal to round-off).
#
# usage (from the repository root):
#
#    pixi run crosscheck
#
# This is not one of fpm's tests: it needs a C compiler and the upstream sources.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

cc -O2 -I"$here" -I"$root/reference/qdldl/include" \
    "$root/reference/qdldl/src/qdldl.c" "$here/crosscheck.c" -o "$work/crosscheck_c"
gfortran -O2 -J"$work" "$root/src/qdldl_kinds.F90" "$root/src/qdldl_core.f90" \
    "$here/crosscheck.f90" -o "$work/crosscheck_f"

cd "$work"
./crosscheck_f
./crosscheck_c
python3 "$here/compare.py"
