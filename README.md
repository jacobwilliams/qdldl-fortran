[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![Build Status](https://github.com/jacobwilliams/qdldl-fortran/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/qdldl-fortran/actions)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/qdldl-fortran)](https://github.com/jacobwilliams/qdldl-fortran/commits/master)
[![Docs](https://img.shields.io/badge/docs-api-blue)](https://jacobwilliams.github.io/qdldl-fortran)

A modern Fortran port of [QDLDL](https://github.com/osqp/qdldl), the sparse
$LDL^T$ solver for quasi-definite matrices used inside the OSQP solver. It
comes with fill-reducing orderings (AMD and reverse Cuthill–McKee) and an
object-oriented interface. You only need [fpm](https://fpm.fortran-lang.org)
to build it.

## What it is, and when to use it

QDLDL factors a sparse symmetric matrix as $PAP^T = LDL^T$ ($L$ unit lower
triangular, $D$ diagonal) **without pivoting**. That makes it small and fast,
and it is the reason for its main limitation. The factorization is guaranteed
stable only for **quasi-definite** matrices,

$$
\begin{bmatrix} H & B^T \cr B & -C \end{bmatrix}, \quad H \succ 0,\ C \succ 0,
$$

which can be factored in *any* symmetric order (Vanderbei, 1995). The order can
therefore be chosen for sparsity alone.

* **Use it for** small to medium quasi-definite systems: regularized KKT
  systems of optimization solvers, and least-squares augmented systems.
  The overhead per call is very low. Factoring and solving never allocate
  memory, the matrix can be of any real kind (single, double, or quadruple
  precision), and the solver objects can be copied.
* **Don't use it for** large 3-D problems or genuinely indefinite matrices.
  It is single-threaded, not supernodal, and doesn't pivot. Use MUMPS, HSL,
  SPRAL, or Pardiso there.

The signs of $D$ give the **inertia** of the matrix for free (Sylvester's law).
This count is exact for a quasi-definite matrix, and only as reliable as the
pivots otherwise.

## Example

```fortran
use qdldl_module
type(qdldl_type) :: ldl
integer :: istat

! the pattern, once: entries of either triangle in coordinate form
! (duplicates are added; missing diagonal entries become explicit zeros)
call ldl%analyze(n, irow, icol, istat)        ! AMD ordering by default

! the values, as often as needed (in the order of irow/icol)
call ldl%factor(val, istat)
call ldl%solve(b, istat)                      ! b is overwritten by x

call ldl%inertia(n_positive, n_negative, n_zero)
call ldl%multiply(x, y)                       ! y = A x, with the last values
```

See [`example/example_kkt.f90`](example/example_kkt.f90) for a complete program
(a KKT system with static regularization and iterative refinement).

To use it from another fpm project:

```toml
[dependencies]
qdldl-fortran = { git = "https://github.com/jacobwilliams/qdldl-fortran" }
```

## The interface

### `qdldl_type`

| Procedure | What it does |
|---|---|
| `analyze(n, irow, icol, istat [, perm] [, ordering])` | Analyzes the pattern, given in coordinate form. Computes the ordering, the elimination tree, and the nonzeros of $L$, and allocates everything else. |
| `analyze_csc(n, Ap, Ai, istat [, perm] [, ordering])` | The same, for a matrix in compressed sparse column (CSC) form (1-based, `Ap(1) = 1`). |
| `set_signs(signs, istat)` | The expected sign of each pivot (`+1` for the $H$ block, `-1` for the $-C$ block, `0` if unknown), used for regularization. |
| `factor(val, istat)` | Factors the matrix with new values. Duplicates are added. Allocates nothing. |
| `solve(b, istat [, refine])` | Solves $Ax = b$ in place, with optional iterative refinement. Allocates nothing. |
| `multiply(x, y [, istat])` | $y = Ax$ with the last values (unregularized). |
| `inertia(n_positive, n_negative, n_zero [, n_regularized])` | The inertia of the last factorization. |
| `get_permutation(perm)` | The permutation: row `perm(k)` of $A$ is row `k` of $PAP^T$. |
| `is_analyzed()`, `is_factored()` | The state. |
| `destroy()` | Frees everything (the options are kept). |

**Input.** An off-diagonal entry $(i,j)$ stands for both $(i,j)$ and
$(j,i)$, and duplicates are added. If you give both triangles, give the values
of only one of them (with zeros for the other), or the off-diagonal values are
doubled. The permutation is applied internally: `b` and `x` are always in the
original order.

**Copies.** All components are allocatable (no pointers), so `ldl2 = ldl1`
makes an independent copy, including the factorization.

**Options** (public components, with their defaults):

| Option | Default | Meaning |
|---|---|---|
| `ordering` | `qdldl_order_default` (AMD) | `qdldl_order_natural`, `qdldl_order_rcm`, `qdldl_order_amd`, or `qdldl_order_user` (with `perm`) |
| `zero_pivot_tol` | `0` | a pivot with $\lvert d_k\rvert \le$ `zero_pivot_tol` $\cdot \max_{ij}\lvert A_{ij}\rvert$ is zero (`0`: an exact zero only, as upstream) |
| `regularize` | `.false.` | dynamic regularization: a zero pivot, or one with $s_k d_k <$ `reg_eps`, is replaced by $s_k$ `reg_delta` instead of stopping the factorization |
| `reg_eps` | $\epsilon^{0.8}$ (3e-13 in double precision) | threshold of dynamic regularization (absolute) |
| `reg_delta` | $\sqrt\epsilon$ (1.5e-8) | magnitude of a regularized pivot (absolute) |
| `static_reg` | `0` | added to each diagonal entry with its expected sign before factoring |
| `max_refine` | `0` | maximum number of steps of iterative refinement in `solve`, against the unregularized matrix |

$s_k$ is the expected sign from `set_signs`. If it is unknown, the sign of
the pivot itself (dynamic regularization) or of the diagonal entry (static
regularization) is used.

**Results** (read-only public components): `n`, `nnz_a`, `nnz_l`,
`n_positive`, `n_negative`, `n_zero`, `n_regularized`, `zero_pivot_column`,
`max_abs_l` and `pivot_ratio` (to detect pivot growth), `refine_steps`, and
`residual`.

**Status codes.** `qdldl_success` (0), and negative named constants for
errors: `qdldl_error_not_upper`, `qdldl_error_overflow`,
`qdldl_error_empty_column`, `qdldl_error_zero_pivot`,
`qdldl_error_invalid_input`, `qdldl_error_not_analyzed`,
`qdldl_error_not_factored`, `qdldl_error_out_of_memory`, and
`qdldl_error_not_finite`. `qdldl_status_message(istat)` describes each one.
Allocations use `stat=`, so running out of memory returns a status code
instead of crashing.

### Regularization and inertia: caveats

* Without regularization (the default), the behaviour is upstream's: the first
  zero pivot stops the factorization (`qdldl_error_zero_pivot`, with the row in
  `zero_pivot_column`).
* **Static regularization** (`static_reg` > 0, with `set_signs`) makes a KKT
  matrix $\left[\matrix{H & J^T \cr J & 0}\right]$ with $H \succ 0$
  quasi-definite. Iterative refinement (`max_refine` > 0) then recovers the
  solution of the *unregularized* system. Choose `static_reg` well above the
  round-off of the real kind ($\sqrt\epsilon$ is a reasonable start).
* **Dynamic regularization** (`regularize`) always completes the factorization.
  `n_regularized` > 0 tells you that the matrix was not quasi-definite in that
  order, or was singular (an optimization solver would then increase its
  Hessian shift).
* The inertia counts the pivots *as computed, before replacement*. It is exact
  for a quasi-definite matrix with `n_regularized = 0`. Otherwise it describes
  the matrix that was actually factored.
* Without pivoting, a matrix that isn't quasi-definite can factor with large
  growth and give a wrong inertia. A large `max_abs_l` or `pivot_ratio` is the
  sign of that.

### Low-level routines

The one-to-one port of the C routines is public too, for callers that manage
their own storage. It uses 1-based indexing, the upper triangle in CSC form,
and no duplicates:

| Fortran | C |
|---|---|
| `qdldl_etree(n, Ap, Ai, work, Lnz, etree)` | `QDLDL_etree` |
| `qdldl_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork [, zero_col])` | `QDLDL_factor` |
| `qdldl_solve(n, Lp, Li, Lx, Dinv, x)` | `QDLDL_solve` |
| `qdldl_lsolve`, `qdldl_ltsolve` | `QDLDL_Lsolve`, `QDLDL_Ltsolve` |
| `qdldl_factor_ext(...)` | (new) `qdldl_factor` with a zero-pivot tolerance, dynamic regularization, and the inertia |

The orderings are also available on their own: `qdldl_amd_order` and
`qdldl_rcm` (on the adjacency structure built by `qdldl_symmetric_pattern`),
and the modernized SuiteSparse routine `amd` (module `qdldl_amd`). As in
SuiteSparse's C AMD, `qdldl_amd_order` removes dense rows (more than
$\max(16, 10\sqrt{n})$ entries) before ordering and places them last;
otherwise a dense row, common in KKT matrices, can make AMD take $O(n^2)$
time. The `amd` routine itself does not do this.

## Kinds

The kinds are chosen by preprocessor macros: `-DREAL32`, `-DREAL64` (the
default), or `-DREAL128` for the reals (`qdldl_wp`), and `-DINT64` for 64-bit
indices (`qdldl_ip`). You need 64-bit indices when the nonzeros of $L$ can
exceed $2^{31}-1$. Otherwise the analysis returns `qdldl_error_overflow`.

```bash
fpm test --profile debug --flag "-DREAL32 -DINT64"
```

## Development

Upstream's C QDLDL is a git submodule (`reference/qdldl`), used only by the
cross-check. Clone with `git clone --recursive`, or run
`git submodule update --init` in an existing clone.

The build environment is [pixi](https://pixi.sh) (gfortran, fpm, FORD,
fortitude, and a C compiler):

```bash
pixi run fpm test                 # all tests (double precision)
pixi run test-real32              # ... in single precision (also test-real128, test-int64)
pixi run fpm run --example benchmark --profile release
pixi run fortitude check          # lint
pixi run stack-check              # no automatic arrays or array temporaries in the library
pixi run crosscheck               # compare with upstream's C QDLDL (in reference/qdldl)
pixi run docs                     # the API documentation, with FORD
```

* **Tests**: `test_core` runs upstream QDLDL's unit tests through the
  low-level routines. `test_type` covers the object-oriented interface.
  `test_orderings` covers the orderings, with the fill of each recorded as a
  regression check. `test_random_qd` checks random quasi-definite matrices,
  with exact inertia. `test_inertia` covers regularization, refinement, and
  zero pivots, and `test_kinds` the kind of the build. Tolerances scale with
  `epsilon(1.0_wp)`, and every test passes in each kind.
* **Cross-check**: `tools/crosscheck` builds the C QDLDL and runs both
  versions on the same random matrices (up to $n = 2000$, 640,000 nonzeros
  in $L$). Their elimination trees and the structure of $L$ are identical,
  and the values of $L$, $D$, and the solutions agree **bit for bit**.
* **No stack use that grows with the problem**: the library has no automatic
  arrays and no array temporaries (`tools/stack_check.sh`). The work arrays
  are kept in the type, and allocated in the analysis.

### Benchmark

`example/benchmark.f90`, release build, on an Apple M-series laptop (times in
seconds):

| Matrix | n | ordering | nnz(L) | analyze | factor | solve |
|---|---:|---|---:|---:|---:|---:|
| 2-D Laplacian 100² | 10,000 | natural | 990,099 | 0.003 | 0.023 | 0.0012 |
| | | RCM | 671,550 | 0.002 | 0.011 | 0.0008 |
| | | AMD | 196,332 | 0.002 | 0.003 | 0.0002 |
| 2-D KKT 100² | 20,000 | natural | 2,019,896 | 0.004 | 0.063 | 0.0020 |
| | | AMD | 659,448 | 0.004 | 0.015 | 0.0007 |
| 3-D Laplacian 20³ | 8,000 | natural | 3,047,619 | 0.006 | 0.206 | 0.0032 |
| | | AMD | 834,282 | 0.004 | 0.057 | 0.0009 |
| 3-D KKT 20³ | 16,000 | natural | 5,996,075 | 0.012 | 0.794 | 0.0059 |
| | | AMD | 3,598,128 | 0.011 | 0.508 | 0.0036 |

### C and Fortran kernel comparison

`pixi run benchmark` compiles both implementations with `-O3` and times their
matching low-level elimination-tree, factorization, and solve routines on the
same 2-D five-point Laplacian matrices using natural ordering. The matrix is
generated once; conversion between the C (0-based) and Fortran (1-based) CSC
indices, allocations, and correctness checks are outside the timed regions.
Each entry is the median per-call time from five samples (20 elimination-tree,
5 factorization, or 100 solve calls per sample). A ratio below 1 means the C
routine was faster.

Results on an Apple M5 (macOS 27.0), using GNU Fortran 15.3.0 and Clang 23.1.1:

| Grid | n | nnz(A) | nnz(L) | Routine | Fortran (us) | C (us) | C / Fortran |
|---|---:|---:|---:|---|---:|---:|---:|
| 20 x 20 | 400 | 1,160 | 7,619 | elimination tree | 9.80 | 8.55 | 0.87 |
| | | | | factorization | 44.40 | 44.40 | 1.00 |
| | | | | solve | 9.81 | 9.69 | 0.99 |
| 50 x 50 | 2,500 | 7,400 | 122,549 | elimination tree | 202.10 | 182.40 | 0.90 |
| | | | | factorization | 1,300.80 | 1,222.60 | 0.94 |
| | | | | solve | 152.31 | 141.00 | 0.93 |
| 100 x 100 | 10,000 | 29,800 | 990,099 | elimination tree | 1,788.45 | 2,126.45 | 1.19 |
| | | | | factorization | 20,319.20 | 18,508.00 | 0.91 |
| | | | | solve | 1,154.61 | 1,102.17 | 0.95 |

These figures compare the low-level kernels, not the complete public Fortran
interface (which also supports ordering, coordinate input, and additional
validation).

## License

Apache License 2.0 ([LICENSE](LICENSE)). This library is a derivative of
QDLDL (Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, Ian
McInerney, The OSQP developers; Apache 2.0). It includes a modernized version
of the Fortran AMD of SuiteSparse (Copyright 1996-2022, Timothy A. Davis,
Patrick R. Amestoy, and Iain S. Duff; BSD 3-clause). See [NOTICE](NOTICE).

## References

* R. J. Vanderbei, *Symmetric quasidefinite matrices*, SIAM J. Optim. 5(1),
  100–113, 1995.
* T. A. Davis, *Algorithm 849: A concise sparse Cholesky factorization
  package*, ACM TOMS 31(4), 587–591, 2005.
* P. R. Amestoy, T. A. Davis, I. S. Duff, *An approximate minimum degree
  ordering algorithm*, SIAM J. Matrix Anal. Appl. 17(4), 886–905, 1996.
* B. Stellato, G. Banjac, P. Goulart, A. Bemporad, S. Boyd, *OSQP: an
  operator splitting solver for quadratic programs*, Math. Prog. Comp. 12,
  637–672, 2020.
