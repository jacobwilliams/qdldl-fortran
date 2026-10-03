# Plan: a Fortran QDLDL package (fpm)

A plan for a standalone fpm package: a modern Fortran port of
[QDLDL](https://github.com/osqp/qdldl), the sparse \(LDL^T\) solver for
quasi-definite matrices inside the OSQP solver, with a fill-reducing
ordering and a small object-oriented interface. It is written to be read on
its own: section 9 says what [sqpopt](https://github.com/jacobwilliams/sqpopt)
needs from it, as its first user, but the package must not depend on sqpopt.

Facts about the upstream C code (function names, arguments, return values)
are as I remember them from QDLDL's repository; check each against the
current upstream source before porting (stage 1 starts with that).

---

## 1. What QDLDL is, and why port it

- **What it does.** Factors a sparse symmetric matrix as
  \(A = LDL^T\) (\(L\) unit lower triangular, \(D\) diagonal), then solves
  \(Ax = b\). It has two stages: a *symbolic* one (the elimination tree and
  the column counts of \(L\), from the pattern alone) and a *numeric* one
  (the values). The symbolic stage is done once per pattern, the numeric one
  for every new set of values.
- **No pivoting.** That is why it is small (a few hundred lines of C), and
  its limitation: it is only guaranteed stable for **quasi-definite**
  matrices, \( \begin{bmatrix} H & B^T \\ B & -C \end{bmatrix} \) with \(H\)
  and \(C\) positive definite (Vanderbei, 1995). Those can be factored in
  *any* symmetric order without pivoting, so the order can be chosen for
  sparsity alone. For other indefinite matrices the factorization may still
  run, but can be unstable (tiny pivots, large growth).
- **Inertia comes free.** By Sylvester's law of inertia, the signs of
  \(D\) are the signs of \(A\)'s eigenvalues. So a successful factorization
  gives the numbers of positive, negative, and (with a tolerance) zero
  eigenvalues, which optimization solvers need for inertia control.
- **No ordering.** QDLDL factors the matrix in the order given. OSQP orders
  it first with AMD (approximate minimum degree, from SuiteSparse). This
  package must provide an ordering, or the fill makes it useless on real
  problems.
- **Why a Fortran port.** Fortran projects (sqpopt first) get a sparse
  symmetric factorization that fpm builds with nothing else: no C compiler
  setup, no MUMPS or HSL, any real kind (single, double, quadruple), and
  types that can be copied (allocatable components, no pointers), unlike a
  MUMPS handle.
- **What it isn't.** A replacement for MUMPS, HSL, SPRAL, or Pardiso on
  large problems, or on genuinely indefinite matrices: it is
  single-threaded, not supernodal, and doesn't pivot. It is a dependable,
  always-available solver for small to medium quasi-definite systems, with
  very little overhead per call.

**Upstream:** <https://github.com/osqp/qdldl> (Paul Goulart, Bartolomeo
Stellato, Goran Banjac; Apache License 2.0). The algorithm is an
"up-looking" \(LDL^T\) driven by the elimination tree, close to Tim Davis's
LDL package (T. A. Davis, *Algorithm 849: A concise sparse Cholesky
factorization package*, ACM TOMS 31(4), 2005).

**Design reference for the extras:** the Rust port of QDLDL used by the
Clarabel solver (also Apache 2.0) adds AMD ordering, refactoring with new
values, static and dynamic regularization with expected pivot signs, and
iterative refinement around it. Look at its interface before designing
stage 4; don't copy its code without keeping its licence notice.

---

## 2. Upstream C interface (to port first, one to one)

As I recall it (verify):

| C function | What it does |
|---|---|
| `QDLDL_etree(n, Ap, Ai, work, Lnz, etree)` | From the **upper triangle** of \(A\) in compressed sparse column (CSC) form (`Ap` column pointers, `Ai` row indices, 0-based), compute the elimination tree `etree` and the nonzero count of each column of \(L\), `Lnz`. Returns the total nonzeros of \(L\), or a negative value if the input isn't upper triangular (an entry below the diagonal) or the count overflows the integer type. |
| `QDLDL_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork)` | The numeric factorization: \(L\) in CSC (`Lp`, `Li`, `Lx`), the diagonal `D` and its inverse `Dinv`. Work arrays: `bwork` (logical, `n`), `iwork` (integer, `3n`), `fwork` (real, `n`). Returns the number of positive entries of `D`, or `-1` if a pivot is exactly zero. |
| `QDLDL_solve(n, Lp, Li, Lx, Dinv, x)` | Solves \(LDL^Tx = b\) in place (`x` holds `b` on entry). |
| `QDLDL_Lsolve`, `QDLDL_Ltsolve` | The two triangular solves, separately. |

Types: `QDLDL_int` (an integer, 32 or 64 bit), `QDLDL_float` (double or
single), `QDLDL_bool`. Markers inside the factorization: `QDLDL_UNKNOWN`
(\(-1\), "no parent" in the elimination tree), `QDLDL_USED` /
`QDLDL_UNUSED` (the flags in `bwork`).

**Inputs the C code assumes:** only the upper triangle, including the
diagonal; CSC with sorted or unsorted row indices within a column (check
which); no duplicates. A missing diagonal entry behaves as a zero pivot,
so callers should always include every diagonal position (an explicit zero
if need be).

---

## 3. Package layout

```
qdldl-fortran/
  fpm.toml
  LICENSE            Apache 2.0 (see section 8)
  NOTICE             credits QDLDL's authors, and AMD's if included
  README.md
  ford.md            FORD settings for the API docs
  src/
    qdldl_kinds.F90        real and integer kinds (preprocessor-selected)
    qdldl_core.f90         the one-to-one port: etree, factor, solve, Lsolve, Ltsolve
    qdldl_ordering.f90     natural, user-given, RCM, and AMD orderings
    qdldl_amd.f            (optional) SuiteSparse's Fortran AMD, unchanged, if its licence allows
    qdldl_module.f90       the object-oriented interface (qdldl_type) and the public API
  test/
    test_core_*.f90        ports of upstream's unit tests
    test_type_*.f90        the high-level interface
    test_random_qd.f90     random quasi-definite matrices
    test_inertia.f90       inertia, zero pivots, regularization
    test_kinds.f90         run in each real kind by CI
  example/
    example_kkt.f90        factor and solve a small KKT system
    benchmark.f90          timings on scalable matrices (section 7)
```

Conventions (from the author's other Fortran libraries):

- `implicit none` everywhere; reals `real(wp)`, literals `1.0_wp`.
- Kind selection like `sqpopt_kinds`: `-DREAL32`, `-DREAL64` (default),
  `-DREAL128` choose `wp`; add `-DINT64` for 64-bit indices (`ip`), since
  the nonzeros of \(L\) can exceed \(2^{31}\) on large problems.
- **No local variable initialized in its declaration** (implicit `save`).
- **Nothing sized by the problem on the stack:** no automatic arrays, no
  array temporaries of size `n` or `nnz` (expressions as actual arguments,
  array constructors, `pack`, ...). Work arrays are allocatable components
  of the type, allocated once in the analysis, so factor and solve allocate
  nothing.
- Allocations use `stat=`, and report failure (out of memory) through the
  status code, never by crashing.
- FORD comments: a `!>` header for each module, `!!` for every procedure,
  type, component, and dummy argument.
- 1-based indexing throughout (see section 4.1).

---

## 4. Stages

### Stage 1: the one-to-one port (`qdldl_core`)

Port `etree`, `factor`, `solve`, `Lsolve`, `Ltsolve` as plain procedures
with explicit-shape or assumed-shape arguments, keeping upstream's
algorithm and argument order so that the two can be compared line by line.

- **Indexing.** Convert to 1-based: column pointers `Ap(1:n+1)` with
  `Ap(1) = 1`; row indices `1..n`. Upstream's "no parent" marker \(-1\)
  becomes `0`. Go through every loop bound and every comparison with a
  marker; this is where a port goes wrong.
- **Return values.** Keep upstream's meaning (total nonzeros of \(L\); number
  of positive pivots; negative codes for errors), but as named integer
  constants (`qdldl_error_not_upper`, `qdldl_error_zero_pivot`,
  `qdldl_error_overflow`, ...).
- **Tests.** Port upstream's test cases (its `tests` directory, as I
  recall: a basic matrix, an identity, a 2x2, a singleton, a matrix with a
  zero on the diagonal, a rank-deficient one, an OSQP KKT matrix, and input
  that isn't upper triangular). Each checks the solution (\(\lVert Ax - b
  \rVert\)), and the error cases their error code.
- **Cross-check (optional).** A script outside fpm's tests that builds the
  C library and runs both on the same matrices, comparing `D`, `L`, and the
  solutions to round-off. Worth doing once, to trust the port.

**Done when:** every upstream test passes, in double and single precision.

### Stage 2: the object-oriented interface (`qdldl_type`)

What most users want: give the pattern once, then factor and solve as
often as needed.

```fortran
type(qdldl_type) :: ldl
call ldl%analyze(n, irow, icol, istat [, perm] [, ordering])   ! the pattern, once
call ldl%factor(val, istat)                                     ! values in the pattern's order
call ldl%solve(b, istat [, refine])                             ! b is overwritten by x
call ldl%inertia(n_positive, n_negative, n_zero)
call ldl%multiply(x, y)                                         ! y = A x (the last values factored)
call ldl%destroy()
```

- **Input in coordinate form**, entries of *either* triangle (or both
  mirrored, if documented which wins), **duplicates added together**,
  missing diagonal positions added as explicit zeros. `analyze` builds the
  upper-triangular CSC of the *permuted* matrix \(PAP^T\), and a map from
  each coordinate entry to its position in that CSC, so that `factor` only
  scatters `val` through the map (no sorting, no allocation).
- **Also accept CSC directly** (upper triangle), for callers that have it.
- **The permutation** \(P\) is stored with its inverse; `solve` applies it
  to `b` and to the solution, so the caller never sees it.
- **Status codes**: success; not analysed / not factored; invalid input
  (index out of range, `n < 0`); zero pivot (with the column); out of
  memory; non-finite values.
- **Copyable**: allocatable components only, so assignment copies a
  factorization (useful for a solver that keeps one while computing
  another). Document it.

**Tests:** the stage-1 matrices through the type; coordinate input with
duplicates, both triangles, and a missing diagonal; refactoring with new
values (same pattern) gives the same result as a fresh analysis; a copied
object solves independently; `multiply` against a dense product.

### Stage 3: orderings (`qdldl_ordering`)

Without a fill-reducing order, \(L\) fills in (a 2-D grid with \(n\)
unknowns has \(O(n^{1.5})\) fill even with a good order, and far more
without one).

1. **Natural** (no permutation) and **user-given** (`perm` argument): for
   testing, and for callers with their own ordering.
2. **Reverse Cuthill–McKee (RCM):** simple to write in Fortran (a
   breadth-first search from a pseudo-peripheral node), reduces bandwidth,
   good on banded and grid-like problems. A safe default if AMD is not
   available.
3. **AMD** (approximate minimum degree; Amestoy, Davis & Duff, SIAM J.
   Matrix Anal. Appl. 17(4), 1996, and ACM TOMS Algorithm 837, 2004): the
   right default, as in OSQP. Options, in order of preference:
   - SuiteSparse's AMD has included Fortran 77 versions (`amd.f`,
     `amdbar.f`) under a BSD licence. If that is still so, include it
     unchanged, with its licence notice, and wrap it (it wants the full
     symmetric pattern, without the diagonal, in its own format: check).
   - Otherwise, port AMD from the C source (larger: about 2,000 lines), or
     ship RCM only at first.
4. **Default:** AMD if built in, else RCM. Report the nonzeros of \(L\) for
   each, so a user can compare.

**Tests:** each ordering gives a valid permutation and the same solution;
on a 2-D grid KKT matrix, AMD's fill is well below the natural order's
(record the counts in the test, as a regression check).

### Stage 4: robustness for optimization solvers

What makes it usable beyond strictly quasi-definite matrices. All optional,
off by default (so the default behaviour is upstream's).

- **Zero-pivot tolerance.** Treat \(|d_k| \le \tau\) (relative to the
  largest diagonal entry of \(A\), say) as zero: report it as
  `n_zero` and either stop (upstream's behaviour) or continue with a
  replacement value (below).
- **Static regularization.** Add \(\epsilon_k\) to each diagonal entry, with
  the sign the caller expects for that row (`+` for the \(H\) block, `-` for
  the \(-C\) block). This makes a KKT matrix \(\begin{bmatrix} H & J^T \\ J &
  0 \end{bmatrix}\) quasi-definite when \(H\) is positive definite; iterative
  refinement against the *unregularized* matrix then removes the effect of
  \(\epsilon\).
- **Dynamic regularization** (as in Clarabel's port): given the expected
  sign \(s_k\) of each pivot, replace a pivot with \(s_k d_k < \epsilon\) by
  \(s_k \delta\), and count the replacements. The factorization always
  completes; the count tells the caller that the matrix was not
  quasi-definite (an optimization solver would then raise its Hessian
  shift).
- **Inertia.** `n_positive`, `n_negative`, `n_zero` from \(D\) (and the
  number of regularized pivots). Document that it is exact for a
  quasi-definite matrix, and only as reliable as the pivots otherwise.
- **Iterative refinement** in `solve`: a few steps of
  \(r = b - Ax\), \(x \leftarrow x + \text{solve}(r)\), against the original
  (unregularized) values, while the residual is above round-off and each
  step at least halves it. Kept work arrays, no allocation.
- **Pivot growth check (optional).** \(\max |L|\) or \(\max |d_k| / \min
  |d_k|\), so a caller can tell an unstable factorization of a matrix that
  isn't quasi-definite.

**Tests:** matrices that are quasi-definite (inertia exact), a KKT matrix
with a singular \(H\) block (zero pivots: reported, or regularized),
regularization plus refinement reaching the unregularized solution to
round-off, an indefinite matrix that is not quasi-definite (factorization
completes with dynamic regularization; the count is non-zero).

### Stage 5: kinds, performance, docs, release

- **Kinds:** CI builds and tests with `REAL32`, `REAL64`, `REAL128`, and
  `INT64`. Tolerances in the tests scale with `epsilon(1.0_wp)`.
- **Stack:** a check that there are no automatic arrays or problem-sized
  temporaries (gfortran's `-Warray-temporaries`, plus a grep for
  explicit-shape locals sized by arguments).
- **Benchmark** (`example/benchmark.f90`): 2-D and 3-D Laplacians, and KKT
  matrices built from them (\(\begin{bmatrix} I & A^T \\ A & -\epsilon I
  \end{bmatrix}\)); analysis, factorization, and solve times, and the
  nonzeros of \(L\), for each ordering. Compare with the C QDLDL (should be
  within a small factor) and, where available, with MUMPS (QDLDL should win
  on small matrices, by its low overhead, and lose on large 3-D ones).
- **Docs:** README (what it is, when to use it and when not, a short
  example, the licence), FORD API docs, and the inertia and regularization
  caveats stated plainly.
- **CI:** GitHub Actions with gfortran (and ifx if convenient), the kind
  matrix above, and the FORD build.
- **Release:** tag `v0.1.0`; usable as an fpm dependency by git URL and
  tag.

---

## 5. Public API (target)

```fortran
module qdldl_module
    ! kinds
    integer, parameter, public :: qdldl_wp   ! real kind
    integer, parameter, public :: qdldl_ip   ! integer kind of indices and counts

    ! orderings
    integer, parameter, public :: qdldl_order_natural = 0
    integer, parameter, public :: qdldl_order_user    = 1
    integer, parameter, public :: qdldl_order_rcm     = 2
    integer, parameter, public :: qdldl_order_amd     = 3
    integer, parameter, public :: qdldl_order_default = -1   ! AMD if built in, else RCM

    ! status codes (0 = success; named constants for every failure)
    integer, parameter, public :: qdldl_success = 0
    ! qdldl_error_not_analyzed, qdldl_error_not_factored, qdldl_error_invalid_input,
    ! qdldl_error_zero_pivot, qdldl_error_out_of_memory, qdldl_error_not_finite, ...

    type, public :: qdldl_type
        ! options (components with defaults)
        integer         :: ordering = qdldl_order_default
        real(qdldl_wp)  :: zero_pivot_tol = 0.0_qdldl_wp   ! 0: exactly zero only (upstream)
        logical         :: regularize = .false.
        real(qdldl_wp)  :: reg_eps = ..., reg_delta = ...
        integer         :: max_refine = 0                   ! 0: no refinement (upstream)
        ! results
        integer(qdldl_ip) :: n_positive, n_negative, n_zero, n_regularized
        integer(qdldl_ip) :: nnz_l
        ! private: the CSC of P A P^T, the map, perm/iperm, etree, Lnz, L, D, Dinv,
        !          the expected signs, the work arrays
    contains
        procedure :: analyze, factor, solve, multiply, inertia, destroy
        procedure :: set_signs        ! the expected pivot signs, for regularization
    end type
end module
```

Low-level procedures (`qdldl_etree`, `qdldl_factor`, `qdldl_solve`,
`qdldl_lsolve`, `qdldl_ltsolve`) stay public too, for callers that manage
their own storage.

---

## 6. Pitfalls to watch

- **Index conversion** (stage 1): every `-1` marker, every `< 0` test,
  every loop that runs to `n-1` or uses `Ap(k+1)`.
- **Upper triangle only**: the C code rejects entries below the diagonal;
  the coordinate front end must swap them, not drop them.
- **Missing diagonals** behave as zero pivots: insert them.
- **Duplicates** must be summed into one CSC entry, and the map must send
  every duplicate to the same position (scatter with `+`, after zeroing).
- **Integer overflow** of the nonzeros of \(L\) with 32-bit indices:
  detect it in `etree` (upstream does) and say to build with `INT64`.
- **Stability**: without pivoting, a matrix that isn't quasi-definite can
  factor with huge growth and give a wrong inertia; say so in the docs, and
  offer the regularization and growth check.
- **Quadruple precision** is slow; fine for correctness, not benchmarks.

---

## 7. Validation matrix (what "done" means)

| Check | Matrices | Pass when |
|---|---|---|
| Upstream tests | upstream's cases | same results and error codes |
| Solution accuracy | random quasi-definite, KKT from Laplacians | \(\lVert Ax-b \rVert / (\lVert A \rVert \lVert x \rVert + \lVert b \rVert) \le c\,\epsilon\) |
| Inertia | quasi-definite with known block sizes | `n_positive`, `n_negative` exact |
| Refactor | same pattern, new values | same as a fresh analysis |
| Orderings | every ordering | same solution; AMD/RCM fill below natural |
| Regularization | singular and non-quasi-definite KKT | completes; refined solution of the unregularized system where it exists |
| Kinds | REAL32, REAL64, REAL128, INT64 | all tests pass |
| Memory | huge `n` | out-of-memory status, no crash |
| Stack | the whole library | no automatic arrays or problem-sized temporaries |

---

## 8. Licensing

- QDLDL is Apache 2.0. A port is a derivative work: keep the licence, add a
  `NOTICE` crediting the original authors, and mark the files as changed
  (a line at the top of each ported file). The simplest course is to
  release the whole package under Apache 2.0.
- If SuiteSparse's Fortran AMD is included, keep its BSD licence text in
  those files and mention it in `NOTICE`.
- An Apache 2.0 dependency can be used by an MIT project such as sqpopt;
  sqpopt then has to pass on the notice if it redistributes the package's
  code. (Not legal advice; check if it matters.)

---

## 9. The first user: sqpopt

sqpopt's only sparse factorization goes through one type,
`sqpopt_symmetric_solver_type` (`src/sqpopt_symmetric_solver_module.F90`),
which today wraps MUMPS and only exists in a build with MUMPS. Its interface
is what this package's `qdldl_type` must be able to sit behind:

| sqpopt needs | qdldl_type provides |
|---|---|
| `initialize(n, irow, icol, ok, threads)`: the pattern of one triangle in coordinate form, duplicates added | `analyze` with coordinate input (stage 2); `threads` ignored |
| `factor(val, ok)`: values in the pattern's order; `n_negative`, `n_null` | `factor`, `inertia` (stages 2 and 4) |
| `solve(b, ok, refine)`: in place, with up to 2 steps of iterative refinement, failing if not finite | `solve` with `max_refine` (stage 4) |
| `multiply(x, y)`: \(y = Ax\) with the last values | `multiply` (stage 2) |
| `out_of_memory` | the out-of-memory status |
| any real kind (sqpopt builds in REAL32/64/128) | the kinds of stage 5; MUMPS only gives double |

How sqpopt would use it, by how quasi-definite each matrix is:

- **Direct least squares** (`sqpopt_least_squares_module`): the matrix
  \(\begin{bmatrix} I & J_S^T \\ J_S & -\epsilon I \end{bmatrix}\) is
  quasi-definite already; QDLDL is exact there.
- **Direct QP** with BFGS or L-BFGS Hessians: \(\begin{bmatrix} H & J_w^T \\
  J_w & 0 \end{bmatrix}\) with \(H\) positive definite; quasi-definite with
  static regularization \(-\epsilon I\) in the second block, plus refinement.
- **Inertia control** with the exact or SR1 Hessian: \(H\) may be
  indefinite. The safe use requires \(H + \delta I\) positive definite on the
  whole space (dynamic regularization tells when it isn't, and sqpopt raises
  \(\delta\)), which is stricter than what inertia control tests with MUMPS
  (positive definite on the null space of the working set), so shifts would
  be larger. This is the case to measure before relying on it.

None of that integration is part of this package; it is listed so the
interface supports it. In sqpopt it would be a roadmap item of its own (a
second backend of `sqpopt_symmetric_solver_type`, available in the default
build), alongside F25 (SPRAL SSIDS) in `plan/ROADMAP.md`.

---

## References

- R. J. Vanderbei, *Symmetric quasidefinite matrices*, SIAM J. Optim.
  5(1), 100–113, 1995.
- T. A. Davis, *Algorithm 849: A concise sparse Cholesky factorization
  package*, ACM TOMS 31(4), 587–591, 2005.
- P. R. Amestoy, T. A. Davis, I. S. Duff, *An approximate minimum degree
  ordering algorithm*, SIAM J. Matrix Anal. Appl. 17(4), 886–905, 1996;
  and *Algorithm 837: AMD*, ACM TOMS 30(3), 381–388, 2004.
- B. Stellato, G. Banjac, P. Goulart, A. Bemporad, S. Boyd, *OSQP: an
  operator splitting solver for quadratic programs*, Math. Prog. Comp. 12,
  637–672, 2020 (QDLDL's use, and its regularization and refinement).
- QDLDL: <https://github.com/osqp/qdldl>. Clarabel (its Rust port of
  QDLDL): <https://github.com/oxfordcontrol/Clarabel.rs>.
