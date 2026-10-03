#include "qdldl.h"

int benchmark_c_etree(const QDLDL_int n, const QDLDL_int *Ap, const QDLDL_int *Ai,
                      QDLDL_int *work, QDLDL_int *Lnz, QDLDL_int *etree) {
    return QDLDL_etree(n, Ap, Ai, work, Lnz, etree);
}

QDLDL_int benchmark_c_factor(const QDLDL_int n, const QDLDL_int *Ap, const QDLDL_int *Ai,
                             const QDLDL_float *Ax, QDLDL_int *Lp, QDLDL_int *Li,
                             QDLDL_float *Lx, QDLDL_float *D, QDLDL_float *Dinv,
                             const QDLDL_int *Lnz, const QDLDL_int *etree,
                             QDLDL_bool *bwork, QDLDL_int *iwork, QDLDL_float *fwork) {
    return QDLDL_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree,
                        bwork, iwork, fwork);
}

void benchmark_c_solve(const QDLDL_int n, const QDLDL_int *Lp, const QDLDL_int *Li,
                       const QDLDL_float *Lx, const QDLDL_float *Dinv,
                       QDLDL_float *x) {
    QDLDL_solve(n, Lp, Li, Lx, Dinv, x);
}
