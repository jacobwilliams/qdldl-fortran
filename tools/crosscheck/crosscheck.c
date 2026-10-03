/* The C side of the cross-check (see run.sh): reads matrices.txt, factors and
 * solves each matrix with upstream QDLDL, and writes c_out.txt. */
#include <stdio.h>
#include <stdlib.h>
#include "qdldl.h"

int main(void) {
    FILE *in = fopen("matrices.txt", "r"), *out = fopen("c_out.txt", "w");
    int nmat, m;
    if (!in || !out || fscanf(in, "%d", &nmat) != 1) return 1;
    for (m = 0; m < nmat; m++) {
        QDLDL_int n, nnz, i, sumLnz, r;
        if (fscanf(in, "%d %d", &n, &nnz) != 2) return 1;
        QDLDL_int *Ap = malloc(sizeof(QDLDL_int) * (n + 1)), *Ai = malloc(sizeof(QDLDL_int) * nnz);
        QDLDL_float *Ax = malloc(sizeof(QDLDL_float) * nnz), *x = malloc(sizeof(QDLDL_float) * n);
        for (i = 0; i <= n; i++) if (fscanf(in, "%d", &Ap[i]) != 1) return 1;
        for (i = 0; i < nnz; i++) if (fscanf(in, "%d", &Ai[i]) != 1) return 1;
        for (i = 0; i < nnz; i++) if (fscanf(in, "%lf", &Ax[i]) != 1) return 1;
        for (i = 0; i < n; i++) if (fscanf(in, "%lf", &x[i]) != 1) return 1;

        QDLDL_int *etree = malloc(sizeof(QDLDL_int) * n), *Lnz = malloc(sizeof(QDLDL_int) * n);
        QDLDL_int *Lp = malloc(sizeof(QDLDL_int) * (n + 1)), *iwork = malloc(sizeof(QDLDL_int) * 3 * n);
        QDLDL_bool *bwork = malloc(sizeof(QDLDL_bool) * n);
        QDLDL_float *D = malloc(sizeof(QDLDL_float) * n), *Dinv = malloc(sizeof(QDLDL_float) * n);
        QDLDL_float *fwork = malloc(sizeof(QDLDL_float) * n);
        sumLnz = QDLDL_etree(n, Ap, Ai, iwork, Lnz, etree);
        QDLDL_int *Li = malloc(sizeof(QDLDL_int) * (sumLnz > 0 ? sumLnz : 1));
        QDLDL_float *Lx = malloc(sizeof(QDLDL_float) * (sumLnz > 0 ? sumLnz : 1));
        r = QDLDL_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork);
        QDLDL_solve(n, Lp, Li, Lx, Dinv, x);

        fprintf(out, "%d %d %d\n", n, sumLnz, r);
        for (i = 0; i < n; i++) fprintf(out, "%d%c", etree[i], i < n - 1 ? ' ' : '\n');
        for (i = 0; i <= n; i++) fprintf(out, "%d%c", Lp[i], i < n ? ' ' : '\n');
        for (i = 0; i < sumLnz; i++) fprintf(out, "%d%c", Li[i], i < sumLnz - 1 ? ' ' : '\n');
        for (i = 0; i < sumLnz; i++) fprintf(out, "%.17e%c", Lx[i], i < sumLnz - 1 ? ' ' : '\n');
        for (i = 0; i < n; i++) fprintf(out, "%.17e%c", D[i], i < n - 1 ? ' ' : '\n');
        for (i = 0; i < n; i++) fprintf(out, "%.17e%c", x[i], i < n - 1 ? ' ' : '\n');
        free(Ap); free(Ai); free(Ax); free(x); free(etree); free(Lnz); free(Lp); free(iwork);
        free(bwork); free(D); free(Dinv); free(fwork); free(Li); free(Lx);
    }
    fclose(in);
    fclose(out);
    return 0;
}
