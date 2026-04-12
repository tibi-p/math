#include "sparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -----------------------------------------------------------------------
 * SparseMatBuilder
 *
 * GMP ownership rules in this file:
 *  - Every Rat stored in an SEntry has its OWN independent mpq_t.
 *  - smat_builder_free calls rat_clear on every stored Rat.
 *  - smat_freeze initialises new Rats for M->vals (does NOT alias builder).
 *  - smat_free calls rat_clear on every M->vals entry.
 *  - smat_builder_add takes `val` by value and treats it as a move:
 *    it rat_set-copies the VALUE then rat_clears the by-value arg.
 * ----------------------------------------------------------------------- */

#define ROW_INIT_CAP 4

SparseMatBuilder *smat_builder_new(int nrows, int ncols)
{
    SparseMatBuilder *B = malloc(sizeof(SparseMatBuilder));
    if (!B) return NULL;
    B->nrows = nrows;
    B->ncols = ncols;
    B->rows  = calloc((size_t)nrows, sizeof(SEntry *));
    B->lens  = calloc((size_t)nrows, sizeof(int));
    B->caps  = calloc((size_t)nrows, sizeof(int));
    if (!B->rows || !B->lens || !B->caps) {
        free(B->rows); free(B->lens); free(B->caps); free(B);
        return NULL;
    }
    return B;
}

void smat_builder_free(SparseMatBuilder *B)
{
    if (!B) return;
    for (int i = 0; i < B->nrows; i++) {
        for (int j = 0; j < B->lens[i]; j++)
            rat_clear(&B->rows[i][j].val);
        free(B->rows[i]);
    }
    free(B->rows);
    free(B->lens);
    free(B->caps);
    free(B);
}

int smat_builder_add(SparseMatBuilder *B, int row, int col, Rat val)
{
    if (row < 0 || row >= B->nrows || col < 0 || col >= B->ncols) {
        rat_clear(&val); return -1;
    }
    if (rat_is_zero(val)) { rat_clear(&val); return 0; }

    SEntry *r = B->rows[row];
    int len = B->lens[row];

    /* Linear scan for existing entry — merge by adding. */
    for (int k = 0; k < len; k++) {
        if (r[k].col == col) {
            Rat sum = rat_zero();
            int rc = rat_add(r[k].val, val, &sum);
            rat_clear(&val);
            if (rc != RAT_OK) { rat_clear(&sum); return RAT_OVERFLOW; }
            rat_set(&r[k].val, sum);
            rat_clear(&sum);
            return 0;
        }
    }

    /* New entry: grow row if needed. */
    if (len == B->caps[row]) {
        int new_cap = B->caps[row] ? B->caps[row] * 2 : ROW_INIT_CAP;
        SEntry *nr = realloc(B->rows[row], (size_t)new_cap * sizeof(SEntry));
        if (!nr) { rat_clear(&val); return -1; }
        B->rows[row] = nr;
        B->caps[row] = new_cap;
        r = nr;
    }
    r[len].col = col;
    r[len].val = rat_zero();           /* init a fresh Rat          */
    rat_set(&r[len].val, val);         /* copy value from arg       */
    rat_clear(&val);                   /* release the by-value copy */
    B->lens[row]++;
    return 0;
}

/* -----------------------------------------------------------------------
 * SparseMat (CSR)
 * ----------------------------------------------------------------------- */

static int sentrycmp(const void *a, const void *b)
{
    return ((const SEntry *)a)->col - ((const SEntry *)b)->col;
}

SparseMat *smat_freeze(const SparseMatBuilder *B)
{
    SparseMat *M = malloc(sizeof(SparseMat));
    if (!M) return NULL;
    M->nrows = B->nrows;
    M->ncols = B->ncols;
    M->row_ptr = malloc((size_t)(B->nrows + 1) * sizeof(int));
    if (!M->row_ptr) { free(M); return NULL; }

    int nnz = 0;
    for (int i = 0; i < B->nrows; i++) nnz += B->lens[i];

    M->col_ind = malloc((size_t)nnz * sizeof(int));
    M->vals    = malloc((size_t)nnz * sizeof(Rat));
    if (!M->col_ind || !M->vals) {
        free(M->row_ptr); free(M->col_ind); free(M->vals); free(M);
        return NULL;
    }

    int pos = 0;
    for (int i = 0; i < B->nrows; i++) {
        M->row_ptr[i] = pos;
        int len = B->lens[i];

        /* Sort a temporary index array by column (avoids aliasing Rat values). */
        int *order = malloc((size_t)len * sizeof(int));
        if (!order && len > 0) { smat_free(M); return NULL; }
        for (int k = 0; k < len; k++) order[k] = k;
        /* Simple insertion sort on order[] by B->rows[i][order[k]].col */
        for (int k = 1; k < len; k++) {
            int tmp = order[k];
            int j = k - 1;
            while (j >= 0 && B->rows[i][order[j]].col > B->rows[i][tmp].col)
                { order[j + 1] = order[j]; j--; }
            order[j + 1] = tmp;
        }
        for (int k = 0; k < len; k++) {
            M->col_ind[pos] = B->rows[i][order[k]].col;
            M->vals[pos]    = rat_zero();                      /* own mpq_t */
            rat_set(&M->vals[pos], B->rows[i][order[k]].val); /* copy value */
            pos++;
        }
        free(order);
    }
    M->row_ptr[B->nrows] = pos;
    return M;
}

void smat_free(SparseMat *M)
{
    if (!M) return;
    int nnz = M->row_ptr ? M->row_ptr[M->nrows] : 0;
    for (int k = 0; k < nnz; k++)
        rat_clear(&M->vals[k]);
    free(M->row_ptr);
    free(M->col_ind);
    free(M->vals);
    free(M);
}

SparseMat *smat_from_dense_int(const int *A, int nrows, int ncols)
{
    SparseMatBuilder *B = smat_builder_new(nrows, ncols);
    if (!B) return NULL;
    for (int i = 0; i < nrows; i++)
        for (int j = 0; j < ncols; j++)
            if (A[i * ncols + j] != 0)
                smat_builder_add(B, i, j, rat_int(A[i * ncols + j]));
    SparseMat *M = smat_freeze(B);
    smat_builder_free(B);
    return M;
}

SparseMat *smat_copy(const SparseMat *M)
{
    SparseMat *C = malloc(sizeof(SparseMat));
    if (!C) return NULL;
    C->nrows = M->nrows;
    C->ncols = M->ncols;
    int nnz = M->row_ptr[M->nrows];
    C->row_ptr = malloc((size_t)(M->nrows + 1) * sizeof(int));
    C->col_ind = malloc((size_t)nnz * sizeof(int));
    C->vals    = malloc((size_t)nnz * sizeof(Rat));
    if (!C->row_ptr || !C->col_ind || !C->vals) { smat_free(C); return NULL; }
    memcpy(C->row_ptr, M->row_ptr, (size_t)(M->nrows + 1) * sizeof(int));
    memcpy(C->col_ind, M->col_ind, (size_t)nnz * sizeof(int));
    for (int k = 0; k < nnz; k++) {
        C->vals[k] = rat_zero();
        rat_set(&C->vals[k], M->vals[k]);
    }
    return C;
}

/* NOTE: smat_get returns a Rat by value.  For builtin/fp this is safe.
 * For GMP it returns a struct-copy of M->vals[mid].q — treat as read-only
 * and do NOT call rat_clear on the result. */
Rat smat_get(const SparseMat *M, int row, int col)
{
    int lo = M->row_ptr[row], hi = M->row_ptr[row + 1];
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if      (M->col_ind[mid] < col) lo = mid + 1;
        else if (M->col_ind[mid] > col) hi = mid;
        else                             return M->vals[mid];
    }
    return rat_zero();
}

void smat_print(const SparseMat *M)
{
    printf("SparseMat %d\xc3\x97%d:\n", M->nrows, M->ncols);
    for (int i = 0; i < M->nrows; i++) {
        printf("  row %d:", i);
        for (int k = M->row_ptr[i]; k < M->row_ptr[i + 1]; k++) {
            printf(" [%d]=", M->col_ind[k]);
            rat_print(M->vals[k]);
        }
        printf("\n");
    }
}
