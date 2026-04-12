#include "sparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -----------------------------------------------------------------------
 * SparseMatBuilder
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
    for (int i = 0; i < B->nrows; i++) free(B->rows[i]);
    free(B->rows);
    free(B->lens);
    free(B->caps);
    free(B);
}

int smat_builder_add(SparseMatBuilder *B, int row, int col, Rat val)
{
    if (row < 0 || row >= B->nrows || col < 0 || col >= B->ncols) return -1;
    if (rat_is_zero(val)) return 0;

    SEntry *r = B->rows[row];
    int len = B->lens[row];

    /* Linear scan for existing entry. */
    for (int k = 0; k < len; k++) {
        if (r[k].col == col) {
            Rat sum;
            if (rat_add(r[k].val, val, &sum) != RAT_OK) return RAT_OVERFLOW;
            r[k].val = sum;
            return 0;
        }
    }
    /* New entry: grow if needed. */
    if (len == B->caps[row]) {
        int new_cap = B->caps[row] ? B->caps[row] * 2 : ROW_INIT_CAP;
        SEntry *nr = realloc(B->rows[row], (size_t)new_cap * sizeof(SEntry));
        if (!nr) return -1;
        B->rows[row] = nr;
        B->caps[row] = new_cap;
        r = nr;
    }
    r[len].col = col;
    r[len].val = val;
    B->lens[row]++;
    return 0;
}

/* -----------------------------------------------------------------------
 * SparseMat (CSR)
 * ----------------------------------------------------------------------- */

/* Compare SEntry by column (for qsort). */
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

    /* Count total nonzeros. */
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
        /* Copy and sort by column. */
        SEntry *tmp = malloc((size_t)len * sizeof(SEntry));
        if (!tmp && len > 0) {
            smat_free(M); return NULL;
        }
        memcpy(tmp, B->rows[i], (size_t)len * sizeof(SEntry));
        qsort(tmp, (size_t)len, sizeof(SEntry), sentrycmp);
        for (int k = 0; k < len; k++) {
            M->col_ind[pos] = tmp[k].col;
            M->vals[pos]    = tmp[k].val;
            pos++;
        }
        free(tmp);
    }
    M->row_ptr[B->nrows] = pos;
    return M;
}

void smat_free(SparseMat *M)
{
    if (!M) return;
    free(M->row_ptr);
    free(M->col_ind);
    free(M->vals);
    free(M);
}

SparseMat *smat_from_dense_int(const int *A, int nrows, int ncols)
{
    SparseMatBuilder *B = smat_builder_new(nrows, ncols);
    if (!B) return NULL;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            int v = A[i * ncols + j];
            if (v != 0) {
                smat_builder_add(B, i, j, rat_int(v));
            }
        }
    }
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
    if (!C->row_ptr || !C->col_ind || !C->vals) {
        smat_free(C); return NULL;
    }
    memcpy(C->row_ptr, M->row_ptr, (size_t)(M->nrows + 1) * sizeof(int));
    memcpy(C->col_ind, M->col_ind, (size_t)nnz * sizeof(int));
    memcpy(C->vals,    M->vals,    (size_t)nnz * sizeof(Rat));
    return C;
}

Rat smat_get(const SparseMat *M, int row, int col)
{
    int lo = M->row_ptr[row], hi = M->row_ptr[row + 1];
    /* Binary search. */
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (M->col_ind[mid] < col) lo = mid + 1;
        else if (M->col_ind[mid] > col) hi = mid;
        else return M->vals[mid];
    }
    return rat_zero();
}

void smat_print(const SparseMat *M)
{
    printf("SparseMat %d×%d:\n", M->nrows, M->ncols);
    for (int i = 0; i < M->nrows; i++) {
        printf("  row %d:", i);
        for (int k = M->row_ptr[i]; k < M->row_ptr[i + 1]; k++) {
            printf(" [%d]=", M->col_ind[k]);
            rat_print(M->vals[k]);
        }
        printf("\n");
    }
}
