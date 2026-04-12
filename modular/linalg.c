#include "linalg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -----------------------------------------------------------------------
 * DenseMat allocation
 *
 * Every cell is initialised via rat_zero() which for GMP calls mpq_init.
 * dmat_free calls rat_clear on every cell (no-op for builtin/fp, mpq_clear
 * for GMP).  No memcpy is ever used on Rat arrays.
 * ----------------------------------------------------------------------- */

DenseMat *dmat_new(int nrows, int ncols)
{
    DenseMat *M = malloc(sizeof(DenseMat));
    if (!M) return NULL;
    M->nrows = nrows;
    M->ncols = ncols;
    M->data  = malloc((size_t)(nrows * ncols) * sizeof(Rat));
    if (!M->data) { free(M); return NULL; }
    for (int i = 0; i < nrows * ncols; i++)
        M->data[i] = rat_zero();   /* also calls mpq_init for GMP */
    return M;
}

void dmat_free(DenseMat *M)
{
    if (!M) return;
    for (int i = 0; i < M->nrows * M->ncols; i++)
        rat_clear(&M->data[i]);
    free(M->data);
    free(M);
}

/* Element-wise copy — never memcpy, to keep GMP mpq_t ownership clean. */
DenseMat *dmat_copy(const DenseMat *M)
{
    DenseMat *C = dmat_new(M->nrows, M->ncols);
    if (!C) return NULL;
    for (int i = 0; i < M->nrows * M->ncols; i++)
        rat_set(&C->data[i], M->data[i]);
    return C;
}

DenseMat *dmat_from_sparse(const SparseMat *S)
{
    DenseMat *M = dmat_new(S->nrows, S->ncols);
    if (!M) return NULL;
    for (int i = 0; i < S->nrows; i++)
        for (int k = S->row_ptr[i]; k < S->row_ptr[i + 1]; k++)
            rat_set(dmat_at(M, i, S->col_ind[k]), S->vals[k]);
    return M;
}

DenseMat *dmat_identity(int n)
{
    DenseMat *I = dmat_new(n, n);
    if (!I) return NULL;
    for (int i = 0; i < n; i++)
        rat_set_one(dmat_at(I, i, i));
    return I;
}

void dmat_print(const DenseMat *M)
{
    printf("DenseMat %d\xc3\x97%d:\n", M->nrows, M->ncols);
    for (int i = 0; i < M->nrows; i++) {
        printf(" [");
        for (int j = 0; j < M->ncols; j++) {
            rat_print(dmat_get(M, i, j));
            if (j + 1 < M->ncols) printf(", ");
        }
        printf("]\n");
    }
}

/* -----------------------------------------------------------------------
 * Swap two matrix rows, columns [col_start, ncols).
 * Uses a temporary Rat that is properly init'd and cleared.
 * ----------------------------------------------------------------------- */
static void swap_rows(DenseMat *M, int r1, int r2, int col_start,
                      DenseMat *also, int also_col_start)
{
    Rat tmp = rat_zero();
    int nc = M->ncols;
    for (int c = col_start; c < nc; c++) {
        rat_set(&tmp, *dmat_at(M, r1, c));
        rat_set(dmat_at(M, r1, c), *dmat_at(M, r2, c));
        rat_set(dmat_at(M, r2, c), tmp);
    }
    if (also) {
        int anc = also->ncols;
        for (int c = also_col_start; c < anc; c++) {
            rat_set(&tmp, *dmat_at(also, r1, c));
            rat_set(dmat_at(also, r1, c), *dmat_at(also, r2, c));
            rat_set(dmat_at(also, r2, c), tmp);
        }
    }
    rat_clear(&tmp);
}

/* -----------------------------------------------------------------------
 * RREF via Gaussian elimination over Q.
 * ----------------------------------------------------------------------- */
int dmat_rref(DenseMat *M, int *pivot_col, DenseMat *apply_to)
{
    int nrows = M->nrows, ncols = M->ncols;
    int pivot_row = 0, rank = 0;

    for (int i = 0; i < nrows; i++) pivot_col[i] = -1;

    Rat piv_inv = rat_zero();
    Rat prod    = rat_zero();
    Rat sub     = rat_zero();
    Rat term    = rat_zero();

    for (int col = 0; col < ncols && pivot_row < nrows; col++) {
        /* Find sparsest non-zero pivot candidate at or below pivot_row. */
        int best = -1, best_nnz = ncols + 1;
        for (int r = pivot_row; r < nrows; r++) {
            if (rat_is_zero(dmat_get(M, r, col))) continue;
            int nnz = 0;
            for (int c2 = 0; c2 < ncols; c2++)
                if (!rat_is_zero(dmat_get(M, r, c2))) nnz++;
            if (best < 0 || nnz < best_nnz) { best = r; best_nnz = nnz; }
        }
        if (best < 0) continue;

        if (best != pivot_row)
            swap_rows(M, pivot_row, best, 0, apply_to, 0);

        /* Scale pivot row so pivot == 1. */
        if (rat_inv(dmat_get(M, pivot_row, col), &piv_inv) != RAT_OK) goto fail;

        for (int c2 = 0; c2 < ncols; c2++) {
            if (rat_is_zero(dmat_get(M, pivot_row, c2))) continue;
            if (rat_mul(dmat_get(M, pivot_row, c2), piv_inv, &prod) != RAT_OK) goto fail;
            rat_set(dmat_at(M, pivot_row, c2), prod);
        }
        if (apply_to) {
            for (int c2 = 0; c2 < apply_to->ncols; c2++) {
                if (rat_mul(dmat_get(apply_to, pivot_row, c2), piv_inv, &prod) != RAT_OK) goto fail;
                rat_set(dmat_at(apply_to, pivot_row, c2), prod);
            }
        }

        /* Eliminate all other rows.
         * factor must be a deep-copied Rat (not an alias of M's cell), because
         * the c2 loop writes back to M[r][col] which would corrupt the alias. */
        Rat factor = rat_zero();
        for (int r = 0; r < nrows; r++) {
            if (r == pivot_row) continue;
            rat_set(&factor, dmat_get(M, r, col)); /* deep copy, not alias */
            if (rat_is_zero(factor)) continue;
            for (int c2 = 0; c2 < ncols; c2++) {
                if (rat_mul(factor, dmat_get(M, pivot_row, c2), &term) != RAT_OK)
                    { rat_clear(&factor); goto fail; }
                if (rat_sub(dmat_get(M, r, c2), term, &sub) != RAT_OK)
                    { rat_clear(&factor); goto fail; }
                rat_set(dmat_at(M, r, c2), sub);
            }
            if (apply_to) {
                for (int c2 = 0; c2 < apply_to->ncols; c2++) {
                    if (rat_mul(factor, dmat_get(apply_to, pivot_row, c2), &term) != RAT_OK)
                        { rat_clear(&factor); goto fail; }
                    if (rat_sub(dmat_get(apply_to, r, c2), term, &sub) != RAT_OK)
                        { rat_clear(&factor); goto fail; }
                    rat_set(dmat_at(apply_to, r, c2), sub);
                }
            }
        }
        rat_clear(&factor);

        pivot_col[pivot_row] = col;
        pivot_row++;
        rank++;
    }

    rat_clear(&piv_inv); rat_clear(&prod); rat_clear(&sub); rat_clear(&term);
    return rank;

fail:
    rat_clear(&piv_inv); rat_clear(&prod); rat_clear(&sub); rat_clear(&term);
    return -1;
}

/* -----------------------------------------------------------------------
 * Kernel basis
 * ----------------------------------------------------------------------- */
DenseMat *dmat_kernel(const DenseMat *M)
{
    DenseMat *W = dmat_copy(M);
    if (!W) return NULL;

    int *pivot_col = malloc((size_t)W->nrows * sizeof(int));
    if (!pivot_col) { dmat_free(W); return NULL; }

    int rank = dmat_rref(W, pivot_col, NULL);
    if (rank < 0) { free(pivot_col); dmat_free(W); return NULL; }

    int ncols   = W->ncols;
    int nullity = ncols - rank;

    char *is_pivot = calloc((size_t)ncols, 1);
    if (!is_pivot) { free(pivot_col); dmat_free(W); return NULL; }
    for (int r = 0; r < W->nrows; r++)
        if (pivot_col[r] >= 0) is_pivot[pivot_col[r]] = 1;

    DenseMat *ker = dmat_new(nullity, ncols);
    if (!ker) { free(is_pivot); free(pivot_col); dmat_free(W); return NULL; }

    Rat neg_val = rat_zero();
    int krow = 0;
    for (int fc = 0; fc < ncols; fc++) {
        if (is_pivot[fc]) continue;
        rat_set_one(dmat_at(ker, krow, fc));
        for (int r = 0; r < rank; r++) {
            int pc = pivot_col[r];
            if (pc < 0) continue;
            if (rat_neg(dmat_get(W, r, fc), &neg_val) != RAT_OK) {
                dmat_free(ker); ker = NULL; goto cleanup;
            }
            rat_set(dmat_at(ker, krow, pc), neg_val);
        }
        krow++;
    }

cleanup:
    rat_clear(&neg_val);
    free(is_pivot);
    free(pivot_col);
    dmat_free(W);
    return ker;
}

/* -----------------------------------------------------------------------
 * Matrix multiply A (m×k) * B (k×n) → C (m×n)
 * ----------------------------------------------------------------------- */
DenseMat *dmat_mul(const DenseMat *A, const DenseMat *B)
{
    if (A->ncols != B->nrows) return NULL;
    int m = A->nrows, k = A->ncols, n = B->ncols;
    DenseMat *C = dmat_new(m, n);
    if (!C) return NULL;

    Rat sum  = rat_zero();
    Rat prod = rat_zero();
    Rat acc  = rat_zero();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            rat_set_zero(&sum);
            for (int l = 0; l < k; l++) {
                if (rat_mul(dmat_get(A, i, l), dmat_get(B, l, j), &prod) != RAT_OK)
                    goto fail;
                if (rat_add(sum, prod, &acc) != RAT_OK)
                    goto fail;
                rat_set(&sum, acc);
            }
            rat_set(dmat_at(C, i, j), sum);
        }
    }

    rat_clear(&sum); rat_clear(&prod); rat_clear(&acc);
    return C;

fail:
    rat_clear(&sum); rat_clear(&prod); rat_clear(&acc);
    dmat_free(C);
    return NULL;
}

/* -----------------------------------------------------------------------
 * Vertical stack [A; B]
 * ----------------------------------------------------------------------- */
DenseMat *dmat_vstack(const DenseMat *A, const DenseMat *B)
{
    if (A->ncols != B->ncols) return NULL;
    DenseMat *C = dmat_new(A->nrows + B->nrows, A->ncols);
    if (!C) return NULL;
    for (int i = 0; i < A->nrows; i++)
        for (int j = 0; j < A->ncols; j++)
            rat_set(dmat_at(C, i, j), dmat_get(A, i, j));
    for (int i = 0; i < B->nrows; i++)
        for (int j = 0; j < B->ncols; j++)
            rat_set(dmat_at(C, A->nrows + i, j), dmat_get(B, i, j));
    return C;
}

/* -----------------------------------------------------------------------
 * Image basis — stub.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_image(const DenseMat *M)
{
    (void)M;
    return NULL;
}
