#include "linalg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -----------------------------------------------------------------------
 * DenseMat allocation
 * ----------------------------------------------------------------------- */

DenseMat *dmat_new(int nrows, int ncols)
{
    DenseMat *M = malloc(sizeof(DenseMat));
    if (!M) return NULL;
    M->nrows = nrows;
    M->ncols = ncols;
    M->data  = calloc((size_t)(nrows * ncols), sizeof(Rat));
    if (!M->data) { free(M); return NULL; }
    /* calloc gives 0-bytes; set every entry to rat_zero() = {0,1}.
     * Since {0,0} is zero-bytes but den==0 is invalid, we must init. */
    for (int i = 0; i < nrows * ncols; i++)
        M->data[i] = rat_zero();
    return M;
}

void dmat_free(DenseMat *M)
{
    if (!M) return;
    free(M->data);
    free(M);
}

DenseMat *dmat_copy(const DenseMat *M)
{
    DenseMat *C = malloc(sizeof(DenseMat));
    if (!C) return NULL;
    C->nrows = M->nrows;
    C->ncols = M->ncols;
    size_t sz = (size_t)(M->nrows * M->ncols) * sizeof(Rat);
    C->data = malloc(sz);
    if (!C->data) { free(C); return NULL; }
    memcpy(C->data, M->data, sz);
    return C;
}

DenseMat *dmat_from_sparse(const SparseMat *S)
{
    DenseMat *M = dmat_new(S->nrows, S->ncols);
    if (!M) return NULL;
    for (int i = 0; i < S->nrows; i++) {
        for (int k = S->row_ptr[i]; k < S->row_ptr[i + 1]; k++) {
            *dmat_at(M, i, S->col_ind[k]) = S->vals[k];
        }
    }
    return M;
}

DenseMat *dmat_identity(int n)
{
    DenseMat *I = dmat_new(n, n);
    if (!I) return NULL;
    for (int i = 0; i < n; i++)
        *dmat_at(I, i, i) = rat_one();
    return I;
}

void dmat_print(const DenseMat *M)
{
    printf("DenseMat %d×%d:\n", M->nrows, M->ncols);
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
 * RREF via Gaussian elimination over Q.
 *
 * For each column left-to-right, find the first non-zero entry in that
 * column at or below the current pivot row (choosing sparsest row first
 * for numerical/fill-in efficiency), swap it up, normalise the pivot row,
 * then eliminate above and below.
 *
 * apply_to receives the same row operations (pass NULL to skip).
 * pivot_col[i] = j means row i has its leading 1 in column j.
 * Rows with no pivot get pivot_col[i] = -1.
 * Returns the rank.
 * ----------------------------------------------------------------------- */
int dmat_rref(DenseMat *M, int *pivot_col, DenseMat *apply_to)
{
    int nrows = M->nrows, ncols = M->ncols;
    int pivot_row = 0;
    int rank = 0;

    for (int i = 0; i < nrows; i++)
        pivot_col[i] = -1;

    for (int col = 0; col < ncols && pivot_row < nrows; col++) {
        /* Find pivot: first non-zero in this column at or below pivot_row.
         * Among candidates, prefer the sparsest row (fewest nonzeros). */
        int best = -1, best_nnz = ncols + 1;
        for (int r = pivot_row; r < nrows; r++) {
            if (rat_is_zero(dmat_get(M, r, col))) continue;
            /* Count nonzeros in this row. */
            int nnz = 0;
            for (int c2 = 0; c2 < ncols; c2++)
                if (!rat_is_zero(dmat_get(M, r, c2))) nnz++;
            if (best < 0 || nnz < best_nnz) {
                best = r; best_nnz = nnz;
            }
        }
        if (best < 0) continue;  /* entire column below pivot_row is zero */

        /* Swap best row with pivot_row. */
        if (best != pivot_row) {
            for (int c2 = 0; c2 < ncols; c2++) {
                Rat tmp = *dmat_at(M, pivot_row, c2);
                *dmat_at(M, pivot_row, c2) = *dmat_at(M, best, c2);
                *dmat_at(M, best, c2) = tmp;
            }
            if (apply_to) {
                for (int c2 = 0; c2 < apply_to->ncols; c2++) {
                    Rat tmp = *dmat_at(apply_to, pivot_row, c2);
                    *dmat_at(apply_to, pivot_row, c2) = *dmat_at(apply_to, best, c2);
                    *dmat_at(apply_to, best, c2) = tmp;
                }
            }
        }

        /* Scale pivot row so pivot entry == 1. */
        Rat piv = dmat_get(M, pivot_row, col);
        Rat piv_inv;
        if (rat_inv(piv, &piv_inv) != RAT_OK) return -1;  /* shouldn't happen */

        for (int c2 = 0; c2 < ncols; c2++) {
            Rat prod;
            if (!rat_is_zero(dmat_get(M, pivot_row, c2))) {
                if (rat_mul(dmat_get(M, pivot_row, c2), piv_inv, &prod) != RAT_OK)
                    return -1;
                *dmat_at(M, pivot_row, c2) = prod;
            }
        }
        if (apply_to) {
            for (int c2 = 0; c2 < apply_to->ncols; c2++) {
                Rat prod;
                if (rat_mul(dmat_get(apply_to, pivot_row, c2), piv_inv, &prod) != RAT_OK)
                    return -1;
                *dmat_at(apply_to, pivot_row, c2) = prod;
            }
        }

        /* Eliminate all other rows in this column. */
        for (int r = 0; r < nrows; r++) {
            if (r == pivot_row) continue;
            Rat factor = dmat_get(M, r, col);
            if (rat_is_zero(factor)) continue;
            for (int c2 = 0; c2 < ncols; c2++) {
                Rat sub, term;
                if (rat_mul(factor, dmat_get(M, pivot_row, c2), &term) != RAT_OK)
                    return -1;
                if (rat_sub(dmat_get(M, r, c2), term, &sub) != RAT_OK)
                    return -1;
                *dmat_at(M, r, c2) = sub;
            }
            if (apply_to) {
                for (int c2 = 0; c2 < apply_to->ncols; c2++) {
                    Rat sub, term;
                    if (rat_mul(factor, dmat_get(apply_to, pivot_row, c2), &term) != RAT_OK)
                        return -1;
                    if (rat_sub(dmat_get(apply_to, r, c2), term, &sub) != RAT_OK)
                        return -1;
                    *dmat_at(apply_to, r, c2) = sub;
                }
            }
        }

        pivot_col[pivot_row] = col;
        pivot_row++;
        rank++;
    }
    return rank;
}

/* -----------------------------------------------------------------------
 * Kernel basis
 *
 * Algorithm: augment M with identity [M | I], reduce M side to RREF,
 * read off the kernel from free-variable columns.
 *
 * Alternatively (cleaner): reduce M^T to RREF and read nullspace.
 * We use the direct method: after RREF of M, for each free column f,
 * the kernel vector has 1 at position f and -pivot_row_value at pivot
 * positions, 0 elsewhere.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_kernel(const DenseMat *M)
{
    DenseMat *W = dmat_copy(M);
    if (!W) return NULL;

    int *pivot_col = malloc((size_t)W->nrows * sizeof(int));
    if (!pivot_col) { dmat_free(W); return NULL; }

    int rank = dmat_rref(W, pivot_col, NULL);
    if (rank < 0) { free(pivot_col); dmat_free(W); return NULL; }

    int ncols = W->ncols;
    int nullity = ncols - rank;

    /* Mark pivot columns. */
    char *is_pivot = calloc((size_t)ncols, 1);
    /* pivot_row -> pivot_col mapping (already in pivot_col[r]). */
    /* Build pivot_col_to_row mapping. */
    int *col_to_pivrow = malloc((size_t)ncols * sizeof(int));
    if (!is_pivot || !col_to_pivrow) {
        free(is_pivot); free(col_to_pivrow);
        free(pivot_col); dmat_free(W); return NULL;
    }
    for (int c = 0; c < ncols; c++) col_to_pivrow[c] = -1;
    for (int r = 0; r < W->nrows; r++) {
        if (pivot_col[r] >= 0) {
            is_pivot[pivot_col[r]] = 1;
            col_to_pivrow[pivot_col[r]] = r;
        }
    }

    /* Kernel has dimension nullity; each free variable gives one basis vector. */
    DenseMat *ker = dmat_new(nullity, ncols);
    if (!ker) {
        free(is_pivot); free(col_to_pivrow);
        free(pivot_col); dmat_free(W); return NULL;
    }

    int krow = 0;
    for (int fc = 0; fc < ncols; fc++) {
        if (is_pivot[fc]) continue;
        /* Kernel basis vector for free column fc:
         * entry at fc = 1.
         * For each pivot column pc (in pivot row r): entry at pc = -W[r][fc]. */
        *dmat_at(ker, krow, fc) = rat_one();
        for (int r = 0; r < rank; r++) {
            int pc = pivot_col[r];
            if (pc < 0) continue;
            Rat neg_val;
            if (rat_neg(dmat_get(W, r, fc), &neg_val) != RAT_OK) {
                /* overflow; shouldn't happen for small N */
                dmat_free(ker); ker = NULL; goto cleanup;
            }
            *dmat_at(ker, krow, pc) = neg_val;
        }
        krow++;
    }

cleanup:
    free(is_pivot);
    free(col_to_pivrow);
    free(pivot_col);
    dmat_free(W);
    return ker;
}

/* -----------------------------------------------------------------------
 * Matrix multiply
 * ----------------------------------------------------------------------- */
DenseMat *dmat_mul(const DenseMat *A, const DenseMat *B)
{
    if (A->ncols != B->nrows) return NULL;
    int m = A->nrows, k = A->ncols, n = B->ncols;
    DenseMat *C = dmat_new(m, n);
    if (!C) return NULL;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Rat sum = rat_zero();
            for (int l = 0; l < k; l++) {
                Rat prod, s2;
                if (rat_mul(dmat_get(A, i, l), dmat_get(B, l, j), &prod) != RAT_OK)
                    { dmat_free(C); return NULL; }
                if (rat_add(sum, prod, &s2) != RAT_OK)
                    { dmat_free(C); return NULL; }
                sum = s2;
            }
            *dmat_at(C, i, j) = sum;
        }
    }
    return C;
}

/* -----------------------------------------------------------------------
 * Vertical stack [A; B]
 * ----------------------------------------------------------------------- */
DenseMat *dmat_vstack(const DenseMat *A, const DenseMat *B)
{
    if (A->ncols != B->ncols) return NULL;
    DenseMat *C = dmat_new(A->nrows + B->nrows, A->ncols);
    if (!C) return NULL;
    size_t row_sz = (size_t)A->ncols * sizeof(Rat);
    for (int i = 0; i < A->nrows; i++)
        memcpy(C->data + i * A->ncols, A->data + i * A->ncols, row_sz);
    for (int i = 0; i < B->nrows; i++)
        memcpy(C->data + (A->nrows + i) * A->ncols,
               B->data + i * B->ncols, row_sz);
    return C;
}

/* -----------------------------------------------------------------------
 * Image basis — not needed yet; stub.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_image(const DenseMat *M)
{
    (void)M;
    return NULL;  /* TODO */
}
