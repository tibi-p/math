#ifndef MODULAR_LINALG_H
#define MODULAR_LINALG_H

#include "sparse.h"
#include "rat.h"

/* -----------------------------------------------------------------------
 * Dense rational matrix — used internally for elimination.
 * Rows are stored contiguously: entry (i,j) is at data[i*ncols + j].
 * ----------------------------------------------------------------------- */
typedef struct {
    int  nrows;
    int  ncols;
    Rat *data;
} DenseMat;

DenseMat *dmat_new(int nrows, int ncols);    /* zero-initialised */
void      dmat_free(DenseMat *M);
DenseMat *dmat_copy(const DenseMat *M);
DenseMat *dmat_from_sparse(const SparseMat *S);
DenseMat *dmat_identity(int n);

static inline Rat *dmat_at(DenseMat *M, int r, int c)
{ return &M->data[r * M->ncols + c]; }

static inline Rat dmat_get(const DenseMat *M, int r, int c)
{ return M->data[r * M->ncols + c]; }

void dmat_print(const DenseMat *M);

/* -----------------------------------------------------------------------
 * Gaussian elimination (RREF) over Q.
 *
 * Reduces M in-place to reduced row echelon form.
 * pivot_col[i] = column index of the pivot in row i (or -1 if zero row).
 * Returns the rank.
 *
 * If apply_to is non-NULL, the same row operations are applied to it
 * (used to track the transformation matrix).
 * ----------------------------------------------------------------------- */
int dmat_rref(DenseMat *M, int *pivot_col, DenseMat *apply_to);

/* -----------------------------------------------------------------------
 * Kernel basis.
 *
 * Given M (nrows × ncols), compute a basis for ker(M) ⊆ Q^ncols.
 * Returns a DenseMat whose ROWS are the kernel basis vectors.
 * Returns NULL on allocation failure.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_kernel(const DenseMat *M);

/* -----------------------------------------------------------------------
 * Image (column space) basis.
 * Returns a DenseMat whose COLUMNS span the image of M.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_image(const DenseMat *M);

/* -----------------------------------------------------------------------
 * Matrix multiply: A (m×k) * B (k×n) → C (m×n).
 * Returns NULL on overflow or allocation failure.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_mul(const DenseMat *A, const DenseMat *B);

/* -----------------------------------------------------------------------
 * Stack two matrices vertically: [A; B] where A and B have same ncols.
 * ----------------------------------------------------------------------- */
DenseMat *dmat_vstack(const DenseMat *A, const DenseMat *B);

#endif /* MODULAR_LINALG_H */
