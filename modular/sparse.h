#ifndef MODULAR_SPARSE_H
#define MODULAR_SPARSE_H

#include "rat.h"

/* -----------------------------------------------------------------------
 * SEntry — one nonzero entry (col, value) in a sparse row.
 * ----------------------------------------------------------------------- */
typedef struct {
    int col;
    Rat val;
} SEntry;

/* -----------------------------------------------------------------------
 * SparseMatBuilder — dynamic row-of-lists format used during construction.
 * ----------------------------------------------------------------------- */
typedef struct {
    int      nrows;
    int      ncols;
    SEntry **rows;   /* rows[i] = malloc'd array of nonzero entries for row i */
    int     *lens;   /* lens[i] = number of entries in rows[i]                */
    int     *caps;   /* caps[i] = allocated capacity of rows[i]               */
} SparseMatBuilder;

SparseMatBuilder *smat_builder_new(int nrows, int ncols);
void              smat_builder_free(SparseMatBuilder *B);

/* Add val to entry (row, col); creates entry if not present.
 * Returns 0 on success, RAT_OVERFLOW on overflow. */
int smat_builder_add(SparseMatBuilder *B, int row, int col, Rat val);

/* -----------------------------------------------------------------------
 * SparseMat — frozen CSR (Compressed Sparse Row) format for elimination.
 *
 * After freezing, rows are sorted by column index.
 * ----------------------------------------------------------------------- */
typedef struct {
    int    nrows;
    int    ncols;
    int   *row_ptr;  /* row_ptr[i]..row_ptr[i+1]-1 = range in col_ind/vals */
    int   *col_ind;  /* column indices                                       */
    Rat   *vals;     /* values, parallel to col_ind                          */
} SparseMat;

/* Convert builder to frozen CSR; caller still owns B (must free separately). */
SparseMat *smat_freeze(const SparseMatBuilder *B);
void       smat_free(SparseMat *M);

/* Build a SparseMat from a dense integer matrix (for small N). */
SparseMat *smat_from_dense_int(const int *A, int nrows, int ncols);

/* Deep copy. */
SparseMat *smat_copy(const SparseMat *M);

/* Return value at (row, col); rat_zero() if not present. */
Rat smat_get(const SparseMat *M, int row, int col);

/* Debug: print to stdout. */
void smat_print(const SparseMat *M);

#endif /* MODULAR_SPARSE_H */
