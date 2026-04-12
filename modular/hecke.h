#ifndef MODULAR_HECKE_H
#define MODULAR_HECKE_H

#include "p1.h"
#include "manin.h"
#include "linalg.h"

/* -----------------------------------------------------------------------
 * CuspidalSpace — the full computational context for S₂(Γ₀(N)).
 *
 * Workflow:
 *   1. cspace_build(N)          — build P1 table, relations, boundary map,
 *                                  and the basis for the cuspidal subspace.
 *   2. cspace_hecke_matrix(cs, p) — compute T_p as a matrix on the basis.
 *   3. cspace_eigenvalues(cs, p)  — return sorted list of eigenvalues.
 *   4. cspace_free(cs)           — tear down.
 * ----------------------------------------------------------------------- */
typedef struct {
    int       N;
    P1Table  *p1;

    /* Cuspidal subspace basis: dim × mu matrix over Q.
     * Row i = i-th basis vector expressed in Manin symbol coordinates. */
    DenseMat *basis;
    int       dim;          /* = basis->nrows = dim S_2(Gamma_0(N))        */

    /* Change-of-basis: Manin-symbol coords → cuspidal basis coords.
     * dim × mu matrix (pseudo-inverse of basis for the relevant subspace). */
    DenseMat *coord;
} CuspidalSpace;

/* Build the cuspidal space for S₂(Γ₀(N)).  Returns NULL on failure. */
CuspidalSpace *cspace_build(int N);

/* Free all resources. */
void cspace_free(CuspidalSpace *cs);

/* Compute the matrix of the Hecke operator T_p on the cuspidal basis.
 * Returns a (dim × dim) DenseMat; caller must dmat_free it.
 * p must be prime (not checked).  Works for any p (including p | N). */
DenseMat *cspace_hecke_matrix(const CuspidalSpace *cs, int p);

/* Compute T_p * x_i (Manin symbol i) and return as a ManinElt.
 * Coefficients are all ±1 integers before applying relations. */
ManinElt *hecke_action_on_symbol(const P1Table *t, int sym_idx, int p);

/* Print a summary of the cuspidal space. */
void cspace_print(const CuspidalSpace *cs);

#endif /* MODULAR_HECKE_H */
