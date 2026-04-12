#ifndef MODULAR_MANIN_H
#define MODULAR_MANIN_H

#include "rat.h"
#include "p1.h"

/* -----------------------------------------------------------------------
 * ManinSym — a single generator of the free module Z[P¹(Z/NZ)].
 * At the code level this is just an index into a P1Table.
 * ----------------------------------------------------------------------- */
typedef struct {
    int idx;   /* index in [0, mu) */
} ManinSym;

/* -----------------------------------------------------------------------
 * ManinElt — a formal Q-linear combination of Manin symbols.
 *
 * Kept sorted by index; no index appears twice.
 * ----------------------------------------------------------------------- */
typedef struct {
    int  len;
    int  cap;
    int *indices;   /* parallel arrays, sorted ascending by index */
    Rat *coeffs;
} ManinElt;

ManinElt *melt_new(void);
void      melt_free(ManinElt *e);
ManinElt *melt_copy(const ManinElt *e);

/* Add c * x_idx to e.  Returns 0 on success, -1 on overflow. */
int  melt_add_term(ManinElt *e, int idx, Rat c);

/* Zero out any terms whose coefficient is zero (after cancellation). */
void melt_compact(ManinElt *e);

/* Scale all coefficients by r in-place.  Returns 0 on success. */
int  melt_scale(ManinElt *e, Rat r);

/* -----------------------------------------------------------------------
 * Relation matrices for M₂(Γ₀(N)) = Z[P¹(Z/NZ)] / (S-rel, T-rel)
 *
 * build_relation_rows writes one row per S-relation and one row per
 * T-relation orbit into the caller-supplied (nrows × mu) integer matrix
 * rel[row][col] (dense, for small N).  Returns the number of rows written.
 *
 * The caller is responsible for allocating rel as int rel[2*mu][mu].
 * ----------------------------------------------------------------------- */
int manin_build_relations(const P1Table *t,
                           int  *rel,       /* flat array [nrows_max * mu] */
                           int   nrows_max,
                           int  *nrows_out);

/* -----------------------------------------------------------------------
 * Boundary map δ : Z[P¹(Z/NZ)] → Z[cusps of Γ₀(N)]
 *
 * δ(x_i) = [cusp_top(i)] - [cusp_bot(i)]  in Z[P¹(Q)/Γ₀(N)].
 *
 * Cusps of Γ₀(N) are enumerated and stored in cusp_list[0..ncusps-1].
 * The boundary matrix bdy[i * ncusps + j] = coefficient of cusp j in δ(x_i).
 * Returns ncusps.
 * ----------------------------------------------------------------------- */
int manin_build_boundary(const P1Table *t,
                          Cusp **cusp_list_out,  /* malloc'd, caller frees */
                          int   *bdy,            /* flat [mu × ncusps]     */
                          int    bdy_ncols_max);

#endif /* MODULAR_MANIN_H */
