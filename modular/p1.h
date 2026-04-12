#ifndef MODULAR_P1_H
#define MODULAR_P1_H

#include <stdint.h>
#include "rat.h"
#include "cusp.h"

/* -----------------------------------------------------------------------
 * P1Elem — a canonical representative of a class in P¹(Z/NZ).
 *
 * (c : d) with 0 ≤ c, d < N, gcd(gcd(c,d), N) == 1, and (c, d) is the
 * lexicographically smallest pair in its equivalence class under scaling
 * by units (Z/NZ)*.
 * ----------------------------------------------------------------------- */
typedef struct {
    int64_t c;
    int64_t d;
} P1Elem;

/* -----------------------------------------------------------------------
 * P1Table — the complete projective line P¹(Z/NZ).
 *
 * Built once for a given conductor N.  Contains μ = [SL₂(Z):Γ₀(N)] elements.
 * ----------------------------------------------------------------------- */
typedef struct {
    int     N;
    int     mu;           /* total element count                            */
    P1Elem *elems;        /* elems[i] = canonical (c,d) for index i        */
    int    *index_table;  /* N×N flat array; index_table[c*N+d] = index    */
                          /* -1 if (c,d) is not a canonical representative  */
} P1Table;

/* Build / free. */
P1Table *p1table_build(int N);
void     p1table_free(P1Table *t);

/* Normalise (c,d) to its canonical representative in P¹(Z/NZ) and return
 * its index.  Returns -1 if gcd(gcd(c,d), N) != 1 (not in P¹). */
int p1_index(const P1Table *t, int64_t c, int64_t d);

/* Lift the element at index idx to a full matrix [[a,b],[c,d]] ∈ SL₂(Z).
 * The lifted (c_out, d_out) may differ from elems[idx] by a multiple of N.
 * Returns 0 on success, -1 on failure. */
int p1_lift_to_sl2z(const P1Table *t, int idx,
                     int64_t *a_out, int64_t *b_out,
                     int64_t *c_out, int64_t *d_out);

/* Convenience: apply the action of S = [[0,-1],[1,0]] to index idx
 * and return the resulting index in P¹(Z/NZ). */
int p1_apply_S(const P1Table *t, int idx);

/* Apply U = [[0,-1],[1,-1]] (order-3 element of PSL₂(Z)) to idx. */
int p1_apply_U(const P1Table *t, int idx);

/* Return the cusp at infinity associated to the matrix lift of idx.
 * This is the cusp a/c (top-left / bottom-left of the lifted matrix). */
Cusp p1_cusp_top(const P1Table *t, int idx);   /* M^{-1}(∞) = a/c  */
Cusp p1_cusp_bot(const P1Table *t, int idx);   /* M^{-1}(0) = b/d  */

/* Debug. */
void p1table_print(const P1Table *t);

#endif /* MODULAR_P1_H */
