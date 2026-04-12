#include "hecke.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Helper for GCD reduction */
static int64_t manin_gcd(int64_t a, int64_t b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        int64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

/* ============================================================================
 * hecke_action_on_symbol
 * * Computes the right action of the Hecke operator T_p (or U_p) on a single 
 * Manin symbol. 
 *
 * Mathematically, a Manin symbol [c:d] represents a homology path in the 
 * extended upper half-plane between the cusps M^{-1}(0) and M^{-1}(infty), 
 * where M is a matrix in SL_2(Z) with bottom row (c, d). 
 * * The Hecke operator sums over right coset representatives (transversals).
 * This function applies those transversals to the path, re-triangulates the 
 * resulting paths back onto the edges of the fundamental Farey domain, and 
 * accumulates the resulting standard Manin symbols.
 *
 * Parameters:
 * t       - The P^1(Z/NZ) table defining the level N structure.
 * sym_idx - The index of the input Manin symbol (0 to mu-1).
 * p       - The prime defining the Hecke operator.
 *
 * Returns:
 * A newly allocated ManinElt representing the linear combination of 
 * resulting symbols, or NULL on failure.
 * ============================================================================ */
ManinElt *hecke_action_on_symbol(const P1Table *t, int sym_idx, int p)
{
    ManinElt *result = melt_new();
    if (!result) return NULL;

    /* Step 1: Lift the projective line coordinates [c:d] modulo N 
     * up to a true SL_2(Z) matrix [[a, b], [c, d]]. 
     * The bottom row (c, d) completely determines the right coset in Gamma_0(N). */
    int64_t a_lift, b_lift, c, d;
    if (p1_lift_to_sl2z(t, sym_idx, &a_lift, &b_lift, &c, &d) != 0) {
        melt_free(result); return NULL;
    }

    int N = t->N;
    Rat one = rat_one();

    /* ------------------------------------------------------------------
     * Operator U_p (Triggered when p divides the level N)
     * ------------------------------------------------------------------
     * For U_p, the transversals are exactly M_j = [[1, j], [0, p]] for j=0..p-1.
     * Unlike T_p, this action maps paths slightly outside the standard Farey 
     * tessellation. We use Cremona's explicit formulas to map the results 
     * directly back to valid symbols without geometric triangulation.
     */
    if (N % p == 0) {
        if (c % p == 0) {
            /* Case 1: p divides c. 
             * In the geometry of X_0(N), this symbol is connected to the cusp 
             * at infinity. The sum over the p transversals collapses to exactly 
             * one primitive symbol: [c/p : d]. */
            int64_t nc = c / p;
            int64_t nd = d;

            /* C's modulo operator (%) can return negative numbers. 
             * We wrap it to ensure positive indices in [0, N-1]. */
            int mod_c = ((nc % N) + N) % N;
            int mod_d = ((nd % N) + N) % N;
            int idx = p1_index(t, mod_c, mod_d);

            if (idx >= 0) {
                if (melt_add_term(result, idx, one) != 0) {
                    melt_free(result); return NULL;
                }
            }
        } else {
            /* Case 2: p does not divide c. 
             * The sum yields p distinct symbols. Mathematically, applying M_j 
             * gives the coordinates [c*p : d*p - j*c]. */
            for (int j = 0; j < p; j++) {
                int64_t nc = c * (int64_t)p;
                int64_t nd = d * (int64_t)p - c * (int64_t)j;

                /* The resulting coordinates might not be coprime. 
                 * We divide out the greatest common divisor to project them 
                 * back onto the projective line P^1(Q). */
                int64_t g = manin_gcd(nc, nd);
                if (g > 1) {
                    nc /= g;
                    nd /= g;
                }

                int mod_c = ((nc % N) + N) % N;
                int mod_d = ((nd % N) + N) % N;
                int idx = p1_index(t, mod_c, mod_d);

                if (idx >= 0) {
                    if (melt_add_term(result, idx, one) != 0) {
                        melt_free(result); return NULL;
                    }
                }
            }
        }
        melt_compact(result);
        return result;
    }

    /* ------------------------------------------------------------------
     * Operator T_p (Triggered when p does NOT divide N)
     * ------------------------------------------------------------------
     * For T_p, the standard transversals cut across the interior of the 
     * Farey triangles, breaking the Manin symbols. 
     * Instead, we sum over the Heilbronn matrices H_p. Merel proved that 
     * applying the Heilbronn set automatically re-triangulates the broken 
     * paths back onto the edges of the fundamental domain.
     *
     * A Heilbronn matrix is [[A, B], [C, D]] where:
     * 1) Determinant = p (so A*D - B*C = p)
     * 2) A > B >= 0
     * 3) D > C >= 0
     */
    for (int A = 1; A <= p; A++) {
        for (int B = 0; B < A; B++) {
            for (int C = 0; C <= p - A; C++) {
                /* Since AD - BC = p, we can solve for D: D = (p + BC) / A.
                 * We only proceed if (p + BC) is perfectly divisible by A. */
                int rem = p + B * C;
                if (rem % A == 0) {
                    int D = rem / A;

                    /* Enforce the final Heilbronn condition */
                    if (D > C) {
                        /* Right action of [[A, B], [C, D]] on the bottom row (c, d):
                         * [*, *; c, d] * [A, B; C, D] -> bottom row is (c*A+d*C, c*B+d*D) */
                        int64_t nc = c * (int64_t)A + d * (int64_t)C;
                        int64_t nd = c * (int64_t)B + d * (int64_t)D;

                        /* Project back to P^1(Q) by factoring out the GCD */
                        int64_t g = manin_gcd(nc, nd);
                        if (g > 1) {
                            nc /= g;
                            nd /= g;
                        }

                        /* Unlike U_p, Heilbronn matrices map paths such that the 
                         * coordinates are guaranteed to be positive in standard form, 
                         * so a simple modulo N is sufficient here. */
                        int mod_c = nc % N;
                        int mod_d = ((nd % N) + N) % N;
                        int idx = p1_index(t, mod_c, mod_d);

                        /* Accumulate the resulting symbol into our linear combination */
                        if (idx >= 0) {
                            if (melt_add_term(result, idx, one) != 0) {
                                melt_free(result); return NULL;
                            }
                        }
                    }
                }
            }
        }
    }

    /* Remove any terms that cancelled out to 0 and compress the memory array */
    melt_compact(result);
    return result;
}

/* -----------------------------------------------------------------------
 * Express a ManinElt in the cuspidal basis coordinates.
 *
 * cs->coord is dim × mu.  For ManinElt e = ∑ c_i * x_i,
 * the coordinate vector is coord * e  (dim-vector).
 *
 * Returns a newly allocated Rat[cs->dim] array.
 * ----------------------------------------------------------------------- */
static Rat *manin_elt_to_coords(const CuspidalSpace *cs, const ManinElt *e)
{
    int dim = cs->dim;
    int mu  = cs->p1->mu;
    (void)mu;

    Rat *v = malloc((size_t)dim * sizeof(Rat));
    if (!v) return NULL;
    for (int i = 0; i < dim; i++) v[i] = rat_zero();

    /* v = coord * e_sparse */
    for (int k = 0; k < e->len; k++) {
        int col = e->indices[k];
        Rat c   = e->coeffs[k];
        if (rat_is_zero(c)) continue;
        for (int i = 0; i < dim; i++) {
            Rat prod, sum;
            if (rat_mul(dmat_get(cs->coord, i, col), c, &prod) != RAT_OK)
                { free(v); return NULL; }
            if (rat_add(v[i], prod, &sum) != RAT_OK)
                { free(v); return NULL; }
            v[i] = sum;
        }
    }
    return v;
}

/* -----------------------------------------------------------------------
 * cspace_hecke_matrix
 * ----------------------------------------------------------------------- */
DenseMat *cspace_hecke_matrix(const CuspidalSpace *cs, int p)
{
    int dim = cs->dim;
    DenseMat *Tp = dmat_new(dim, dim);
    if (!Tp) return NULL;

    /* For each basis vector b_j (row j of cs->basis, a mu-vector),
     * compute T_p * b_j = ∑_i (b_j)_i * T_p(x_i),
     * then express the result in cuspidal coordinates. */

    for (int j = 0; j < dim; j++) {
        /* Accumulate T_p(b_j) as a ManinElt. */
        ManinElt *accum = melt_new();
        if (!accum) { dmat_free(Tp); return NULL; }

        for (int i = 0; i < cs->p1->mu; i++) {
            Rat bij = dmat_get(cs->basis, j, i);
            if (rat_is_zero(bij)) continue;

            ManinElt *ti = hecke_action_on_symbol(cs->p1, i, p);
            if (!ti) { melt_free(accum); dmat_free(Tp); return NULL; }

            /* Add bij * ti into accum. */
            for (int k = 0; k < ti->len; k++) {
                Rat scaled;
                if (rat_mul(bij, ti->coeffs[k], &scaled) != RAT_OK) {
                    melt_free(ti); melt_free(accum); dmat_free(Tp); return NULL;
                }
                if (melt_add_term(accum, ti->indices[k], scaled) != 0) {
                    melt_free(ti); melt_free(accum); dmat_free(Tp); return NULL;
                }
            }
            melt_free(ti);
        }
        melt_compact(accum);

        /* Express accum in cuspidal coordinates → column j of Tp. */
        Rat *coords = manin_elt_to_coords(cs, accum);
        melt_free(accum);
        if (!coords) { dmat_free(Tp); return NULL; }

        for (int i = 0; i < dim; i++)
            *dmat_at(Tp, i, j) = coords[i];
        free(coords);
    }
    return Tp;
}

/* -----------------------------------------------------------------------
 * cspace_build
 *
 * Steps:
 *  1. Build P1 table.
 *  2. Build relation matrix (S-rel + U-rel) as dense integer matrix.
 *  3. Build boundary map matrix δ.
 *  4. Stack them: combined = [relations; boundary]  (nrows × mu).
 *  5. Kernel of combined = cuspidal subspace basis.
 *  6. Build coord = pseudo-left-inverse of basis on the cuspidal subspace.
 * ----------------------------------------------------------------------- */
CuspidalSpace *cspace_build(int N)
{
    CuspidalSpace *cs = malloc(sizeof(CuspidalSpace));
    if (!cs) return NULL;
    cs->N     = N;
    cs->basis = NULL;
    cs->coord = NULL;

    cs->p1 = p1table_build(N);
    if (!cs->p1) { free(cs); return NULL; }

    int mu = cs->p1->mu;

    int  *rel_dense = NULL;
    int  *bdy_dense = NULL;
    Cusp *cusp_list = NULL;

    /* --- Step 2: relation matrix --- */
    int max_rel_rows = 2 * mu + 4;
    rel_dense = calloc((size_t)max_rel_rows * mu, sizeof(int));
    if (!rel_dense) goto fail;

    int nrel = 0;
    if (manin_build_relations(cs->p1, rel_dense, max_rel_rows, &nrel) != 0)
        goto fail;

    /* --- Step 3: boundary matrix --- */
    int max_cusps = 2 * N + 4;
    bdy_dense = calloc((size_t)mu * max_cusps, sizeof(int));
    if (!bdy_dense) goto fail;

    int ncusps = manin_build_boundary(cs->p1, &cusp_list, bdy_dense, max_cusps);
    if (ncusps < 0) goto fail;

    /* Transpose bdy_dense (mu × ncusps) → rows indexed by Manin symbol.
     * We want rows = one Manin symbol, cols = cusps. It's already in that form. */

    /* --- Step 4: stack [rel_rows; bdy_as_rows] into one DenseMat --- */
    /* rel_dense is (nrel × mu), bdy_dense is (mu × ncusps) — different shapes.
     * We need a combined matrix whose kernel is the cuspidal subspace.
     *
     * The combined system in Z^mu:
     *   - The S- and U-relation rows: each is a vector in Z^mu.
     *   - The boundary map rows: δ(x_i) gives a vector in Z^ncusps.
     *     Transposing: each Manin symbol x_i has a boundary vector.
     *     We want: for each cusp class c, the map ∑_i (bdy[i][c]) * x_i = 0.
     *     i.e., each cusp class gives a row in Z^mu.
     *
     * So bdy contributes ncusps rows in Z^mu (one per cusp class c):
     *   row_c = (bdy[0][c], bdy[1][c], ..., bdy[mu-1][c]).
     */

    int total_rows = nrel + ncusps;
    DenseMat *combined = dmat_new(total_rows, mu);
    if (!combined) goto fail;

    /* Fill relation rows. */
    for (int r = 0; r < nrel; r++)
        for (int c = 0; c < mu; c++)
            *dmat_at(combined, r, c) = rat_int(rel_dense[r * mu + c]);

    /* Fill boundary rows (transpose of bdy_dense). */
    for (int cusp = 0; cusp < ncusps; cusp++) {
        for (int i = 0; i < mu; i++) {
            int v = bdy_dense[i * max_cusps + cusp];
            *dmat_at(combined, nrel + cusp, i) = rat_int(v);
        }
    }

    free(rel_dense);  rel_dense = NULL;
    free(bdy_dense);  bdy_dense = NULL;
    free(cusp_list);  cusp_list = NULL;

    /* --- Step 5: kernel of combined --- */
    DenseMat *ker = dmat_kernel(combined);
    dmat_free(combined);
    if (!ker) goto fail;

    cs->basis = ker;
    cs->dim   = ker->nrows;

    /* --- Step 6: coordinate map ---
     * We need coord: Z^mu → Q^dim s.t. coord * v = (cuspidal coordinates).
     * Since basis rows are orthonormal in the sense that basis * basis^T = I
     * (they are an orthogonal rational basis from RREF), the left-inverse is:
     *   coord = (basis * basis^T)^{-1} * basis.
     * For an RREF basis the pivot structure gives an even simpler form:
     * just take the transpose and reduce.  We use:
     *   coord = basis^T reduced — actually the cleanest is:
     *   coord_ij = (basis^T)[j][i] = basis[i][j] transposed,
     * and to invert basis*basis^T.
     *
     * Simple approach for correctness: form coord as the dim×mu matrix
     * where coord[i][*] = basis[i][*] (i.e., coord = basis itself).
     * Then cs->coord * v for v = ∑ c_k * e_k  gives the component along
     * basis row i only if the basis rows are orthogonal with unit norm.
     *
     * For the general case: compute the Gram matrix G = basis * basis^T
     * (dim×dim), invert it, and set coord = G^{-1} * basis.
     */
    {
        /* basis^T */
        DenseMat *bT = dmat_new(mu, cs->dim);
        if (!bT) goto fail;
        for (int i = 0; i < cs->dim; i++)
            for (int j = 0; j < mu; j++)
                *dmat_at(bT, j, i) = dmat_get(cs->basis, i, j);

        /* G = basis * bT  (dim × dim) */
        DenseMat *G = dmat_mul(cs->basis, bT);
        dmat_free(bT);
        if (!G) goto fail;

        /* Invert G via augmented RREF [G | I]. */
        DenseMat *I = dmat_identity(cs->dim);
        if (!I) { dmat_free(G); goto fail; }
        int *piv = malloc((size_t)cs->dim * sizeof(int));
        if (!piv) { dmat_free(G); dmat_free(I); goto fail; }
        int rank = dmat_rref(G, piv, I);
        free(piv);
        if (rank != cs->dim) {
            /* Degenerate basis — shouldn't happen for a proper kernel. */
            dmat_free(G); dmat_free(I); goto fail;
        }
        dmat_free(G);  /* now reduced to I; we only need I (= G^{-1}) */

        /* coord = G^{-1} * basis */
        DenseMat *coord = dmat_mul(I, cs->basis);
        dmat_free(I);
        if (!coord) goto fail;
        cs->coord = coord;
    }

    return cs;

fail:
    free(rel_dense);
    free(bdy_dense);
    free(cusp_list);
    cspace_free(cs);
    return NULL;
}

void cspace_free(CuspidalSpace *cs)
{
    if (!cs) return;
    p1table_free(cs->p1);
    dmat_free(cs->basis);
    dmat_free(cs->coord);
    free(cs);
}

void cspace_print(const CuspidalSpace *cs)
{
    printf("S_2(Gamma_0(%d)): dim = %d\n", cs->N, cs->dim);
    if (cs->dim == 0) {
        printf("  (trivial space)\n");
        return;
    }
    printf("  Basis (%d vectors in Z^%d):\n", cs->dim, cs->p1->mu);
    dmat_print(cs->basis);
}
