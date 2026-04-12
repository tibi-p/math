#include "manin.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -----------------------------------------------------------------------
 * ManinElt
 * ----------------------------------------------------------------------- */

#define MELT_INIT_CAP 8

ManinElt *melt_new(void)
{
    ManinElt *e = malloc(sizeof(ManinElt));
    if (!e) return NULL;
    e->len = 0;
    e->cap = MELT_INIT_CAP;
    e->indices = malloc(MELT_INIT_CAP * sizeof(int));
    e->coeffs  = malloc(MELT_INIT_CAP * sizeof(Rat));
    if (!e->indices || !e->coeffs) {
        free(e->indices); free(e->coeffs); free(e);
        return NULL;
    }
    return e;
}

void melt_free(ManinElt *e)
{
    if (!e) return;
    for (int i = 0; i < e->len; i++)
        rat_clear(&e->coeffs[i]);
    free(e->indices);
    free(e->coeffs);
    free(e);
}

ManinElt *melt_copy(const ManinElt *e)
{
    ManinElt *c = malloc(sizeof(ManinElt));
    if (!c) return NULL;
    c->len = e->len;
    c->cap = e->cap;
    c->indices = malloc(e->cap * sizeof(int));
    c->coeffs  = malloc(e->cap * sizeof(Rat));
    if (!c->indices || !c->coeffs) {
        free(c->indices); free(c->coeffs); free(c);
        return NULL;
    }
    memcpy(c->indices, e->indices, (size_t)e->len * sizeof(int));
    for (int i = 0; i < e->len; i++) {
        c->coeffs[i] = rat_zero();
        rat_set(&c->coeffs[i], e->coeffs[i]);
    }
    return c;
}

/* Grow capacity to at least new_cap. */
static int melt_grow(ManinElt *e, int new_cap)
{
    if (new_cap <= e->cap) return 0;
    int *ni = realloc(e->indices, (size_t)new_cap * sizeof(int));
    Rat *nc = realloc(e->coeffs,  (size_t)new_cap * sizeof(Rat));
    if (!ni || !nc) return -1;
    e->indices = ni;
    e->coeffs  = nc;
    e->cap = new_cap;
    return 0;
}

int melt_add_term(ManinElt *e, int idx, Rat c)
{
    if (rat_is_zero(c)) return 0;

    /* Binary search for idx. */
    int lo = 0, hi = e->len;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (e->indices[mid] < idx) lo = mid + 1;
        else                        hi = mid;
    }
    /* lo == position where idx should be. */
    if (lo < e->len && e->indices[lo] == idx) {
        /* Merge: add coefficients. */
        Rat sum = rat_zero();
        if (rat_add(e->coeffs[lo], c, &sum) != RAT_OK) { rat_clear(&sum); return -1; }
        rat_set(&e->coeffs[lo], sum);
        rat_clear(&sum);
        return 0;
    }
    /* Insert new term at position lo. */
    if (e->len == e->cap) {
        if (melt_grow(e, e->cap * 2) != 0) return -1;
    }
    memmove(e->indices + lo + 1, e->indices + lo,
            (size_t)(e->len - lo) * sizeof(int));
    memmove(e->coeffs  + lo + 1, e->coeffs  + lo,
            (size_t)(e->len - lo) * sizeof(Rat));
    /* lo+1..len are valid (memmoved from lo..len-1); init lo fresh. */
    e->indices[lo]  = idx;
    e->coeffs[lo]   = rat_zero();
    rat_set(&e->coeffs[lo], c);
    e->len++;
    return 0;
}

void melt_compact(ManinElt *e)
{
    int w = 0;
    for (int r = 0; r < e->len; r++) {
        if (!rat_is_zero(e->coeffs[r])) {
            if (w < r) {
                e->indices[w] = e->indices[r];
                rat_set(&e->coeffs[w], e->coeffs[r]);
            }
            w++;
        }
    }
    /* Clear positions w..len-1 (zeros or values that were moved). */
    for (int r = w; r < e->len; r++)
        rat_clear(&e->coeffs[r]);
    e->len = w;
}

int melt_scale(ManinElt *e, Rat r)
{
    Rat prod = rat_zero();
    for (int i = 0; i < e->len; i++) {
        if (rat_mul(e->coeffs[i], r, &prod) != RAT_OK) { rat_clear(&prod); return -1; }
        rat_set(&e->coeffs[i], prod);
    }
    rat_clear(&prod);
    melt_compact(e);
    return 0;
}

/* -----------------------------------------------------------------------
 * Relations
 *
 * We produce a set of integer row vectors in Z^mu, each encoding a
 * relation among Manin symbols.
 *
 * S-relation for symbol i:  e_i + e_{S(i)} = 0
 *   → row: +1 at col i, +1 at col S(i).  (If i == S(i), row is 2*e_i.)
 *
 * U-relation for orbit {i, U(i), U²(i)}: e_i + e_{U(i)} + e_{U²(i)} = 0
 *   → row: +1 at col i, +1 at U(i), +1 at U²(i).
 *   We emit one row per orbit (avoid duplicates with a visited[] array).
 * ----------------------------------------------------------------------- */

int manin_build_relations(const P1Table *t,
                           int  *rel,
                           int   nrows_max,
                           int  *nrows_out)
{
    int mu = t->mu;
    int nrows = 0;
    memset(rel, 0, (size_t)nrows_max * mu * sizeof(int));

    /* --- S-relations --- */
    /* Track which pairs we've already emitted (avoid i and S(i) both emitting). */
    char *s_done = calloc((size_t)mu, 1);
    if (!s_done) return -1;

    for (int i = 0; i < mu; i++) {
        if (s_done[i]) continue;
        int si = p1_apply_S(t, i);
        if (si < 0) { free(s_done); return -1; }

        if (nrows >= nrows_max) { free(s_done); return -1; }
        int *row = rel + nrows * mu;
        row[i]  += 1;
        row[si] += 1;
        nrows++;

        s_done[i]  = 1;
        s_done[si] = 1;
    }
    free(s_done);

    /* --- U-relations (3-term) --- */
    char *u_done = calloc((size_t)mu, 1);
    if (!u_done) return -1;

    for (int i = 0; i < mu; i++) {
        if (u_done[i]) continue;
        int ui  = p1_apply_U(t, i);
        int uui = p1_apply_U(t, ui);
        if (ui < 0 || uui < 0) { free(u_done); return -1; }

        if (nrows >= nrows_max) { free(u_done); return -1; }
        int *row = rel + nrows * mu;
        row[i]   += 1;
        row[ui]  += 1;
        row[uui] += 1;
        nrows++;

        u_done[i]   = 1;
        u_done[ui]  = 1;
        u_done[uui] = 1;
    }
    free(u_done);

    *nrows_out = nrows;
    return 0;
}

/* -----------------------------------------------------------------------
 * Boundary map
 *
 * For each Manin symbol i, lift to [[a,b],[c,d]] ∈ SL₂(Z).
 * δ(x_i) = [M^{-1}(∞)] - [M^{-1}(0)] = [(-d/c)] - [(-b/a)]  (as cusps).
 *
 * We collect the distinct Γ₀(N)-cusp classes encountered, assign them
 * indices, and fill bdy[i * ncusps + cusp_idx] with ±1.
 * ----------------------------------------------------------------------- */

int manin_build_boundary(const P1Table *t,
                          Cusp **cusp_list_out,
                          int   *bdy,
                          int    bdy_ncols_max)
{
    int mu = t->N;  /* upper bound on cusp count is small */
    (void)mu;

    int max_cusps = bdy_ncols_max;
    Cusp *cusps = malloc((size_t)max_cusps * sizeof(Cusp));
    if (!cusps) return -1;
    int ncusps = 0;

    int N = t->N;
    memset(bdy, 0, (size_t)t->mu * max_cusps * sizeof(int));

    /* Find or insert a Γ₀(N)-cusp class; return its index. */
    /* (Linear scan — cusp counts are tiny, O(N) at most.) */
    #define FIND_OR_INSERT(cusp_val, cusp_idx_out) do { \
        int found = 0; \
        for (int k = 0; k < ncusps; k++) { \
            if (cusp_gamma0_eq(cusps[k], (cusp_val), N)) { \
                (cusp_idx_out) = k; found = 1; break; \
            } \
        } \
        if (!found) { \
            if (ncusps >= max_cusps) { free(cusps); return -1; } \
            cusps[ncusps] = (cusp_val); \
            (cusp_idx_out) = ncusps++; \
        } \
    } while(0)

    for (int i = 0; i < t->mu; i++) {
        Cusp ctop = p1_cusp_top(t, i);
        Cusp cbot = p1_cusp_bot(t, i);

        int itop, ibot;
        FIND_OR_INSERT(ctop, itop);
        FIND_OR_INSERT(cbot, ibot);

        /* δ(x_i) = [ctop] - [cbot] */
        bdy[i * max_cusps + itop] += 1;
        bdy[i * max_cusps + ibot] -= 1;
    }
    #undef FIND_OR_INSERT

    *cusp_list_out = cusps;
    return ncusps;
}
