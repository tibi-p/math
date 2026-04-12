#include "p1.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

/* Extended Euclidean algorithm: returns gcd(a, b) and sets *x, *y s.t.
 * a*x + b*y = gcd(a, b). */
static int64_t ext_gcd(int64_t a, int64_t b, int64_t *x, int64_t *y)
{
    if (b == 0) { *x = 1; *y = 0; return a; }
    int64_t x1, y1;
    int64_t g = ext_gcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

/* Return ((a % m) + m) % m — a proper non-negative remainder. */
static int64_t mod_pos(int64_t a, int64_t m)
{
    return ((a % m) + m) % m;
}

/* Compute the canonical (lex-min) representative of (c, d) in P¹(Z/NZ).
 * Iterates over all λ ∈ (Z/NZ)* and picks the smallest (λc, λd).
 * Returns 1 if (c, d) is primitive (belongs to P¹), 0 otherwise.
 * Sets *cc, *dd to the canonical pair. */
static int p1_canonicalise(int N, int64_t c, int64_t d,
                            int64_t *cc, int64_t *dd)
{
    c = mod_pos(c, N);
    d = mod_pos(d, N);

    /* Primitivity: gcd(gcd(c,d), N) must be 1. */
    int64_t g = i64_gcd(i64_gcd(c, d), (int64_t)N);
    if (g != 1) return 0;

    int64_t best_c = c, best_d = d;
    for (int64_t lam = 1; lam < (int64_t)N; lam++) {
        if (i64_gcd(lam, (int64_t)N) != 1) continue;
        int64_t lc = (lam * c) % N;
        int64_t ld = (lam * d) % N;
        if (lc < best_c || (lc == best_c && ld < best_d)) {
            best_c = lc;
            best_d = ld;
        }
    }
    *cc = best_c;
    *dd = best_d;
    return 1;
}

/* -----------------------------------------------------------------------
 * Build / free
 * ----------------------------------------------------------------------- */

P1Table *p1table_build(int N)
{
    P1Table *t = malloc(sizeof(P1Table));
    if (!t) return NULL;

    t->N  = N;
    t->mu = 0;

    /* Upper bound on mu: N * prod(1 + 1/p) <= 2*N for any N. */
    int max_mu = 3 * N + 4;
    t->elems = malloc((size_t)max_mu * sizeof(P1Elem));

    t->index_table = malloc((size_t)(N * N) * sizeof(int));
    if (!t->elems || !t->index_table) {
        free(t->elems); free(t->index_table); free(t);
        return NULL;
    }
    memset(t->index_table, -1, (size_t)(N * N) * sizeof(int));

    for (int64_t c = 0; c < N; c++) {
        for (int64_t d = 0; d < N; d++) {
            int64_t cc, dd;
            if (!p1_canonicalise(N, c, d, &cc, &dd)) continue;
            /* Only process each canonical representative once. */
            if (cc != c || dd != d) continue;
            if (t->index_table[c * N + d] != -1) continue;

            int idx = t->mu++;
            t->elems[idx].c = c;
            t->elems[idx].d = d;
            t->index_table[c * N + d] = idx;
        }
    }

    return t;
}

void p1table_free(P1Table *t)
{
    if (!t) return;
    free(t->elems);
    free(t->index_table);
    free(t);
}

/* -----------------------------------------------------------------------
 * Index lookup
 * ----------------------------------------------------------------------- */

int p1_index(const P1Table *t, int64_t c, int64_t d)
{
    int64_t cc, dd;
    if (!p1_canonicalise(t->N, c, d, &cc, &dd)) return -1;
    return t->index_table[cc * t->N + dd];
}

/* -----------------------------------------------------------------------
 * Lift to SL₂(Z)
 *
 * Given the canonical representative (c, d) at index idx, find a, b ∈ Z
 * with a*d - b*c = 1 (i.e., [[a,b],[c,d]] ∈ SL₂(Z)).
 *
 * Strategy: if gcd(c, d) == 1 over Z, use ext_gcd directly.
 * Otherwise, replace c with c + k*N for the smallest k ≥ 0 s.t.
 * gcd(c + k*N, d) == 1 (guaranteed to exist since gcd(gcd(c,d),N)==1
 * implies gcd(c,d) is coprime to N, so by CRT we can always find such k).
 * ----------------------------------------------------------------------- */
int p1_lift_to_sl2z(const P1Table *t, int idx,
                     int64_t *a_out, int64_t *b_out,
                     int64_t *c_out, int64_t *d_out)
{
    if (idx < 0 || idx >= t->mu) return -1;

    int64_t c = t->elems[idx].c;
    int64_t d = t->elems[idx].d;
    int N = t->N;

    /* Find c' = c + k*N with gcd(c', d) == 1. */
    int64_t c_prime = c;
    int found = 0;
    for (int k = 0; k < N + 1; k++) {
        if (i64_gcd(c_prime, d) == 1) { found = 1; break; }
        c_prime += N;
    }
    if (!found) return -1;  /* Should never happen for valid P1 elements. */

    /* Solve a*d - b*c' = 1 via ext_gcd on (d, c'): d*x + c'*y = 1. */
    int64_t x, y;
    int64_t g = ext_gcd(d, c_prime, &x, &y);
    if (g != 1) return -1;

    /* a = x, b = -y gives a*d - b*c' = x*d + y*c' = 1. */
    *a_out = x;
    *b_out = -y;
    *c_out = c_prime;
    *d_out = d;
    return 0;
}

/* -----------------------------------------------------------------------
 * Matrix actions of S and U on P¹(Z/NZ)
 *
 * S = [[0,-1],[1,0]]:  (c:d) -> (-d:c) mod N
 * U = [[0,-1],[1,-1]]: (c:d) -> (-d:c-d) mod N
 * ----------------------------------------------------------------------- */

int p1_apply_S(const P1Table *t, int idx)
{
    if (idx < 0 || idx >= t->mu) return -1;
    int64_t c = t->elems[idx].c;
    int64_t d = t->elems[idx].d;
    /* S*(c:d) = (-d : c) */
    return p1_index(t, -d, c);
}

int p1_apply_U(const P1Table *t, int idx)
{
    if (idx < 0 || idx >= t->mu) return -1;
    int64_t c = t->elems[idx].c;
    int64_t d = t->elems[idx].d;
    /* U*(c:d) = (-d : c-d) */
    return p1_index(t, -d, c - d);
}

/* -----------------------------------------------------------------------
 * Cusp extraction
 *
 * For Manin symbol idx with lift [[a,b],[c,d]] ∈ SL₂(Z):
 *   The corresponding path runs from M^{-1}(0) to M^{-1}(∞).
 *   M^{-1} = [[d,-b],[-c,a]].
 *   M^{-1}(0) = (-b)/d  ... wait: Möbius action on 0 gives -b/d.
 *   Hmm — as a Möbius: M^{-1}(z) = (dz - b)/(-cz + a).
 *   M^{-1}(0) = -b/a   ... no: (d*0 - b)/(-c*0 + a) = -b/a.
 *   M^{-1}(∞) = d/(-c) = -d/c.
 *
 * So the path goes from cusp (-b/a) to cusp (-d/c).
 * We follow the convention: top = M^{-1}(∞) = -d/c,
 *                           bot = M^{-1}(0)  = -b/a.
 * ----------------------------------------------------------------------- */

Cusp p1_cusp_top(const P1Table *t, int idx)
{
    int64_t a, b, c, d;
    if (p1_lift_to_sl2z(t, idx, &a, &b, &c, &d) != 0)
        return cusp_inf();
    if (c == 0) return cusp_inf();
    return cusp_from_frac(-d, c);
}

Cusp p1_cusp_bot(const P1Table *t, int idx)
{
    int64_t a, b, c, d;
    if (p1_lift_to_sl2z(t, idx, &a, &b, &c, &d) != 0)
        return cusp_inf();
    if (a == 0) return cusp_inf();
    return cusp_from_frac(-b, a);
}

/* -----------------------------------------------------------------------
 * Debug
 * ----------------------------------------------------------------------- */

void p1table_print(const P1Table *t)
{
    printf("P1(Z/%dZ): mu = %d\n", t->N, t->mu);
    for (int i = 0; i < t->mu; i++) {
        printf("  [%d] (%" PRId64 " : %" PRId64 ")\n",
               i, t->elems[i].c, t->elems[i].d);
    }
}
