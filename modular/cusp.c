#include "cusp.h"

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

/* -----------------------------------------------------------------------
 * Constructors
 * ----------------------------------------------------------------------- */

Cusp cusp_inf(void)
{
    return (Cusp){ CUSP_INF, 1, 0 };
}

Cusp cusp_from_frac(int64_t p, int64_t q)
{
    if (q == 0) {
        /* Only (1:0) is valid for infinity; anything else is an error. */
        return cusp_inf();
    }
    /* Normalise: ensure q > 0 and gcd(|p|,q) == 1. */
    if (q < 0) { p = -p; q = -q; }
    int64_t g = i64_gcd(p < 0 ? -p : p, q);
    if (g > 1) { p /= g; q /= g; }
    return (Cusp){ CUSP_FINITE, p, q };
}

Cusp cusp_from_rat(Rat r)
{
    return cusp_from_frac(r.num, r.den);
}

/* -----------------------------------------------------------------------
 * Predicates
 * ----------------------------------------------------------------------- */

int cusp_is_inf(Cusp c)
{
    return c.kind == CUSP_INF;
}

int cusp_eq(Cusp a, Cusp b)
{
    if (a.kind != b.kind) return 0;
    if (a.kind == CUSP_INF) return 1;
    return a.p == b.p && a.q == b.q;
}

/* -----------------------------------------------------------------------
 * Gamma_0(N) equivalence of cusps
 *
 * Two cusps a/c and b/d are Gamma_0(N)-equivalent iff there exists
 * [[alpha, beta], [gamma, delta]] in Gamma_0(N) (N | gamma) mapping one
 * to the other.  The standard criterion (Shimura / Diamond–Shurman):
 *
 *   a/c ~_{Γ₀(N)} b/d
 *   iff  gcd(c, N) == gcd(d, N)
 *        AND a ≡ b (mod gcd(c,N) · gcd(d,N) / gcd(c,d,N)^2 · ...)
 *
 * For a quick exact test we use the explicit characterisation from
 * Cremona §2.2: a/c ~ b/d mod Γ₀(N) iff gcd(c,N) | d·gcd(c,N)/(something).
 *
 * Practical implementation: represent each cusp by its Γ₀(N)-class
 * identifier (c₁, g₁) where c₁ = c / gcd(c, N) and g₁ = gcd(c, N),
 * then two cusps are equivalent iff they share the same class.
 *
 * Theorem (Cremona 2.1): The cusps of Γ₀(N) are in bijection with pairs
 * (e, f) where ef | N, gcd(e, f) = 1, e > 0, f > 0, and f | gcd(c, N)
 * with e = N / gcd(c, N).  Two rationals a/c and b/d are equivalent iff
 * gcd(c, N) == gcd(d, N)  AND  (a/gcd(a,c)) ≡ (b/gcd(b,d)) mod gcd(c,N).
 *
 * We implement the simpler sufficient-and-necessary condition:
 *   Let g_c = gcd(c, N),  g_d = gcd(d, N).
 *   Equivalent iff g_c == g_d AND a·(c/g_c)^{-1} ≡ b·(d/g_d)^{-1} (mod g_c).
 * ----------------------------------------------------------------------- */
int cusp_gamma0_eq(Cusp a, Cusp b, int N)
{
    if (cusp_is_inf(a) && cusp_is_inf(b)) return 1;
    if (cusp_is_inf(a) || cusp_is_inf(b)) {
        /* ∞ = 1/0; treat denominator as 0, gcd(0, N) = N. */
        int64_t g_inf = N;
        int64_t g_other = cusp_is_inf(a) ? i64_gcd(b.q, (int64_t)N)
                                          : i64_gcd(a.q, (int64_t)N);
        return g_inf == g_other; /* both must have gcd(den,N) = N */
    }

    int64_t g_a = i64_gcd(a.q, (int64_t)N);
    int64_t g_b = i64_gcd(b.q, (int64_t)N);
    if (g_a != g_b) return 0;

    /* Check: a.p * (a.q / g_a) ≡ b.p * (b.q / g_b) (mod g_a).
     * Both sides are computed mod g_a. */
    int64_t g = g_a;
    int64_t lhs = ((a.p % g) * ((a.q / g) % g)) % g;
    int64_t rhs = ((b.p % g) * ((b.q / g) % g)) % g;
    /* Normalise to [0, g). */
    lhs = ((lhs % g) + g) % g;
    rhs = ((rhs % g) + g) % g;
    return lhs == rhs;
}

/* -----------------------------------------------------------------------
 * Debug
 * ----------------------------------------------------------------------- */

void cusp_print(Cusp c)
{
    if (c.kind == CUSP_INF)
        printf("inf");
    else
        printf("%" PRId64 "/%" PRId64, c.p, c.q);
}
