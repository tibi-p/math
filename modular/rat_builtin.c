#include "rat.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

void rat_init(Rat *r)  { (void)r; }
void rat_clear(Rat *r) { (void)r; }

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

/* Absolute value of int64, safe for INT64_MIN via unsigned cast. */
static inline uint64_t u64_abs(int64_t x)
{
    return (x < 0) ? (uint64_t)(-(x + 1)) + 1u : (uint64_t)x;
}

/* Checked multiplication: returns 0 and sets *out = a*b on success,
 * returns RAT_OVERFLOW on overflow. */
static int checked_mul(int64_t a, int64_t b, int64_t *out)
{
    if (a == 0 || b == 0) { *out = 0; return RAT_OK; }
    uint64_t ua = u64_abs(a), ub = u64_abs(b);
    /* ua * ub overflows int64 iff ua > INT64_MAX / ub */
    if (ua > (uint64_t)INT64_MAX / ub) return RAT_OVERFLOW;
    int64_t product = a * b;
    /* Sign check as a sanity guard. */
    int neg = (a < 0) ^ (b < 0);
    if (neg && product > 0) return RAT_OVERFLOW;
    if (!neg && product < 0) return RAT_OVERFLOW;
    *out = product;
    return RAT_OK;
}

/* Checked addition. */
static int checked_add(int64_t a, int64_t b, int64_t *out)
{
    if (b > 0 && a > INT64_MAX - b) return RAT_OVERFLOW;
    if (b < 0 && a < INT64_MIN - b) return RAT_OVERFLOW;
    *out = a + b;
    return RAT_OK;
}

/* -----------------------------------------------------------------------
 * Public API
 * ----------------------------------------------------------------------- */

int64_t i64_gcd(int64_t a, int64_t b)
{
    uint64_t ua = u64_abs(a), ub = u64_abs(b);
    while (ub) {
        uint64_t t = ua % ub;
        ua = ub;
        ub = t;
    }
    return (int64_t)ua;
}

int rat_reduce(Rat *r)
{
    if (r->den == 0) return RAT_DIVZERO;
    if (r->num == 0) { r->den = 1; return RAT_OK; }
    int64_t g = i64_gcd(r->num, r->den);
    r->num /= g;
    r->den /= g;
    /* Ensure den > 0. */
    if (r->den < 0) { r->num = -r->num; r->den = -r->den; }
    return RAT_OK;
}

Rat rat_zero(void)               { return (Rat){0, 1}; }
Rat rat_one(void)                { return (Rat){1, 1}; }
Rat rat_int(int64_t n)           { return (Rat){n, 1}; }

Rat rat_make(int64_t num, int64_t den)
{
    Rat r = {num, den};
    if (rat_reduce(&r) != RAT_OK) {
        /* Caller passed den==0; return a sentinel that will propagate. */
        r.num = 0; r.den = 0;
    }
    return r;
}

int rat_is_zero(Rat a) { return a.num == 0; }

int rat_eq(Rat a, Rat b)
{
    /* Both are in canonical form, so component-wise equality suffices. */
    return a.num == b.num && a.den == b.den;
}

int rat_neg(Rat a, Rat *out)
{
    /* -INT64_MIN is undefined; but num cannot be INT64_MIN in canonical
     * form because den > 0 and gcd=1 implies num is unrestricted in
     * magnitude, so we must guard. */
    if (a.num == INT64_MIN) return RAT_OVERFLOW;
    *out = (Rat){-a.num, a.den};
    return RAT_OK;
}

int rat_inv(Rat a, Rat *out)
{
    if (a.num == 0) return RAT_DIVZERO;
    Rat r;
    if (a.num > 0) {
        r.num = a.den;
        r.den = a.num;
    } else {
        /* a.num < 0: flip sign to keep den positive. */
        if (a.num == INT64_MIN) return RAT_OVERFLOW;
        r.num = -a.den;
        r.den = -a.num;
    }
    *out = r;
    return RAT_OK;
}

/* a/b + c/d = (a*d + c*b) / (b*d), then reduce. */
int rat_add(Rat a, Rat b, Rat *out)
{
    /* Reduce common factors between cross-denominators first (Knuth trick). */
    int64_t g = i64_gcd(a.den, b.den);
    int64_t bden_g = b.den / g;   /* b.den / gcd */
    int64_t aden_g = a.den / g;   /* a.den / gcd */

    int64_t num_part1, num_part2, num, den;
    if (checked_mul(a.num, bden_g, &num_part1) != RAT_OK) return RAT_OVERFLOW;
    if (checked_mul(b.num, aden_g, &num_part2) != RAT_OK) return RAT_OVERFLOW;
    if (checked_add(num_part1, num_part2, &num) != RAT_OK) return RAT_OVERFLOW;
    if (checked_mul(a.den, bden_g, &den) != RAT_OK) return RAT_OVERFLOW;

    Rat r = {num, den};
    if (rat_reduce(&r) != RAT_OK) return RAT_OVERFLOW;
    *out = r;
    return RAT_OK;
}

int rat_sub(Rat a, Rat b, Rat *out)
{
    Rat nb;
    if (rat_neg(b, &nb) != RAT_OK) return RAT_OVERFLOW;
    return rat_add(a, nb, out);
}

/* (a/b) * (c/d) = (a*c) / (b*d); cross-reduce first. */
int rat_mul(Rat a, Rat b, Rat *out)
{
    /* Cross-cancel gcd(a.num, b.den) and gcd(b.num, a.den). */
    int64_t g1 = i64_gcd(u64_abs(a.num), (uint64_t)b.den);
    int64_t g2 = i64_gcd(u64_abs(b.num), (uint64_t)a.den);

    int64_t an = a.num / g1;
    int64_t bd = b.den / g1;
    int64_t bn = b.num / g2;
    int64_t ad = a.den / g2;

    int64_t num, den;
    if (checked_mul(an, bn, &num) != RAT_OK) return RAT_OVERFLOW;
    if (checked_mul(ad, bd, &den) != RAT_OK) return RAT_OVERFLOW;

    Rat r = {num, den};
    if (rat_reduce(&r) != RAT_OK) return RAT_OVERFLOW;
    *out = r;
    return RAT_OK;
}

int rat_div(Rat a, Rat b, Rat *out)
{
    Rat ib;
    int rc = rat_inv(b, &ib);
    if (rc != RAT_OK) return rc;
    return rat_mul(a, ib, out);
}

int rat_scale(Rat *r, int64_t n)
{
    Rat s = rat_int(n);
    Rat result;
    int rc = rat_mul(*r, s, &result);
    if (rc != RAT_OK) return rc;
    *r = result;
    return RAT_OK;
}

void rat_set(Rat *dst, Rat src)    { *dst = src; }
void rat_set_zero(Rat *r)          { r->num = 0; r->den = 1; }
void rat_set_one(Rat *r)           { r->num = 1; r->den = 1; }
void rat_set_si(Rat *r, int64_t n) { r->num = n; r->den = 1; }

void rat_print(Rat r)
{
    if (r.den == 1)
        printf("%" PRId64, r.num);
    else
        printf("%" PRId64 "/%" PRId64, r.num, r.den);
}
