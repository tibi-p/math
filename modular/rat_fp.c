#include "rat.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

/* All arithmetic is mod P. */
#define P RAT_FP_PRIME

/* -----------------------------------------------------------------------
 * Lifecycle — no-ops
 * ----------------------------------------------------------------------- */
void rat_init(Rat *r)  { (void)r; }
void rat_clear(Rat *r) { (void)r; }

/* -----------------------------------------------------------------------
 * GCD — still needed by p1.c and cusp.c which use i64_gcd directly.
 * ----------------------------------------------------------------------- */
int64_t i64_gcd(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = b; b = a % b; a = t; }
    return a;
}

/* -----------------------------------------------------------------------
 * Fast modular exponentiation: base^exp mod P.
 * ----------------------------------------------------------------------- */
static int64_t powmod(int64_t base, int64_t exp, int64_t mod)
{
    int64_t result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) result = (__int128)result * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

/* Modular inverse via Fermat's little theorem: a^{P-2} mod P.
 * Valid because P is prime. Returns 0 if a == 0. */
static int64_t modinv(int64_t a)
{
    a = ((a % P) + P) % P;
    if (a == 0) return 0;
    return powmod(a, P - 2, P);
}

/* Normalise any int64 into [0, P). */
static inline int64_t norm(int64_t x)
{
    return ((x % P) + P) % P;
}

/* -----------------------------------------------------------------------
 * Constructors
 * ----------------------------------------------------------------------- */
Rat rat_zero(void)         { return (Rat){0}; }
Rat rat_one(void)          { return (Rat){1}; }
Rat rat_int(int64_t n)     { return (Rat){norm(n)}; }

Rat rat_make(int64_t num, int64_t den)
{
    int64_t d = norm(den);
    if (d == 0) return rat_zero();   /* caller passed 0 denominator */
    int64_t n = norm(num);
    /* n/d = n * d^{-1} mod P */
    return (Rat){ (__int128)n * modinv(d) % P };
}

int rat_reduce(Rat *r) { r->val = norm(r->val); return RAT_OK; }

/* -----------------------------------------------------------------------
 * Predicates
 * ----------------------------------------------------------------------- */
int rat_is_zero(Rat a) { return a.val == 0; }
int rat_eq(Rat a, Rat b) { return a.val == b.val; }

/* -----------------------------------------------------------------------
 * Arithmetic
 * ----------------------------------------------------------------------- */
int rat_neg(Rat a, Rat *out)
{
    out->val = a.val == 0 ? 0 : P - a.val;
    return RAT_OK;
}

int rat_add(Rat a, Rat b, Rat *out)
{
    out->val = (a.val + b.val) % P;
    return RAT_OK;
}

int rat_sub(Rat a, Rat b, Rat *out)
{
    out->val = (a.val + P - b.val) % P;
    return RAT_OK;
}

int rat_mul(Rat a, Rat b, Rat *out)
{
    out->val = (__int128)a.val * b.val % P;
    return RAT_OK;
}

int rat_inv(Rat a, Rat *out)
{
    if (a.val == 0) return RAT_DIVZERO;
    out->val = modinv(a.val);
    return RAT_OK;
}

int rat_div(Rat a, Rat b, Rat *out)
{
    if (b.val == 0) return RAT_DIVZERO;
    out->val = (__int128)a.val * modinv(b.val) % P;
    return RAT_OK;
}

int rat_scale(Rat *r, int64_t n)
{
    r->val = (__int128)r->val * norm(n) % P;
    return RAT_OK;
}

void rat_set(Rat *dst, Rat src)    { *dst = src; }
void rat_set_zero(Rat *r)          { r->val = 0; }
void rat_set_one(Rat *r)           { r->val = 1; }
void rat_set_si(Rat *r, int64_t n) { r->val = norm(n); }

/* -----------------------------------------------------------------------
 * Debug — print as signed integer in (-P/2, P/2] for readability.
 * ----------------------------------------------------------------------- */
void rat_print(Rat r)
{
    int64_t v = r.val;
    if (v > P / 2) v -= P;   /* map to nearest integer */
    printf("%" PRId64, v);
}
