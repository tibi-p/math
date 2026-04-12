#include "rat.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* -----------------------------------------------------------------------
 * Lifecycle
 * ----------------------------------------------------------------------- */
void rat_init(Rat *r)  { mpq_init(r->q); }
void rat_clear(Rat *r) { mpq_clear(r->q); }

/* -----------------------------------------------------------------------
 * GCD — delegates to GMP.
 * ----------------------------------------------------------------------- */
int64_t i64_gcd(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = b; b = a % b; a = t; }
    return a;
}

/* -----------------------------------------------------------------------
 * Constructors
 * ----------------------------------------------------------------------- */

/* Helper: build a Rat by value (for functions that return Rat).
 * The caller must have already called rat_init on the destination,
 * OR this is a by-value return where the caller owns the result.
 * For GMP, returning mpq_t by value is impossible (it's an array type),
 * so we return a struct wrapping it — the struct copy copies the pointer,
 * which is only safe as a temporary before the caller writes it into an
 * already-initialised slot via rat_add/rat_mul etc.
 *
 * Convention: rat_zero/rat_one/rat_int/rat_make return a NEWLY initialised
 * Rat that the caller must eventually rat_clear.  This is the one case
 * where we allocate inside a constructor — acceptable because these are
 * only called during matrix init (via dmat_new → rat_init + assignment). */

Rat rat_zero(void)
{
    Rat r; mpq_init(r.q); mpq_set_si(r.q, 0, 1); return r;
}

Rat rat_one(void)
{
    Rat r; mpq_init(r.q); mpq_set_si(r.q, 1, 1); return r;
}

Rat rat_int(int64_t n)
{
    Rat r; mpq_init(r.q); mpq_set_si(r.q, (long)n, 1); return r;
}

Rat rat_make(int64_t num, int64_t den)
{
    Rat r; mpq_init(r.q);
    if (den == 0) { mpq_set_si(r.q, 0, 1); return r; }
    mpq_set_si(r.q, (long)num, (unsigned long)(den < 0 ? -den : den));
    if (den < 0) mpq_neg(r.q, r.q);
    mpq_canonicalize(r.q);
    return r;
}

int rat_reduce(Rat *r) { mpq_canonicalize(r->q); return RAT_OK; }

/* -----------------------------------------------------------------------
 * Predicates
 * ----------------------------------------------------------------------- */
int rat_is_zero(Rat a) { return mpq_sgn(a.q) == 0; }

int rat_eq(Rat a, Rat b) { return mpq_equal(a.q, b.q); }

/* -----------------------------------------------------------------------
 * Arithmetic
 *
 * All functions write into an already-initialised *out.
 * ----------------------------------------------------------------------- */
int rat_neg(Rat a, Rat *out)  { mpq_neg(out->q, a.q);       return RAT_OK; }
int rat_add(Rat a, Rat b, Rat *out) { mpq_add(out->q, a.q, b.q); return RAT_OK; }
int rat_sub(Rat a, Rat b, Rat *out) { mpq_sub(out->q, a.q, b.q); return RAT_OK; }
int rat_mul(Rat a, Rat b, Rat *out) { mpq_mul(out->q, a.q, b.q); return RAT_OK; }

int rat_inv(Rat a, Rat *out)
{
    if (mpq_sgn(a.q) == 0) return RAT_DIVZERO;
    mpq_inv(out->q, a.q);
    return RAT_OK;
}

int rat_div(Rat a, Rat b, Rat *out)
{
    if (mpq_sgn(b.q) == 0) return RAT_DIVZERO;
    mpq_div(out->q, a.q, b.q);
    return RAT_OK;
}

int rat_scale(Rat *r, int64_t n)
{
    /* Multiply r->q by integer n in-place. */
    mpq_t tmp; mpq_init(tmp);
    mpq_set_si(tmp, (long)n, 1);
    mpq_mul(r->q, r->q, tmp);
    mpq_clear(tmp);
    return RAT_OK;
}

/* -----------------------------------------------------------------------
 * In-place setters (dst must already be mpq_init'd)
 * ----------------------------------------------------------------------- */
void rat_set(Rat *dst, Rat src)    { mpq_set(dst->q, src.q); }
void rat_set_zero(Rat *r)          { mpq_set_si(r->q, 0, 1); }
void rat_set_one(Rat *r)           { mpq_set_si(r->q, 1, 1); }
void rat_set_si(Rat *r, int64_t n) { mpq_set_si(r->q, (long)n, 1); }

/* -----------------------------------------------------------------------
 * Debug
 * ----------------------------------------------------------------------- */
void rat_print(Rat r)
{
    gmp_printf("%Qd", r.q);
}
