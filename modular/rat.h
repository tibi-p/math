#ifndef MODULAR_RAT_H
#define MODULAR_RAT_H

#include <stdint.h>

/* -----------------------------------------------------------------------
 * Rat — an element of Q in canonical form.
 *
 * Invariants (enforced by rat_reduce; must hold at all times):
 *   - den > 0
 *   - gcd(|num|, den) == 1
 *   - Zero is represented as {0, 1}; {0, 0} is a hard error.
 * ----------------------------------------------------------------------- */
typedef struct {
    int64_t num;
    int64_t den;
} Rat;

/* Return codes used by all rat_* operations. */
#define RAT_OK        0
#define RAT_OVERFLOW -1
#define RAT_DIVZERO  -2

/* Euclidean GCD; result is always >= 0. gcd(0,0) == 0. */
int64_t i64_gcd(int64_t a, int64_t b);

/* Reduce r in-place to canonical form.
 * Returns RAT_DIVZERO if r->den == 0.  */
int rat_reduce(Rat *r);

/* Constructors. */
Rat rat_zero(void);                 /* 0/1  */
Rat rat_one(void);                  /* 1/1  */
Rat rat_int(int64_t n);             /* n/1  */
Rat rat_make(int64_t num, int64_t den); /* reduces automatically */

/* Predicates. */
int rat_is_zero(Rat a);
int rat_eq(Rat a, Rat b);

/* Arithmetic — return RAT_OK or RAT_OVERFLOW.
 * On overflow *out is left unchanged. */
int rat_add(Rat a, Rat b, Rat *out);
int rat_sub(Rat a, Rat b, Rat *out);
int rat_mul(Rat a, Rat b, Rat *out);
int rat_inv(Rat a, Rat *out);       /* RAT_DIVZERO if a.num == 0 */
int rat_div(Rat a, Rat b, Rat *out);
int rat_neg(Rat a, Rat *out);

/* Scale r by integer n in-place. */
int rat_scale(Rat *r, int64_t n);

/* Debug: print to stdout as "num/den". */
void rat_print(Rat r);

#endif /* MODULAR_RAT_H */
