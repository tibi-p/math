#ifndef MODULAR_RAT_H
#define MODULAR_RAT_H

#include <stdint.h>

/* -----------------------------------------------------------------------
 * Backend selection — compile with one of:
 *   -DRAT_BACKEND_BUILTIN   int64_t num/den (default)
 *   -DRAT_BACKEND_GMP       GMP mpq_t
 *   -DRAT_BACKEND_FP        F_p modular arithmetic
 * ----------------------------------------------------------------------- */

#if defined(RAT_BACKEND_GMP)
/* ---- GMP backend ---- */
#  include <gmp.h>
   typedef struct { mpq_t q; } Rat;

#elif defined(RAT_BACKEND_FP)
/* ---- F_p backend ---- */
#  ifndef RAT_FP_PRIME
#    define RAT_FP_PRIME 1000000007LL
#  endif
   typedef struct { int64_t val; } Rat;   /* val in [0, RAT_FP_PRIME) */

#else
/* ---- Builtin backend (default) ---- */
#  if !defined(RAT_BACKEND_BUILTIN)
#    define RAT_BACKEND_BUILTIN
#  endif
   /* Invariants: den > 0, gcd(|num|, den) == 1, zero == {0,1}. */
   typedef struct { int64_t num; int64_t den; } Rat;

#endif  /* backend selection */

/* -----------------------------------------------------------------------
 * Return codes (all backends).
 * ----------------------------------------------------------------------- */
#define RAT_OK        0
#define RAT_OVERFLOW -1
#define RAT_DIVZERO  -2

/* -----------------------------------------------------------------------
 * Lifecycle — MUST be called for every Rat that lives in heap memory.
 * No-ops for builtin and fp; calls mpq_init/mpq_clear for GMP.
 * ----------------------------------------------------------------------- */
void rat_init(Rat *r);
void rat_clear(Rat *r);

/* -----------------------------------------------------------------------
 * GCD (builtin and fp expose this; GMP backend provides a stub).
 * ----------------------------------------------------------------------- */
int64_t i64_gcd(int64_t a, int64_t b);

/* -----------------------------------------------------------------------
 * Constructors.
 * ----------------------------------------------------------------------- */
Rat rat_zero(void);
Rat rat_one(void);
Rat rat_int(int64_t n);
Rat rat_make(int64_t num, int64_t den);

/* Reduce r in-place to canonical form. RAT_DIVZERO if den == 0. */
int rat_reduce(Rat *r);

/* -----------------------------------------------------------------------
 * Predicates.
 * ----------------------------------------------------------------------- */
int rat_is_zero(Rat a);
int rat_eq(Rat a, Rat b);

/* -----------------------------------------------------------------------
 * Arithmetic — return RAT_OK, RAT_OVERFLOW, or RAT_DIVZERO.
 * GMP and Fp backends never return RAT_OVERFLOW.
 * ----------------------------------------------------------------------- */
int rat_add(Rat a, Rat b, Rat *out);
int rat_sub(Rat a, Rat b, Rat *out);
int rat_mul(Rat a, Rat b, Rat *out);
int rat_inv(Rat a, Rat *out);
int rat_div(Rat a, Rat b, Rat *out);
int rat_neg(Rat a, Rat *out);
int rat_scale(Rat *r, int64_t n);

/* -----------------------------------------------------------------------
 * In-place setters — operate on an ALREADY INITIALISED Rat.
 * For GMP these call mpq_set* without allocating a new mpq_t.
 * For builtin/fp they are plain assignments.
 * NEVER call these on an uninitialised Rat.
 * ----------------------------------------------------------------------- */
void rat_set(Rat *dst, Rat src);    /* copy value of src into dst */
void rat_set_zero(Rat *r);          /* set r to 0 */
void rat_set_one(Rat *r);           /* set r to 1 */
void rat_set_si(Rat *r, int64_t n); /* set r to integer n */

/* -----------------------------------------------------------------------
 * Debug.
 * ----------------------------------------------------------------------- */
void rat_print(Rat r);

#endif /* MODULAR_RAT_H */
