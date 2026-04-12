#ifndef MODULAR_CUSP_H
#define MODULAR_CUSP_H

#include <stdint.h>
#include "rat.h"

/* -----------------------------------------------------------------------
 * Cusp — an element of P¹(Q) = Q ∪ {∞}.
 *
 * Invariants:
 *   - CUSP_INF:    p == 1, q == 0.  Represents the cusp at infinity.
 *   - CUSP_FINITE: q > 0, gcd(|p|, q) == 1.  Represents the cusp p/q.
 * ----------------------------------------------------------------------- */
typedef enum { CUSP_FINITE = 0, CUSP_INF = 1 } CuspKind;

typedef struct {
    CuspKind kind;
    int64_t  p;   /* numerator   (defined when CUSP_FINITE) */
    int64_t  q;   /* denominator (defined when CUSP_FINITE, always > 0) */
} Cusp;

/* Constructors. */
Cusp cusp_inf(void);
Cusp cusp_from_frac(int64_t p, int64_t q); /* reduces; aborts if q == 0 */
Cusp cusp_from_rat(Rat r);

/* Predicates. */
int cusp_is_inf(Cusp c);
int cusp_eq(Cusp a, Cusp b);

/* Action of Gamma_0(N) equivalence: returns 1 if a ~_{Gamma_0(N)} b. */
int cusp_gamma0_eq(Cusp a, Cusp b, int N);

/* Debug. */
void cusp_print(Cusp c);

#endif /* MODULAR_CUSP_H */
