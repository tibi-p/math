#include <stdio.h>
#include <stdlib.h>

#include "hecke.h"

/* Simple primality test for small p. */
static int is_prime(int n)
{
    if (n < 2) return 0;
    for (int i = 2; (long long)i * i <= n; i++)
        if (n % i == 0) return 0;
    return 1;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s N [p_max]\n", argv[0]);
        fprintf(stderr, "  Compute Hecke eigenvalues for S_2(Gamma_0(N)).\n");
        return 1;
    }

    int N = atoi(argv[1]);
    int p_max = (argc >= 3) ? atoi(argv[2]) : 20;

    if (N < 1 || N > 300) {
        fprintf(stderr, "N must be in [1, 300] for the int64 prototype.\n");
        return 1;
    }

    printf("Building S_2(Gamma_0(%d))...\n", N);
    CuspidalSpace *cs = cspace_build(N);
    if (!cs) {
        fprintf(stderr, "Failed to build cuspidal space.\n");
        return 1;
    }

    cspace_print(cs);

    if (cs->dim == 0) {
        printf("dim S_2(Gamma_0(%d)) = 0. No eigenforms.\n", N);
        cspace_free(cs);
        return 0;
    }

    printf("\nHecke eigenvalue matrices T_p for p <= %d:\n\n", p_max);
    for (int p = 2; p <= p_max; p++) {
        if (!is_prime(p)) continue;

        DenseMat *Tp = cspace_hecke_matrix(cs, p);
        if (!Tp) {
            fprintf(stderr, "Failed to compute T_%d.\n", p);
            continue;
        }

        printf("T_%d =\n", p);
        dmat_print(Tp);
        printf("\n");
        dmat_free(Tp);
    }

    cspace_free(cs);
    return 0;
}
