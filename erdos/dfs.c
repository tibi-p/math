/*
 * dfs.c — DFS search for an odd distinct covering system
 *
 * We look for a finite set of congruences {n ≡ a_i (mod m_i)} using distinct
 * odd moduli drawn from divisors of an odd abundant number L, such that every
 * integer is covered by at least one congruence.
 *
 * Working in Z/LZ: every residue in {0, ..., L-1} must be hit.
 *
 * Pruning:
 *   (1) Density:  suffix_budget[idx] >= |uncovered|
 *       (If remaining moduli collectively can't cover enough residues, abort.)
 *   (2) Tight:    Σ_{i≥idx} max_a |mask(i,a) ∩ uncov| >= |uncovered|
 *       (Sum of best individual coverages must still reach the target.)
 *
 * Residues are tried in decreasing order of coverage of the current uncovered
 * set, which both finds solutions quickly and drives the tight pruning hard.
 *
 * Usage: ./dfs [L]   (default L=945; try also L=315 for the trivial case)
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* ---- Runtime configuration ---- */
static int L;
static int NW;          /* ceil(L / 64) — number of 64-bit words in a bitset */

static int  moduli[64]; /* divisors of L > 1, sorted ascending = descending coverage */
static int  nmod;

/* ---- Precomputed data ---- */
static uint64_t *mask_pool;  /* flat storage: mask_pool[(i*max_m + a)*NW .. +NW-1] */
static int       max_m;      /* largest modulus + 1, for indexing */
static int       suffix_budget[65]; /* suffix_budget[i] = Σ_{j≥i} L/moduli[j] */

/* ---- Search state ---- */
static int      chosen_a[64];  /* chosen_a[i] = residue for moduli[i], or -1=skip */
static long long node_count;
static time_t   t_start;

/* ---- Bitset helpers (NW-word dynamic width) ---- */

static inline int bs_empty(const uint64_t *b) {
    for (int i = 0; i < NW; i++) if (b[i]) return 0;
    return 1;
}

static inline int bs_popcount(const uint64_t *b) {
    int c = 0;
    for (int i = 0; i < NW; i++) c += __builtin_popcountll(b[i]);
    return c;
}

static inline void bs_andnot(uint64_t *dst, const uint64_t *a, const uint64_t *b) {
    /* dst = a & ~b */
    for (int i = 0; i < NW; i++) dst[i] = a[i] & ~b[i];
}

static inline int bs_and_popcount(const uint64_t *a, const uint64_t *b) {
    int c = 0;
    for (int i = 0; i < NW; i++) c += __builtin_popcountll(a[i] & b[i]);
    return c;
}

/* ---- Mask accessor ---- */
static inline uint64_t *get_mask(int mod_idx, int a) {
    return mask_pool + ((size_t)mod_idx * max_m + a) * NW;
}

/* ---- Initialisation ---- */

static void compute_divisors(void) {
    nmod = 0;
    for (int d = 2; d <= L; d++)
        if (L % d == 0)
            moduli[nmod++] = d;
    /* Result is sorted ascending = coverage descending (L/3 > L/5 > ...) */
}

static void precompute(void) {
    NW    = (L + 63) / 64;
    max_m = moduli[nmod - 1] + 1;

    mask_pool = calloc((size_t)nmod * max_m * NW, sizeof(uint64_t));
    if (!mask_pool) { perror("calloc"); exit(1); }

    for (int i = 0; i < nmod; i++) {
        int m = moduli[i];
        for (int r = 0; r < L; r++) {
            int a = r % m;
            uint64_t *mask = get_mask(i, a);
            mask[r / 64] |= (1ULL << (r % 64));
        }
    }

    suffix_budget[nmod] = 0;
    for (int i = nmod - 1; i >= 0; i--)
        suffix_budget[i] = suffix_budget[i + 1] + L / moduli[i];
}

/* ---- DFS ---- */

/*
 * uncov  : bitset of currently uncovered residues (lives on caller's stack)
 * idx    : index into moduli[] of the next modulus to decide on
 * returns: 1 if a cover was found (chosen_a[] holds the solution), 0 otherwise
 */
static int dfs(uint64_t *uncov, int idx) {
    ++node_count;

    /* Periodic progress report (~every 128M nodes) */
    if ((node_count & ((1LL << 27) - 1)) == 0) {
        printf("  [%.3e nodes | %lds] depth=%d uncov=%d\n",
               (double)node_count, (long)(time(NULL) - t_start),
               idx, bs_popcount(uncov));
        fflush(stdout);
    }

    if (bs_empty(uncov)) return 1;   /* complete cover — success */
    if (idx == nmod)     return 0;   /* ran out of moduli — fail */

    int uc = bs_popcount(uncov);

    /* (1) Density pruning */
    if (suffix_budget[idx] < uc) return 0;

    /* (2) Tighter bound: sum of best single-modulus coverage for each remaining
           modulus.  If this sum < uc, even the most optimistic assignment fails. */
    {
        int cap_sum = 0;
        for (int i = idx; i < nmod; i++) {
            int m = moduli[i], best = 0, cap = L / m;
            for (int a = 0; a < m; a++) {
                int c = bs_and_popcount(uncov, get_mask(i, a));
                if (c > best) { best = c; if (best >= cap) break; }
            }
            cap_sum += best;
            if (cap_sum >= uc) break; /* upper bound already sufficient */
        }
        if (cap_sum < uc) return 0;
    }

    /* Compute coverage of each residue for moduli[idx] */
    int m = moduli[idx];
    int ord[945], cov[945];
    for (int a = 0; a < m; a++) {
        ord[a] = a;
        cov[a] = bs_and_popcount(uncov, get_mask(idx, a));
    }

    /* Sort residues by coverage descending (insertion sort; m ≤ 945, depth ≤ 15) */
    for (int i = 1; i < m; i++) {
        int oi = ord[i], ci = cov[i], j = i - 1;
        while (j >= 0 && cov[j] < ci) { ord[j+1] = ord[j]; cov[j+1] = cov[j]; --j; }
        ord[j+1] = oi; cov[j+1] = ci;
    }

    uint64_t new_uc[NW]; /* VLA — NW set at startup from L */

    /* Branch: try each residue for moduli[idx] */
    for (int i = 0; i < m; i++) {
        int a = ord[i];
        bs_andnot(new_uc, uncov, get_mask(idx, a));
        chosen_a[idx] = a;
        if (dfs(new_uc, idx + 1)) return 1;
    }

    /* Branch: skip moduli[idx] entirely (only if density still feasible) */
    if (suffix_budget[idx + 1] >= uc) {
        chosen_a[idx] = -1;
        if (dfs(uncov, idx + 1)) return 1;
    }

    return 0;
}

/* ---- Entry point ---- */

int main(int argc, char *argv[]) {
    L = (argc > 1) ? atoi(argv[1]) : 945;

    compute_divisors();
    precompute();

    /* Summary */
    printf("=== Odd Covering System Search ===\n");
    printf("L = %d\n", L);
    printf("Moduli  (%d):", nmod);
    for (int i = 0; i < nmod; i++) printf(" %d", moduli[i]);
    printf("\n");
    printf("Coverage   :");
    for (int i = 0; i < nmod; i++) printf(" %d", L / moduli[i]);
    printf("\n");
    printf("Budget: %d  |  L: %d  |  Slack: %d\n\n",
           suffix_budget[0], L, suffix_budget[0] - L);

    if (suffix_budget[0] < L) {
        printf("TRIVIALLY IMPOSSIBLE: total density %.4f < 1.\n",
               (double)suffix_budget[0] / L);
        free(mask_pool);
        return 0;
    }

    /* Initial state: all L residues uncovered */
    uint64_t uncov[NW];
    memset(uncov, 0, NW * sizeof(uint64_t));
    for (int r = 0; r < L; r++)
        uncov[r / 64] |= (1ULL << (r % 64));

    node_count = 0;
    t_start    = time(NULL);
    printf("Starting DFS...\n");
    fflush(stdout);

    int found = dfs(uncov, 0);

    long elapsed = (long)(time(NULL) - t_start);
    printf("\n--- %lld nodes explored in %lds ---\n\n", node_count, elapsed);

    if (found) {
        printf("SOLUTION FOUND! Congruences:\n");
        for (int i = 0; i < nmod; i++)
            if (chosen_a[i] >= 0)
                printf("  n ≡ %3d  (mod %3d)\n", chosen_a[i], moduli[i]);
        printf("\nRun: python3 verify.py <solution above> to confirm.\n");
    } else {
        printf("PROVED: No odd covering system exists for L = %d.\n", L);
        printf("(Exhaustive DFS found no valid assignment.)\n");
    }

    free(mask_pool);
    return 0;
}
