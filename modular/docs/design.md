# Design Document: Modular Symbols Engine for S₂(Γ₀(N))

## Architectural Philosophy

Every struct is a **formal algebraic object**. Fields are never accessed directly outside their own module. All arithmetic operations check for overflow and return a status code. The dependency graph is strictly layered:

```
rat → cusp → p1 → manin → sparse → hecke → linalg
```

No layer may include a header from a layer above it.

---

## Phase 1: Data Structures & Memory

### 1. Exact Rational Arithmetic — `Rat`

**Mathematical object:** An element of ℚ, stored in canonical form.

**Invariants (enforced after every operation):**
- `den > 0` always.
- `gcd(|num|, den) = 1` (fully reduced).
- Zero is `{0, 1}`. It is a hard error to construct `{0, 0}`.

**Struct:**
```c
typedef struct {
    int64_t num;   /* signed numerator             */
    int64_t den;   /* strictly positive denominator */
} Rat;
```

**Functions (signatures only, no bodies yet):**
```c
/* Euclidean algorithm; result is always >= 0 */
int64_t i64_gcd(int64_t a, int64_t b);

/* Reduce r in-place. Returns 0 on success, -1 on den==0. */
int rat_reduce(Rat *r);

/* Arithmetic; return -1 on overflow, 0 on success. */
int rat_add(Rat a, Rat b, Rat *out);
int rat_sub(Rat a, Rat b, Rat *out);
int rat_mul(Rat a, Rat b, Rat *out);
int rat_inv(Rat a, Rat *out);   /* fails if a.num == 0 */

int rat_eq(Rat a, Rat b);       /* 1 if equal, 0 otherwise */
int rat_is_zero(Rat a);
```

**Overflow strategy:** Before computing `a*b` for `int64_t`, check `|a| <= INT64_MAX / |b|`. If not, abort with a hard error — the caller must use a smaller conductor. No silent truncation, ever.

**Memory:** All `Rat` values are stack-allocated. No heap allocation in this layer.

---

### 2. The Projective Line P¹(ℚ) — Cusps

**Mathematical object:** An element of ℙ¹(ℚ) = ℚ ∪ {∞}, representing a cusp in the extended upper half-plane.

**Invariants:**
- A *finite* cusp p/q is stored with `q > 0` and `gcd(|p|, q) = 1`.
- The cusp ∞ = (1:0) is the unique case where `q == 0` and `p == 1`.
- No other `q == 0` is valid.

**Struct:**
```c
typedef enum { CUSP_FINITE = 0, CUSP_INF = 1 } CuspKind;

typedef struct {
    CuspKind kind;
    int64_t  p;    /* numerator   (meaningful when kind == CUSP_FINITE) */
    int64_t  q;    /* denominator (positive, meaningful when FINITE)    */
} Cusp;
```

**Key operations:**
```c
Cusp cusp_inf(void);
Cusp cusp_from_frac(int64_t p, int64_t q);  /* reduces and normalises */
int  cusp_eq(Cusp a, Cusp b);
```

**Role:** Used only in the boundary map δ : M₂(Γ₀(N)) → ℤ[cusps]. The cusps of Γ₀(N) are equivalence classes of P¹(ℚ) under the Γ₀(N) action; their count equals Σ_{d|N} φ(gcd(d, N/d)).

---

### 3. P¹(ℤ/Nℤ) — The Manin Symbol Index Table

**Mathematical object:** The projective line over ℤ/Nℤ.

An element (c : d) ∈ ℙ¹(ℤ/Nℤ) satisfies:
- `0 ≤ c, d < N`
- **Primitivity condition:** For every prime p | N, NOT (p | c AND p | d). Equivalently, `gcd(gcd(c,d), N) == 1`.
- Two pairs (c, d) and (c', d') are equivalent iff ∃ λ ∈ (ℤ/Nℤ)* s.t. c' ≡ λc, d' ≡ λd (mod N).

The total count is μ = [SL₂(ℤ) : Γ₀(N)] = N · ∏_{p | N} (1 + 1/p). This is the **degree of the modular curve X₀(N)** and equals the number of Manin symbols.

**Canonical representative:** Among all (λc mod N, λd mod N) for λ ∈ (ℤ/Nℤ)*, choose the lexicographically smallest pair. This is the *normalised* representative.

**Struct (a single element):**
```c
typedef struct {
    int64_t c;
    int64_t d;
} P1Elem;
```

**The P1 Table** (allocated once per conductor N, used throughout):
```c
typedef struct {
    int     N;
    int     mu;           /* total count = [SL_2(Z) : Gamma_0(N)]        */
    P1Elem *elems;        /* elems[i] = canonical (c,d) for index i      */
    int    *index_table;  /* flat N×N array; index_table[c*N+d] = index  */
                          /* -1 if (c,d) is not canonical                 */
} P1Table;
```

**Memory:** `elems` is a single `malloc(mu * sizeof(P1Elem))`. `index_table` is `malloc(N * N * sizeof(int))`. Both are freed together via `p1table_free`. No other dynamic allocation in this layer.

**Key functions:**
```c
P1Table *p1table_build(int N);
void     p1table_free(P1Table *t);

/* Returns index of (c,d) after normalisation, or -1 if not in P1. */
int      p1_index(const P1Table *t, int64_t c, int64_t d);

/* Lift canonical (c,d) to a matrix [[a,b],[c,d]] in SL_2(Z).      */
/* Uses extended Euclidean + adjustment of c by multiples of N.     */
int      p1_lift_to_sl2z(const P1Table *t, int idx,
                          int64_t *a, int64_t *b,
                          int64_t *c, int64_t *d);
```

**Lift algorithm sketch:** Given canonical (c, d) with gcd(c, d) possibly > 1, find c' = c + k·N such that gcd(c', d) = 1 (always possible when gcd(gcd(c,d), N) = 1). Then apply extended Euclidean to get a, b with ad - bc' = 1. Store c' internally for lift purposes only; the P1 index is still keyed on (c mod N, d mod N).

---

### 4. Manin Symbols — The Free Module Generators

**Mathematical object:** A Manin symbol for Γ₀(N) is an element of ℙ¹(ℤ/Nℤ), understood as a generator of the free ℤ-module M₂(N) = ℤ[ℙ¹(ℤ/Nℤ)] before quotienting by relations.

At the code level, a Manin symbol is simply its P1 index. The struct exists for clarity at call sites:

```c
typedef struct {
    int idx;   /* index into P1Table::elems, range [0, mu) */
} ManinSym;
```

A **formal linear combination** of Manin symbols (a module element):

```c
typedef struct {
    int      len;
    int     *indices;  /* ManinSym indices                    */
    Rat     *coeffs;   /* corresponding rational coefficients */
    int      cap;      /* allocated capacity                  */
} ManinElt;
```

`ManinElt` is heap-allocated. `indices` and `coeffs` are parallel arrays, always kept sorted by index for fast merging. An invariant is that no index appears twice; they are merged with coefficient addition on insert.

---

### 5. Sparse Matrix over ℚ — `SparseMat`

**Mathematical object:** A morphism M : ℚⁿ → ℚᵐ, represented with minimal memory use.

This is the workhorse for: encoding the relation matrix, the boundary map, and the Hecke operator matrices.

**Format:** Compressed Sparse Row (CSR) for read-only passes (Gaussian elimination). During construction, use a dynamic **row-of-lists** format, then convert.

**Construction-phase struct:**
```c
typedef struct {
    int   col;
    Rat   val;
} SEntry;

typedef struct {
    int      nrows, ncols;
    SEntry **rows;     /* rows[i] = malloc'd array of nonzero entries */
    int     *lens;     /* lens[i] = current length of rows[i]         */
    int     *caps;     /* caps[i] = allocated capacity of rows[i]     */
} SparseMatBuilder;
```

**Frozen (elimination-phase) struct — CSR format:**
```c
typedef struct {
    int   nrows, ncols;
    int  *row_ptr;   /* row_ptr[i]..row_ptr[i+1]-1 are indices for row i */
    int  *col_ind;   /* column indices of nonzero entries                 */
    Rat  *vals;      /* values, parallel to col_ind                       */
} SparseMat;
```

**Key operations:**
```c
/* Construction */
SparseMatBuilder *smat_builder_new(int nrows, int ncols);
int               smat_add(SparseMatBuilder *B, int row, int col, Rat val);
SparseMat        *smat_freeze(SparseMatBuilder *B);
void              smat_builder_free(SparseMatBuilder *B);

/* Elimination */
void  smat_free(SparseMat *M);
int   smat_gaussian_elim(SparseMat *M,  /* in-place, over Q           */
                          int **pivot_cols,
                          int  *rank);
void  smat_kernel_basis(SparseMat *M, SparseMat **ker); /* after elim */
```

**Pivoting strategy:** Partial pivoting by sparsity (choose the row with the fewest nonzeros as pivot — this is essential for performance with modular-symbol matrices which have a very specific block structure). All row operations are exact rational arithmetic.

---

## Phase 2: Core Algorithms

### Algorithm 1: Generating ℙ¹(ℤ/Nℤ)

**Input:** conductor N ∈ ℤ₊.
**Output:** `P1Table` with μ = [SL₂(ℤ) : Γ₀(N)] elements.

**Steps:**
1. Factorise N into distinct prime factors p₁, …, pₖ (trial division; N is small).
2. Allocate `index_table[N][N]` initialised to -1.
3. For each (c, d) with 0 ≤ c, d < N:
   a. Check primitivity: compute g = gcd(gcd(c, d), N). If g ≠ 1, skip.
   b. Compute canonical representative: iterate λ over (ℤ/Nℤ)* (all λ with gcd(λ, N) = 1), compute (λc mod N, λd mod N), find lex-min pair (c₀, d₀).
   c. If `index_table[c₀][d₀] == -1`: assign next index i, store (c₀, d₀) in `elems[i]`, set `index_table[c₀][d₀] = i`.
4. Assert final count equals N · ∏_{p | N}(1 + 1/p).

**Complexity:** O(N² · φ(N)) — acceptable for N ≤ 200 (μ ≤ ~400 for typical conductors).

---

### Algorithm 2: The Relation Matrix — Quotient by S and T Relations

The free module ℤ[ℙ¹(ℤ/Nℤ)] must be quotiented by two families of relations to yield the space of modular symbols M₂(Γ₀(N)). The **cuspidal subspace** C₂(Γ₀(N)) (which equals S₂(Γ₀(N)) as a ℚ-vector space via integration) is then the kernel of the boundary map δ within this quotient.

**The two generators of PSL₂(ℤ) used:**

Let S = [[0, -1], [1, 0]] and U = [[0, -1], [1, -1]] (order 3 in PSL₂(ℤ)).

Their left actions on P¹(ℤ/Nℤ) are:
```
S · (c : d)  =  (-d : c)  mod N
U · (c : d)  =  (-d : c-d) mod N
U² · (c : d) =  (d-c : -c) mod N
```

**S-relation (2-term):** For each index i with element (c, d):
```
x_i + x_{S·i} = 0
```
where `S·i = p1_index(t, -d mod N, c mod N)`.
Special case: if `S·i == i`, then `2·x_i = 0` (torsion; over ℚ this forces x_i = 0).

**T-relation (3-term):** For each orbit {i, U·i, U²·i}:
```
x_i + x_{U·i} + x_{U²·i} = 0
```
Orbits of size 1 or 2 (fixed points of U or U²) yield lower-order torsion relations; over ℚ, fixed points are forced to zero.

**Implementation:** Build a `SparseMatBuilder` with μ columns and one row per independent relation. The kernel of this matrix (via `smat_gaussian_elim` + `smat_kernel_basis`) is a ℚ-basis for the full modular symbols space M₂(Γ₀(N)).

**Boundary map δ:**

For Manin symbol with index i (element (c, d), lifted to [[a, b], [c, d]] ∈ SL₂(ℤ)):
```
δ(x_i) = [a/c]_{Γ₀(N)} - [b/d]_{Γ₀(N)}
```
where a/c and b/d are cusps of Γ₀(N) (with ∞ = 1/0). The cusps of Γ₀(N) are indexed by divisors d | N with gcd(d, N/d) = 1. Build δ as a `SparseMat` with rows = cusp indices, columns = Manin symbol indices.

The **cuspidal subspace** is ker(δ) ∩ (quotient by S and T relations). This is computed by chaining the two row reduction steps.

---

### Algorithm 3: Hecke Operator T_p Action on Manin Symbols

**Input:** prime p, Manin symbol index i (with bottom row (c, d)).
**Output:** `ManinElt` representing T_p · x_i (or U_p · x_i when p | N).

The two cases — p ∤ N and p | N — require fundamentally different formulas.

---

#### Case 1: T_p, p ∤ N — Heilbronn matrices

The naïve approach of picking coset representatives `[[1,j],[0,p]]` and left-multiplying the SL₂(ℤ) lift γ fails: the bottom rows `(pc, pd)` of `M·γ` are identical for all j, producing p copies of the same symbol. The correct tool is the **Heilbronn matrix set**.

A **Heilbronn matrix** for p is an integer matrix `H = [[A,B],[C,D]]` satisfying:
```
det(H) = A*D - B*C = p,   A > B ≥ 0,   D > C ≥ 0.
```

Merel proved that summing over all Heilbronn matrices for p gives the correct Hecke action on modular symbols: the set automatically re-triangulates broken Farey paths back onto the fundamental domain edges without requiring explicit coset reduction.

**Right action on the bottom row:** For Manin symbol [c:d] and Heilbronn matrix H:
```
(c, d) · [[A, B], [C, D]]  =  (c*A + d*C,  c*B + d*D)
```
Reduce by GCD, then look up in the P1 table:
```
nc = c*A + d*C;   nd = c*B + d*D
g  = gcd(nc, nd);  nc /= g;  nd /= g
idx = p1_index(t, nc % N, nd % N)
```

**Enumeration of Heilbronn matrices:** Iterate over A ∈ [1,p], B ∈ [0,A), C ∈ [0,p−A]. Compute `rem = p + B*C`; emit `[[A,B],[C,D]]` with `D = rem/A` only when `A | rem` and `D > C`.

---

#### Case 2: U_p, p | N — explicit formula

When p | N, the operator is U_p (not T_p) and the Heilbronn approach does not apply. The action splits on whether p divides the bottom-left coordinate c.

**Subcase 2a — p | c:** The p transversals collapse to a single image:
```
[c : d]  →  [c/p : d]
```
Look up `p1_index(t, (c/p) % N, d % N)`.

**Subcase 2b — p ∤ c:** The action produces p distinct images:
```
for j = 0..p-1:  [c : d]  →  [c*p : d*p - j*c]   (then reduce by GCD)
```
For each j: `nc = c*p`, `nd = d*p - j*c`, reduce by `gcd(nc, nd)`, then look up.

---

**Building the T_p matrix:** For each basis vector of the cuspidal subspace (a row of `cs->basis`, a μ-vector), compute the Hecke image as a `ManinElt` by applying the action to each nonzero Manin symbol and accumulating. Then project the result into cuspidal coordinates via `cs->coord` to get one column of the dim × dim Hecke matrix.

---

### Algorithm 4: Exact Linear Algebra — Kernel & Eigenvectors

**Gaussian elimination over ℚ (fraction-free preferred):**

Use the **Bareiss algorithm** for intermediate steps to keep numerator/denominator sizes bounded. If overflow is detected, abort and report that N is too large for int64_t arithmetic.

**Kernel computation:** After row-reducing the relation+boundary matrix to reduced row echelon form (RREF), the kernel is read off from the free (non-pivot) columns via back-substitution. Each kernel basis vector is a column of rational numbers.

**Eigenvector decomposition:**

The Hecke operators T_p (for p prime, p ≤ some bound B) are simultaneously diagonalisable over ℚ (Atkin–Lehner theory). The algorithm is:

1. Start with the full cuspidal subspace V.
2. For each prime p in order: compute the T_p matrix, compute its minimal polynomial (via the Cayley-Hamilton sequence `I, T_p, T_p², …` using the current basis), and split V into T_p-eigenspaces.
3. Repeat until every subspace is 1-dimensional. At that point, each basis vector is a Hecke eigenform.
4. Normalise each eigenvector: scale so the leading coefficient is 1 (corresponding to the q-expansion coefficient a₁ = 1).

**Overflow guard:** At each matrix multiply step, check if any intermediate `Rat` value has `|num|` or `den` exceeding a threshold (e.g., 10¹⁵). If so, report a graceful overflow error.

---

## Proposed File Layout

```
modular/
├── docs/
│   └── design.md      # this document
├── rat.h              # Rat struct, i64_gcd, arithmetic ops
├── rat.c
├── cusp.h             # Cusp struct, P¹(Q) representation
├── cusp.c
├── p1.h               # P1Elem, P1Table, p1_index, p1_lift_to_sl2z
├── p1.c
├── manin.h            # ManinSym, ManinElt, relation builder
├── manin.c
├── sparse.h           # SparseMatBuilder, SparseMat, CSR format
├── sparse.c
├── linalg.h           # Gaussian elim, kernel, RREF
├── linalg.c
├── hecke.h            # T_p action, Hecke matrix builder, eigenvectors
├── hecke.c
└── main.c             # Driver: read N from argv, print eigenvalues
```

---

## Summary of Dependencies and Ownership

| Layer | Owns | Depends on |
|---|---|---|
| `rat` | `Rat`, `i64_gcd` | `<stdint.h>` only |
| `cusp` | `Cusp`, cusp arithmetic | `rat` |
| `p1` | `P1Table`, `P1Elem`, lift | `cusp`, `rat` |
| `manin` | `ManinSym`, `ManinElt`, relations, boundary map | `p1`, `rat` |
| `sparse` | `SparseMat`, builder, CSR | `rat` |
| `linalg` | elimination, kernel, RREF | `sparse`, `rat` |
| `hecke` | T_p action, Hecke matrix, eigenvectors | all above |
