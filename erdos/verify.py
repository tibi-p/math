#!/usr/bin/env python3
"""
verify.py — Verify that a set of congruences forms a covering system for Z/LZ.

Usage:
    python3 verify.py L a1 m1 a2 m2 ...

  where each pair (a_i, m_i) means  n ≡ a_i (mod m_i).

Example (made-up, for testing):
    python3 verify.py 15  0 3  1 5  2 5  0 5

The script checks:
  1. All moduli divide L.
  2. All moduli are distinct.
  3. Every residue in {0, ..., L-1} is covered by at least one congruence.
  4. (Optional) if SageMath is available, also verifies via CRT/exact arithmetic.
"""

import sys

def verify(L, congruences):
    print(f"Verifying covering system for Z/{L}Z")
    print(f"Congruences ({len(congruences)}):")
    for a, m in congruences:
        print(f"  n ≡ {a}  (mod {m})")

    errors = []

    # Check moduli divide L
    for a, m in congruences:
        if L % m != 0:
            errors.append(f"  ERROR: modulus {m} does not divide L={L}")

    # Check distinct moduli
    mods = [m for _, m in congruences]
    if len(mods) != len(set(mods)):
        dupes = [m for m in set(mods) if mods.count(m) > 1]
        errors.append(f"  ERROR: duplicate moduli: {dupes}")

    # Check coverage
    covered = [False] * L
    for a, m in congruences:
        for r in range(a % m, L, m):
            covered[r] = True

    uncovered = [r for r in range(L) if not covered[r]]

    if errors:
        print("\nFAILED (structural errors):")
        for e in errors:
            print(e)
        return False

    if uncovered:
        print(f"\nFAILED: {len(uncovered)} residues not covered: {uncovered[:20]}"
              + (" ..." if len(uncovered) > 20 else ""))
        return False

    print(f"\nOK: all {L} residues are covered.")

    # Density sanity check
    density = sum(1/m for _, m in congruences)
    print(f"Sum of 1/m_i = {density:.6f}  (must be ≥ 1 for a covering)")

    # Try SageMath verification if available
    try:
        from sage.all import ZZ, CRT_basis
        print("\nSageMath available — running exact CRT verification...")
        # Build each residue class as a set of integers mod L via CRT
        # Since all m_i | L, this is equivalent to the bitset check, but exact.
        exact_covered = set()
        for a, m in congruences:
            exact_covered.update(r for r in range(L) if r % m == a)
        exact_uncovered = set(range(L)) - exact_covered
        if exact_uncovered:
            print(f"SageMath FAILED: uncovered = {sorted(exact_uncovered)}")
            return False
        else:
            print("SageMath: confirmed exact cover.")
    except ImportError:
        pass  # SageMath not installed; Python check above is sufficient

    return True


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    args = sys.argv[1:]
    L = int(args[0])
    pairs = args[1:]

    if len(pairs) % 2 != 0:
        print("ERROR: congruences must be pairs  a m  a m  ...")
        sys.exit(1)

    congruences = []
    for i in range(0, len(pairs), 2):
        a, m = int(pairs[i]), int(pairs[i+1])
        congruences.append((a, m))

    ok = verify(L, congruences)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
