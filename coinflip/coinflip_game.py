from collections import defaultdict
import math
import sys


def binomial(n, k):
    if n < 0 and k == n:
        return 1
    elif n < 0 or k < 0:
        return 0
    else:
        return math.comb(n, k)


def inclusive_range(l, r):
    return range(l, r + 1)


def flip_combinations(n, a, b):
    u = binomial(a + b - 1, b - 1) * binomial(n - a - b, b)
    v = binomial(a + b, b) * binomial(n - a - b - 1, b)
    return u + v


def flip_combinations_h(n, a, b):
    u = binomial(a + b - 1, b - 1) * binomial(n - a - b - 1, b - 1)
    v = binomial(a + b, b) * binomial(n - a - b - 2, b - 1)
    return u + v


def flip_combinations_t(n, a, b):
    return flip_combinations(n - 1, a, b)


def sum_for_diff(flip_callback, n, diff):
    s = 0
    for b in inclusive_range(0, n // 2):
        a = b + diff
        if a >= 0 and a <= n - 2 * b:
            s += flip_callback(n, a, b)
    return s


def simulate(n):
    results = [0, 0, 0]
    scores = defaultdict(int)
    for i in range(1 << n):
        v = []
        k = i
        for j in range(n):
            v.append(k & 1)
            k >>= 1

        left = 0
        right = 0
        for j in range(n - 1):
            if v[j] == 0:
                if v[j + 1] == 0:
                    left += 1
                else:
                    right += 1

        score = left - right
        if score > 0:
            idx = 0
        elif score == 0:
            idx = 1
        else:
            idx = 2
        results[idx] += 1
        scores[(left, right)] += 1

    diff = results[2] - results[0]
    print(f"n={n}:", results, diff, len(scores))

    d = 0
    for b in inclusive_range(0, n // 2):
        for a in inclusive_range(0, n - 2 * b):
            f = flip_combinations(n, a, b)
            g = scores.get((a, b), 0)
            if f != g:
                print(f"a={a} b={b} f={f} g={g}")
                print()
            d += 1
    print(len(scores), d)

    sa0 = sum_for_diff(flip_combinations_h, n - 1, 0)
    sa1 = sum_for_diff(flip_combinations_h, n - 1, -1)
    sb0 = sum_for_diff(flip_combinations_t, n - 1, 0)
    sb1 = sum_for_diff(flip_combinations_t, n - 1, 1)
    print("A gain", sa0, sa1)
    print("B gain", sb0, sb1)

    return results


def main(n):
    left = n
    for i in range(left, n + 1):
        simulate(i)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    else:
        n = 3
    main(n)
