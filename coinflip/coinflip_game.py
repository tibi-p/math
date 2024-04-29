import sys


def simulate(n):
    results = [0, 0, 0]
    scores = []
    for i in range(1 << n):
        v = []
        k = i
        for j in range(n):
            v.append(k & 1)
            k >>= 1

        score = 0
        for j in range(n - 1):
            if v[j] == 0:
                if v[j + 1] == 0:
                    score += 1
                else:
                    score -= 1

        if score > 0:
            idx = 0
        elif score == 0:
            idx = 1
        else:
            idx = 2
        results[idx] += 1
        #scores.append(score)

    #scores.sort()
    diff = results[2] - results[0]
    print(f"n={n}:", results, diff)
    return results


def main(n):
    for i in range(1, n + 1):
        simulate(i)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    else:
        n = 3
    main(n)
