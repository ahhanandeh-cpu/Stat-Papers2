"""Exact small-mean expansions for the Poisson-Gold family PGD(k, theta).

Uses exact truncated power-series arithmetic over the rationals.
Prints (1) the dispersion curves D_k(mu); (2) pmf differences between two
members at the same mean; (3) the leading term of half the chi-square
divergence, which equals the KL divergence at matched mean to leading order.
"""
from fractions import Fraction as Fr
from math import factorial

N = 12  # keep terms up to theta^N / mu^N


def mul(a, b):
    c = [Fr(0)] * (N + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(N + 1 - i):
            c[i + j] += ai * b[j]
    return c


def inv(a):  # 1/a, a[0] != 0
    b = [Fr(0)] * (N + 1)
    b[0] = 1 / a[0]
    for n in range(1, N + 1):
        b[n] = -sum(a[i] * b[n - i] for i in range(1, n + 1)) / a[0]
    return b


def shift(a, s):  # multiply by x^s
    return [Fr(0)] * s + a[: N + 1 - s] if s <= N else [Fr(0)] * (N + 1)


def compose(a, b):  # a(b(x)), b[0] = 0
    out = [Fr(0)] * (N + 1)
    p = [Fr(1)] + [Fr(0)] * N
    for i in range(N + 1):
        out = [o + a[i] * pi for o, pi in zip(out, p)]
        p = mul(p, b)
    return out


def poly(coefs):
    a = [Fr(0)] * (N + 1)
    for e, c in coefs.items():
        if e <= N:
            a[e] += Fr(c)
    return a


def pgd_series(k):
    """Series in theta of: mean, second factorial moment, and pmf(x)."""
    S0 = poly({j - 1: factorial(j - 1) for j in range(1, k + 1)})
    m1 = mul(poly({j: factorial(j) for j in range(1, k + 1)}), inv(S0))
    m2 = mul(poly({j + 1: factorial(j + 1) for j in range(1, k + 1)}), inv(S0))
    S1 = poly({j: factorial(j - 1) for j in range(1, k + 1)})
    # (1+theta)^(-m) as series
    def binom_neg(m):
        c = [Fr(1)]
        for n in range(1, N + 1):
            c.append(c[-1] * Fr(-(m + n - 1), n))
        return c
    def pmf(x):
        tot = [Fr(0)] * (N + 1)
        for j in range(1, k + 1):
            term = shift(binom_neg(x + j), x + j)
            tot = [a + Fr(factorial(x + j - 1), factorial(x)) * b for a, b in zip(tot, term)]
        return mul(tot, inv(S1)) if S1[0] != 0 else None
    return m1, m2, pmf


def pmf_series_theta(k, x):
    # S1 has zero constant term (starts at theta^1); factor theta out.
    S1 = poly({j - 1: factorial(j - 1) for j in range(1, k + 1)})  # S1/theta
    def binom_neg(m):
        c = [Fr(1)]
        for n in range(1, N + 1):
            c.append(c[-1] * Fr(-(m + n - 1), n))
        return c
    tot = [Fr(0)] * (N + 1)
    for j in range(1, k + 1):
        term = shift(binom_neg(x + j), x + j - 1)  # divided by theta
        tot = [a + Fr(factorial(x + j - 1), factorial(x)) * b for a, b in zip(tot, term)]
    return mul(tot, inv(S1))


def reverse(m1):
    """theta as a series in mu, given mu = m1(theta), m1[0]=0, m1[1]!=0."""
    g = [Fr(0), 1 / m1[1]] + [Fr(0)] * (N - 1)
    for _ in range(N + 1):
        comp = compose(m1, g)
        g = [gi - (ci - (1 if i == 1 else 0)) / m1[1] for i, (gi, ci) in enumerate(zip(g, comp))]
    return g


def fmt(a, var="mu", maxterms=6):
    out = []
    for i, c in enumerate(a):
        if c != 0:
            out.append(f"{c}*{var}^{i}")
        if len(out) >= maxterms:
            break
    return " + ".join(out) if out else "0"


th = {}
pm = {}
for k in range(1, 6):
    m1, m2, _ = pgd_series(k)
    th[k] = reverse(m1)
    m2mu = compose(m2, th[k])
    # D = 1 + (m2 - mu^2)/mu
    num = [a - (1 if i == 2 else 0) for i, a in enumerate(m2mu)]
    D = [Fr(1) + num[1]] + num[2:] + [Fr(0)]
    D[0] = Fr(1) + num[1]
    print(f"D_{k}(mu) = 1 + " + fmt([Fr(0)] + num[2:] + [Fr(0)], maxterms=6))
    for x in range(0, 9):
        pm[(k, x)] = compose(pmf_series_theta(k, x), th[k])

print()
for a, b in [(1, 2), (2, 3), (3, 4), (4, 5), (2, 4), (1, 3)]:
    diffs = [[u - v for u, v in zip(pm[(a, x)], pm[(b, x)])] for x in range(0, 9)]
    lead = []
    for d in diffs:
        nz = [(i, c) for i, c in enumerate(d) if c != 0]
        lead.append(f"{nz[0][1]}*mu^{nz[0][0]}" if nz else "0")
    # half chi-square: sum_x d^2 / p_a ; p_a(x) = mu^x * (series with nonzero const)
    half = [Fr(0)] * (N + 1)
    for x, d in enumerate(diffs):
        pa = pm[(a, x)]
        # divide by mu^x: shift down
        lowest = next(i for i, c in enumerate(pa) if c != 0)
        pa_red = pa[lowest:] + [Fr(0)] * lowest
        d2 = mul(d, d)
        q = mul(d2, inv(pa_red))
        q = q[lowest:] + [Fr(0)] * lowest
        half = [h + Fr(1, 2) * v for h, v in zip(half, q)]
    nz = [(i, c) for i, c in enumerate(half) if c != 0]
    print(f"k={a} vs k={b}: pmf differences x=0..4: {lead[:5]}")
    print(f"   half chi-square at matched mean = {nz[0][1]}*mu^{nz[0][0]} + {nz[1][1]}*mu^{nz[1][0]} + ...")
