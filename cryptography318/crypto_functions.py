import operator
from math import gcd, isqrt
from .prime import IsPrime, NextPrime
from .deprecated import deprecated


def apply(proc, lst):
    lst_proc = {'*': operator.mul, '+': operator.add, '-': operator.sub, '/': operator.truediv}
    if proc in lst_proc:
        acc = lst[0]
        for e in lst[1:]:
            acc = lst_proc[proc](acc, e)
        return acc
    return proc(lst[:])


def toBase(n, base):
    # finds highest exponent of base
    exp = 0
    while True:
        if n < pow(base, exp + 1):
            break
        exp += 1

    # iterates starting at highest exponent down to 0, creates list of 'base' base
    list_coeff = [0 for _ in range(exp + 1)]
    ind = 0
    for i in range(exp, -1, -1):
        k = n // pow(base, i)
        list_coeff[ind] = k
        n -= k * pow(base, i)
        ind += 1

    return list_coeff


def fromBase(lst, base):
    acc = 0
    l = len(lst)
    for i in range(l):
        acc += lst[i] * pow(base, l - i - 1)
    return acc


def StringToNum(s, base=128):
    result, i = 0, 0
    for c in s[::-1]:
        result += (base ** i) * ord(c)
        i += 1
    return result


def NumToString(n, base=128):
    string = ''
    chars = toBase(n, base)
    for c in chars:
        string += chr(c)
    return string


@deprecated
def ApplyMult(lst):
    acc = 1
    for e in lst:
        acc *= e
    return acc


def ExtendedGCD(a, b):
    if a == 0:
        return b, 0, 1

    g, x, y = ExtendedGCD(b % a, a)
    return g, y - (b // a) * x, x


@deprecated
def GCD(a, b):
    return gcd(a, b)


@deprecated
def ModularInverse(x, m):
    return pow(x, -1, m)


def PrimesLT(p):
    if p < 2:
        raise ValueError("Must enter a number greater than the smallest prime (2)")
    primes = [2]
    for n in range(3, ((p - 1) | 1) + 2, 2):
        if IsPrime(n):
            primes.append(n)
    return primes


def PrimePi(p):
    """Returns number of primes <= given number"""

    return len(PrimesLT(p))


def BSmoothQ(n, B):
    """Returns True if all prime factors of given number are <= given B"""

    factor_base = [p := 2]
    while len(factor_base) < PrimePi(B):
        factor_base.append(p := NextPrime(p))

    for p in factor_base:
        while n % p == 0:
            n //= p

    return True if n == 1 else False


def makeChineseRemainder():
    # gets lists of solutions and moduli from user for chinese remainder
    nums, mods = [], []
    equations = int(input("How many equations: "))
    for i in range(equations):
        print("x = ", end="")
        x = int(input())
        print("mod ", end="")
        m = int(input())
        print()
        nums.append(x)
        mods.append(m)
    return nums, mods


def ChineseRemainder(nums, mods):
    # initializes lists of moduli, M = product of all moduli
    M = ApplyMult(mods)

    # maps list of moduli and their inverses to x and y respectively
    x, y = [], []
    for m in mods:
        mi = M // m
        x.append(mi)
        y.append(pow(mi, -1, m))

    # accumulates product of numbers and moduli and their inverses
    acc = 0
    for i in range(len(nums)):
        acc = (acc + nums[i] * x[i] * y[i]) % M

    return acc


def BSGS(g, h, p, prog=False, N=None):
    if N is None:
        N = p - 1
    n = isqrt(N) + 1

    if prog:
        increment = n // 25
        print("Starting Part 1: Creating Lists")
        A, B = {}, {}
        count = 0
        for i in range(n):
            if count >= increment:
                count = 0
                print("-", end="")
            count += 1
            j = i * n
            A[pow(g, i, p)] = i
            B[(h * pow(g, -1 * j, p)) % p] = j

        print("\nDone With Part 1. Starting Part 2: Finding Matches")

        count = 0
        for e in A:
            if count >= increment:
                count = 0
                print("-", end="")
            count += 1
            if e in B:
                print()
                return A[e] + B[e]
        print()

    else:
        A, B = {}, {}
        for i in range(n):
            j = i * n
            A[pow(g, i, p)] = i
            B[(h * pow(g, -1 * j, p)) % p] = j

        for e in A:
            if e in B:
                return A[e] + B[e]
    return None


def PohligHellman(g, h, p, q, exp, prog=False):
    X = []
    if prog:
        print("Starting process for X0")

    r = pow(g, pow(q, exp - 1), p)
    X0 = BSGS(r, pow(h, pow(q, exp - 1), p), p, prog, q)
    X.append(X0)

    if prog:
        print(f"Found X0 = {X0}\n")

    for i in range(1, exp):
        if prog:
            print(f"Starting process for X{i}")
        exp_term = fromBase(X[::-1], q)
        h_term = pow(h * pow(pow(g, exp_term, p), -1, p), pow(q, exp - i - 1), p)
        Xi = BSGS(r, h_term, p, prog, q)
        X.append(Xi)
        if prog:
            print(f"Found X{i} = {Xi}\n")

    return fromBase(X[::-1], q)


def DSA(D, S1, S2, g, p, q, A):
    S2_inv = pow(S2, -1, q)
    V1 = (D * S2_inv) % q
    V2 = (S1 * S2_inv) % q
    if ((pow(g, V1, p) * pow(A, V2, p)) % p) % q == S1:
        return True
    return False


def PollardP1(n, limit=pow(10, 6), first_n=4):
    """Pollard's p - 1 algorithm for factoring large composites.
    Returns one non-trivial factor if factor-able, False if otherwise."""

    if IsPrime(n):
        raise ValueError("Make sure to enter a composite number")

    if not 1 < first_n < 8:
        first_n = 8

    for a in [2, 3, 5, 7, 11, 13, 17, 19][:first_n]:
        m = a
        for j in range(2, limit):
            m = pow(m, j, n)
            k = gcd(m - 1, n)
            if 1 < k < n:
                return k

    return False


def FactorInt(n):
    """
    Function that checks if number has small prime factors, then attempts Pollards p-1 algorithm for
    factoring large composite numbers in the form N = p * q, returning one non-trivial factor of N. If neither
    of these methods factor N, sympy.factorint function is used to further factor N, if possible.

    Returns a Python dictionary with each key being a prime factor and the associated value being the power of
    that prime factor.
    """

    if IsPrime(n):
        return {n: 1}

    known_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
                    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
                    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443,
                    449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
                    587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
                    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
                    853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983,
                    991, 997]

    factors = {}
    for p in known_primes:
        while n % p == 0:
            if IsPrime(n):
                if n not in factors:
                    factors[n] = 1
                else:
                    factors[n] += 1
                return factors
            if p not in factors:
                factors[p] = 1
            else:
                factors[p] += 1
            n //= p

    while not IsPrime(n):
        k = PollardP1(n)
        # if Pollard p-1 returns False, try using sympy.factorint
        if not k:
            from sympy import factorint
            sy_factors = factorint(n)
            for e in sy_factors:
                if e not in factors:
                    factors[e] = sy_factors[e]
                else:
                    factors[e] += sy_factors[e]
            return factors

        n //= k
        if k not in factors:
            factors[k] = 1
        else:
            factors[k] += 1

    if n != 1 and IsPrime(n):
        if n not in factors:
            factors[n] = 1
        else:
            factors[n] += 1
        return factors

    return QuadraticSieve(n)


def _factorWithKnown(p, q, N):
    factors = {}
    factor_p = FactorInt(p)
    for n in factor_p:
        factors[n] = factor_p[n]
    factor_q = FactorInt(q)
    for n in factor_q:
        factors[n] = factor_q[n]

    while N % p != 0:
        factors[p] += 1
        N //= p
    while N % q != 0:
        factors[q] += 1
        N //= q
    if N == 1:
        return factors
    if IsPrime(N):
        factors[N] = 1
        return factors
    more_factors = FactorInt(N)
    for f in more_factors:
        factors[f] = more_factors[f]
    return factors


def _factorPerfectSquare(N, B=7):
    from itertools import combinations_with_replacement as _all

    m = PrimePi(B)
    b_smooth_nums = []
    squared_nums = {}
    a = isqrt(N) - 1
    while len(b_smooth_nums) < m:
        ci = pow(a, 2, N)
        if BSmoothQ(ci, B):
            b_smooth_nums.append(ci)
            squared_nums[ci] = a
        a += 1

    ci_factors = {}
    factor_base = PrimesLT(B)
    mat = []
    for num in b_smooth_nums:
        ci = num
        exp = []
        for p in factor_base:
            count = 0
            while num % p == 0:
                num //= p
                count += 1
            exp.append(count)
        ci_factors[ci] = exp
        mat.append(exp)

    for c in range(2, len(mat)):
        possible = list(map(list, list(_all(mat, c))))
        for comb in possible:
            acc = comb[0][:]
            perfect_2 = True

            for i in range(1, len(comb)):
                for j in range(len(acc)):
                    acc[j] += comb[i][j]
            for j in acc:
                if j % 2 != 0:
                    perfect_2 = False

            if not perfect_2:
                continue
            a, b = 1, 1
            for num in b_smooth_nums:
                choice = ci_factors[num]
                if choice in comb:
                    a *= squared_nums[num]
            for k in range(len(acc)):
                b *= pow(factor_base[k], acc[k] // 2)

            p, q = gcd(a - b, N), gcd(a + b, N)
            if 1 < p < N and 1 < q < N:
                print(f"p: {p}, q: {q}")
                if p * q == N:
                    return {p: 1, q: 1}
                return _factorWithKnown(p, q, N)
            if 1 < p < N:
                q = N // p
                print(f"p: {p}, q: {q}")
                if IsPrime(q) and N == p * q:
                    return {p: 1, q: 1}
                if IsPrime(q):
                    return _factorWithKnown(p, q, N)
                q_factors = FactorInt(q)
                if p in q_factors:
                    q_factors[p] += 1
                else:
                    q_factors[p] = 1
                return q_factors
            if 1 < q < N:
                p = N // q
                print(f"p: {p}, q: {q}")
                if IsPrime(p) and N == p * q:
                    return {p: 1, q: 1}
                if IsPrime(p):
                    return _factorWithKnown(p, q, N)
                p_factors = FactorInt(p)
                if q in p_factors:
                    p_factors[q] += 1
                else:
                    p_factors[q] = 1
                return p_factors

    return False


def QuadraticSieve(N, B=None):
    """Performs Quadratic Sieve Algorithm with a given Smoothness value B on a given composite N"""

    from itertools import combinations_with_replacement as _all
    from math import e, log, sqrt

    if B is None:
        L = pow(e, sqrt(log(N) * log(log(N))))
        B = int(pow(L, 1 / sqrt(2)))
        print(f"L: {L}\n"
              f"B: {B}")

    m = PrimePi(B)
    sq = isqrt(N)
    while pow(sq, 2) < N:
        sq += 1

    b_smooth_nums = []
    squared_nums = {}
    factor_base = PrimesLT(B)
    while True:
        c_i = pow(sq, 2, N)
        sq += 1
        if BSmoothQ(c_i, B):
            exp = [0 for _ in range(m)]
            for i in range(m):
                prime = factor_base[i]
                while c_i % prime == 0:
                    exp[i] += 1
                    c_i //= prime
            squared_nums[sq] = exp
            b_smooth_nums.append(exp)
            if len(b_smooth_nums) < 3:
                continue

            for c in range(2, len(b_smooth_nums)):
                choices = list(map(list, list(_all(b_smooth_nums, c))))
                for choice in choices:
                    if exp not in choice:
                        continue
                    exp_sum = [0 for _ in range(m)]
                    b = 1
                    valid = True
                    for power in choice:
                        for i in range(m):
                            exp_sum[i] += power[i]
                    for i in range(m):
                        n = exp_sum[i]
                        if n > 2 and n % 2 != 0:
                            valid = False
                            break
                        b *= pow(factor_base[i], n // 2)
                    if not valid:
                        continue

                    a = 1
                    for sq in squared_nums:
                        if squared_nums[sq] in choice:
                            a *= sq

                    # checks to see if a + b or a - b is a multiple of a factor of N
                    p = gcd(N, a + b)
                    q = gcd(N, a - b)

                    if 1 < p < N and 1 < q < N:
                        return _factorWithKnown(p, q, N)
                    if 1 < p < N:
                        q = N // p
                        if IsPrime(q) and IsPrime(p):
                            return {p: 1, q: 1}
                        return _factorWithKnown(p, q, N)
                    if 1 < q < N:
                        p = N // q
                        if IsPrime(p) and IsPrime(q):
                            return {p: 1, q: 1}
                        return _factorWithKnown(p, q, N)
