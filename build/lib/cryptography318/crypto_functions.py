from math import gcd, isqrt
from prime import IsPrime, NextPrime


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


def GCD(a, b):
    return gcd(a, b)


def ModularInverse(x, m):
    return pow(x, -1, m)


def PrimePi(p):
    """Returns number of primes <= given number"""

    if p < 2:
        raise ValueError("Must enter a number greater than the smallest prime (2)")
    primes = 1
    for n in range(3, ((p - 1) | 1) + 2, 2):
        if IsPrime(n):
            primes += 1
    return primes


def BSmoothQ(n, B):
    """Returns True if all prime factors of given number are <= given B"""

    factor_base = [p := 2]
    while len(factor_base) < PrimePi(B):
        factor_base.append(p := NextPrime(p))

    for p in factor_base:
        while n % p == 0:
            n //= p

    return True if n == 1 else False


def MakeChineseRemainder():
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
