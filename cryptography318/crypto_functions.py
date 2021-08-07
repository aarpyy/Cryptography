import operator
from math import gcd, isqrt, sqrt
from .prime import IsPrime, NextPrime
from .tools import deprecated


class EllipticCurve:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.mod = p

    def point(self, x=None, y=None):
        return Elliptic(self, x, y)

    def string(self, s):
        return StringToElliptic(self, s)


class Elliptic:
    def __init__(self, E, x=None, y=None):
        self.E = E
        if x is not None and y is not None:
            self.point = (x, y)
        else:
            self.point = None

    def __neg__(self):
        return Elliptic(self.E, self.point[0], -self.point[1])

    def __hash__(self):
        return hash(self.point)

    def __getitem__(self, item):
        return self.point[item]

    def __add__(self, other):
        if other.point is None:
            return self
        if self.point is None:
            return other
        x1, y1 = self.point
        x2, y2 = other.point
        if (x1 - x2) % self.E.mod == 0 and (y1 + y2) % self.E.mod == 0:
            return None
        if (x1 - x2) % self.E.mod == 0 and (y1 - y2) % self.E.mod == 0:
            slope = ((3 * pow(x1, 2) + self.E.a) * pow(2 * y1, -1, self.E.mod)) % self.E.mod
        else:
            slope = ((y2 - y1) * pow(x2 - x1, -1, self.E.mod)) % self.E.mod
        x3 = (pow(slope, 2) - x1 - x2) % self.E.mod
        y3 = (slope * (x1 - x3) - y1) % self.E.mod
        return Elliptic(self.E, x3, y3)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        P = -self
        return P.__add__(other)

    def __mul__(self, other):
        return self.__rmul__(other)

    def __rmul__(self, other):
        if isinstance(other, float):
            other = int(other)
        P = self
        if other < 0:
            other = abs(other)
            P = -self
        Q = P
        R = Elliptic(self.E)
        while other > 0:
            if other & 1:
                R += Q
            Q += Q
            other //= 2
        return R

    def __str__(self):
        return str(self.point)

    def __eq__(self, other):
        if isinstance(other, Elliptic):
            return self.point == other.point
        return self.point == other

    def to_string(self):
        return EllipticToString(self)


def apply(proc, lst):
    lst_proc = {'*': operator.mul, '+': operator.add, '-': operator.sub, '/': operator.truediv}
    if proc in lst_proc:
        acc = lst[0]
        for e in lst[1:]:
            acc = lst_proc[proc](acc, e)
        return acc
    return proc(*lst)


def sqrt_safe(n):
    try:
        result = sqrt(n)
    except OverflowError:
        result = isqrt(n)
    finally:
        return result


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
        k, n = divmod(n, pow(base, i))
        list_coeff[ind] = k
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


def StringToElliptic(E, s):
    for i in range(100):
        x = StringToNum(s) * 100 + i
        y = sqrt_safe(pow(x, 3, E.mod) + E.a * x + E.b) % E.mod
        if isinstance(y, int):
            return E.point(x, y)
    return None


def EllipticToString(n):
    return NumToString(n[0] // 100)


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


def BSmoothQ(n, B=None, factors=None):
    """Returns True if all prime factors of given number are <= given B"""
    if factors is None:
        factor_base = [p := 2]
        while len(factor_base) < PrimePi(B):
            factor_base.append(p := NextPrime(p))
    else:
        factor_base = factors

    for p in factor_base:
        while n % p == 0:
            n //= p

    return n == 1


def factor_base_exp(n, factors):
    """Function that takes in number and list of factors and returns a list of each of the powers of each factor."""

    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while n % f == 0:
            n //= f
            exp[i] += 1
    return exp


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


def baby_step_giant_step(g, h, p, prog=False, N=None):
    """Function attempts to solve DLP using classic baby-step-giant-step algorithm."""
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


def elliptic_bsgs(P, Q, N=None):
    if N is None:
        N = P.E.mod + 1 + 2 * isqrt(P.E.mod)

    n = isqrt(N)

    # creates baby-step table, P * i for i from 0 to n
    a = list(map(lambda i: P * i, range(n)))
    # creates giant-step table, Q - P * jn for j from 0 to n
    b = list(map(lambda j: Q + (P * (-j * n)), range(n)))

    # creates sets out of lists
    A, B = set(a), set(b)

    U = A.intersection(B)
    if not U:
        return None
    point = U.pop()

    # uses lists to find index of intersection
    index_a, index_b = -1, -1
    for i, pair in enumerate(zip(a, b)):
        if pair[0] == point:
            index_a = i
        if pair[1] == point:
            index_b = i
        if -1 not in (index_a, index_b):
            return index_a + n * index_b
    return None


def index_calculus_dlp(g, h, p):
    """Function attempts to solve DLP through index calculus algorithm, computing a series of smaller dlp's
    used to solve the larger."""

    from math import e, sqrt, log
    from functools import reduce

    L = pow(e, sqrt(log(p) * log(log(p))))
    B = int(pow(L, 1 / sqrt(2)))

    primes = PrimesLT(B)

    # currently brute forces solutions to log_g_x with x for each prime <= B
    logs = {}
    for n in primes:
        logs[n] = baby_step_giant_step(g, n, p)

    k = 0
    while True:
        k += 1
        x = (h * pow(g, -k, p)) % p
        if BSmoothQ(x, B):
            exponents = factor_base_exp(x, primes)
            # reduce sums list returned by map, which returns list of products of exponents and logs
            return reduce(lambda a, b: a + b, list(map(lambda i, n: i * logs[n], exponents, logs))) + k


def SolveDLP(g, h, p, q=None):
    """Uses Pollard's Rho algorithm for logarithms to solve given discrete log problem."""

    xstate = (1, 0, 0)
    ystate = (1, 0, 0)

    if q is None:
        q = p - 1

    while True:
        xstate = calculate_state(xstate, g, h, p)
        ystate = calculate_state(calculate_state(ystate, g, h, p), g, h, p)

        if xstate[0] == ystate[0]:
            try:

                # try to return result right away, fails if beta value is not invertible mod q
                result = (xstate[1] - ystate[1]) * pow(ystate[2] - xstate[2], -1, q)
            except ValueError:

                # try to reduce entire equation by gcd, then try to invert again
                alpha = xstate[1] - ystate[1]
                beta = ystate[2] - xstate[2]
                e = gcd(alpha, beta, q)
                if e > 1:
                    alpha //= e
                    beta //= e
                    mod = q // e

                    # if reduced until invertible, find solution
                    if gcd(beta, mod) == 1:
                        log_g_h = alpha * pow(beta, -1, mod)

                        # current solution is mod q // e, but real solution could be mod any increment of
                        # q up until its original value, find real solution by checking
                        for i in range(e):
                            if pow(g, log_g_h, p) == h % p:
                                return log_g_h
                            log_g_h += mod
                continue
            else:
                if pow(g, result, p) == h % p:
                    return result


def calculate_state(state, g, h, p):
    """This function is a helper method for Pollard's Rho algorithm for solving DLP's."""

    x, alpha, beta = state[0], state[1], state[2]

    if 0 <= x < p//3:
        x *= g
        if x >= p:
            x %= p
        alpha += 1
        if alpha >= p - 1:
            alpha %= (p - 1)
        return x, alpha, beta
    elif p//3 <= x < 2 * p//3:
        x = pow(x, 2, p)
        alpha *= 2
        beta *= 2
        if x >= p:
            x %= p
        if alpha >= p - 1:
            alpha %= (p - 1)
        if beta >= p - 1:
            beta %= (p - 1)
        return x, alpha, beta
    elif 2 * p//3 <= x < p:
        x *= h
        if x >= p:
            x %= p
        beta += 1
        if beta >= p - 1:
            beta %= (p - 1)
        return x, alpha, beta


def pohlig_hellman(g, h, p, q, exp, prog=False):
    """Function attempts to solve DLP mod p where order of g, q, is some prime raised to a power."""

    X = []
    if prog:
        print("Starting process for X0")

    r = pow(g, pow(q, exp - 1), p)
    X0 = baby_step_giant_step(r, pow(h, pow(q, exp - 1), p), p, prog, q)
    X.append(X0)

    if prog:
        print(f"Found X0 = {X0}\n")

    for i in range(1, exp):
        if prog:
            print(f"Starting process for X{i}")
        exp_term = fromBase(X[::-1], q)
        h_term = pow(h * pow(pow(g, exp_term, p), -1, p), pow(q, exp - i - 1), p)
        Xi = baby_step_giant_step(r, h_term, p, prog, q)
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


def pollard_p1(n, limit=pow(10, 6), first_n=4):
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
        k = pollard_p1(n)
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
                if p * q == N:
                    return {p: 1, q: 1}
                return _factorWithKnown(p, q, N)
            if 1 < p < N:
                q = N // p
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


def QuadraticSieve1(N, B=None):
    """Performs Quadratic Sieve Algorithm with a given Smoothness value B on a given composite N"""

    from itertools import combinations_with_replacement as _all
    from math import e, log, sqrt

    if B is None:
        L = pow(e, sqrt(log(N) * log(log(N))))
        B = int(pow(L, 1 / sqrt(2)))

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


def QuadraticSieve(N, B=None):
    """Performs Quadratic Sieve Algorithm with a given Smoothness value B on a given composite N"""

    from itertools import combinations as _all
    from math import e, log, sqrt

    if B is None:
        L = pow(e, sqrt(log(N) * log(log(N))))
        B = int(pow(L, 1 / sqrt(2)))

    m = PrimePi(B)
    sq = isqrt(N)
    while pow(sq, 2) < N:
        sq += 1

    b_smooth_nums = []
    squared_nums = {}
    factor_base = PrimesLT(B)
    nums_found = 0
    while nums_found < m + 1:
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
            nums_found += 1

    choices = list(map(list, list(_all(b_smooth_nums, 2))))
    for c in range(2, len(b_smooth_nums)):
        for choice in choices:
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

        choices = combinations_cumulative(choices, b_smooth_nums)

    return {}


def combinations_cumulative(combinations, source):
    """Function takes in set of combinations of k-choices without replacement and returns a new
    set of combinations of k + 1 choices without replacement."""

    new_combination = []
    for choice in combinations:
        for e in source:
            if e not in choice:
                temp = choice[:]
                temp.append(e)
                add = True
                for c in new_combination:
                    different = False
                    i = 0
                    while i < len(temp) and not different:
                        t = temp[i]
                        i += 1
                        if t not in c:
                            different = True
                            continue
                    if not different:
                        add = False
                        break
                if add:
                    new_combination.append(temp)

    return new_combination
