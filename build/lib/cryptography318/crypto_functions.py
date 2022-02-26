from random import randrange
from math import gcd, isqrt, sqrt, prod
from .prime import IsPrime, NextPrime, PrimesLT, PrimePi
from .tools import deprecated, join_dict
from .quadratic_sieve import quadratic_sieve, _factor_with_known


class EllipticCurve:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.mod = p

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f'EllipticCurve({self.a}, {self.b}, {self.mod})'

    @classmethod
    def safe_curve(cls, p):
        """Constructs a random curve mod p that is not singular."""

        while True:
            a = randrange(p)
            b = randrange(p)
            if 4 * pow(a, 3) + 27 * pow(b, 2) != 0:
                return cls(a, b, p)

    @classmethod
    def safe_curve_and_point(cls, p):
        """Constructs a random curve mod p that is not singular, and a point that exists
        on that curve."""

        while True:
            a = randrange(p)
            x = randrange(p)
            y = randrange(p)
            b = (pow(y, 2, p) - pow(x, 3, p) - a * x) % p
            if 4 * pow(a, 3) + 27 * pow(b, 2) == 0:
                continue
            curve = cls(a, b, p)
            return curve, curve.point(x, y)

    def point(self, x=None, y=None):

        # identity element
        if x is None and y is None:
            return Elliptic(self)

        # checks if point given is on curve
        if pow(y, 2) % self.mod == (pow(x, 3) + self.a * x + self.b) % self.mod:
            return Elliptic(self, x, y)

        # if neither identity or point on curve, this point doesn't exist
        raise ValueError("Argument values incompatible with this curve")

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

        # if points add to identity, return identity (Elliptic with point = None)
        if (x1 - x2) % self.E.mod == 0 and (y1 + y2) % self.E.mod == 0:
            return Elliptic(self.E)

        # two try/catch clauses allow for information to be obtained if crash occurs, mostly useful in
        # lenstra's elliptic curve factorization algorithm
        if (x1 - x2) % self.E.mod == 0 and (y1 - y2) % self.E.mod == 0:
            try:
                slope = ((3 * pow(x1, 2) + self.E.a) * pow(2 * y1, -1, self.E.mod)) % self.E.mod
            except ValueError as e:
                raise ValueError(str(e) + f" base: {2 * y1}")
        else:
            try:
                slope = ((y2 - y1) * pow(x2 - x1, -1, self.E.mod)) % self.E.mod
            except ValueError as e:
                raise ValueError(str(e) + f" base: {x2 - x1}")

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

        # if point is identity, return point
        if self.point is None:
            return self

        # only work with integers in field
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

    def copy(self):
        return Elliptic(self.E, self.point[0], self.point[1])


def sqrt_safe(n):
    try:
        return sqrt(n)
    except OverflowError:
        pass
    finally:
        return isqrt(n)


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
        x_term = pow(x, 3, E.mod) + E.a * x + E.b
        if quadratic_residue(x_term, E.mod):
            y = 1
            while pow(y, 2, E.mod) != x_term:
                y += 1
            return E.point(x, y)
    return None


def EllipticToString(point):
    return NumToString(point[0] // 100)


def ExtendedGCD(a, b):
    if a == 0:
        return b, 0, 1

    g, x, y = ExtendedGCD(b % a, a)
    return g, y - (b // a) * x, x


def quadratic_residue(a, p):
    """Returns True if n is a quadratic residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == 1


def BSmoothQ(n, B=None, factors=None):
    """Returns True if all prime factors of given number are <= given B"""

    if factors is None:
        factors = [p := 2]
        while len(factors) < PrimePi(B):
            factors.append(p := NextPrime(p))

    for p in factors:
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


def ChineseRemainder(nums, mods):
    # initializes lists of moduli, M = product of all moduli
    M = prod(mods)

    # maps list of moduli and their inverses to x and y respectively
    x, y = [], []
    for m in mods:
        mi = M // m
        x.append(mi)
        y.append(pow(mi, -1, m))

    # accumulates product of number and moduli and their inverses
    acc = 0
    for i in range(len(nums)):
        acc = (acc + nums[i] * x[i] * y[i]) % M

    return acc


def baby_step_giant_step(g, h, p, order=None):
    """Function attempts to solve DLP using classic baby-step-giant-step algorithm."""
    if order is None:
        order = p - 1

    n = isqrt(order) + 1

    # find lists A = g^i B = h * g^-jn s.t. A[i] == B[j] for some indices i, j, this collision allows us to solve
    A = list(map(lambda e: pow(g, e, p), range(n)))
    B = list(map(lambda e: (h * pow(g, -e * n, p)) % p, range(n)))

    # convert to set for intersection calculation
    a, b = set(A), set(B)
    U = a.intersection(b)

    # if empty set, no collisions found
    if not U:
        return None

    # otherwise, find first indices of match and use to solve
    match = U.pop()
    i, j = A.index(match), B.index(match)
    return i + (j * n)


def elliptic_bsgs(P, Q, order=None):
    if order is None:
        order = P.E.mod + 1 + 2 * isqrt(P.E.mod)

    n = isqrt(order)

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
    i, j = a.index(point), b.index(point)
    return i + j * n


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


def pollard_rho_dlp(g, h, p, q=None):
    """Uses Pollard's Rho algorithm for logarithms to solve given discrete log problem. Function will run
    indefinitely until a solution is found.

    :param g: integer base
    :param h: integer solution to g^x for some x
    :param p: integer prime modulus
    :param q: integer order of g (smallest integer s.t. pow(g, q, p) == 1
    :return: solution to g^x = h for integer x"""

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
                return (xstate[1] - ystate[1]) * pow(ystate[2] - xstate[2], -1, q)
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
                        for i in range(e - 1):
                            if pow(g, log_g_h, p) == h % p:
                                return log_g_h
                            log_g_h += mod
                continue


def calculate_state(state, g, h, p):
    """This function is a helper method for Pollard's Rho algorithm for solving DLP's."""

    x, alpha, beta = state[0], state[1], state[2]

    if 0 <= x < p // 3:
        x *= g
        if x >= p:
            x %= p
        alpha += 1
        if alpha >= p - 1:
            alpha %= (p - 1)
        return x, alpha, beta
    elif p // 3 <= x < 2 * p // 3:
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
    elif 2 * p // 3 <= x < p:
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


def lenstra_elliptic(n, limit=pow(10, 6)):
    """Performs Lenstra's Elliptic Curve Factorization algorithm on integer n."""

    # choose random point and elliptic curve, unless curve meets conditions not suitable, use it
    while True:
        a = randrange(n)
        x = randrange(n)
        y = randrange(n)
        b = (pow(y, 2, n) - pow(x, 3, n) - a * x) % n
        if 4 * pow(a, 3) + 27 * pow(b, 2) == 0:
            continue
        break

    E = EllipticCurve(a, b, n)
    P = E.point(x, y)

    try:
        for j in range(1, limit):
            P *= j

    # if there was a crash, it is because pow(x, -1, n) not invertible, so gcd(x, n) is a factor
    except ValueError as e:

        # gets the number not invertible from the value error thrown
        k = int(str(e).split('base: ')[1])
        return gcd(k, n)

    return None


def FactorInt(n):
    """
    Attempts to factor given integer with four methods, returning None if un-factorable.
    Function first checks if number is prime, then iterates through all primes < 1000 attempting
    to divide n. Function then Lenstra's Elliptic Curve factorization algorithm to try finding small
    factors, then tries Pollard's P-1 algorithm to find one non-trivial factor of n.
    If it succeeds, adds factor to solution set. If n is still factorable, tries quadratic sieve method
    to return all remaining factors. If this returns None, uses sympy's factorint() method and returns result.

    :param n: int number to be factored
    :return: dictionary of all primes factors and their powers, or None if not factorable
    """

    if n == 1:
        return {}

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
            if p not in factors:
                factors[p] = 1
            else:
                factors[p] += 1
            n //= p

    if n == 1:
        return factors

    if IsPrime(n):
        return join_dict(factors, {n: 1})

    k = pollard_p1(n)
    if k is not None:
        factors_k = {k: 1} if IsPrime(k) else FactorInt(k)
        factors = join_dict(factors, factors_k)
        n //= k

    k = lenstra_elliptic(n)
    if k is not None:
        factors_k = {k: 1} if IsPrime(k) else FactorInt(k)
        factors = join_dict(factors, factors_k)
        n //= k

    factors_qs = quadratic_sieve(n)
    if factors_qs is None:
        return join_dict(factors, {n: 1})

    return join_dict(factors, factors_qs)


def pollard_p1(n, limit=pow(10, 6)):
    """Pollard's p - 1 algorithm for factoring large composites.
    Returns one non-trivial factor if factor-able, False if otherwise."""

    if IsPrime(n):
        return n

    for a in [2, 3]:
        m = a
        for j in range(2, limit):
            m = pow(m, j, n)
            k = gcd(m - 1, n)
            if 1 < k < n:
                return k

    return None


def combinations_cumulative(found, source):
    """Function takes in set of combinations of k-choices without replacement and returns a new
    set of combinations of k + 1 choices without replacement.

    :param found: list of already found combinations, cannot be empty
    :param source: list of source for combinations"""

    new_combination = []
    for choice in found:
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


def gen_choice(indices, matrix):
    """Generator for choosing elements from matrix."""

    for i in indices:
        yield matrix[i]


@deprecated
def __factor_perfect_square(n, B=7):
    """Attempts a similar attack to quadratic sieve, finding B-smooth perfect squares in an attempt
    to find a multiple of n. Function written for learning purposes, if trying to factor integer, using
    FactorInt."""

    from itertools import combinations_with_replacement as _all

    m = PrimePi(B)
    b_smooth_nums = []
    squared_nums = {}
    a = isqrt(n) - 1
    while len(b_smooth_nums) < m:
        ci = pow(a, 2, n)
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

            p, q = gcd(a - b, n), gcd(a + b, n)
            if 1 < p < n or 1 < q < n:
                return _factor_with_known(p, q, n)

    return False
