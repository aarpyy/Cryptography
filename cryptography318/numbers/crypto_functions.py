from random import randrange
from math import gcd, isqrt, sqrt, prod

from .prime import primes_gen
from cryptography318.linalg.array import ArrayMod, Array


# classes and methods for working with elliptic curve cryptography
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
        """Creates Elliptic object that is point on instance EllipticCurve. If neither x or y are given,
        the identity element is returned. If just x is given, the closest z value >= x is found s.t.
        z^3 + az + b has a square root mod p, and the point is returned with x value = z, y value =
        square root of z mod p. If both are given, the point is verified to be on the curve then returned.
        If just y is given, or point given is not on the curve, error is thrown."""

        # identity element
        if x is None and y is None:
            return Elliptic(self)

        if y is None:
            x_term = pow(x, 3) + self.a * x + self.b
            while not quadratic_residue(x_term, self.mod):
                x += 1
                x_term = pow(x, 3) + self.a * x + self.b
            y = find_roots(x_term, self.mod)

        # checks if point given is on curve
        if pow(y, 2) % self.mod == (pow(x, 3) + self.a * x + self.b) % self.mod:
            return Elliptic(self, x, y)

        # if neither identity or point on curve, this point doesn't exist
        raise ValueError("Argument values incompatible with this curve")

    def on_curve(self, point):
        return pow(point[1], 2, self.mod) == (pow(point[0], 3) + self.a * point[0] + self.b) % self.mod

    def string(self, s):
        return string_to_elliptic(self, s)


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
        if (x1 - x2) % self.E._mod == 0 and (y1 + y2) % self.E._mod == 0:
            return Elliptic(self.E)

        # two try/catch clauses allow for information to be obtained if crash occurs, mostly useful in
        # lenstra's elliptic curve factorization algorithm
        if (x1 - x2) % self.E._mod == 0 and (y1 - y2) % self.E._mod == 0:
            try:
                slope = ((3 * pow(x1, 2) + self.E.a) * pow(2 * y1, -1, self.E._mod)) % self.E._mod
            except ValueError as e:
                raise ValueError(str(e) + f" base: {2 * y1}")
        else:
            try:
                slope = ((y2 - y1) * pow(x2 - x1, -1, self.E._mod)) % self.E._mod
            except ValueError as e:
                raise ValueError(str(e) + f" base: {x2 - x1}")

        x3 = (pow(slope, 2) - x1 - x2) % self.E._mod
        y3 = (slope * (x1 - x3) - y1) % self.E._mod
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
        return elliptic_to_string(self)

    def copy(self):
        return Elliptic(self.E, self.point[0], self.point[1])


def string_to_elliptic(curve, s):
    for i in range(100):
        x = string_to_int(s) * 100 + i
        x_term = pow(x, 3, curve._mod) + curve.a * x + curve.b
        if quadratic_residue(x_term, curve._mod):
            y = find_roots(x_term, curve._mod)
            return curve.point(x, y)
    return None


def elliptic_to_string(point):
    return int_to_string(point[0] // 100)


# basic encryption systems
def vigenere_encrypt(key, plaintext):
    key = key.lower()
    plaintext = plaintext.lower()
    full_key = [plaintext[i - len(key)] if i > len(key) - 1 else key[i] for i in range(len(plaintext))]

    int_key = ArrayMod(list(map(lambda c: ord(c) - 97, full_key)), 26)
    int_text = ArrayMod(list(map(lambda c: ord(c) - 97, plaintext)), 26)
    return ''.join(map(lambda c: chr(c + 97), int_key + int_text))


def caeser_shift(plaintext):
    frequencies = [.082, .015, .028, .043, .127, .022, .020, .061, .070, .002, .008, .040, .024, .067, .075,
                   .019, .001, .060, .063, .091, .028, .010, .023, .001, .020, .001]
    plaintext = plaintext.lower()
    plain_freq = Array([0] * 26)
    for c in plaintext:
        index = ord(c) - 97
        plain_freq[index] += 1
    plain_freq /= len(plaintext)

    dot = lambda a, b: sum(map(lambda e1, e2: e1 * e2, a, b))

    products = []
    for k in range(26):
        products.append(dot(plain_freq, frequencies))
        plain_freq.shift_elements(shift=1)
    products = abs(Array(products) - dot(frequencies, frequencies))
    index = products.index(min(products)) + 1
    return ''.join(map(lambda c: chr(((ord(c) - 97) + index) % 26 + 96), plaintext))


# helper methods for working with cryptographic functions
def to_base(n, base):
    # finds highest exponent of base
    exp = 0
    limit = base
    while n > limit:
        exp += 1
        limit *= base

    # iterates starting at highest exponent down to 0, creates list of 'base' base
    list_coeff = [0] * (exp + 1)
    for i in range(exp, -1, -1):
        list_coeff[i], n = divmod(n, pow(base, i))
    return list_coeff  # returns list starting at lowest base with increasing powers


def from_base(lst, base):
    return sum(map(lambda i, n: n * pow(base, i), range(len(lst)), lst))


def string_to_int(s, base=128):

    # string s is reversed on purpose, first char is highest power base
    return sum(map(lambda i, c: ord(c) * pow(base, i), range(len(s)), s[::-1]))


def int_to_string(n, base=128):
    return ''.join(map(lambda c: chr(c), to_base(n, base)[::-1]))


def sqrt_safe(n):
    try:
        return sqrt(n)
    except OverflowError:
        pass
    finally:
        return isqrt(n)


def extended_gcd(*args):

    # actual extended gcd function
    def ext_gcd(a, b):
        if a == 0:
            return b, 0, 1

        g, x, y = ext_gcd(b % a, a)
        return g, y - (b // a) * x, x

    # if just two arguments, return normal extended gcd
    if len(args) == 2:
        return ext_gcd(args[0], args[1])

    g, u, v = ext_gcd(args[0], args[1])  # gcd of first two args; u, v s.t. args sum to gcd
    values = Array([u, v])
    for e in args[2:]:
        g, u1, v1 = ext_gcd(g, e)
        values *= u1  # u1 s.t. u1 * gcd prev two values + v1 * curr val = gcd, so u1 applied to prev values
        values.append(v1)
    return g, *values


# methods for working within finite field of integers
def find_roots(a, p):
    """Finds a solution for x to equation x^2 = a (mod p). If a solution is returned, a second
    solution s2 will also exist where s2 = -x (mod p)."""

    if not quadratic_residue(a, p):
        return None

    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    q = p - 1
    s = 0
    while not q & 1:
        q >>= 1
        s += 1

    z = randrange(2, p)
    while not quadratic_non_residue(z, p):
        z = randrange(2, p)

    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)

    while True:
        if t == 0:
            return 0
        if t == 1:
            return r

        i = 0
        x = t
        while x != 1:
            x = pow(x, 2, p)
            i += 1

            if i == m:
                return None

        b = pow(c, pow(2, m - i - 1), p)
        c = pow(b, 2, p)
        m = i

        t *= c
        if t >= p:
            t %= p
        r *= b
        if r >= p:
            r %= p


def quadratic_residue(a, p):
    """Returns True if n is a quadratic residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == 1


def quadratic_non_residue(a, p):
    """Returns True if n is a quadratic non-residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == p - 1


def chinese_remainder(values, moduli):
    # initializes lists of moduli, mod = product of all moduli
    mod = prod(moduli)

    # maps list of moduli and their inverses to x and y respectively
    x, y = [], []
    for m in moduli:
        mi = mod // m
        x.append(mi)
        y.append(pow(mi, -1, m))

    # accumulates product of numbers and moduli and their inverses
    acc = 0
    for i in range(len(values)):
        acc = (acc + values[i] * x[i] * y[i]) % mod

    return acc


# methods for solving dlp's
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


def elliptic_bsgs(point, result, order=None):
    if order is None:
        order = point.E._mod + 1 + 2 * isqrt(point.E._mod)

    n = isqrt(order)

    # creates baby-step table, P * i for i from 0 to n
    a = list(map(lambda x: point * x, range(n)))
    # creates giant-step table, Q - P * jn for j from 0 to n
    b = list(map(lambda x: result + (point * (-x * n)), range(n)))

    # creates sets out of lists
    set_a, set_b = set(a), set(b)

    matches = set_a.intersection(set_b)

    if not matches:
        return None
    point = matches.pop()

    # uses lists to find index of intersection
    i, j = a.index(point), b.index(point)
    return i + j * n


def pollard_rho_dlp(g, h, p, order=None):
    """Uses Pollard's Rho algorithm for logarithms to solve given discrete log problem. Function will run
    indefinitely until a solution is found.

    :param g: integer base
    :param h: integer solution to g^x for some x
    :param p: integer prime modulus
    :param order: integer order of g (smallest integer s.t. pow(g, order, p) == 1
    :return: solution to g^x = h for integer x"""

    xstate = (1, 0, 0)
    ystate = (1, 0, 0)

    if order is None:
        order = p - 1

    while True:
        xstate = calculate_state(xstate, g, h, p)
        ystate = calculate_state(calculate_state(ystate, g, h, p), g, h, p)

        if xstate[0] == ystate[0]:
            try:

                # try to return result right away, fails if beta value is not invertible mod q
                return (xstate[1] - ystate[1]) * pow(ystate[2] - xstate[2], -1, order)
            except ValueError:

                # try to reduce entire equation by gcd, then try to invert again
                alpha = xstate[1] - ystate[1]
                beta = ystate[2] - xstate[2]
                e = gcd(alpha, beta, order)
                if e > 1:
                    alpha //= e
                    beta //= e
                    mod = order // e

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


def index_calculus_dlp(g, h, p):
    """Function attempts to solve DLP through index calculus algorithm, computing a series of smaller dlp's
    used to solve the larger."""

    from math import e, sqrt, log
    from functools import reduce

    L = pow(e, sqrt(log(p) * log(log(p))))
    B = int(pow(L, 1 / sqrt(2)))

    primes = primes(B)

    # currently brute forces solutions to log_g_x with x for each prime <= B
    logs = {}
    for n in primes:
        logs[n] = baby_step_giant_step(g, n, p)

    k = 0
    while True:
        k += 1
        x = (h * pow(g, -k, p)) % p
        if b_smooth(x, B):
            exponents = factor_with_base(x, primes)
            # reduce sums list returned by map, which returns list of products of exponents and logs
            return reduce(lambda a, b: a + b, list(map(lambda i, n: i * logs[n], exponents, logs))) + k


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
    X0 = baby_step_giant_step(r, pow(h, pow(q, exp - 1), p), p, q)
    X.append(X0)

    if prog:
        print(f"Found X0 = {X0}\n")

    for i in range(1, exp):
        if prog:
            print(f"Starting process for X{i}")
        exp_term = from_base(X[::-1], q)
        h_term = pow(h * pow(pow(g, exp_term, p), -1, p), pow(q, exp - i - 1), p)
        Xi = baby_step_giant_step(r, h_term, p, q)
        X.append(Xi)
        if prog:
            print(f"Found X{i} = {Xi}\n")

    return from_base(X[::-1], q)


def lenstra_elliptic(n, limit=pow(10, 7)):
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
        q = gcd(k, n)
        return q

    return None


def b_smooth(n, b=None, factors=None):
    """Returns True if all prime factors of given number are <= B"""

    if factors is None:
        factors = primes_gen(b)

    for p in factors:
        while n % p == 0:
            n //= p

    return n == 1
