from math import floor, gcd, isqrt, log, sqrt
from random import Random

from cryptography318.prime.primesieve import primesieve
from cryptography318.utils.misc import as_int


rndm = Random()


class Identity:
    """ Simple additive identity element to enable operations on a Weierstrass curve """

    def __eq__(self, other):
        return isinstance(other, Identity)

    def __neg__(self):
        return self

    def __add__(self, other):
        return other

    def __mul__(self, other):
        return self


class WeierstrassPoint:
    __slots__ = 'x', 'n', 'y', 'a'

    infinity = Identity()

    def __new__(cls, x, y, a, n):
        self = super().__new__(cls)
        self.x = x
        self.y = y
        self.n = n
        self.a = a
        return self

    def __eq__(self, other):
        if other is WeierstrassPoint.infinity:
            return False
        return (self.x - other.x) % self.n == 0 and (self.y - other.y) % self.n == 0

    def __neg__(self):
        return WeierstrassPoint(self.x, -self.y % self.n, self.a, self.n)

    def __add__(self, other):

        # If identity element, return itself
        if other is WeierstrassPoint.infinity:
            return self

        if (self.x - other.x) % self.n == 0 and (self.y + other.y) % self.n == 0:
            return WeierstrassPoint.infinity

        try:

            # If they are the same point, use the doubling algorithm
            if self == other:
                slope = ((3 * self.x * self.x + self.a) * pow(2 * self.y, -1, self.n)) % self.n
            else:
                slope = ((other.y - self.y) * pow(other.x - self.x, -1, self.n)) % self.n

        # This except is the basis of lenstra's ecm. raising an error with the specific non-invertible value
        # allows for gcd to be taken to find factor of modulus
        except ValueError:

            # Raise ValueError with message as the value that does not have an inverse
            raise ValueError(str((2 * self.y) if self == other else (other.x - self.x)))

        x3 = (slope * slope - self.x - other.x) % self.n
        y3 = (slope * (self.x - x3) - self.y) % self.n
        return WeierstrassPoint(x3, y3, self.a, self.n)

    def __mul__(self, other):
        P = self
        if other < 0:
            other = -other
            P = -P

        Q = P
        R = WeierstrassPoint.infinity
        while other > 0:
            if other & 1:
                R += Q
            Q += Q
            other >>= 1
        return R


def gen_weierstrass_point(n):
    while True:
        a = rndm.randrange(n)
        x = rndm.randrange(n)
        y = rndm.randrange(n)

        # generate b based on point
        b = (pow(y, 2, n) - pow(x, 3, n) - a * x) % n

        # if singular curve, try again
        if 4 * pow(a, 3) + 27 * pow(b, 2) != 0:
            return WeierstrassPoint(x, y, a, n)


class MontgomeryPoint:
    """
    Point on Elliptic curve in Montgomery form over finite field with field operations '+' and '*'.
    For detailed explanation of operations refer to Gaj.pdf.
    """

    __slots__ = 'z', 'x', 'n', 'a_24'

    def __new__(cls, x, z, a_24, n):
        self = super().__new__(cls)
        self.x = x
        self.z = z
        self.a_24 = a_24
        self.n = n
        return self

    def __eq__(self, other):
        return (self.x - other.x) % self.n == 0 and (self.z - other.z) % self.n == 0

    def add(self, Q, difference):
        """
        :param self: point 1
        :param Q: point 2
        :param difference: P - Q
        Cost:
            *: 6;
            +: 6;
        """
        u, v = (self.x - self.z) * (Q.x + Q.z), (self.x + self.z) * (Q.x - Q.z)
        add, sub = u + v, u - v

        x3 = (difference.z * add * add) % self.n
        z3 = (difference.x * sub * sub) % self.n
        return MontgomeryPoint(x3, z3, self.a_24, self.n)

    def double(self):
        """
        Computes addition P + Q on E over FF(m) when P == Q

        Cost:
            *: 5;
            +: 4;

        Refer to Gaj p.4, algorithm 3
        """
        u, v = self.x - self.z, self.x + self.z
        u, v = u * u, v * v
        _4xz = v - u
        x2 = (v * u) % self.n
        z2 = (_4xz * (u + self.a_24 * _4xz)) % self.n
        return MontgomeryPoint(x2, z2, self.a_24, self.n)

    def ladder(self, k):
        """
        Refer to Gaj p.4, algorithm 2
        """

        # Here, the difference between P and Q is self, since P = 2 * Q, so we reuse self as the difference
        P = self.double()
        Q = self
        for bit in bin(k)[3:]:

            # compare via strings instead of convert to int is faster
            if bit == '1':
                Q = Q.add(P, self)
                P = P.double()
            else:
                P = P.add(Q, self)
                Q = Q.double()
        return Q


class MontgomeryPointZ1(MontgomeryPoint):

    def add(self, Q, difference):
        """
        Computes addition P + Q on E over FF(m) when z_P == 1 (z coordinate of initial point)

        Refer to Gaj p.4, algorithm 3. For explanation of z coordinate == 1, refer to just below
        algorithm 3 p.4

        Cost:
            *: 5;
            +: 6;
        """
        u, v = (self.x - self.z) * (Q.x + Q.z), (self.x + self.z) * (Q.x - Q.z)
        add, sub = u + v, u - v
        x3 = (add * add) % self.n   # Here we skip a multiply
        z3 = (difference.x * sub * sub) % self.n
        return MontgomeryPoint(x3, z3, self.a_24, self.n)


def make_point(x, z, a_24, n):
    if z == 1:
        return MontgomeryPointZ1(x, z, a_24, n)
    else:
        return MontgomeryPoint(x, z, a_24, n)


def ecm_weierstrass(N, B=None, retry=50, verbose=False):
    """
    Lenstra's Elliptic Curve Factorization method over a short Weierstrass curve
    with parameters a, b s.t. 4a^3 + 27b^2 != 0
    """

    if B is None:
        from math import e
        L = pow(e, isqrt(int(log(N) * log(log(N)))))
        B = int(pow(L, 1 / pow(2, .5)))

    if verbose:
        print(f"Weierstrass curve factorization with B = {B}")

    for _ in range(retry):

        P = gen_weierstrass_point(N)

        if verbose:
            b = (pow(P.y, 2, N) - pow(P.x, 3, N) - P.a * P.x) % N
            print(f"Trying curve y^2 = x^3 + {P.a}x + {b} (mod {N})")

        try:
            # Keep multiplying until we perform an operation that is not invertible
            for j in range(1, B):
                P *= j
        except ValueError as exc:
            # If an error was thrown, a number was inverted that is not invertible, take the gcd
            try:
                k = int(str(exc).strip())
            except Exception:
                raise exc

            g = gcd(k, N)
            if 1 < g < N:
                if verbose:
                    print(f"Found factor {g} from gcd of {k} and {N}")
                return g

    return None


def ecm_mont_phase1(N, use_z1=False):
    """
    Helper function for elliptic curve factorization over a Montgomery curve.

    :return: initial point for elliptic curve factorization or non-trivial factor of N
    """

    sigma = rndm.randint(6, N - 1)
    u = pow(sigma, 2, N) - 5 % N
    u_3 = pow(u, 3, N)
    v = 4 * sigma % N
    try:
        _4u3v_inv = pow(4 * u_3 * v, -1, N)
        C = pow(v - u, 3, N) * (3 * u + v) * _4u3v_inv - 2 % N
    except ValueError:
        # If there is no mod-inverse for 4 * u_3 * v (mod N) then they share a non-trivial factor
        return gcd(4 * u_3 * v, N)

    a24 = (C + 2) * pow(4, -1, N) % N  # a24 = ((v-u)^3 * (3uv) / u^3v

    # Initial point is [u^3 : v^3] but we can divide by v^3 if 'use_z1' is True
    if use_z1:
        return make_point(pow(4 * u_3 * u * _4u3v_inv, 3, N), 1, a24, N)
    else:
        return make_point(u_3, pow(v, 3, N), a24, N)


def ecm_mont(N, B1=None, B2=None, retry=50, verbose=False, use_z1=False):
    """
    Performs Lenstra's elliptic curve factorization method with elliptic curve E
    in Montgomery form over FF(N) with Suyama's parameterization.

    References
    ----------
    [1] Gaj K. et al. Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware
    https://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
    [2]
    https://cloudflare-ipfs.com/ipfs/bafykbzacedds3ejoo5mlrx6zh2npeyxu7vqqrrrq5mzrjjwajlgm5zbbies4u

    :param N: number to be factored
    :param B1: limit for phase 1; start for phase 2
    :param B2: limit for phase 2
    :param retry: number of trials to be run
    :param verbose: print information about each trial
    :param use_z1: use z = 1 for initial point to reduce cost of addition
    :return: one non-trivial factor of N or None
    """

    if B1 is None:
        from math import e
        g = isqrt(N)
        B1 = int(pow(e, sqrt(log(g) * log(log(g)) * 0.5)))
    elif B1 & 1 == 1:
        B1 += 1

    if B2 is None:
        B2 = B1 * 100
    elif B2 & 1 == 1:
        B2 += 1

    primesieve.extend(B2)

    if verbose:
        print(f"Montgomery curve factorization with B1 = {B1}, B2 = {B2}, and {retry} curves")

    # with initial point P, get k to find initial Q = kP
    k = 1
    for p in primesieve.primerange(B1 + 1):
        k *= pow(p, floor(log(B1, p)))

    D = isqrt(B2)
    S: list[int | MontgomeryPoint] = [0] * D
    beta = [0] * D

    for trial in range(retry):

        P = ecm_mont_phase1(N, use_z1=use_z1)
        if isinstance(P, int):
            return P

        Q = P.ladder(k)

        # If found factor already, return, otherwise continue to phase 2
        g = gcd(Q.z, N)
        if 1 < g < N:
            if verbose:
                print(f"Found factor {g} from gcd of {Q.z} and {N}")
            return g

        if g == N:
            continue

        """ Begin Phase 2 """

        S[0] = Q.double()  # S_1 = double(Q)
        S[1] = S[0].double()  # S_2 = double(S_1)
        beta[0] = (S[0].x * S[0].z) % N
        beta[1] = (S[1].x * S[1].z) % N
        for d in range(2, D):
            S[d] = S[d - 1].add(S[0], S[d - 2])
            beta[d] = (S[d].x * S[d].z) % N

        D2 = 2 * D
        g = 1
        B = B1 - 1
        T = Q.ladder(B - D2)
        R = Q.ladder(B)

        for r in range(B, B2, D2):
            alpha = (R.x * R.z) % N
            for q in primesieve.primerange(r + 2, D2 + 1):  # Loop over primes
                delta = ((q - r) // 2) - 1  # Distance to next prime - minus one bc 0-based index

                # g = g((X(R) − X(Sδ))(Z(R) + Z(Sδ)) − α + βδ) mod n;
                g = (g * ((R.x - S[delta].x) * (R.z + S[delta].z)) - alpha + beta[delta]) % N

            R, T = R.add(S[D - 1], T), R

        g = gcd(g, N)
        if 1 < g < N:
            if verbose:
                print(f"Found factor {g} from gcd of {g} and {N}")
            return g

    return None


def ecm_mont_basic(N, B1=None, B2=None, retry=50):
    """
    A more basic version of the below elliptic factorization over a Montgomery form elliptic curve.
    This method opts for a direct multiplication approach for phase 2. Phase 1 is computed
    identically to find Q = kP. If kP missed the order of E(FF(q)), for some q factor of N, by a small
    integer, then we compute some multiple of kP to find the factor. In this phase 2, the integer
    multiplied against kP is the z coord of the product of Q with each prime p, B1 < p < B2. This results
    in a sometimes quicker calculation if N is small enough, but in general the below function should
    be used for elliptic curve factorization.
    """
    if B1 is None:
        from math import e
        B1 = (int(pow(e, isqrt(int((log(N) * log(log(N))) / 2)))) | 1) + 1  # | 1 + 1 ensures even

    if B2 is None:
        B2 = B1 * 100

    primesieve.extend(B2)

    # With initial point P, get k to find initial Q = kP
    k = 1
    for p in primesieve.primerange(B1 + 1):
        k *= floor(log(B1, p))

    for _ in range(retry):

        P = ecm_mont_phase1(N)
        if isinstance(P, int):
            return P

        Q = P.ladder(k)

        # If found factor already, return, otherwise continue to phase 2
        q = gcd(Q.z, N)
        if 1 < q < N:
            return q

        if q == N:
            continue

        """ Begin Phase 2 """

        d = 1
        for p in primesieve.primerange(B1, B2):
            d = (d * Q.ladder(p).z) % N

        q = gcd(d, N)
        if 1 < q < N:
            return q

    return None


def ecm(N, B1=None, B2=None, retry=50, verbose=False, use_z1=False, use_weierstrass=False, use_basic=False):
    """
    General purpose function to compute Lenstra's elliptic curve factorization method
    on a composite integer N. If number is a 32-bit integer, use a Weierstrass curve
    to find a factor. If this succeeds, return the factor. If this fails, or number is
    greater than a 32-bit integer, use a Montgomery curve to find a factor.

    :param N: number to be factored
    :param B1: limit for phase 1; start for phase 2
    :param B2: limit for phase 2
    :param retry: number of trials to be run
    :param verbose: print information about each trial
    :param use_z1: use z = 1 for initial point to reduce cost of addition
    :param use_weierstrass: use a Weierstrass curve for factorization
    :param use_basic: use a basic Montgomery curve factorization
    :return: one non-trivial factor of N or None
    """

    N = as_int(N)
    if N < 0:
        return -1 * ecm(-N, B1=B1, B2=B2, retry=retry, verbose=verbose, use_z1=use_z1)

    if use_weierstrass:
        return ecm_weierstrass(N, B=B1, retry=retry, verbose=verbose)

    if use_basic:
        return ecm_mont_basic(N, B1=B1, B2=B2, retry=retry)

    return ecm_mont(N, B1=B1, B2=B2, retry=retry, verbose=verbose, use_z1=use_z1)
