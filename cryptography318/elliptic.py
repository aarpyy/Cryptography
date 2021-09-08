from random import randrange, Random
from abc import ABCMeta, abstractmethod
from typing import overload, Union
from math import gcd, log, isqrt, floor, ceil

from .crypto_functions import find_roots
from .prime import primes_gen, isprime


rndm = Random()


class Curve(metaclass=ABCMeta):

    # parameters of each curve are dependent on the form of each curve, but all operate over FF(m)
    __slots__ = 'modulus'

    # __new__ is required for all curves but is not an abstract method since each curve
    # relies on calling super().__new__ which should inherit __new__ from ABCMeta
    # each curve relying on __new__ instead of __init__ also indicates that all curves are immutables

    @abstractmethod
    def __eq__(self, other: 'Curve') -> bool:
        """
        Determines if two instances of Elliptic Curves are the same.
        """
        raise NotImplementedError

    @overload
    @abstractmethod
    def point(self, x: int) -> 'Point':
        """
        Given integer x in field, returns Point object with coordinates as integer given and
        integer solution to the specific curve function.
        """
        raise NotImplementedError

    @overload
    @abstractmethod
    def point(self, x: int, param: int) -> 'Point':
        """
        point() method when curve does not store all parameters as instance attributes.
        """
        raise NotImplementedError

    @abstractmethod
    def point(self, *args) -> 'Point':
        raise NotImplementedError


class Weierstrass(Curve):

    __slots__ = 'a', 'b',

    def __new__(cls, a: int, b: int, mod: int):
        self = super().__new__(cls)
        self.a, self.b, self.modulus = a, b, mod
        return self

    def __repr__(self):
        return f"Weierstrass({self.a}, {self.b}, {self.modulus})"

    def __eq__(self, other: 'Weierstrass'):
        return other.modulus == self.modulus and other.a == self.a and other.b == self.b

    def point(self, x: int) -> 'WeierstrassPoint':
        return WeierstrassPoint(x, find_roots((pow(x, 3, self.modulus) + self.a * x + self.b) % self.modulus, self.modulus), self)

    @classmethod
    def safe_curve_and_point(cls, p: int) -> tuple['Weierstrass', 'WeierstrassPoint']:
        """Computes a non-singular Weierstrass elliptic curve safe for elliptic
        curve cryptography, and a point (x, y) that lies on the curve."""

        while 1:
            a = randrange(p)
            x = randrange(p)
            y = randrange(p)

            # generate b based on point
            b = (pow(y, 2, p) - pow(x, 3, p) - a * x) % p

            # if singular curve, try again
            if 4 * pow(a, 3) + 27 * pow(b, 2) == 0:
                continue
            curve = cls(a, b, p)
            return curve, WeierstrassPoint(x, y, curve)


class Montgomery(Curve):
    """
    Elliptic curve in Montgomery form over finite field. For all explanations of parameterization,
    operations, omission of parameters, etc. refer to
    https://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf (hereon referred to as Gaj.pdf or Gaf)
    """

    __slots__ = 'A_24'

    def __new__(cls, A: int, mod: int):
        self = super().__new__(cls)
        self.A_24, self.modulus = A, mod
        return self

    def __repr__(self):
        return f"Montgomery({self.A_24}, {self.modulus})"

    def __eq__(self, other: 'Montgomery') -> bool:
        """
        Determines if two instances of Montgomery curves are equal. Returns False iff two curves
        are not equal, True if equal for all parameters other than B since B is not
        saved as an instanced attribute.
        """
        return self.modulus == other.modulus and self.A_24 == other.A_24

    def point(self, x: int, param: int) -> 'MontgomeryPoint':
        """
        Currently, I am unsure of how to determine a point on the curve, starting with parameters
        (a + 2) / 4, x, and z, thus this method will remain unimplemented as it is unimportant
        in its utility in Lenstra's ECM.
        """
        raise NotImplementedError

    @classmethod
    def safe_curve_and_point(cls, p: int) -> Union[tuple['Montgomery', 'MontgomeryPoint'], int]:
        """
        Implementation of Suyama's parameterization of Montgomery curves.

        Curve constructed s.t. E over FF(q) for some prime q has order = 12 * k for some
        positive integer k. This helps in finding some q factor of N s.t. q is a product of small
        primes, in the computational process of Lenstra's ECM.

        Refer to section 2.3 p.4 of Gaj.pdf
        """
        sigma = rndm.randrange(6, p - 1)
        sig_sq = sigma * sigma
        u = (sig_sq - 5) % p
        v = (4 * sigma) % p
        u_sub_v = v - u
        u_3 = u * u * u
        try:
            a_sig = ((u_sub_v * u_sub_v * u_sub_v) * (3 * u + v) * pow(4 * u_3 * v, -1, p) - 2) % p
        except ValueError:

            # if composite m for FF(m) then this curve is likely intended for Lenstra's ECM and gcd is factor
            return gcd(4 * u_3 * v, p)

        x, z = u_3 % p, (v * v * v) % p
        E = Montgomery(a_sig, p)
        P = MontgomeryPoint(x, z, E)
        return E, P


class Point(metaclass=ABCMeta):

    __slots__ = 'x', 'curve', 'modulus'

    # __new__ is required for all points but is not an abstract method since each point
    # relies on calling super().__new__ which should inherit __new__ from ABCMeta
    # each point relying on __new__ instead of __init__ also indicates that all points are immutables

    @abstractmethod
    def __repr__(self): ...

    @abstractmethod
    def __eq__(self, other: 'Point') -> bool: ...

    @abstractmethod
    def __neg__(self) -> 'Point': ...

    @abstractmethod
    def __add__(self, other: 'Point') -> 'Point': ...

    @abstractmethod
    def __mul__(self, other: int) -> 'Point': ...

    @property
    @abstractmethod
    def points(self) -> tuple: ...


class Identity(Point):

    def __new__(cls):
        self = super().__new__(cls)
        self.x, self.curve, self.modulus = None, None, None
        return self

    def __repr__(self):
        return f'<{__class__.__name__} object>'

    def __eq__(self, other: 'Point') -> bool:
        return False

    def __neg__(self) -> 'Point':
        return self

    def __add__(self, other: 'Point') -> 'Point':
        return other

    __radd__ = __add__

    def __mul__(self, other: int) -> 'Point':
        return self

    __rmul__ = __mul__

    @property
    def points(self) -> tuple:
        return None, None


class WeierstrassPoint(Point):

    __slots__ = 'y',

    def __new__(cls, x: int, y: int, curve: Weierstrass):
        self = super().__new__(cls)
        self.x, self.y, self.curve = x % curve.modulus, y % curve.modulus, curve
        self.modulus = curve.modulus
        return self

    def __repr__(self):
        return f'<{__class__.__name__} object; point={self.points}; ' \
               f'curve={self.curve}>'

    def __eq__(self, other: 'WeierstrassPoint'):
        return not (self.x - other.x) % self.modulus and not (self.y - other.y) % self.modulus

    def __neg__(self):
        return WeierstrassPoint(self.x, -self.y, self.curve)

    def __add__(self, other: 'WeierstrassPoint'):

        # if identity element, return itself
        if isinstance(other, Identity):
            return self

        if not (self.x - other.x) % self.modulus and not (self.y + other.y) % self.modulus:
            return Identity()

        try:
            if self == other:
                slope = ((3 * self.x * self.x + self.curve.a) * pow(2 * self.y, -1, self.curve.modulus)) \
                        % self.curve.modulus
            else:
                slope = ((other.y - self.y) * pow(other.x - self.x, -1, self.curve.modulus)) % self.curve.modulus

        # this except is the basis of lenstra's ecm. raising an error with the specific non-invertible value
        # allows for gcd to be taken to find factor of modulus
        except ValueError as exc:
            raise ValueError(str(exc) + f"base: {(2 * self.y) if self == other else (other.x - self.x)}")

        x3 = (slope * slope - self.x - other.x) % self.curve.modulus
        y3 = (slope * (self.x - x3) - self.y) % self.curve.modulus
        return WeierstrassPoint(x3, y3, self.curve)

    __radd__ = __add__

    def __mul__(self, other: int) -> Union['WeierstrassPoint', 'Identity']:
        P = self
        if other < 0:
            other = -other
            P = -P

        Q = P
        R = Identity()
        while other > 0:
            if other & 1:
                R += Q
            Q += Q
            other //= 2
        return R

    __rmul__ = __mul__

    @property
    def points(self) -> tuple:
        return self.x, self.y


class MontgomeryPoint(Point):
    """
    Point on Elliptic curve in Montgomery form over finite field with field operations '+' and '*'
    as defined by the _add and __mul__ methods below. For detailed explanation of operations
    refer to Gaj.pdf.
    """

    __slots__ = 'z',

    def __new__(cls, x: int, z: int, curve: Montgomery):
        self = super().__new__(cls)
        self.x, self.z, self.curve = x % curve.modulus, z % curve.modulus, curve
        self.modulus = curve.modulus
        return self

    def __repr__(self):
        return f'<{__class__.__name__} object; point=({self.x}, {self.z}); ' \
               f'curve={repr(self.curve)}>'

    def __eq__(self, other: 'MontgomeryPoint') -> bool:
        """
        __eq__ assumes that both points are from the same curve over the same field FF(m). It
        is easy to check if the curves are different (simply compare P.curve == Q.curve) but this
        step is omitted since it is easily done outside __eq__ and would almost always evaluate as True.
        """
        return not (self.x - other.x) % self.modulus and not (self.z - other.z) % self.modulus

    def __neg__(self) -> 'MontgomeryPoint':
        return MontgomeryPoint(self.x, -self.z, self.curve)

    def add_general(self, Q: 'MontgomeryPoint', initial: 'MontgomeryPoint'):
        """
        Computes addition on E over FF(m). For specifically Montgomery curves, operation involving
        symbol '+' is unsupported, due to the necessity of the initial point in addition. For
        all methods of addition this function should be used. If it is known that the points are different and
        z == 1, use _add_z1, otherwise use _add if z != 1. If pints are the same, use _double.

        Refer to Gaj p.4, algorithm 3
        """
        # if same point, double
        if self == Q:
            return self.double()

        x0, z0 = initial.x, initial.z

        # if z0 == 1, then multiplication is simplified (this is common)
        if z0 == 1:
            return self.add_z1(Q, initial)

        return self.add(Q, initial)

    def add(self, Q: 'MontgomeryPoint', initial: 'MontgomeryPoint'):
        """
        Cost:
            *: 6;
            +: 6;
        """
        u, v = self.x - self.z * Q.x + Q.z, self.x + self.z * Q.x - Q.z
        add, sub = u + v, u - v

        x3 = (initial.z * add * add) % self.modulus
        z3 = (initial.x * sub * sub) % self.modulus
        return MontgomeryPoint(x3, z3, self.curve)

    def add_z1(self, Q: 'MontgomeryPoint', initial: 'MontgomeryPoint'):
        """
        Computes addition P + Q on E over FF(m) when z_P == 1 (z coordinate of initial point)

        Refer to Gaj p.4, algorithm 3. For explanation of z coordinate == 1, refer to just below
        algorithm 3 p.4

        Cost:
            *: 5;
            +: 6;

        Note
        ----
        While this cost is small, _add_z1 has one less % operation than _add in addition to the reduction in *
        """
        selfdiff, Qdiff = self.x - self.z, Q.x - Q.z
        selfsum, Qsum = self.x + self.z, Q.x + Q.z
        u, v = selfdiff * Qsum, selfsum * Qdiff
        add, sub = u + v, u - v
        x3 = (add * add) % self.modulus
        z3 = (initial.x * sub * sub) % self.modulus
        return MontgomeryPoint(x3, z3, self.curve)

    def double(self):
        """
        Computes addition P + Q on E over FF(m) when P == Q

        Refer to Gaj p.4, algorithm 3
        """
        selfdiff, selfsum = self.x - self.z, self.x + self.z
        selfdiff, selfsum = selfdiff * selfdiff, selfsum * selfsum
        self4xz = selfsum - selfdiff
        x2 = (selfsum * selfdiff) % self.modulus
        z2 = (self4xz * (selfdiff + self.curve.A_24 * self4xz)) % self.modulus
        return MontgomeryPoint(x2, z2, self.curve)

    def __add__(self, other):
        """
        __add__ is required to be implemented since it is an abstract method but _add requires two
        additional parameters so this is left as null.
        """
        return NotImplemented

    def ladder(self, k: int) -> 'MontgomeryPoint':
        """
        Refer to Gaj p.4, algorithm 2
        """
        P = self.double()
        Q = self
        for bit in bin(k)[3:]:

            # compare via strings instead of convert to int is faster
            if bit == '1':
                Q = Q.add(P, self)
                P = P.double()
            else:
                Q = Q.double()
                P = P.add(Q, self)
        return Q

    def ladder_z1(self, k: int) -> 'MontgomeryPoint':
        """
        Refer to Gaj p.4, algorithm 2
        """
        P = self.double()
        Q = self
        for bit in bin(k)[3:]:
            if bit == '1':
                Q = Q.add_z1(P, self)
                P = P.double()
            else:
                Q = Q.double()
                P = P.add_z1(Q, self)
        return Q

    def __mul__(self, k: int) -> 'MontgomeryPoint':
        if k < 0:
            return (-self).__mul__(-k)

        if self.z == 1:
            return self.ladder_z1(k)

        return self.ladder(k)

    @property
    def points(self) -> tuple:
        return self.x, self.z


def ecm_mont(N: int, B1: int = None, B2: int = None, _retry: int = 50):
    """
    Performs Lenstra's elliptic curve factorization method with elliptic curve E
    in Montgomery form over FF(N) with Suyama's parameterization.

    References
    ----------
    [1] Gaj K. et al. Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware
    https://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf

    :param N: integer to be factored
    :param B1: integer limit for phase 1; start for phase 2
    :param B2: integer limit for phase 2
    :param _retry: integer number of trials to be run
    :return: one integer factor of N or None
    """
    if B1 is None:
        from math import e
        B1 = int(pow(e, isqrt(int((log(N) * log(log(N))) / 2))))

    if B2 is None:
        B2 = B1 * 100

    # with initial point P, get k to find initial Q = kP
    k = 1
    for p in primes_gen(2, B1):
        k *= floor(log(B1, p))

    trials = 0
    while trials <= _retry:
        trials += 1

        # get montgomery curve for phase 1
        res = Montgomery.safe_curve_and_point(N)
        if isinstance(res, int):
            return res

        E, P = res
        Q = P.ladder(k)

        # if found factor already, return, otherwise continue to phase 2
        q = gcd(Q.z, N)
        if 1 < q < N:
            return q

        d = 1

        for p in primes_gen(B1, B2):
            d = (d * Q.ladder(p).z) % N

        q = gcd(d, N)
        if 1 < q < N:
            return q

    return None


def ecm_mont_large(N: int, B1: int = None, B2: int = None, _retry: int = 50):
    """
    A modification of the above ecm_mont() function that relies on a pre-processing
    step for phase 2 to find m, j s.t. for prime p B1 < p < B2, p = m*D +/- j for
    D = sqrt(B2), p*P = Identity meaning (m*D +/- j) * P = O, thus mD * P = +/- jP
    for mD != +/- j which allows for d = prod(X_mDP * Z_jP - X_jP * Z_mDP) to be
    found s.t. q = gcd(d, N), 1 < q < N, q is a factor of N (not necessarily prime).

    Notes
    -----
    While the exact value is currently unknown, this method should be used for larger
    composites, where the cost of the pre-processing steps are absorbed by the
    time complexity of factoring. Since the time of Lenstra's ECM is dominated by the
    smallest factor of composite N, this method should be reserved for when N's smallest
    prime factor is larger (until this pre-processing is optimized).
    """
    if B1 is None:
        from math import e
        B1 = int(pow(e, isqrt(int((log(N) * log(log(N))) / 2))))

    if B2 is None:
        B2 = B1 * 100

    # with initial point P, get k to find initial Q = kP
    k = 1
    for p in primes_gen(2, B1):
        k *= floor(log(B1, p))

    D = isqrt(B2)

    # faster to do both separate computations than to do one and copy s
    S = [0] * (D + 1)
    gcd_table = [0] * (D + 1)
    set_j = set()
    for j in range(1, D//2):
        if gcd(j, D) == 1:
            gcd_table[j] = 1
            set_j.add(j)

    prime_table = []
    M_min = floor((B1 + D/2)/D)
    M_max = ceil((B2 - D/2)/D)
    M_range = range(M_max - M_min)

    for m in M_range:
        prime_table.append([0] * (D + 1))
        for j in set_j:
            if isprime(m * D + j) or isprime(m * D - j):
                prime_table[m][j] = 1

    trials = 0
    while trials <= _retry:
        trials += 1

        # get montgomery curve for phase 1
        res = Montgomery.safe_curve_and_point(N)
        if isinstance(res, int):
            return res

        E, P = res
        Q_0 = P.ladder(k)

        # if found factor already, return, otherwise continue to phase 2
        q = gcd(Q_0.z, N)
        if 1 < q < N:
            return q

        Q = Q_0
        R = Q.ladder(M_min)

        for j in gcd_table:
            if gcd_table[j]:
                S[j] = Q
            Q = Q.add(Q_0.double(), Q_0)

        Q = Q_0.ladder(D)
        d = 1
        for m in M_range:
            for j in prime_table[m]:
                if prime_table[m][j]:
                    j_Q = S[j]
                    d = (d * (R.x * j_Q.z - R.z * j_Q.x)) % N
            R = R.add(Q, Q_0)

        q = gcd(d, N)
        if 1 < q < N:
            return q

    return None


def ecm_weierstrass(N: int, B: int = None, _retry: int = 50):

    if B is None:
        from math import e
        L = pow(e, isqrt(int(log(N) * log(log(N)))))
        B = int(pow(L, 1 / pow(2, .5)))

    trials = 0
    while trials <= _retry:
        trials += 1

        E, P = Weierstrass.safe_curve_and_point(N)

        try:
            for j in range(1, B):
                P = j * P
        except ValueError as exc:
            k = int(str(exc).split('base: ')[1])
            q = gcd(k, N)
            if 1 < q < N:
                return q

    return None
