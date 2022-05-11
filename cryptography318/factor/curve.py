from abc import ABCMeta, abstractmethod
from random import Random
from cryptography318.prime.prime import sqrt_mod


rndm = Random()


class Curve(metaclass=ABCMeta):

    # parameters of each curve are dependent on the form of each curve, but all operate over FF(m)
    __slots__ = 'modulus'

    # __new__ is required for all curves but is not an abstract method since each curve
    # relies on calling super().__new__ which should inherit __new__ from ABCMeta
    # each curve relying on __new__ instead of __init__ also indicates that all curves are immutables

    @abstractmethod
    def __eq__(self, other):
        """
        Determines if two instances of Elliptic Curves are the same.
        """
        raise NotImplementedError

    @abstractmethod
    def point(self, *args):
        raise NotImplementedError


class Weierstrass(Curve):

    __slots__ = 'a', 'b',

    def __new__(cls, a, b, mod):
        self = super().__new__(cls)
        self.a, self.b, self.modulus = a, b, mod
        return self

    def __repr__(self):
        return f"Weierstrass({self.a}, {self.b}, {self.modulus})"

    def __eq__(self, other):
        return other.modulus == self.modulus and other.a == self.a and other.b == self.b

    def point(self, x):
        return WeierstrassPoint(
            x, sqrt_mod((pow(x, 3, self.modulus) + self.a * x + self.b) % self.modulus, self.modulus), self
        )


class Montgomery(Curve):
    """
    Elliptic curve in Montgomery form over finite field. For all explanations of parameterization,
    operations, omission of parameters, etc. refer to
    https://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf (hereon referred to as Gaj.pdf or Gaf)
    """

    __slots__ = 'A_24'

    def __new__(cls, A, mod):
        self = super().__new__(cls)
        self.A_24, self.modulus = A, mod
        return self

    def __repr__(self):
        return f"Montgomery({self.A_24}, {self.modulus})"

    def __eq__(self, other):
        """
        Determines if two instances of Montgomery curves are equal. Returns False iff two curves
        are not equal, True if equal for all parameters other than B since B is not
        saved as an instanced attribute.
        """
        return self.modulus == other.modulus and self.A_24 == other.A_24

    def point(self, *args):
        """
        Currently, I am unsure of how to determine a point on the curve, starting with parameters
        (a + 2) / 4, x, and z, thus this method will remain unimplemented as it is unimportant
        in its utility in Lenstra's ECM.
        """
        raise NotImplementedError


class Point(metaclass=ABCMeta):

    __slots__ = 'x', 'curve', 'modulus'

    # __new__ is required for all points but is not an abstract method since each point
    # relies on calling super().__new__ which should inherit __new__ from ABCMeta
    # each point relying on __new__ instead of __init__ also indicates that all points are immutables

    @abstractmethod
    def __repr__(self): ...

    @abstractmethod
    def __eq__(self, other): ...

    @abstractmethod
    def __neg__(self): ...

    @abstractmethod
    def __add__(self, other): ...

    @abstractmethod
    def __mul__(self, other): ...

    @property
    @abstractmethod
    def points(self): ...


class Identity(Point):

    def __new__(cls):
        self = super().__new__(cls)
        self.x, self.curve, self.modulus = None, None, None
        return self

    def __repr__(self):
        return f'<{__class__.__name__} object>'

    def __eq__(self, other):
        return False

    def __neg__(self):
        return self

    def __add__(self, other):
        return other

    def __mul__(self, other):
        return self

    @property
    def points(self):
        return None, None


class WeierstrassPoint(Point):

    __slots__ = 'y',

    def __new__(cls, x, y, curve):
        self = super().__new__(cls)
        self.x, self.y, self.curve = x % curve.modulus, y % curve.modulus, curve
        self.modulus = curve.modulus
        return self

    def __repr__(self):
        return f'<{__class__.__name__} object; point={self.points}; ' \
               f'curve={self.curve}>'

    def __eq__(self, other):
        return not (self.x - other.x) % self.modulus and not (self.y - other.y) % self.modulus

    def __neg__(self):
        return WeierstrassPoint(self.x, -self.y, self.curve)

    def __add__(self, other):

        # If identity element, return itself
        if isinstance(other, Identity):
            return self

        if not (self.x - other.x) % self.modulus and not (self.y + other.y) % self.modulus:
            return Identity()

        try:

            # If they are the same point, use the doubling algorithm
            if self == other:
                slope = ((3 * self.x * self.x + self.curve.a) * pow(2 * self.y, -1, self.modulus)) % self.modulus
            else:
                slope = ((other.y - self.y) * pow(other.x - self.x, -1, self.modulus)) % self.modulus

        # This except is the basis of lenstra's ecm. raising an error with the specific non-invertible value
        # allows for gcd to be taken to find factor of modulus
        except ValueError:

            # Raise ValueError with message as the value that does not have an inverse
            raise ValueError(str((2 * self.y) if self == other else (other.x - self.x)))

        x3 = (slope * slope - self.x - other.x) % self.modulus
        y3 = (slope * (self.x - x3) - self.y) % self.modulus
        return WeierstrassPoint(x3, y3, self.curve)

    def __mul__(self, other):
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

    @property
    def points(self):
        return self.x, self.y


class MontgomeryPoint(Point):
    """
    Point on Elliptic curve in Montgomery form over finite field with field operations '+' and '*'
    as defined by the _add and __mul__ methods below. For detailed explanation of operations
    refer to Gaj.pdf.
    """

    __slots__ = 'z',

    def __new__(cls, x, z, curve):
        self = super().__new__(cls)
        self.x, self.z, self.curve = x % curve.modulus, z % curve.modulus, curve
        self.modulus = curve.modulus
        return self

    def __repr__(self):
        return f'<{__class__.__name__} object; point={self.points}; ' \
               f'curve={repr(self.curve)}>'

    def __eq__(self, other):
        """
        __eq__ assumes that both points are from the same curve over the same field FF(m). It
        is easy to check if the curves are different (simply compare P.curve == Q.curve) but this
        step is omitted since it is easily done outside __eq__ and would almost always evaluate as True.
        """
        return not (self.x - other.x) % self.modulus and not (self.z - other.z) % self.modulus

    def __neg__(self):
        return MontgomeryPoint(self.x, -self.z, self.curve)

    def add_general(self, Q, difference):
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

        x0, z0 = difference.x, difference.z

        # if z0 == 1, then multiplication is simplified (this is common)
        if z0 == 1:
            return self.add_z1(Q, difference)

        return self.add(Q, difference)

    def add(self, Q, difference):
        """
        Cost:
            *: 6;
            +: 6;
        """
        u, v = (self.x - self.z) * (Q.x + Q.z), (self.x + self.z) * (Q.x - Q.z)
        add, sub = u + v, u - v

        x3 = (difference.z * add * add) % self.modulus
        z3 = (difference.x * sub * sub) % self.modulus
        return MontgomeryPoint(x3, z3, self.curve)

    def add_z1(self, Q, difference):
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
        x3 = (add * add) % self.modulus
        z3 = (difference.x * sub * sub) % self.modulus
        return MontgomeryPoint(x3, z3, self.curve)

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
        x2 = (v * u) % self.modulus
        z2 = (_4xz * (u + self.curve.A_24 * _4xz)) % self.modulus
        return MontgomeryPoint(x2, z2, self.curve)

    def __add__(self, other):
        """
        __add__ is required to be implemented since it is an abstract method but _add requires two
        additional parameters so this is left as null.
        """
        return NotImplemented

    def ladder(self, k):
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

    def ladder_z1(self, k):
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

    def __mul__(self, k):
        if k < 0:
            return (-self).__mul__(-k)

        if self.z == 1:
            return self.ladder_z1(k)

        return self.ladder(k)

    @property
    def points(self):
        return self.x, self.z
