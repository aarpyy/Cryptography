import re

from numbers import *

from math import sqrt, isqrt

from cryptography318.numbers.factor import factor
from .expr import Add, Mul

# val inside radical cannot be decimal with leading 0., even though this would be valid square root
# this is not parsed because it doesn't ever create a valid square root object, would only return float
_SQUARE_ROOT_FORMAT = re.compile(r"""
    \A\s*                       # any leading whitespace
    (?P<coeff>\d*)              # leading coefficient, optional, not necessary for coeff decimal
    (?:\.(?P<coeff_dec>\d*))?   # optional decimal value if coefficient is float
    \s*?                        # optional leading whitespace on multiplication operator
    \*?                         # optional '*' symbol for multiplication operator
    \s*?                        # optional trailing whitespace on multiplication operator
    (?:sqrt\(                   # prefix of sqrt( w/ any capitalization
        (?P<val>\d+))?          # get integer value of float if float else full value, required
    (?:\.(?P<decimal>\d*))?     # get decimal value if float
    \)                          # match to termination of sqrt(), if not matched will return None, invalid sqrt obj
""", re.VERBOSE | re.IGNORECASE)


def normalize(operand, s=1):
    """
    Given an integer radicand x and an integer coefficient y, attempts to reduce x, y
    such that y is the largest integer possible while ensuring that x remains an integer.

    Example
    -------
    >>> normalize(1350, 3)  # x = 1350, y = 3; x = 2 * 3**3 * 5**2; sqrt(x) = 3 * 5 * sqrt(2 * 3)
    6, 45

    above, since x has two primes with powers > 1, they can be taken out of the square root, thus
    2 * 3**3 * 5**2 can be reduced to 2 * 3, removing all of the even powers of primes
    """
    factors = factor(operand)

    if factors == {operand: 1} or factors is None:
        return operand, s

    for f in factors:
        if factors[f] > 1:  # if there is a square or greater in the factors, push it to the coeff
            power = pow(f, factors[f] // 2)
            s *= power
            operand //= (power * power)
    return operand, s


class Sqrt(Real):

    __slots__ = '_radicand'

    def __new__(cls, radicand, *, _normalize=True):
        """
        Constructs new immutable instance of SquareRoot object.

        Parsing of input is done such that the most common instances require the least
        if-statements. A table of the costs are listed below for each combination of types.

        Notes
        -----
        Cost Table [radicand, coefficient] - the value of each cost is the number of boolean statements
        to be evaluated prior to returning. type/int refers to a float or rational that has a value
        of an integer (i.e. float.is_integer() or Rational.denominator == 1). Costs are not given in
        a visually logical order, not in order of cost (but are in general order of frequency as integer
        operands are expected to be most common).

        [int, None]: 4;
        [int, int]: 5;
        [int, float/int]: 8;
        [int, Rational/int]: 10;
        [float/int, None]: 6;
        [float/int, int]: 9;
        [float/int, float/int]: 10;
        [float/int, Rational/int]: 12;
        [Rational/int, None]: 8;
        [Rational/int, int]: 11;
        [Rational/int, float/int]: 13;
        [Rational/int, Rational/int]: 15;
        [str, None]: 11;

        Note
        ----
        The largest cost of creating the SquareRoot object is the normalization, taking out all prime powers > 2
        inside radicand. This step requires a complete factorization of the operand, then several
        multiplications to normalize the expression. Simply put, if efficiency is at question, provide
        prime integers or composite integers with all prime powers < 2 as the radicand, and _normalize=False.
        """

        self = super().__new__(cls)

        if isinstance(radicand, int):
            coefficient = 1
        elif (isinstance(radicand, float) and radicand.is_integer()) or \
                (isinstance(radicand, Rational) and radicand.denominator == 1):
            radicand = int(radicand)
            coefficient = 1
        elif isinstance(radicand, Sqrt):
            return radicand
        elif isinstance(radicand, str):
            match = _SQUARE_ROOT_FORMAT.match(radicand)
            if not match:
                raise ValueError(f"invalid string literal for {self.__class__.__name__}: {radicand}")

            # get groups for coefficient, and coefficient decimal
            co, co_d = match.group('coeff'), match.group('coeff_dec')

            # if coefficient decimal is valued, default integer value of coeff to 0, otherwise if coeff-decimal
            # if not None, default coeff to 1
            scalar = float(co or '1') if not co_d else (float((co or '0') + '.' + co_d))
            value = float(match.group('val') + '.' + (match.group('decimal') or '0'))

            if scalar.is_integer() and value.is_integer():
                radicand, coefficient = int(value), int(scalar)
            else:
                raise ValueError(f"invalid string literal for {self.__class__.__name__}: {radicand}")
        else:
            raise TypeError(f"invalid radicand for {self.__class__.__name__}: {radicand}")

        if _normalize:
            radicand, coefficient = normalize(radicand, coefficient)
        if radicand == 1:
            return coefficient
        elif coefficient == 1:
            self._radicand = radicand
            return self
        else:
            self._radicand = radicand
            return Mul(coefficient, self)

    # all methods below that containt type hints: 'SquareRoot' throw unresolved reference to _ attribute warnings
    # if no hint is given (because they reference private variables _value or _scalar)
    # the hint must also be given in the .py file, not .pyi
    @property
    def eval(self: 'Sqrt'):
        return sqrt(self._radicand)

    @property
    def radicand(self: 'Sqrt'):
        return self._radicand

    # these two properties allow Fraction to easily perform operations with SquareRoot even though its not Rational
    @property
    def numerator(self):
        return self

    @property
    def denominator(self):
        return 1

    def __str__(self: 'Sqrt'):
        return f"sqrt({self._radicand})"

    def __repr__(self: 'Sqrt'):
        return f'Sqrt({self._radicand})'

    def __add__(self, other):
        if isinstance(other, int):
            return Add(self, other)
        elif isinstance(other, float):
            return self.eval + other
        elif isinstance(other, Sqrt):
            if (ar := self.radicand) == other.radicand:
                return Mul(Sqrt(ar, _normalize=False), 2)
            else:
                return Add(self, other)
        elif isinstance(other, complex):
            return complex(self) + other
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, int):
            return Add(other, self)
        elif isinstance(other, float):
            return other + self.eval
        elif isinstance(other, Sqrt):
            if (ar := self.radicand) == other.radicand:
                return Mul(2, Sqrt(ar, _normalize=False))
            else:
                return Add(other, self)
        elif isinstance(other, Real):
            return Add(other, self)
        elif isinstance(other, Complex):
            return complex(other) + complex(self)
        else:
            return NotImplemented

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self: 'Sqrt', other):
        if isinstance(other, int):
            if other == 1:
                return self
            else:
                return Mul(self, other)
        elif isinstance(other, float):
            if other == 1:
                return self
            elif other.is_integer():
                return Mul(self, int(other))
            else:
                return float(self) * other
        elif isinstance(other, Sqrt):
            ar, br = self._radicand, other.radicand
            if ar == br:
                return ar
            else:
                return Sqrt(ar * br)
        elif isinstance(other, complex):
            return complex(self) * other
        else:
            return NotImplemented

    def __rmul__(self: 'Sqrt', other):
        if isinstance(other, int):
            if other == 1:
                return self
            else:
                return Mul(other, self)
        elif isinstance(other, Sqrt):
            ar, br = self._radicand, other.radicand
            if ar == br:
                return ar
            else:
                return Sqrt(ar * br)
        elif isinstance(other, float):
            if other == 1:
                return self
            elif other.is_integer():
                return Mul(int(other), self)
            else:
                return other * float(self)
        elif isinstance(other, Real):
            return Mul(other, self)
        elif isinstance(other, Complex):
            return complex(other) * complex(self)
        else:
            return NotImplemented

    def __floordiv__(self: 'Sqrt', other):
        return self.eval // other

    def __rfloordiv__(self: 'Sqrt', other):
        return other // self.eval

    def __truediv__(self: 'Sqrt', other):
        if isinstance(other, Sqrt):
            rad, r = divmod(self._radicand, other.radicand)
            if r == 0:
                if rad == 1:
                    return 1
                else:
                    return Sqrt(rad)
            else:
                return self.eval / other.eval
        elif other == 1:
            return self
        elif other == -1:
            return -self
        else:
            return self.eval / other

    def __rtruediv__(self: 'Sqrt', other):
        if other % self._radicand == 0:
            return (other // self._radicand) * self
        else:
            return other / self.eval

    def __mod__(self: 'Sqrt', other):
        if isinstance(other, Sqrt) and other.radicand == self._radicand:
            return 0
        else:
            return self.eval % other

    def __rmod__(self, other):
        return other % self.eval

    def __pow__(self: 'Sqrt', power):
        if isinstance(power, int):
            if power > 0:
                d, r = divmod(power, 2)
                if r == 0:
                    return pow(self._radicand, d)
                else:
                    return Mul(pow(self._radicand, d), Sqrt(self._radicand))
            else:
                return pow(self.eval, power)
        else:
            return NotImplemented

    def __rpow__(self: 'Sqrt', other):
        return pow(other, self.eval)

    def __int__(self: 'Sqrt'):
        return isqrt(self._radicand)

    def __float__(self: 'Sqrt'):
        return sqrt(self._radicand)

    def __trunc__(self: 'Sqrt'):
        return int(self.eval // 1)

    def __floor__(self: 'Sqrt'):
        return int(self.eval // 1)

    def __ceil__(self: 'Sqrt'):
        return -int(-self.eval // 1)

    def __round__(self, ndigits=None):
        if ndigits is None:
            return self.__trunc__()
        elif isinstance(ndigits, int):
            shift = pow(10, ndigits)
            return int(self.eval * shift) / shift
        else:
            raise TypeError(f"{ndigits} cannot be interpreted as integer")

    def __abs__(self: 'Sqrt'):
        return +self

    def __pos__(self: 'Sqrt'):
        return Sqrt(self._radicand, _normalize=False)

    def __neg__(self: 'Sqrt'):
        return Mul(-1, Sqrt(self._radicand, _normalize=False))

    def __hash__(self: 'Sqrt'):
        return hash(self.eval)

    def __eq__(self: 'Sqrt', other):
        if isinstance(other, Sqrt):
            return other.radicand == self._radicand
        else:
            return self.eval == other

    def __lt__(self: 'Sqrt', other):
        if isinstance(other, Sqrt):
            return self._radicand < other.radicand
        else:
            return self.eval < other

    def __le__(self: 'Sqrt', other):
        if isinstance(other, Sqrt):
            return self._radicand <= other.radicand
        else:
            return self.eval <= other

    def __gt__(self: 'Sqrt', other):
        if isinstance(other, Sqrt):
            return self._radicand > other.radicand
        else:
            return self.eval > other

    def __ge__(self: 'Sqrt', other):
        if isinstance(other, Sqrt):
            return self._radicand >= other.radicand
        else:
            return self.eval >= other
