import re
import operator

from typing import *
from numbers import *

from abc import abstractmethod
from math import sqrt, gcd
from functools import reduce
from typing import Any, overload
from sympy.ntheory.primetest import is_square
from math import nan
from decimal import Decimal

from .crypto_functions import factor
from .tools import lcm, number_to_integer
from .prime import isprime


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


def sqrt_gcd(*args):
    def _gcd(a, b):
        if isinstance(a, SquareRoot):
            a = a._scalar
        if isinstance(b, SquareRoot):
            b = b._scalar
        return gcd(a, b)

    return reduce(lambda i, c: _gcd(i, c), args)


def reduce_sqrt(operand, s=1):
    factors = factor(operand)

    if factors is None or factors == {operand: 1}:
        return operand, s

    for f in factors:
        if factors[f] > 1:  # if there is a square or greater in the factors, push it to the coeff
            power = pow(f, factors[f] // 2)
            s *= power
            operand //= pow(power, 2)
    return operand, s


class SquareRoot(Real):

    __slots__ = ('_scalar', '_value')  # sets data stored in RAM to just these two vals (not using python dict)

    def __new__(cls, operand, coefficient=None, _normalize=True):

        self = super().__new__(cls)

        if coefficient is None:
            if isinstance(operand, int):
                self._value, self._scalar = reduce_sqrt(operand) if _normalize else (operand, 1)
                return self._scalar if self._value == 1 else self
            elif isinstance(operand, float):
                if operand.is_integer():
                    self._value, self._scalar = reduce_sqrt(int(operand)) if _normalize else (int(operand), 1)
                    return self._scalar if self._value == 1 else self
                return sqrt(operand)
            elif isinstance(operand, str):
                match = _SQUARE_ROOT_FORMAT.match(operand)
                if not match:
                    raise ValueError(f"invalid string literal for {__class__.__name__}: {operand}")

                # get groups for coefficient, and coefficient decimal
                co, co_d = match.group('coeff'), match.group('coeff_dec')

                # if coefficient decimal is valued, default integer value of coeff to 0, otherwise if coeff-decimal
                # if not None, default coeff to 1
                scalar = float(co or '1') if not co_d else (float((co or '0') + '.' + co_d))
                value = float(match.group('val') + '.' + (match.group('decimal') or '0'))

                if scalar.is_integer() and value.is_integer():
                    self._value, self._scalar = reduce_sqrt(int(value), int(scalar)) if _normalize \
                        else (int(value), int(scalar))
                    return self._scalar if self._value == 1 else self
                else:
                    res = float(scalar * sqrt(value))
                    return int(res) if res.is_integer() else res
            # if rational but denominator can be removed or is 1
            elif isinstance(operand, Rational) and operand.denominator == 1:
                self._value, self._scalar = reduce_sqrt(operand.numerator) if _normalize else (operand.numerator, 1)
                return self._scalar if self._value == 1 else self
            elif isinstance(operand, Complex) and operand.imag == 0:
                operand = operand.real
                if isinstance(operand, int) or (isinstance(operand, float) and operand.is_integer()):
                    self._value, self._scalar = reduce_sqrt(int(operand)) if _normalize else (int(operand), 1)
                    return self._scalar if self._value == 1 else self
                raise ValueError(f"complex numbers must have imaginary value of 0")
            else:
                try:
                    return sqrt(operand)
                except TypeError:
                    raise ValueError(f"invalid argument(s) for {__class__.__name__}: {operand}")
        elif isinstance(operand, int) and isinstance(coefficient, int):
            pass
        elif isinstance(operand, (float, Complex, Rational)) or isinstance(coefficient, (float, Complex, Rational)):
            operand = number_to_integer(operand)
            coefficient = number_to_integer(coefficient)
            if isinstance(operand, int) and isinstance(coefficient, int):
                pass

            # only instance of when operand can't be converted into integer but still is valid is if denominator
            # can be pulled out of square root
            elif isinstance(coefficient, int) and isinstance(operand, Rational) and is_square(operand.denominator):
                if coefficient % sqrt(operand.denominator) == 0:
                    operand = int(operand.numerator)
                    coefficient //= int(sqrt(operand.denominator))
            else:
                try:
                    return coefficient * sqrt(operand)
                except TypeError:
                    raise ValueError(f"invalid argument(s) for {__class__.__name__}: "
                                     f"operand={operand}, coefficient={coefficient}")
        else:
            try:
                return coefficient * sqrt(operand)
            except TypeError:
                raise ValueError(f"invalid argument(s) for {__class__.__name__}: "
                                 f"operand={operand}, coefficient={coefficient}")

        # at this point, every input has either been returned, an error raised, or converted into at int, meaning
        # values here are guaranteed to be integers

        self._value, self._scalar = reduce_sqrt(operand, coefficient) if _normalize else (operand, coefficient)
        return self._scalar if self._value == 1 else self

    def _operator_fallbacks(square_root_operator: Callable[['SquareRoot', Number], Number],
                            standard_operator: Callable[[Any, Any], Any]):
        def forward(a, b):
            if isinstance(b, (int, float, SquareRoot)):
                return square_root_operator(a, b)
            elif isinstance(b, complex):
                return standard_operator(complex(a), b)
            else:

                # if can't match to a known type, convert sqrt to float (nearest type) and try again
                return standard_operator(float(a), b)

        forward.__name__ = '__' + standard_operator.__name__ + '__'
        forward.__doc__ = square_root_operator.__doc__

        def reverse(b, a):
            if isinstance(a, SquareRoot):
                return square_root_operator(a, b)
            elif isinstance(a, (int, float)):
                return standard_operator(a, float(b))
            elif isinstance(a, Real):
                return standard_operator(float(a), float(b))
            elif isinstance(a, Complex):
                return standard_operator(complex(a), complex(b))
            else:

                # if can't match to a known type, convert sqrt to float (nearest type) and try again
                return standard_operator(a, float(b))

        reverse.__name__ = '__r' + standard_operator.__name__ + '__'
        reverse.__doc__ = square_root_operator.__doc__

        return forward, reverse

    def __str__(self):
        if self._value == 1:
            return str(self._scalar)
        return f'({self._scalar}*sqrt({self._value}))' if self._scalar != 1 else f'sqrt({self._value})'

    def __repr__(self):
        return f'SquareRoot({self._value}, {self._scalar})'

    def __int__(self):
        return int(self.value)

    def __float__(self):
        return float(self.value)

    def __complex__(self):
        return complex(self.value)

    def __abs__(self):
        return SquareRoot(self._value, abs(self._scalar))

    def __pos__(self):
        return SquareRoot(self._value, self._scalar)

    def __neg__(self):
        return SquareRoot(self._value, -self._scalar)

    def __eq__(self, other):
        if isinstance(other, SquareRoot):

            # second part of or should never be true if first part is false, there as a fail safe
            return (other._value == self._value and other._scalar == self._scalar) or (self.value == other.value)
        return self.value == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def _add(a, b):
        if isinstance(b, SquareRoot) and a._value == b._value:
            return SquareRoot(a._value, a._scalar + b._scalar)
        return a.value + b

    __add__, __radd__ = _operator_fallbacks(_add, operator.add)

    def _sub(a, b):
        if isinstance(b, SquareRoot) and a._value == b._value:
            return SquareRoot(a._value, a._scalar - b._scalar)
        return a.value - b

    __sub__, __rsub__ = _operator_fallbacks(_sub, operator.sub)

    def _mul(a, b):
        if isinstance(b, SquareRoot):
            if a._value == b._value:
                return a._scalar * b._scalar * a._value
            return SquareRoot(a._value * b._value, a._scalar * b._scalar)
        elif isinstance(b, int):
            return SquareRoot(a._value, a._scalar * b)
        return a.value * b

    __mul__, __rmul__ = _operator_fallbacks(_mul, operator.mul)

    def _floordiv(a, b):
        if isinstance(b, SquareRoot):
            if a._value == b._value:
                return a._scalar // b._scalar
            elif a._value % b._value == 0:
                return SquareRoot(a._value // b._value, a._scalar // b._scalar)
            return a.value // b.value
        elif isinstance(b, int) and a._scalar % b == 0:
            return SquareRoot(a._value, a._scalar // b)
        return a.value // b

    __floordiv__, __rfloordiv__ = _operator_fallbacks(_floordiv, operator.floordiv)

    def _truediv(a, b):
        if isinstance(b, SquareRoot):
            if a._value == b._value:
                return a._scalar / b._scalar
            elif a._value % b._value == 0 and a._scalar % b._scalar == 0:
                return SquareRoot(a._value // b._value, a._scalar / b._scalar)
            return a.value / b.value
        elif isinstance(b, int) and a._scalar % b == 0:
            return SquareRoot(a._value, a._scalar / b)
        return a.value / b

    __truediv__, __rtruediv__ = _operator_fallbacks(_truediv, operator.truediv)

    def __hash__(self) -> int:
        return float(self.value).__hash__()

    def __trunc__(self) -> int:
        return float(self.value).__trunc__()

    def __floor__(self) -> int:
        return float(self.value).__floor__()

    def __ceil__(self) -> int:
        return float(self.value).__ceil__()

    @overload
    def __round__(self, ndigits: None = ...) -> Integral: ...

    @overload
    def __round__(self, ndigits: int = ...) -> Real: ...

    def __round__(self, ndigits: int = None) -> Any:
        return float(self.value).__round__(ndigits=ndigits)

    def __mod__(self, other: Any) -> Any:
        return float(self.value).__mod__(other)

    def __rmod__(self, other: Any) -> Any:
        return float(self.value).__rmod__(other)

    def __lt__(self, other: Any) -> bool:
        return float(self.value).__lt__(other)

    def __le__(self, other: Any) -> bool:
        return float(self.value).__le__(other)

    def __pow__(self, power):
        if power % 2 == 0:
            return pow(self._value, power // 2) * pow(self._scalar, power)  # power // 2 for value since inside sqrt
        return pow(self.value, power)

    def __rpow__(self, other):
        return float(self.value).__rpow__(other)

    @property
    def value(self):
        if self._value == 1:
            return self._scalar
        return sqrt(self._value) * self._scalar


class Fraction:
    def __init__(self, num, denom=1.0):

        # attempts to check if number if is a square root, creating SR obj if it is
        if isinstance(denom, float) and not denom.is_integer():
            if abs(int(k := pow(denom, 2)) / k) < 1 + pow(10, -6):  # if less than .000001% away from sqrt, it is sqrt
                denom = SquareRoot(int(k))
            else:
                raise ValueError(f"denominator must be rational or square root value")
        if isinstance(num, float) and not num.is_integer():
            if abs(int(k := pow(num, 2)) / k) < 1 + pow(10, -6):  # if less than .000001% away from sqrt, it is sqrt
                num = SquareRoot(int(k))
            else:
                raise ValueError(f"denominator must be rational or square root value")

        # square root in denominator, move to numerator
        if isinstance(denom, SquareRoot):
            f = SquareRoot(denom._value)
            denom *= f
            num *= f

        # make sure integers, if float then guaranteed to be integer-float
        self.numerator = num if isinstance(num, SquareRoot) else int(num)
        self.denominator = int(denom)

        g = sqrt_gcd(self.numerator, self.denominator)  # if fraction can be reduced, do it, use sqrt_gcd for sqrt objs
        self.numerator //= g
        self.denominator //= g

    def __str__(self):
        return f'{self.numerator}/{self.denominator}'

    def __repr__(self):
        return f'Fraction({self.numerator}, {self.denominator})'

    def __eq__(self, other):
        if isinstance(other, Fraction):
            return self.numerator == other.numerator and self.denominator == other.denominator
        return self.eval() == other

    def __add__(self, other):
        if isinstance(other, Fraction):
            f = lcm(self.denominator, other.denominator)
            numer = self.numerator * (f // self.denominator) + other.numerator * (f // other.denominator)
            return Fraction(numer, f)

        return Fraction(self.numerator + other * self.denominator, self.denominator)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Fraction):
            f = lcm(self.denominator, other.denominator)
            numer = self.numerator * (f // self.denominator) - other.numerator * (f // other.denominator)
            return Fraction(numer, f)

        return Fraction(self.numerator - other * self.denominator, self.denominator)

    def __rsub__(self, other):
        if isinstance(other, Fraction):
            f = lcm(self.denominator, other.denominator)
            numer = other.numerator * (f // other.denominator) - self.numerator * (f // self.denominator)
            return Fraction(numer, f)

        return Fraction(self.numerator - other * self.denominator, self.denominator)

    def __mul__(self, other):
        if isinstance(other, Fraction):
            return Fraction(self.numerator * other.numerator, self.denominator * other.denominator)
        return Fraction(self.numerator * other, self.denominator)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Fraction):
            return Fraction(self.numerator * other.denominator, self.denominator * other.numerator)
        return Fraction(self.numerator, self.denominator * other)

    def __rtruediv__(self, other):
        pass

    def eval(self):
        return self.numerator / self.denominator
