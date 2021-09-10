import re
import operator

from typing import *
from numbers import *

from abc import abstractmethod
from math import sqrt, gcd
from functools import reduce
from typing import Any, overload
from decimal import Decimal

from .factor import factor


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
    def _gcd(a: 'SquareRoot', b: 'SquareRoot'):
        return gcd(a.radicand, b.radicand)
    return reduce(lambda i, c: _gcd(i, c), args)


def reduce_sqrt(operand, s=1):
    factors = factor(operand)

    if factors == {operand: 1} or factors is None:
        return operand, s

    for f in factors:
        if factors[f] > 1:  # if there is a square or greater in the factors, push it to the coeff
            power = pow(f, factors[f] // 2)
            s *= power
            operand //= (power * power)
    return operand, s


def try_sqrt(v, c=None):
    """
    Helper function for SquareRoot.__new__() that attempts to return square root of parameters
    given, throwing a more helpful error if a problem occurs.
    """
    try:
        return sqrt(v) if c is None else (c * sqrt(v))
    except TypeError:
        raise TypeError(f"all input(s) must be rational not type(s): {type(v)}" +
                        f", {type(c)}" if c is not None else "")
    except OverflowError:
        raise ValueError(f"{v} too large to be converted to SquareRoot")


class SquareRoot(Real):

    __slots__ = ('_value', '_scalar')

    def __new__(cls, operand, coefficient=None, _normalize=True):
        """
        Constructs new immutable instance of SquareRoot object.

        Parsing of input is done such that the most common instances require the least
        if-statements. A table of the costs are listed below for each combination of types.

        Notes
        -----
        Cost Table [operand, coefficient] - the value of each cost is the number of boolean statements
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
        inside operand. This step requires a complete factorization of the operand, then several
        multiplications to normalize the expression. Simply put, if efficiency is at question, provide
        prime integers or composite integers with all prime powers < 2 as the operand, and _normalize=False.
        """

        self = super().__new__(cls)

        if coefficient is None:
            if isinstance(operand, int):
                self._value, self._scalar = reduce_sqrt(operand) if _normalize else (operand, 1)
                return self._scalar if self._value == 1 else self
            elif (isinstance(operand, float) and operand.is_integer()) or \
                    (isinstance(operand, Rational) and operand.denominator == 1):
                operand = int(operand)
                self._value, self._scalar = reduce_sqrt(operand) if _normalize else (operand, 1)
                return self._scalar if self._value == 1 else self
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
                    return try_sqrt(value, scalar)
            else:
                return try_sqrt(operand)

        # if statement in this order allows for quicker evaluation, since it is more common that coefficient
        # will not be int than operand
        elif isinstance(coefficient, int) and isinstance(operand, int):
            pass
        elif (isinstance(operand, int) or ((isinstance(operand, float) and operand.is_integer()) or
              (isinstance(operand, Rational) and operand.denominator == 1))) and \
                (isinstance(coefficient, int) or ((isinstance(coefficient, float) and coefficient.is_integer()) or
                 (isinstance(coefficient, Rational) and coefficient.denominator == 1))):
            operand = int(operand)
            coefficient = int(coefficient)
        else:
            return try_sqrt(operand, coefficient)

        # at this point, every input has either been returned, an error raised, or converted into at int, meaning
        # values here are guaranteed to be integers

        self._value, self._scalar = reduce_sqrt(operand, coefficient) if _normalize else (operand, coefficient)
        return self._scalar if self._value == 1 else self

    @property
    def eval(self):
        return sqrt(self._value) * self._scalar

    @property
    def radicand(self):
        return self._value

    @property
    def coefficient(self):
        return self._scalar

    def _operator_fallbacks(sr_operator: Callable[['SquareRoot', 'SquareRoot'], Real], std_operator):
        """
        Takes in two callable operators and returns two callable operators analogous to __op__, and
        __rop__

        :param sr_operator: square root operator: _op (ex. _add())
        :param std_operator: standard operator __op__ (ex. __add__())
        """
        def forward(a, b):
            if isinstance(b, SquareRoot):
                return sr_operator(a, b)
            elif isinstance(b, (int, float)):
                return std_operator(float(a), b)
            elif isinstance(b, complex):
                return std_operator(complex(a), b)
            else:
                return NotImplemented

        forward.__name__ = '__' + std_operator.__name__ + '__'
        forward.__doc__ = sr_operator.__doc__

        def reverse(b, a):
            if isinstance(a, SquareRoot):
                return sr_operator(a, b)
            elif isinstance(a, Real):
                return std_operator(float(a), float(b))
            elif isinstance(a, Complex):
                return std_operator(complex(a), complex(b))
            else:

                # since there is nothing to fall back on, instead of returning not implemented, operation
                # is attempted with square root being converted to nearest type: float, incase the other
                # type knows how to deal with floats
                return std_operator(a, float(b))

        reverse.__name__ = '__r' + std_operator.__name__ + '__'
        reverse.__doc__ = sr_operator.__doc__

        return forward, reverse

    def _add(a, b: 'SquareRoot'):
        ar = a.radicand
        if ar == b.radicand:
            return SquareRoot(ar, a.coefficient + b.coefficient)
        return a.eval + b.eval

    __add__, __radd__ = _operator_fallbacks(_add, operator.add)

    def _sub(a, b: 'SquareRoot'):
        ar = a.radicand
        if ar == b.radicand:
            return SquareRoot(ar, a.coefficient - b.coefficient)
        return a.eval - b.eval

    __sub__, __rsub__ = _operator_fallbacks(_sub, operator.sub)

    def _mul(a, b: 'SquareRoot'):
        ar, br = a.radicand, b.radicand
        if ar == br:
            return a.coefficient * b.coefficient * ar
        return SquareRoot(ar * br, a.coefficient * b.coefficient)

    __mul__, __rmul__ = _operator_fallbacks(_mul, operator.mul)

    def _floordiv(a, b: 'SquareRoot'):
        ar, br = a.radicand, b.radicand
        ac, bc = a.coefficient, b.coefficient
        if ar == br:
            return ac // bc
        elif not ar % br and not ac % bc:
            return SquareRoot(ar // br, ac // bc)
        return a.eval // b.eval

    __floordiv__, __rfloordiv__ = _operator_fallbacks(_floordiv, operator.floordiv)

    def _truediv(a, b: 'SquareRoot'):
        ar, br = a.radicand, b.radicand
        ac, bc = a.coefficient, b.coefficient
        if ar == br:
            return ac / bc
        elif not ar % br and not ac % bc:
            return SquareRoot(ar // br, ac // bc)
        return a.eval / b.eval

    __truediv__, __rtruediv__ = _operator_fallbacks(_truediv, operator.truediv)

    def __mod__(self, other):
        if isinstance(other, int):
            r = self.radicand
            c = self.coefficient % other
            if c * sqrt(self.radicand) < other:
                return SquareRoot(r, c)
        return self.eval % other

    def __rmod__(self, other):
        return other % self.eval

    def __pow__(self, power):
        if not power % 2:
            return pow(self.radicand, power // 2) * pow(self.coefficient, power)  # power // 2 for value inside sqrt
        return pow(self.eval, power)

    def __rpow__(self, other):
        return self.eval.__rpow__(other)

    def __float__(self):
        return self.eval

    def __int__(self):
        return int(self.eval)

    def __complex__(self):
        return complex(self.eval)

    def __trunc__(self) -> int:
        return self.eval.__trunc__()

    def __floor__(self) -> int:
        return self.eval.__floor__()

    def __ceil__(self) -> int:
        return self.eval.__ceil__()

    def __abs__(self):
        return SquareRoot(self.radicand, abs(self.coefficient), _normalize=False)

    def __neg__(self):
        return SquareRoot(self.radicand, -self.coefficient, _normalize=False)

    def __pos__(self):
        return SquareRoot(self._value, self._scalar, _normalize=False)

    def __hash__(self):
        return self.eval.__hash__()

    def __round__(self, ndigits=None):
        return self.eval.__round__(ndigits=ndigits)

    def __str__(self):
        r, c = self.radicand, self.coefficient
        return f'sqrt({r})' if c == 1 else f'({c}*sqrt({r}))'

    def __repr__(self):
        return f'SquareRoot({self._value}, {self._scalar})'

    def __eq__(self, other):
        if isinstance(other, SquareRoot):
            return other.radicand == self.radicand and other.coefficient == self.coefficient
        return self.eval == other

    def __ne__(self, other):
        if isinstance(other, SquareRoot):
            return not (other.radicand == self.radicand or other.coefficient == self.coefficient)
        return self.eval != other

    def _compare(self, other, op):
        if isinstance(other, SquareRoot):
            return op(self.eval, other.eval)
        return op(self.eval, other)

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)


class Fraction(Rational):
    __slots__ = '_numerator', '_denominator'

    # the , *, allows for numerator and denominator to be passed as positional arguments, but anything past
    # that (i.e. _normalize) must be passed explicitly as a kwarg
    def __new__(cls, numerator, denominator=None, *, _normalize=True):

        self = super().__new__(cls)

        if denominator is None:
            if type(numerator) is int:
                self._numerator = numerator
                self._denominator = 1
                return self

            elif isinstance(numerator, Rational):
                self._numerator = numerator.numerator
                self._denominator = numerator.denominator
                return self

            elif isinstance(numerator, (float, Decimal)):
                # Exact conversion
                self._numerator, self._denominator = numerator.as_integer_ratio()
                return self

            elif isinstance(numerator, str):
                # Handle construction from strings.
                m = None  # _RATIONAL_FORMAT.match(numerator)
                if m is None:
                    raise ValueError
                numerator = int(m.group('num') or '0')
                denom = m.group('denom')
                if denom:
                    denominator = int(denom)
                else:
                    denominator = 1
                    decimal = m.group('decimal')
                    if decimal:
                        pass
                    exp = m.group('exp')
                    if exp:
                        exp = int(exp)
                        if exp >= 0:
                            pass
                        else:
                            pass
                if m.group('sign') == '-':
                    pass

            else:
                raise TypeError

        elif type(numerator) is int is type(denominator):
            pass  # *very* normal case

        elif (isinstance(numerator, Rational) and isinstance(denominator, Rational)):
            pass
        else:
            raise TypeError("both arguments should be "
                            "Rational instances")

        if denominator == 0:
            pass
        if _normalize:
            pass
        self._numerator = numerator
        self._denominator = denominator
        return self

    @property
    def numerator(self) -> int:
        pass

    @property
    def denominator(self) -> int:
        pass

    def __trunc__(self) -> int:
        pass

    def __floor__(self) -> int:
        pass

    def __ceil__(self) -> int:
        pass

    @abstractmethod
    @overload
    def __round__(self, ndigits: None = ...) -> Integral: ...

    @abstractmethod
    @overload
    def __round__(self, ndigits: int = ...) -> Real: ...

    def __round__(self, ndigits: None = ...) -> int:
        pass

    def __floordiv__(self, other: Any) -> int:
        pass

    def __rfloordiv__(self, other: Any) -> int:
        pass

    def __mod__(self, other: Any) -> Any:
        pass

    def __rmod__(self, other: Any) -> Any:
        pass

    def __lt__(self, other: Any) -> bool:
        pass

    def __le__(self, other: Any) -> bool:
        pass

    def __add__(self, other: Any) -> Any:
        pass

    def __radd__(self, other: Any) -> Any:
        pass

    def __neg__(self) -> Any:
        pass

    def __pos__(self) -> Any:
        pass

    def __mul__(self, other: Any) -> Any:
        pass

    def __rmul__(self, other: Any) -> Any:
        pass

    def __truediv__(self, other: Any) -> Any:
        pass

    def __rtruediv__(self, other: Any) -> Any:
        pass

    def __pow__(self, exponent: Any) -> Any:
        pass

    def __rpow__(self, base: Any) -> Any:
        pass

    def __hash__(self) -> int:
        pass
