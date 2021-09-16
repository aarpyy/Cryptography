# much of this file is based off of fractions.py including the regex expressions, parsing of fraction and
# square root inputs, operator overriders, and wherever else noted

import re
import operator

from typing import *
from numbers import *

from math import gcd, isnan, isinf
from decimal import Decimal

from fractions import Fraction as PyFraction

from .expr import Expr, Add, Sqrt, Mul

__all__ = ['Fraction']

# credit to fractions for the base of this regex
_FRACTION_SQRT_FORMAT = re.compile(r"""
        \A\s*                       # optional whitespace at the start
        (?P<sign>[-+]?)             # an optional sign, then
        (?P<coeff>\d*)              # sqrt coefficient, optional
        (?:sqrt\((?P<num>\d+))?     # sqrt expression, required
        \)                          # termination of sqrt expression
        \s*                         # optional whitespace between numerator and '/'
        (?:/\s*(?P<denom>\d+))?     # optional denominator (and optional whitespace between '/' and denom)
    """, re.VERBOSE | re.IGNORECASE)

# credit to fractions for this regex
_FRACTION_FORMAT = re.compile(r"""
    \A\s*                           # optional whitespace at the start, then
    (?P<sign>[-+]?)                 # an optional sign, then
    (?=\d|\.\d)                     # lookahead for digit or .digit
    (?P<num>\d*)                    # numerator (possibly empty)
    \s*                             # optional whitespace between numerator and '/'
    (?:                             # followed by
       (?:/\s*(?P<denom>\d+))?      # an optional denominator (and optional whitespace between '/' and denom)
    |                               # or
       (?:\.(?P<decimal>\d*))?      # an optional fractional part 
       (?:E(?P<exp>[-+]?\d+))?      # and optional exponent
    )
    \s*\Z                           # and optional whitespace to finish
""", re.VERBOSE | re.IGNORECASE)


class Fraction(Rational):

    __slots__ = '_numerator', '_denominator'

    # the , *, allows for numerator and denominator to be passed as positional arguments, but anything past
    # that (i.e. _normalize) must be passed explicitly as a kwarg
    def __new__(cls, numerator, denominator=None, *, _normalize=True):

        self = super().__new__(cls)

        if denominator is None:
            if isinstance(numerator, (int, Expr)):
                self._numerator, self._denominator = numerator, 1
                return self
            elif isinstance(numerator, Rational):
                self._numerator, self._denominator = numerator.numerator, numerator.denominator
                return self
            elif isinstance(numerator, (float, Decimal)):
                self._numerator, self._denominator = numerator.as_integer_ratio()  # built in conversion
                return self
            elif isinstance(numerator, str):
                if 'sqrt' in numerator:
                    m = _FRACTION_SQRT_FORMAT.match(numerator)
                    if m is None:
                        raise ValueError(f"invalid string literal for {__class__.__name__}: {numerator}")
                    coeff, radicand = int(m.group('coeff') or '1'), int(m.group('num'))
                    denominator = int(m.group('denom') or '1')
                    numerator = Add(Sqrt(radicand, coeff))
                else:
                    m = _FRACTION_FORMAT.match(numerator)
                    if m is None:
                        raise ValueError(f"invalid string literal for {__class__.__name__}: {numerator}")
                    numerator = int(m.group('num') or '0')
                    denom = m.group('denom')
                    if denom:
                        denominator = int(denom)
                    else:
                        denominator = 1
                        decimal = m.group('decimal')
                        if decimal:
                            # shift numerator by length of decimal then add decimal (ex. 1.5, n=1, dec=5, n=n*10+5=15
                            shift = pow(10, len(decimal))
                            numerator = numerator * shift + int(decimal)
                            denominator *= shift
                        exp = m.group('exp')
                        if exp:
                            exp = int(exp)
                            if exp >= 0:
                                numerator *= pow(10, exp)
                            else:
                                denominator *= pow(10, -exp)
                if m.group('sign') == '-':
                    numerator = -numerator
            else:
                raise TypeError(f"all input(s) must be rational not type: {type(numerator)}")
        elif isinstance(numerator, (int, Expr)) and isinstance(denominator, int):
            pass
        elif isinstance(numerator, Rational) and isinstance(denominator, Rational):
            self._numerator = numerator.numerator * denominator.denominator
            self._denominator = numerator.denominator * denominator.numerator
            return self
        else:
            try:
                return numerator / denominator
            except TypeError:
                raise TypeError(f"all input(s) must be rational not type(s): {type(numerator)}, {type(denominator)}")

        if denominator == 0:
            raise ZeroDivisionError(f"fraction invalid with denominator equal to zero")
        if _normalize:
            if isinstance(numerator, Add):
                g = 1
            elif isinstance(numerator, Mul):
                g = gcd(numerator.integer, denominator)
            else:
                g = gcd(numerator, denominator)
            if denominator < 0:
                g = -g
            numerator //= g
            denominator //= g
        self._numerator = numerator
        self._denominator = denominator
        return self

    # all methods below that containt type hints: 'Fraction' throw unresolved reference to _ attribute warnings
    # if no hint is given (because they reference private variables _numerator or _denominator)
    # the hint must also be given in the .py file, not .pyi
    @property
    def numerator(self: 'Fraction') -> int:
        return self._numerator

    @property
    def denominator(self: 'Fraction') -> int:
        return self._denominator

    @property
    def reciprocal(self: 'Fraction') -> 'Fraction':
        if self._numerator < 0:
            return Fraction(-self._denominator, -self._numerator, _normalize=False)
        else:
            return Fraction(self._denominator, self._numerator, _normalize=False)

    def __str__(self: 'Fraction'):
        if self._denominator == 1:
            return str(self._numerator)
        elif isinstance(self._numerator, Expr):
            return f"({self._numerator})/{self._denominator}"
        else:
            return f"{self._numerator}/{self._denominator}"

    def __repr__(self: 'Fraction'):
        return f'{self.__class__.__name__}({self._numerator}, {self._denominator})'

    def _operator_fallbakcs(fraction_operator: Callable, std_operator):

        def forward(a, b):
            if isinstance(b, (int, Fraction, PyFraction, Expr)):
                return fraction_operator(a, b)
            elif isinstance(b, float):
                if b.is_integer():
                    return fraction_operator(a, int(b))
                return std_operator(float(a), b)
            elif isinstance(b, complex):
                return std_operator(complex(a), b)
            else:
                return NotImplemented

        forward.__name__ = '__' + fraction_operator.__name__ + '__'
        forward.__doc__ = std_operator.__doc__

        def reverse(b, a):
            if isinstance(b, (Rational, Expr)):
                return fraction_operator(a, b)
            elif isinstance(b, float):
                if b.is_integer():
                    return fraction_operator(a, int(b))
                return std_operator(float(a), b)
            elif isinstance(b, Real):
                return std_operator(float(a), float(b))
            elif isinstance(b, Complex):
                return std_operator(complex(a), complex(b))
            else:
                return std_operator(float(a), b)

        reverse.__name__ = '__r' + fraction_operator.__name__ + '__'
        reverse.__doc__ = std_operator.__doc__

        return forward, reverse

    def _add(a, b):
        da, db = a.denominator, b.denominator
        return Fraction(a.numerator * db + b.numerator * da, da * db)

    __add__, __radd__ = _operator_fallbakcs(_add, operator.add)

    def _sub(a, b):
        da, db = a.denominator, b.denominator
        return Fraction(a.numerator * db - b.numerator * da, da * db)

    __sub__, __rsub__ = _operator_fallbakcs(_sub, operator.sub)

    def _mul(a, b):
        return Fraction(a.numerator * b.numerator, a.denominator * b.denominator)

    __mul__, __rmul__ = _operator_fallbakcs(_mul, operator.mul)

    def _floordiv(a, b):
        return (a.numerator * b.denominator) // (a.denominator * b.numerator)

    __floordiv__, __rfloordiv__ = _operator_fallbakcs(_floordiv, operator.floordiv)

    def _truediv(a, b):
        return Fraction(a.numerator * b.denominator, a.denominator * b.numerator)

    __truediv__, __rtruediv__ = _operator_fallbakcs(_truediv, operator.truediv)

    def _mod(a, b):
        da, db = a.denominator, b.denominator
        return Fraction((a.numerator * db) % (da * b.numerator), da * db)

    __mod__, __rmod__ = _operator_fallbakcs(_mod, operator.mod)

    def _divmod(a, b):
        da, db = a.denominator, b.denominator
        d, r = divmod(a.numerator * db, da * b.numerator)
        return d, Fraction(r, da * db)

    __divmod__, __rdivmod__ = _operator_fallbakcs(_divmod, divmod)

    def _pow(a, b):
        if b.denominator == 1 and isinstance(b.numerator, int):
            power = b.numerator
            if power >= 0:
                return Fraction(pow(a.numerator, power), pow(a.denominator, power), _normalize=False)
            elif a.numerator >= 0:
                return Fraction(pow(a.denominator, -power), pow(a.numerator, -power), _normalize=False)
            else:
                return Fraction(pow(-a.denominator, -power), pow(-a.numerator, -power), _normalize=False)
        else:
            return pow(float(a), float(b))

    def __pow__(self: 'Fraction', power):
        # if power is integer
        if isinstance(power, Rational):
            return self._pow(power)
        elif isinstance(power, Real):
            return pow(float(self), float(power))
        else:
            return NotImplemented

    def __rpow__(self, base):
        if self.denominator == 1 and isinstance(self.numerator, int) and self.numerator >= 0:
            return pow(base, self.numerator)
        elif isinstance(base, Rational):
            return Fraction(base.numerator, base.denominator, _normalize=False)._pow(self)
        else:
            return pow(float(base), float(self))

    def __int__(self: 'Fraction'):
        if isinstance(self._numerator, Expr):
            return int(self._numerator.eval / self._denominator)
        else:
            return int(self._numerator / self._denominator)

    def __float__(self: 'Fraction'):
        if isinstance(self._numerator, Expr):
            return float(self._numerator.eval / self._denominator)
        else:
            return float(self._numerator / self._denominator)

    def __complex__(self: 'Fraction'):
        if isinstance(self._numerator, Expr):
            return complex(self._numerator.eval / self._denominator)
        else:
            return complex(self._numerator / self._denominator)

    def __bool__(self: 'Fraction'):
        if isinstance(self._numerator, Expr):
            return bool(self._numerator.eval / self._denominator)
        else:
            return bool(self._numerator / self._denominator)

    def __trunc__(self: 'Fraction'):
        if self._numerator < 0:
            return -(-self._numerator // self._denominator)
        return self._numerator // self._denominator

    def __floor__(self: 'Fraction'):
        return self._numerator // self._denominator

    def __ceil__(self: 'Fraction'):
        return -(-self._numerator // self._denominator)

    def __round__(self, ndigits=None):
        if ndigits is None:
            return round(float(self))
        shift = pow(10, abs(ndigits))
        if ndigits < 0:
            return Fraction(round(self / shift) * shift)
        return Fraction(round(self * shift), shift)

    def __abs__(self: 'Fraction'):
        return Fraction(abs(self._numerator), self._denominator, _normalize=False)

    def __pos__(self: 'Fraction'):
        return Fraction(self._numerator, self._denominator, _normalize=False)

    def __neg__(self: 'Fraction'):
        return Fraction(-self._numerator, self._denominator, _normalize=False)

    def __hash__(self):
        return hash(float(self))

    def __eq__(self: 'Fraction', other):
        if isinstance(other, int):
            return self._numerator == other and self._denominator == 1
        elif isinstance(other, (Rational, Add)):
            return self._numerator == other.numerator and self._denominator == other.denominator
        if isinstance(other, Complex) and other.imag == 0:
            other = other.real
        if isinstance(other, float):
            if isnan(other) or isinf(other):

                # following recommendations from fractions that all finite values should compare to inf/nan as float(0)
                return 0.0 == other
            return float(self) == other
        else:
            return NotImplemented

    def __ne__(self: 'Fraction', other):
        if isinstance(other, int):
            return not (self._numerator == other or self._denominator == 1)
        elif isinstance(other, (Rational, Add)):
            return not (self._numerator == other.numerator or self._denominator == other.denominator)
        if isinstance(other, Complex) and other.imag == 0:
            other = other.real
        if isinstance(other, float):
            if isnan(other) or isinf(other):

                # following recommendations from fractions that all finite values should compare to inf/nan as float(0)
                return 0.0 != other
            else:
                return self != Fraction(*other.as_integer_ratio())
        else:
            return NotImplemented

    def _compare(a, b, op):
        if isinstance(b, (Rational, Add)):
            return op(a.numerator * b.denominator, a.denominator * b.numerator)
        elif isinstance(b, Complex):
            b = b.real
        if isinstance(b, float):
            if isnan(b) or isinf(b):
                return op(0.0, b)
            else:
                return op(a, Fraction(*b.as_integer_ratio()))
        else:
            return NotImplemented

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)
