import re

from numbers import *

from math import sqrt, isqrt, prod, gcd
from abc import ABCMeta, abstractmethod
from functools import reduce

from cryptography318.numbers.factor import factor
from cryptography318.core.expr2 import Mul, Add


class Expr(metaclass=ABCMeta):

    # integer value, Real (not float) value(s), types of Real instances
    __slots__ = '_int', '_lazy', '_types'

    @abstractmethod
    def __init__(self): ...

    @property
    @abstractmethod
    def eval(self): ...

    @abstractmethod
    def simplify(self, *, _update=False): ...

    @abstractmethod
    def __str__(self): ...

    @abstractmethod
    def __repr__(self): ...

    @abstractmethod
    def __int__(self): ...

    @abstractmethod
    def __float__(self): ...

    @abstractmethod
    def __complex__(self): ...

    @abstractmethod
    def __abs__(self): ...

    @abstractmethod
    def __pos__(self): ...

    @abstractmethod
    def __neg__(self): ...

    @abstractmethod
    def __eq__(self, other): ...

    @abstractmethod
    def __le__(self, other): ...

    @abstractmethod
    def __lt__(self, other): ...

    @abstractmethod
    def __ge__(self, other): ...

    @abstractmethod
    def __gt__(self, other): ...

    @property
    def numerator(self):
        return +self

    @property
    def denominator(self):
        return 1

    @property
    def integer(self: 'Expr'):
        return self._int

    @property
    def reals(self: 'Expr'):
        return self._lazy

    @property
    def types(self: 'Expr'):
        return self._types


class Mul2(Expr):
    """
    Algebraic expression containing integer value and a list of objects that directly
    or indirectly inherit from Real and must meet the following requirements.

    All Reals must:
    1. not be floats
    2. be able to be squared (i.e. type(x) is type(x * x))
    3. return something when multiplied by an object of same type (i.e. type1 * type1 != NotImplemented
    and does not raise Error)

    Mul objects are communative and associative, so Mul(x, Mul(y, z)) == Mul(x, y, z) == Mul(y, z, x)
    """

    def __init__(self: 'Mul', *exprs, _int=None, _lazy=None):
        """
        Mul(x, y, ..., z) creates Expr object with a list of Real instances and an integer value.
        Mul(_int=x, _lazy=[Reals]) should not be used, as the input is not validated, instead
        an equivalent expression would be Mul(x, *Reals)
        """
        if not exprs:
            self._int = 1 if _int is None else _int
            self._lazy = [] if _lazy is None else _lazy[:]
            self._types = set(type(e) for e in self._lazy)
        elif all(isinstance(e, (Real, Mul)) and not isinstance(e, float) for e in exprs):
            self._int = 1
            self._lazy = []
            self._types = set()
            for e in exprs:
                if isinstance(e, int):
                    self._int *= e
                elif isinstance(e, Mul):
                    self._types.update(e.types)
                    self._lazy += e.reals[:]
                    self._int *= e.integer
                else:
                    self._types.add(type(e))
                    self._lazy.append(e)
        else:
            raise TypeError(f"{self.__class__.__name__} object must be constructed from object inheriting from Real")

    @property
    def eval(self: 'Mul'):
        return self._int * prod(map(lambda e: float(e), self._lazy))

    def equivalent(self: 'Mul', other):
        return sorted(self._lazy) == sorted(other.reals)

    def simplify(self: 'Mul', *, _update=False):
        if len(self._lazy) < 2:
            return None if _update else +self
        elif len(self._lazy) == len(self._types):
            return None if _update else +self
        elif _update:
            res = reduce(lambda i, c: i * Mul(c), self._lazy, Mul(self._int))
            self._int, self._lazy = res.integer, res.reals
        else:
            return reduce(lambda i, c: i * Mul(c), self._lazy, Mul(self._int))

    def __str__(self: 'Mul'):
        self.simplify(_update=True)
        if len(self._lazy) == 0:
            return str(self._int)
        elif self._int == 1:
            if len(self._lazy) == 1:
                return str(self._lazy[0])
            elif len(self._lazy) == 2:
                return f"{self._lazy[0]}*{self._lazy[1]}"
            else:
                printable = f"{self._lazy[0]}"
                for m in self._lazy[1:]:
                    printable += f"*{m}"
                return printable
        elif self._int == -1:
            if len(self._lazy) == 1:
                return f"-{self._lazy[0]}"
            elif len(self._lazy) == 2:
                return f"-{self._lazy[0]}*{self._lazy[1]}"
            else:
                printable = f"-{self._lazy[0]}"
                for m in self._lazy[1:]:
                    printable += f"*{m}"
                return printable
        else:
            if len(self._lazy) == 1:
                return f"{self._int}*{self._lazy[0]}"
            elif len(self._lazy) == 2:
                return f"{self._int}*{self._lazy[0]}*{self._lazy[1]}"
            else:
                printable = f"{self._int}*{self._lazy[0]}"
                for m in self._lazy[1:]:
                    printable += f"*{m}"
                return printable

    def __repr__(self: 'Mul'):
        return f"{self.__class__.__name__}(_int={self._int}, _lazy={self._lazy})"

    def __add__(self: 'Mul', other):
        if isinstance(other, Mul):
            if sorted(self._lazy) == sorted(other.reals):
                return Mul(_int=self._int + other.integer, _lazy=self._lazy)
            else:
                return Add(_lazy=[self, other])
        elif isinstance(other, Add):
            reals = other.reals[:]
            if len(self._lazy) == 0:
                return Add(_int=other.integer + self._int, _lazy=reals)
            elif len(self._lazy) == 1:
                value = self._lazy[0]
                for i, r in enumerate(reals):
                    if isinstance(r, Mul):
                        if r.reals == [value]:
                            _try = r.reals[0] + value
                            if isinstance(_try, (type(value), Mul)):
                                reals[i] = _try
                            else:
                                reals[i] = Mul(2, value)
                            return Add(_int=other.integer, _lazy=reals)
                    elif r == value:
                        reals[i] = Mul(_int=self._int + 1, _lazy=self._lazy)
                        return Add(_int=other.integer, _lazy=reals)
            else:
                for i, r in enumerate(reals):
                    if isinstance(r, Mul) and self.equivalent(r):
                        _try = r + self
                        if isinstance(_try, Add):
                            # if result is an Add, then reals[i] is inside _try.reals so don't overlap
                            del reals[i]
                            return Add(_int=other.integer + _try.integer, _lazy=reals + _try.reals)
                        else:
                            reals[i] = _try
                            return Add(_int=other.integer, _lazy=reals)
            return Add(_int=other.integer, _lazy=reals + self._lazy)
        elif isinstance(other, Real) and not isinstance(other, float):
            return Mul(other) + self
        else:
            return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        return -other + self

    def __rsub__(self, other):
        return other + -self

    def __mul__(self: 'Mul', other):
        if len(self._lazy) == 0:
            if isinstance(other, int):
                return Mul(_int=self._int * other)
            elif isinstance(other, float) and other.is_integer():
                return Mul(_int=self._int * int(other))
            elif isinstance(other, Rational) and other.denominator == 1:
                other = other.numerator
                return Mul(_int=self._int * other)
            elif isinstance(other, Mul):
                return Mul(_int=self._int * other.integer, _lazy=other.reals)
            elif isinstance(other, Add):
                return other * self._int
            elif isinstance(other, Real) and not isinstance(other, float):
                return Mul(_int=self._int, _lazy=[other])
            else:
                return NotImplemented
        else:
            if isinstance(other, int):
                return Mul(_int=self._int * other, _lazy=self._lazy)
            elif isinstance(other, float) and other.is_integer():
                other = int(other)
                return Mul(_int=self._int * other, _lazy=self._lazy)
            elif isinstance(other, Rational) and other.denominator == 1:
                other = other.numerator
                return Mul(_int=self._int * other, _lazy=self._lazy)
            elif isinstance(other, Mul):
                if len(self._lazy) == 1 == len(other.reals):
                    r1, r2 = self._lazy[0], other.reals[0]
                    if (t := type(r1)) is type(r2):
                        _try = r1 * r2
                        if isinstance(_try, int):
                            return Mul(_int=self._int * other.integer * _try)
                        elif isinstance(_try, float) and _try.is_integer():
                            return Mul(_int=self._int * other.integer * int(_try))
                        elif isinstance(_try, Rational) and _try.denominator == 1 and isinstance(_try.numerator, int):
                            return Mul(_int=self._int * other.integer * _try.numerator)
                        # if simplification succeeded but in form of Mul
                        # ex. Sqrt(2) * Sqrt(6) = Sqrt(12) = Sqrt(4) * Sqrt(3) = Mul(2, Sqrt(3))
                        elif isinstance(_try, Mul) and len(_try.reals) == 1:
                            return Mul(_int=self._int * other.integer * _try.integer, _lazy=[*_try.reals])
                        elif isinstance(_try, t):
                            return Mul(_int=self._int * other.integer, _lazy=[_try])
                        # if simplification succeeded but in form of Mul with more than one Real, then
                        # its not actually simplified so there is no case for that
                    return Mul(_int=self._int * other.integer, _lazy=[r1, r2])
                elif len(other.reals) == 1:
                    value = other.reals[0]
                    _type = type(value)
                    _lazy = self._lazy[:]
                    for i, r in enumerate(_lazy):
                        if isinstance(r, _type):
                            _try = r * value

                            # if multiplication resulted in something that can be combined with integer
                            # delete the entry from current lazy and return
                            if isinstance(_try, int):
                                del _lazy[i]
                                return Mul(_int=self._int * other.integer * _try, _lazy=_lazy)
                            elif isinstance(_try, Rational) and \
                                    _try.denominator == 1 and isinstance(_try.numerator, int):
                                del _lazy[i]
                                return Mul(_int=self._int * other.integer * _try.numerator, _lazy=_lazy)
                            elif isinstance(_try, float) and _try.is_integer():
                                del _lazy[i]
                                return Mul(_int=self._int * other.integer * int(_try), _lazy=_lazy)

                            # if multiplication returned same type, add that to list of multiplicands and return
                            elif isinstance(_try, _type):
                                _lazy[i] = _try
                                return Mul(_int=self._int * other.integer, _lazy=_lazy)

                    # if at no point during the iteration could it multiply to something useful, add it as new product
                    return Mul(_int=self._int * other.integer, _lazy=_lazy + [value])
                else:
                    # breaks down multiplication into smaller mults each with len(_lazy) == 1
                    # i * Mul(c) allows __mul__ to be called with len(other.reals) == 1
                    return reduce(
                        lambda i, c: i * Mul(c), other.reals, Mul(_int=self._int * other.integer, _lazy=self._lazy)
                    )
            elif isinstance(other, Add):
                _int = 0
                _lazy = []
                for e in map(lambda e: self * e, other.reals):
                    if isinstance(e, int) or (isinstance(e, Mul) and not e.reals) or\
                            (isinstance(e, Rational) and e.denominator == 1 and isinstance(e.numerator, int)) \
                            or (isinstance(e, float) and e.is_integer()):
                        _int += int(e)
                    else:
                        _lazy.append(e)
                return Add(_int=_int, _lazy=_lazy)
            elif isinstance(other, Real) and not isinstance(other, float):
                return self * Mul(other)
            else:
                return NotImplemented

    def __rmul__(self: 'Mul', other):
        if isinstance(other, int):
            return Mul(_int=self._int * other, _lazy=self._lazy)
        elif isinstance(other, float) and other.is_integer():
            other = int(other)
            return Mul(_int=self._int * other, _lazy=self._lazy)
        elif isinstance(other, Rational) and other.denominator == 1:
            other = other.numerator
            return Mul(_int=self._int * other, _lazy=self._lazy)
        elif isinstance(other, Real) and not isinstance(other, float):

            # switch to __mul__ which is actually implemented for Mul x Mul
            return self * Mul(other)
        else:
            return NotImplemented

    def __floordiv__(self: 'Mul', other):
        if isinstance(other, int):
            if self._int % other == 0:
                return Mul(_int=self._int // other, _lazy=self._lazy)
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Mul):
            if self.equivalent(other):
                d, r = divmod(self._int, other.integer)
                if r == 0:
                    return d
            raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Add):
            if other.integer == 0 and len(other.reals) == 1 and isinstance(other.reals[0], Mul):
                return self // other.reals[0]
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        else:
            return NotImplemented

    def __truediv__(self, other):
        raise TypeError(f"{self.__class__.__name__}.__truediv__() is an unsupported operation. use "
                        f"{self.__class__.__name__}.__floordiv__() instead")

    def __int__(self):
        return int(self.eval)

    def __float__(self):
        return float(self.eval)

    def __complex__(self):
        return complex(self.eval)

    def __abs__(self: 'Mul'):
        if self._int < 0:
            return -self
        else:
            return +self

    def __pos__(self: 'Mul'):
        return Mul(_int=self._int, _lazy=self._lazy)

    def __neg__(self: 'Mul'):
        return Mul(_int=-self._int, _lazy=self._lazy)

    def __hash__(self):
        return hash(self.eval)

    def __eq__(self: 'Mul', other):
        if isinstance(other, Mul):
            return self._int == other.integer and len(self._lazy) == len(other.reals) and \
                   all(e in other.reals for e in self._lazy)
        else:
            return self.eval == other

    def __le__(self, other):
        return self.eval <= other

    def __lt__(self, other):
        return self.eval < other

    def __ge__(self, other):
        return self.eval <= other

    def __gt__(self, other):
        return self.eval < other


class Add2(Expr):

    def __init__(self: 'Add', *exprs, _int=None, _lazy=None):
        if not exprs:
            self._lazy = [] if _lazy is None else _lazy[:]
            self._int = 0 if _int is None else _int
            if self._lazy:
                self._types = set(type(e) for e in self._lazy)
            else:
                self._types = set()
        elif all(isinstance(e, (Real, Expr)) and not isinstance(e, float) for e in exprs):
            self._lazy = []
            self._types = set()
            self._int = 0
            for e in exprs:
                if isinstance(e, int):
                    self._int += e
                elif isinstance(e, Rational) and e.denominator == 1 and isinstance(e.numerator, int):
                    self._int += e.numerator
                elif isinstance(e, float) and e.is_integer():
                    self._int += int(e)
                else:
                    self._lazy.append(e)
                    self._types.add(type(e))
        else:
            raise ValueError(f"invalid arguments for {self.__class__.__name__}")

    @property
    def eval(self: 'Add'):
        return sum(map(lambda sq: float(sq), self._lazy)) + self._int

    def simplify(self: 'Add', *, _update=False):
        if len(self._lazy) == len(self._types):
            return None if _update else +self
        elif _update:
            new = reduce(lambda i, c: i + c.simplify() if isinstance(c, Mul) else c, self._lazy, Add(_int=self._int))
            self._int, self._lazy = new.integer, new.reals
        else:
            return reduce(lambda i, c: i + c.simplify() if isinstance(c, Mul) else c, self._lazy, Add(_int=self._int))

    def __str__(self: 'Add'):
        # self.simplify(_update=True)
        if len(self._lazy) == 0:
            return str(self._int)
        elif self._int > 0:
            printable = f"{self._int}"
            if len(self._lazy) == 1:
                val = self._lazy[0]
                if val < 0:
                    return printable + str(val)
                else:
                    return printable + "+" + str(val)
            else:
                printable += "+" + str(self._lazy[0])
                for val in self._lazy[1:]:
                    if val < 0:
                        printable += str(val)
                    else:
                        printable += "+" + str(val)
                return printable
        else:
            if len(self._lazy) == 1:
                val = self._lazy[0]
                if self._int == 0:
                    return str(val)
                else:
                    return str(val) + str(self._int)
            else:
                printable = str(self._lazy[0])
                for val in self._lazy[1:]:
                    if val < 0:
                        printable += str(val)
                    else:
                        printable += "+" + str(val)
                if self._int == 0:
                    return printable
                else:
                    return printable + str(self._int)

    def __repr__(self: 'Add'):
        return f"{self.__class__.__name__}(_int={self._int}, _lazy={self._lazy})"

    def __add__(self: 'Add', other):
        """
        Combine Add expression with integer or Add or Mul expressions.

        For Mul expressions, the idea is to iterate over the list of Reals in the Add, if at any point
        one of the Reals is equivalent to the Mul (if Mul.reals is just one Real, then this looks like
        direct equivalence, otherwise equivalence is only achieved between two Muls) then increment the
        coefficient value of the Mul expression inside the Add and return it.
        """
        if len(self._lazy) == 0:
            if isinstance(other, int):
                return Add(_int=self._int + other)
            elif isinstance(other, float) and other.is_integer():
                return Add(_int=self._int + int(other))
            elif isinstance(other, Mul) or (isinstance(other, Real) and not isinstance(other, float)):
                return Add(_int=self._int, _lazy=[other])
            elif isinstance(other, Add):
                return Add(_int=self._int + other.integer, _lazy=other.reals)
            elif isinstance(other, Rational) and other.denominator == 1 and isinstance(other.numerator, int):
                return Add(_int=self._int + other.numerator)
            else:
                return NotImplemented
        else:
            if isinstance(other, int):
                return Add(_int=self._int + other, _lazy=self._lazy)
            elif isinstance(other, float) and other.is_integer():
                return Add(_int=self._int + int(other), _lazy=self._lazy)
            elif isinstance(other, Mul):
                reals = self._lazy[:]
                if len(other.reals) == 1:
                    value = other.reals[0]
                    for i, r in enumerate(reals):
                        # if another Mul expr with the same term, increment its coefficient
                        if isinstance(r, Mul) and r.reals == [value]:
                            reals[i] = Mul(_int=r.integer + other.integer, _lazy=r.reals)
                            return Add(_int=self._int, _lazy=reals)
                        # if another Real instance equivalent, increment coeff by 1
                        elif r == value:
                            reals[i] = Mul(_int=1 + other.integer, _lazy=[r])
                            return Add(_int=self._int, _lazy=reals)
                    # base case, if no similar terms, this term just gets added as new term
                    return Add(_int=self._int, _lazy=reals + [other])
                else:
                    for i, r in enumerate(reals):
                        if isinstance(r, Mul) and other.equivalent(r):
                            reals[i] = Mul(_int=r.integer + other.integer, _lazy=[r])
                            return Add(_int=self._int, _lazy=reals)
                    return Add(_int=self._int, _lazy=reals + other.reals)
            elif isinstance(other, Add):
                _expr = Add(self._int, other.integer)
                for e in other.reals:
                    _expr += Mul(e)
                return _expr
            elif isinstance(other, Rational) and other.denominator == 1 and isinstance(other.numerator, int):
                return Add(_int=self._int + other.numerator, _lazy=self._lazy)
            elif isinstance(other, Real) and not isinstance(other, float):
                return self + Mul(other)
            else:
                return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self - other

    def __rsub__(self: 'Add', other):
        return -self + other

    def __mul__(self: 'Add', other):
        if isinstance(other, int):
            _int = self._int * other
            other = Mul(other)
            return Add(_int=_int, _lazy=list(map(lambda e: other * e, self._lazy)))
        elif isinstance(other, Add):
            _int, _lazy = 0, []
            for e in reduce(lambda i, c: i + list(map(lambda e: e * c, other)), self, []):
                if isinstance(e, int):
                    _int += e
                else:
                    _lazy.append(e)
            return Add(_int=_int, _lazy=_lazy)
        elif isinstance(other, Mul):
            return other.__mul__(self)
        elif isinstance(other, float) and other.is_integer():
            _int = self._int * int(other)
            other = Mul(int(other))
            return Add(_int=_int, _lazy=list(map(lambda e: other * e, self._lazy)))
        elif isinstance(other, Rational) and other.denominator == 1 and isinstance(other.numerator, int):
            _int = self._int * other.numerator
            other = Mul(other.numerator)
            return Add(_int=_int, _lazy=list(map(lambda e: other * e, self._lazy)))
        elif isinstance(other, Real) and not isinstance(other, float):
            other = Mul(other)
            _int, _lazy = self._int, []
            for e in map(lambda e: other * e if isinstance(e, Mul) else Mul(e), self._lazy):
                if isinstance(e, Mul) and not e.reals:
                    _int += e.integer
                else:
                    _lazy.append(e)
            return Add(_int=_int, _lazy=_lazy)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __floordiv__(self: 'Add', other):
        if isinstance(other, int):
            if all(isinstance(e, Mul) for e in self._lazy):
                g = gcd(*(e.integer for e in self._lazy), self._int)
            else:
                g = 1
            if g % other == 0:
                return Add(_int=self._int // other, _lazy=list(map(lambda e: e // other, self._lazy)))
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Add):
            quotient = None
            for a, b in zip(sorted(self), sorted(other)):
                e = a // b
                if isinstance(e, int) or (isinstance(e, Mul) and not e.reals):
                    if quotient is not None and quotient != int(e):
                        raise ValueError(f"{self} is not divisible by {other}")
                    else:
                        quotient = int(e)
                else:
                    raise ValueError(f"{self} is not divisible by {other}")
                return quotient
        elif isinstance(other, Mul):
            if len(self.reals) == 1 and isinstance(self.reals[0], Mul):
                return self.reals[0] // other
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        else:
            return NotImplemented

    def __truediv__(self, other):
        return NotImplemented

    def __int__(self):
        return int(self.eval)

    def __float__(self):
        return float(self.eval)

    def __complex__(self):
        return complex(self.eval)

    def __iter__(self: 'Add'):
        if not self._int:
            return iter(self._lazy)
        else:
            return iter([self._int] + self._lazy)

    def __abs__(self):
        if self.eval < 0:
            return -self
        else:
            return +self

    def __pos__(self: 'Add'):
        return Add(_int=self._int, _lazy=self._lazy)

    def __neg__(self: 'Add'):
        return Add(_int=-self._int, _lazy=list(map(lambda e: -e, self._lazy)))

    def __eq__(self, other):
        if isinstance(other, Add):
            return self.eval == other.eval
        else:
            return self.eval == other

    def __lt__(self, other):
        return self.eval < other

    def __le__(self, other):
        return self.eval <= other

    def __gt__(self, other):
        return self.eval > other

    def __ge__(self, other):
        return self.eval >= other
