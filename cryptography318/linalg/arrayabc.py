from abc import abstractmethod, ABCMeta
from collections import Sequence, MutableSequence
from functools import reduce
from numbers import *
from fractions import Fraction
from typing import Union, Callable
from sympy.core.expr import Expr

import operator


def get_shape(o):
    if hasattr(o, '__len__') and hasattr(o, '__getitem__'):
        return len(o), *get_shape(o[0])
    else:
        return ()


class DimensionError(Exception):

    def __init__(self, a, *b, op=None):
        if isinstance(op, Callable):
            self._msg = f"{op.__name__} incompatible with objects of shape(s): {self.get_shape(a)}"
        else:
            self._msg = f"operation incompatible with objects of shape(s): {self.get_shape(a)}"
        if b:
            self._msg += ", " + ", ".join(str(self.get_shape(e)) for e in b)

    def __str__(self):
        return self._msg

    @staticmethod
    def get_shape(o):
        if hasattr(o, '__len__') and hasattr(o, '__getitem__'):
            return len(o), *get_shape(o[0])
        else:
            return ()


# TypeResolutionOrder
TRO = {Number: 0, Complex: 1, Real: 2, Rational: 3, Integral: 4}


def get_dtype(a: Union[type, tuple], *b):
    """
    Determines the Type Resolution Order (TRO) of array objects.
    Takes in two (or more) types and returns one type that should be the dtype of the array.
    This TRO function ensures that an array can be re-constructed using a function similar to
    map(lambda e: dtype(e), array), where dtype is the returned type. This is guaranteed since if any complex
    values in array, dtype will be complex, otherwise dtype will be a type that can take in basically any
    type of Number (float accepts all Reals, Fraction accepts all Rationals, and int accepts
    basically everything)

    :param a: either tuple of two (or more) types or one type (if b is second type)
    :param b: all remaining types
    """
    if b:
        if isinstance(a, tuple):
            a = (*a, *b)
        else:
            a = (a, *b)
    # if b not provided and a is a type (not a tuple) then only one type given, return it
    elif isinstance(a, type):
        return a

    return reduce(lambda r, c: r if issubclass(c, Expr) or TRO[get_ntype(c)] >= TRO[get_ntype(r)] else c, a)


def get_ntype(n):
    res = reduce(lambda r, c: r if r or not issubclass(n, c) else c, (Integral, Rational, Real, Complex, Number), 0)
    if res:
        return res
    else:
        raise TypeError(f"{n} is not a Number")


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


class ABCArray(MutableSequence, metaclass=ABCMeta):
    __slots__ = '_array', '_dtype'

    @abstractmethod
    def __init__(self): ...

    @abstractmethod
    def __str__(self): ...

    @abstractmethod
    def __repr__(self): ...

    @abstractmethod
    def __iter__(self): ...

    @abstractmethod
    def __setitem__(self, key, value): ...

    @abstractmethod
    def __add__(self, other): ...

    @abstractmethod
    def __radd__(self, other): ...

    @abstractmethod
    def __sub__(self, other): ...

    @abstractmethod
    def __rsub__(self, other): ...

    # __mul__, __rmul__ for arrays are dot products when isinstance(other, Sequence)
    @abstractmethod
    def __mul__(self, other): ...

    @abstractmethod
    def __rmul__(self, other): ...

    @abstractmethod
    def __eq__(self, other): ...

    @abstractmethod
    def _compare(a, b, op): ...

    @abstractmethod
    def __pos__(self): ...

    @abstractmethod
    def __neg__(self): ...

    @abstractmethod
    def append(self, value): ...

    @abstractmethod
    def reverse(self): ...

    @abstractmethod
    def extend(self, values): ...

    @abstractmethod
    def pop(self, index=None): ...

    @abstractmethod
    def remove(self, value): ...

    @abstractmethod
    def copy(self): ...

    @abstractmethod
    def dot(self, other): ...

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)

    @property
    def array(self):
        return self._array

    @property
    def dtype(self):
        return self._dtype
