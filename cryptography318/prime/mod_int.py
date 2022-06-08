from numbers import Integral
from typing import Any
from cryptography318.prime.mod_div_error import ModularDivisionError


class mod_int(Integral):
    """
    Integer where all operations are performed over FF(modulus).
    """

    __slots__ = '_value', '_mod'

    def __new__(cls, value, modulus):
        self = super().__new__(cls)
        self._mod = abs(int(modulus))
        self._value = int(value) % self._mod
        return self

    @property
    def modulus(self):
        """
        Modulus of integer.

        :return: modulus
        """
        return self._mod

    @property
    def value(self):
        """
        Integer value.

        :return: int value
        """
        return self._value

    def inverse(self):
        """
        Modular inverse.

        :raises ValueError: if mod inverse does not exist
        :return: mod inverse if one exists
        """
        return self.__class__(pow(self._value, -1, self._mod), self._mod)

    def __str__(self):
        return f"{self._value}mod{self._mod}"

    def __repr__(self):
        return f"{self.__class__}({self._value}, {self._mod})"

    def __int__(self):
        return self._value

    def __pow__(self, exponent, modulus=None):
        if modulus is None:
            return self.__class__(pow(self._value, exponent, self._mod), self._mod)
        else:
            return self.__class__(pow(self._value, exponent, modulus), self._mod)

    def __lshift__(self, other):
        return self.__class__(self._value << other, self._mod)

    def __rlshift__(self, other):
        return other << self._value

    def __rshift__(self, other):
        return self.__class__(self._value >> other, self._mod)

    def __rrshift__(self, other: Any) -> Any:
        return other >> self._value

    def __and__(self, other):
        return self.__class__(self._value & other, self._mod)

    def __rand__(self, other):
        return self.__class__(other & self._value, self._mod)

    def __xor__(self, other):
        return self.__class__(self._value ^ other, self._mod)

    def __rxor__(self, other):
        return self.__class__(other ^ self._value, self._mod)

    def __or__(self, other):
        return self.__class__(self._value | other, self._mod)

    def __ror__(self, other):
        return self.__class__(other | self._value, self._mod)

    def __invert__(self):
        return self.__class__(~self._value, self._mod)

    def __trunc__(self):
        return int(self)

    def __floor__(self):
        return int(self)

    def __ceil__(self):
        return int(self)

    def __round__(self, ndigits=None):
        return round(self._value, ndigits)

    def __floordiv__(self, other):
        if self._value % other == 0:
            return self.__class__(self._value // other, self._mod)
        else:
            raise ModularDivisionError(self._value, other)

    def __rfloordiv__(self, other):
        return other // self._value

    def __mod__(self, other):
        return self.__class__(self._value % other, self._mod)

    def __rmod__(self, other):
        return other % self._value

    def __lt__(self, other):
        return self._value < other

    def __le__(self, other):
        return self._value <= other

    def __add__(self, other):
        return self.__class__(self._value + other, self._mod)

    def __radd__(self, other):
        return self.__class__(other + self._value, self._mod)

    def __neg__(self):
        return self.__class__(-self._value, self._mod)

    def __pos__(self):
        return self

    def __mul__(self, other):
        return self.__class__(self._value * other, self._mod)

    def __rmul__(self, other):
        return self.__class__(other * self._value, self._mod)

    def __truediv__(self, other):
        return self._value / other

    def __rtruediv__(self, other):
        return other / self._value

    def __rpow__(self, base):
        return pow(base, self._value)

    def __hash__(self):
        return hash(self._value + self._mod)
