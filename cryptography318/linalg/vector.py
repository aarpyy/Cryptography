from collections import UserList
from numbers import Real
from typing import Iterable


class Vector(UserList[Real]):
    """
    Class representing a list of real numbers whose binary operations are
    vector operations.
    """

    def __init__(self, initlist=None):

        super().__init__(initlist)

        # Make sure we are receiving real numbers, do this after constructor so iterable is not consumed
        if not all(isinstance(x, Real) for x in self.data):
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

    def __add__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(a + b for a, b in zip(self.data, other, strict=True))
        else:
            return self.__class__(e + other for e in self.data)

    def __radd__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(b + a for a, b in zip(other, self.data, strict=True))
        else:
            return self.__class__(other + e for e in self.data)

    def __iadd__(self, other):
        if isinstance(other, Iterable):
            for i, e in enumerate(other):
                self.data[i] += e
        else:
            for i in range(len(self.data)):
                self.data[i] += other
        return self

    def __sub__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(a - b for a, b in zip(self.data, other, strict=True))
        else:
            return self.__class__(e - other for e in self.data)

    def __rsub__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(b - a for a, b in zip(other, self.data, strict=True))
        else:
            return self.__class__(other - e for e in self.data)

    def __isub__(self, other):
        if isinstance(other, Iterable):
            for i, e in enumerate(other):
                self.data[i] -= e
        else:
            for i in range(len(self.data)):
                self.data[i] -= other
        return self

    def __mul__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(a * b for a, b in zip(self.data, other, strict=True))
        else:
            return self.__class__(e * other for e in self.data)

    def __rmul__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(b * a for a, b in zip(other, self.data, strict=True))
        else:
            return self.__class__(other * e for e in self.data)

    def __imul__(self, other):
        if isinstance(other, Iterable):
            for i, e in enumerate(other):
                self.data[i] *= e
        else:
            for i in range(len(self.data)):
                self.data[i] *= other
        return self

    def __truediv__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(a / b for a, b in zip(self.data, other, strict=True))
        else:
            return self.__class__(e / other for e in self.data)

    def __floordiv__(self, other):
        if isinstance(other, Iterable):
            return self.__class__(a // b for a, b in zip(self.data, other, strict=True))
        else:
            return self.__class__(e // other for e in self.data)

    def __mod__(self, other):
        return self.__class__(e % other for e in self.data)

    def __matmul__(self, other):
        """
        Computes the dot product with another iterable.

        :param other: iterable of numbers
        :return: dot product of self and other
        """

        # If real-valued list, return dot product
        if all(isinstance(e, Real) for e in other):
            return sum(a * b for a, b in zip(self.data, other, strict=True))
        else:

            # Otherwise, we don't know how to work with it. If it's a Matrix instance, it's __rmatmul__ will be called
            return NotImplemented

    # Since all that our matmul is able to do it dot product is it commutative
    __rmatmul__ = __matmul__

    # For all methods that involve adding elements to our data list, we need to make sure that they are
    # also Real instances before letting them add. Included are append, extend, insert, and __setitem__
    def __setitem__(self, key, value):
        if isinstance(value, Real):
            self.data[key] = value
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

    def append(self, item):
        if isinstance(item, Real):
            self.data.append(item)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

    def extend(self, other):
        if all(isinstance(e, Real) for e in other):
            self.data.extend(other)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

    def insert(self, i, item):
        if isinstance(item, Real):
            self.data.insert(i, item)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")
