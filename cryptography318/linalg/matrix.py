from collections import UserList
from numbers import Real
from typing import Iterable
from .vector import Vector


class Matrix(UserList):
    """
    Class representing a list of Vectors whose binary operations are
    matrix operations.
    """

    def __init__(self, initlist=None):
        super().__init__(initlist)

        length = None
        # Make sure we are receiving lists of real numbers
        for i, row in enumerate(self.data):
            if not isinstance(row, Iterable):
                raise TypeError(f"{self.__class__.__name__} only accepts lists of real numbers as list items!")

            # If it is an iterable that isn't already a Vector instance, make it
            if not isinstance(row, Vector):
                self.data[i] = Vector(row)

            # Get length if not set
            if length is None:
                length = len(self.data[i])
            else:
                # Ensure that we aren't accepting a jagged matrix
                if len(self.data[i]) != length:
                    raise ValueError(f"{self.__class__.__name__} cannot be a jagged matrix!")

        self._rows = len(self.data)
        self._cols = length     # This could be None and that's okay

    @property
    def transpose(self):
        return self.__class__(zip(*self.data))

    @property
    def shape(self):
        """
        Shape of matrix as (m, n) where n is optional length of rows

        :return: shape of matrix as tuple
        """
        if self._cols is not None:
            return self._rows, self._cols
        else:
            return self._rows,

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
        if isinstance(other, Matrix):
            if self._cols == other.shape[0]:
                return self.__class__(a * b for a, b in zip(self.data, other.transpose))
            else:
                raise ValueError(f"Matrix multiplication undefined for matrices of shapes "
                                 f"{self.shape} and {other.shape}")
        else:

            # This is NotImplemented since we need to easily check shape and it's more work to just
            # let the other class handle it if it's implemented
            return NotImplemented

    def __rmatmul__(self, other):

        # If vector, do vector-matrix multiplication
        if all(isinstance(e, Real) for e in other):
            return self.__class__(other * e for e in self.transpose)
        else:

            # Otherwise, do matrix multiplication with transpose of self
            return self.__class__(a * b for a, b in zip(other, self.transpose, strict=True))

    # For all methods that involve adding elements to our data list, we need to make sure that they are
    # also lists of Real instances before letting them add. Included are append, extend, insert, and __setitem__
    def __setitem__(self, key, value):
        if isinstance(value, Vector):
            pass
        elif all(isinstance(e, Real) for e in value):
            value = Vector(value)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

        # Make sure our new row matches matrix shape
        if self._cols != len(value):
            raise ValueError(f"List length of {len(value)} cannot be set in {self.__class__.__name__} "
                             f"of shape {self.shape}")
        self.data[key] = value

    def append(self, item):
        if isinstance(item, Vector):
            pass
        elif all(isinstance(e, Real) for e in item):
            item = Vector(item)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

        if self._cols is None:
            self._cols = len(item)
        elif self._cols != len(item):
            raise ValueError(f"List length of {len(item)} cannot be appended to {self.__class__.__name__} "
                             f"of shape {self.shape}")
        self.data.append(item)
        self._rows += 1

    def extend(self, other):

        # Append to values until we have confirmed that all are valid inputs, since extend should either
        # throw an error or have all values extended we don't append unless we can do all
        values = []
        length = self._cols
        n = 0
        for e in other:
            v = Vector(e)
            if length is None:
                length = len(v)
            elif length != len(v):
                raise ValueError(f"List length of {len(v)} cannot be appended to {self.__class__.__name__} "
                                 f"of shape {self.shape}")
            values.append(e)
            n += 1

        self._cols = length
        # Now we can append all values
        self.data.extend(values)
        self._rows += n

    def insert(self, i, item):
        if isinstance(item, Vector):
            pass
        elif all(isinstance(e, Real) for e in item):
            item = Vector(item)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

        if self._cols != len(item):
            raise ValueError(f"List length of {len(item)} cannot be inserted into {self.__class__.__name__} "
                             f"of shape {self.shape}")

        self.data.insert(i, item)
        self._rows += 1
