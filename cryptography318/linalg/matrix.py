from collections import UserList
from numbers import Real
from typing import Iterable
from .vector import Vector


class Matrix(UserList[Vector]):
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

        self._m = len(self.data)
        self._n = length     # This could be None and that's okay

    @classmethod
    def identity(cls, size):
        """
        Returns an n by n identity matrix.

        :param size: dimension of matrix
        :return: identity matrix
        """
        m = []
        for i in range(size):
            row = [0] * size
            row[i] = 1
            m.append(row)
        return cls(m)

    @property
    def transpose(self):
        return self.__class__(zip(*self.data))

    @property
    def shape(self):
        """
        Shape of matrix as (m, n) where n is optional length of rows

        :return: shape of matrix as tuple
        """
        if self._n is not None:
            return self._m, self._n
        else:
            return self._m,

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
            if self._n == other.shape[0]:

                # Here, we could just use a @ other since __rmatmul__ handles vector-matrix multiplication
                # but that would mean we would be taking the transpose of other once per row, so instead
                # we defined the transpose ahead of time, and run through the multiplication ourselves
                T = other.transpose
                return self.__class__(Vector(a @ e for e in T) for a in self.data)
            else:
                raise ValueError(f"Matrix multiplication undefined for matrices of shapes "
                                 f"{self.shape} and {other.shape}")
        else:

            # This is NotImplemented since we need to easily check shape, and it's more work than just
            # letting the other class handle it if it's implemented
            return NotImplemented

    def __rmatmul__(self, other):

        # If vector, do vector-matrix multiplication
        if all(isinstance(e, Real) for e in other):
            return self.__class__(e.__rmul__(other) for e in self.transpose)
        else:

            # Otherwise, do matrix multiplication with transpose of self
            T = self.transpose
            return self.__class__(Vector(a @ e for e in T) for a in other)

    # For all methods that involve adding elements to our data list, we need to make sure that they are
    # also lists of Real instances before letting them add. Included are append, extend, insert, and __setitem__
    def __setitem__(self, key, value):
        if isinstance(key, slice):

            # If it's a slice, vet all inputs before adding to index
            values = []
            for v in value:
                if not isinstance(v, Iterable) or not all(isinstance(e, Real) for e in v):
                    raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")
                elif not isinstance(v, Vector):
                    v = Vector(v)
                if len(v) != self._n:
                    raise ValueError(f"List length of {len(v)} cannot be set in {self.__class__.__name__} "
                                     f"of shape {self.shape}")
                values.append(v)
            self.data[key] = values
        else:
            if isinstance(value, Vector):
                pass
            elif all(isinstance(e, Real) for e in value):
                value = Vector(value)
            else:
                raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

            # Make sure our new row matches matrix shape
            if self._n != len(value):
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

        if self._n is None:
            self._n = len(item)
        elif self._n != len(item):
            raise ValueError(f"List length of {len(item)} cannot be appended to {self.__class__.__name__} "
                             f"of shape {self.shape}")
        self.data.append(item)
        self._m += 1

    def extend(self, other):

        # Append to values until we have confirmed that all are valid inputs, since extend should either
        # throw an error or have all values extended we don't append unless we can do all
        values = []
        length = self._n
        n = 0
        for e in other:
            v = e if isinstance(e, Vector) else Vector(e)
            if length is None:
                length = len(v)
            elif length != len(v):
                raise ValueError(f"List length of {len(v)} cannot be appended to {self.__class__.__name__} "
                                 f"of shape {self.shape}")
            values.append(v)
            n += 1

        self._n = length
        # Now we can append all values
        self.data.extend(values)
        self._m += n

    def insert(self, i, item):
        if isinstance(item, Vector):
            pass
        elif all(isinstance(e, Real) for e in item):
            item = Vector(item)
        else:
            raise TypeError(f"{self.__class__.__name__} only accepts real numbers as list items!")

        if self._n != len(item):
            raise ValueError(f"List length of {len(item)} cannot be inserted into {self.__class__.__name__} "
                             f"of shape {self.shape}")

        self.data.insert(i, item)
        self._m += 1

    # Now we can get into the actual matrix-specific methods
    def slice(self, s):
        """
        Slices matrix vertically, applying an integer or slice index to all rows.

        :param s: index
        :return: matrix with index applied to each row
        """
        if isinstance(s, slice):
            return self.__class__(r[s] for r in self.data)
        else:
            return self.__class__(Vector([r[s]]) for r in self.data)

    def split(self, i):
        """
        Splits matrix by an index, returning the partitioned matrix as a tuple.

        :param i: index
        :return: left and right partitions
        """
        lm, rm = [], []
        for row in self.data:
            lm.append(row[:i])
            rm.append(row[i:])
        return self.__class__(lm), self.__class__(rm)

    def rref(self, offset=0):
        """
        Computes the reduced row-echelon form of the matrix. If offset is provided,
        computes the RREF ignoring the last offset columns.

        :param self: matrix
        :param offset: column index offset from last column
        :return: matrix in reduced row-echelon form
        """

        # If there are no rows, return new empty matrix
        if self._n is None:
            return Matrix()

        pivot_row = 0  # First pivot belongs in first row

        array = self.copy()

        for j in range(self._n - offset):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, self._m):

                # If non-zero element, this row can become pivot row
                if array[i][j] != 0:

                    # Make j'th element the pivot, reducing rest of row as well
                    array[i] /= array[i][j]
                    if i > pivot_row:  # If pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row], array[i]

                    # Row reduce everything else
                    for k in range(self._m):
                        if k != pivot_row:
                            array[k] -= array[pivot_row] * array[k][j]

                    # Increment pivot row and stop looking in this column for pivot, move next
                    pivot_row += 1
                    break

        return array
