from .arrayabc import DimensionError, get_dtype
from .array import binarray, rarray, marray
from abc import ABCMeta, abstractmethod
from typing import overload, Sequence, Union, Optional, Iterable
from collections import MutableSequence
from numbers import *
from numpy.linalg import inv as np_inv
from numpy.linalg import det as np_det
from numpy import asarray
from sympy import solve, Symbol
from sympy.core.expr import Expr
from fractions import Fraction
from math import gcd, inf
from functools import reduce
import operator


def where(iterable):
    return tuple(map(lambda arr: list(map(lambda n: int(n), arr)), asarray(iterable).nonzero()))


class ABCMatrix(MutableSequence, metaclass=ABCMeta):
    __slots__ = '_array', '_dtype'

    @abstractmethod
    def __init__(self):
        ...

    def __str__(self):
        str_array = []
        max_len = 0
        for i in range(len(self)):
            str_array.append([])
            for j in range(len(self[0])):
                n = self[i][j]
                if isinstance(n, Integral) or (isinstance(n, float) and n.is_integer()):
                    string = str(int(n))
                elif isinstance(n, Fraction):
                    string = str(n)
                elif isinstance(n, Real):
                    string = str(round(float(n), 3))
                else:
                    string = str(n)
                if len(string) > max_len:
                    max_len = len(string)
                str_array[i].append(string)

        padding = (max_len + 1) | 1
        formatted = "["
        for i in range(l := len(str_array)):
            if i == 0:
                formatted += "["
            else:
                formatted += " ["
            for j in range(len(str_array[0])):
                e = str_array[i][j]
                pad_left = (padding - len(e)) // 2
                pad_right = padding - len(e) - pad_left

                if j == 0 and e[0] != '-':  # sets numbers back from left [ to not squish, doesn't with negatives
                    formatted += max(pad_left, 1) * " " + f"{e}" + " " * pad_right
                else:
                    formatted += pad_left * " " + f"{e}" + " " * pad_right
            if i == l - 1:
                formatted += "]"
            else:
                formatted += "]\n"
        return formatted + "]"

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array})"

    def __iter__(self):
        return iter(self._array)

    @abstractmethod
    def __setitem__(self, key, value):
        ...

    def __contains__(self, item):
        return item in self.flat

    @abstractmethod
    def _add(a, b):
        ...

    @abstractmethod
    def _sub(a, b):
        ...

    @abstractmethod
    def _mul(a, b):
        ...  # _mul not static since c*a == a*c for scalar c

    @staticmethod
    @overload
    @abstractmethod
    def _matmul(a, b):
        ...  # _matmul static since a*b != b*a

    @staticmethod
    @overload
    @abstractmethod
    def _matmul(a, b, m):
        ...  # overloaded _matmul w/ mod

    @staticmethod
    @abstractmethod
    def _matmul(*args):
        ...

    @staticmethod
    @overload
    @abstractmethod
    def _dot(a, b):
        ...

    @staticmethod
    @overload
    @abstractmethod
    def _dot(a, b, m):
        ...

    @staticmethod
    @abstractmethod
    def _dot(*args):
        ...

    # matrix addition is commutative
    def __add__(self, other):
        return self._add(other)

    __radd__ = __add__

    @abstractmethod
    def __iadd__(self, other):
        ...

    # subtraction is (w/ negation) commutative
    def __sub__(self, other):
        return self._sub(other)

    def __rsub__(self, other):
        return -self + other

    # multiplication against a scalar is commutative, matrix-multiplication is not
    def __mul__(self, other):
        if isinstance(other, Sequence):
            return self._matmul(self, other)
        else:
            return self._mul(other)

    def __rmul__(self, other):
        if isinstance(other, Sequence):
            return self._matmul(other, self)
        else:
            return self._mul(other)

    def __eq__(self, other):
        if isinstance(other, Number):
            return bmatrix([r == other for r in self])
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return bmatrix([x == y for x, y in zip(self, other)])
            else:
                raise DimensionError(self, other, op=operator.eq)
        else:
            return NotImplemented

    @abstractmethod
    def __pos__(self):
        ...

    @abstractmethod
    def __neg__(self):
        ...

    @abstractmethod
    def append(self, value):
        ...

    def reverse(self):
        self._array.reverse()

    def __reversed__(self):
        i = len(self) - 1
        while i > -1:
            yield self._array[i]
            i -= 1

    @abstractmethod
    def extend(self, values):
        ...

    def pop(self, index=None):
        if index is None:
            return self._array.pop()
        else:
            return self._array.pop(index)

    def remove(self, value):
        self._array.remove(value)

    def clear(self):
        self._array.clear()

    def count(self, value):
        return self._array.count(value)

    def index(self, value, start=None, stop=None):
        for i, row in self._array:
            if value in row:
                return i, row.index(value)
        raise ValueError(f"{value} is not in matrix")

    @abstractmethod
    def copy(self):
        ...

    @abstractmethod
    def transpose(self):
        ...

    @abstractmethod
    def invert(self):
        ...

    @abstractmethod
    def _compare(a, b, op):
        ...

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)

    def __len__(self):
        return len(self._array)

    def __delitem__(self, key):
        del self._array[key]

    @property
    def dtype(self):
        return self._dtype

    @property
    def array(self):
        return self._array

    @property
    def aslist(self):
        return list([float(e) for e in r] for r in self._array)

    @property
    def T(self):
        return tuple(zip(*self._array))

    @property
    def flat(self):
        return sum(self, start=[])

    @staticmethod
    @abstractmethod
    def identity(size): ...


class bmatrix(ABCMatrix):

    def __init__(self, array):
        super().__init__()
        if isinstance(array, Sequence):
            if all(isinstance(r, binarray) for r in array):
                self._array = array
            elif all(isinstance(r, Sequence) and all(isinstance(e, Integral) for e in r) for r in array):
                self._array = list(binarray(list(r)) for r in array)
            else:
                raise TypeError(f"{self.__class__.__name__} must be constructed from Sequence of Integral values")
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from Sequence")

    @staticmethod
    def identity(size):
        array = []
        row = [0] * size
        for i in range(size):
            id_row = binarray(row[:])
            id_row[i] = 1
            array.append(id_row)
        return bmatrix(array)

    def _add(a, b):
        if isinstance(b, Integral):
            if b & 1:
                return a.complement
            else:
                width = len(a[0])
                return bmatrix([binarray([0] * width) for _ in range(len(a))])
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bmatrix(x + y for x, y in zip(a, b))
            else:
                raise DimensionError(a, b, op=operator.add)
        else:
            return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, Integral):
            if other & 1:
                return self.complement
            else:
                return self
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                for i in range(len(self)):
                    self._array[i] += other[i]
                return self
            else:
                raise DimensionError(self, other, op=operator.iadd)
        else:
            return NotImplemented

    _sub = _add

    def _mul(a, b):
        if isinstance(b, Integral):
            if b & 1:
                return a.copy()
            else:
                width = len(a[0])
                return bmatrix([binarray([0] * width) for _ in range(len(a))])
        else:
            return NotImplemented

    @staticmethod
    def _matmul(a, b):
        T = tuple(zip(*b))
        return bmatrix([binarray([bmatrix._dot(x, y) for y in T]) for x in a])

    @staticmethod
    def _dot(a, b):
        return sum(y for x, y in zip(a, b) if x)

    def __setitem__(self, key, value):
        if isinstance(value, Integral):
            self._array[key] = value & 1
        else:
            raise TypeError(f"value must be Integral not {type(value)}")

    def __pos__(self):
        return bmatrix(+r for r in self)

    def __neg__(self):
        return bmatrix(-r for r in self)

    def append(self, value):
        if isinstance(value, Integral):
            self._array.append(value & 1)
        else:
            raise NotImplementedError

    def extend(self, values):
        if all(isinstance(e, Integral) for e in values):
            self._array.extend(e & 1 for e in values)
        else:
            raise NotImplementedError

    def copy(self):
        return bmatrix(r[:] for r in self)

    def transpose(self):
        return bmatrix(list(map(lambda *args: binarray([*args]), *self._array)))

    def invert(self):
        return bmatrix(np_inv(self._array).tolist())

    def _compare(a, b, op):
        if isinstance(b, Integral):
            return bmatrix([op(r, b) for r in a])
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bmatrix([op(x, y) for x, y in zip(a, b)])
            else:
                raise DimensionError(a, b, op=op)
        else:
            return NotImplemented

    def insert(self, index, value):
        if isinstance(value, Integral):
            self._array.insert(index, value & 1)
        else:
            raise NotImplementedError

    def __getitem__(self, item):
        if isinstance(item, slice):
            return bmatrix(self._array[item])
        else:
            return self._array[item]

    @property
    def complement(self):
        return bmatrix([e.complement for e in self._array])


class rmatrix(ABCMatrix):

    __slots__ = '_augmented',

    def __init__(self, array, *, aug=False):
        super().__init__()
        if isinstance(array, rmatrix):
            self._array = array.array
        elif isinstance(array, Iterable):
            array = list(map(rarray, array))
            if not reduce(lambda r, c: r if r and len(r) == len(c) else 0, array):
                raise ValueError(f"all rows must be of same length")
        elif isinstance(array, Real):
            self._array = [rarray([float(array)])]
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from 2-dimensional Sequence of Reals")

        _types = set(v.dtype for v in self._array)
        self._dtype = get_dtype(*_types) if self._array else int
        for r in self._array:
            r.astype(self._dtype, _update=True)
        self._augmented = aug

    @property
    def augmented(self):
        return self._augmented

    @staticmethod
    def identity(size):
        array = []
        row = [0] * size
        for i in range(size):
            id_row = rarray(row[:])
            id_row[i] = 1
            array.append(id_row)
        return rmatrix(array)

    def astype(self, dtype, *, _update=False):
        if _update:
            for r in self._array:
                r.astype(self._dtype, _update=True)
        else:
            return rmatrix([r.astype(dtype) for r in self])

    def __setitem__(self, key, value):
        if isinstance(value, rarray):
            if value.dtype is self._dtype:
                self._array[key] = value
            else:
                self[key] = list(value)
        elif isinstance(value, Sequence):
            self._array[key] = rarray(list(value)).astype(self._dtype)
        else:
            raise NotImplementedError

    def _add(a, b):
        if isinstance(b, Sequence):
            if len(a) == len(b):
                return rmatrix([x + y for x, y in zip(a, b)])
            else:
                raise DimensionError(a, b, op=operator.add)
        elif isinstance(b, Real):
            return rmatrix([r + b for r in a])
        else:
            return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, Sequence):
            if len(self) == len(other):
                for i in range(len(self)):
                    self._array[i] += other[i]
                return self
            else:
                raise DimensionError(self, other, op=operator.iadd)
        elif isinstance(other, Real):
            for r in self._array:
                r += other
            return self
        else:
            return NotImplemented

    def _sub(a, b):
        if isinstance(b, Sequence):
            if len(a) == len(b):
                return rmatrix([x - y for x, y in zip(a, b)])
            else:
                raise DimensionError(a, b, op=operator.sub)
        elif isinstance(b, Real):
            return rmatrix([r - b for r in a])
        else:
            return NotImplemented

    def _mul(a, b):
        if isinstance(b, Real):
            return rmatrix([[r * b] for r in a])
        else:
            return NotImplemented

    @staticmethod
    def _matmul(a, b):
        if len(a[0]) == len(b):
            T = tuple(zip(*b))
            return rmatrix([rarray([rmatrix._dot(x, y) for y in T]) for x in a])
        else:
            raise DimensionError(a, b, op=operator.matmul)

    @staticmethod
    def _dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def __truediv__(self, other):
        return rmatrix([r / other for r in self])

    def __floordiv__(self, other):
        return rmatrix([r // other for r in self])

    def __mod__(self, other):
        if isinstance(other, Real):
            return rmatrix(r % other for r in self)
        else:
            return NotImplemented

    def __pos__(self):
        return rmatrix([+r for r in self])

    def __neg__(self):
        return rmatrix([-r for r in self])

    def append(self, value):
        if isinstance(value, rarray):
            if not self._array:
                self._array.append(value)
                self._dtype = value.dtype
            elif len(value) == len(self[0]):
                self._array.append(value)
                self._dtype = get_dtype(self._dtype, value.dtype)
                self.astype(self._dtype, _update=True)
            else:
                raise DimensionError(self, value, op=self.append)
        elif isinstance(value, Sequence):
            try:
                self.append(rarray(value))
            except TypeError:
                raise TypeError(f"must append matrix with array")
        else:
            raise NotImplementedError

    def extend(self, values):
        if isinstance(values, rmatrix):
            if len(values[0]) == len(self[0]):
                self._array.extend(values.array)
                self._dtype = get_dtype(self._dtype, values.dtype)
                self.astype(self._dtype, _update=True)
            else:
                raise DimensionError(self, values, op=self.extend)
        elif isinstance(values, Sequence):
            try:
                self.extend(rmatrix(values))
            except TypeError:
                raise TypeError(f"must extend matrix with another matrix")
        else:
            raise NotImplementedError

    def copy(self):
        return rmatrix([r.copy() for r in self])

    def transpose(self):
        return rmatrix(list(map(lambda *a: rarray([*a]), *self._array)))

    def transpose2(self):
        itermatrix = tuple(map(iter, self._array))
        return rmatrix([rarray(r.__next__() for r in itermatrix) for _ in range(len(self[0]))])

    def invert(self):
        if self.dtype in (int, Fraction):
            matrix = []
            det = self.det()
            for j in range(len(self[0])):
                matrix.append(rarray([]))
                for i in range(len(self)):
                    matrix[j].append(self._cofactor(i, j) / det)
            return rmatrix(matrix)
        else:
            return rmatrix([rarray([float(e) for e in r]) for r in np_inv(self._array)])

    def _compare(a, b, op):
        if isinstance(b, Real):
            return bmatrix(list(op(r, b) for r in a))
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bmatrix(list(op(x, y) for x, y in zip(a, b)))
            else:
                raise DimensionError(a, b, op=op)
        else:
            return NotImplemented

    def insert(self, index, value):
        if isinstance(value, rarray):
            if len(self[0]) == len(value):
                if self._dtype is value.dtype:
                    self._array.insert(index, value)
                else:
                    self._dtype = get_dtype(self._dtype, value.dtype)
                    self._array.insert(index, value)
                    self.astype(self._dtype, _update=True)
            else:
                raise DimensionError(self, value, op=self.insert)
        elif isinstance(value, Sequence):
            self.insert(index, rarray(value))
        else:
            return NotImplemented

    def __getitem__(self, item):
        if isinstance(item, slice):
            return rmatrix(self._array[item])
        else:
            return self._array[item]

    def rm_row(self, value, *, _update=False):
        """
        Removes row at index: value. If value is None, removes all null rows (rows with just 0)
        """
        if _update:
            if value is None:
                self._array = list(r for r in self if any(r))
            else:
                del self._array[value]
        elif value is None:
            return rmatrix(list(r for r in self if any(r)))
        else:
            return rmatrix(list(self[i] for i in range(len(self)) if i != value))

    def rm_col(self, value, *, _update=False):
        if _update:
            if value is None:
                self._array = list(map(lambda *a: rarray([*a]), *(r for r in self.T if any(r))))
            else:
                for r in self._array:
                    del r[value]
        elif value is None:
            return rmatrix(list(map(lambda *a: rarray([*a]), *(r for r in self.T if any(r)))))
        else:
            matrix = self.copy()
            matrix.rm_col(value, _update=True)
            return matrix

    def separate(self, index=-1):
        if index < -len(self[0]) or index >= len(self[0]):
            raise IndexError(f"index of {index} out of range")
        else:
            index %= len(self[0])
            left, right = [], []
            for i in range(len(self)):
                left.append(self[i][:index])
                right.append(self[i][index:])
            return rmatrix(left), rmatrix(right)

    def gauss(self, row, column):
        pivot = self[row]
        self._array = list(map(lambda j: self[j] if j == row else self[j] - (pivot * self[j][column]),
                               range(len(self))))

    def rref(self):
        """Uses Gaussian elimination to row reduce matrix, returning a new instance of a matrix in reduced row
        echelon form. If instance matrix is constructed over the finite field of integers (matrix.mod is an int)
        then a modified version of Gaussian elimination is used to row reduce the matrix to construct rref
        (see examples and notes for more details).

        Examples
        --------
        [1]
        add non-mod example here

        [2]

        >>> a = Matrix([[2, 4, 5], [5, 1, -2], [5, 6, 8]], mod=14, aug=True)
        >>> repr(a.rref())
        Matrix(array=[[1, 0, 2], [0, 1, 2], [0, 0, 7]], aug=True, mod=14)


        Notes
        -----
        Example 2: In this example, the matrix given is augmented and is constructed over Z14. The rref process
        is as follows (for specifics see below paragraph):
        First, the matrix is searched iteratively for a row at index 0 with a value that can be converted to 1. This
        is the second row, which has value 5 at index 0 which is invertible mod 14, with an inverse of 3. The
        entire second row is multiplied by 3, resulting in [1, 3, 8]. This row is then subtracted from all other rows
        such that the value at index 0 for all other rows is 0. This process is then repeated, this time with index
        1 and row 3, which now has value 5 at index 1. Row is multiplied by inverse and subtracted from all other
        rows. At this point, since the matrix is augmented, rref computation is complete.

        At each instance of compuation for example 2, each value for the given column index is checked
        not only to be invertible, but more broadly if it cna be converted to 1, through a combination of
        floor division and multiplication. If this is not possible for any value at the index, rref is still
        possible, but not currently computed. It is possible because, through finding pivots in other rows,
        the value at that index for all rows may become invertible given the modulus, but this requires
        returning to the unfinished column, and this computation is currently not supported. Thus, rref
        over the finite field of integers is not guaranteed to be rref.

        Specifics for example 2:

        after second row is identified as having an invertible value at index 0, second row becomes [1, 3, 8]

        >>> Matrix(array=[[2, 4, 5], [1, 3, 8], [5, 6, 8]], aug=True, mod=14)

        rows are swapped

        >>> Matrix(array=[[1, 3, 8], [2, 4, 5], [5, 6, 8]], aug=True, mod=14)

        rows are reduced

        >>> Matrix(array=[[1, 3, 8], [0, 12, 3], [0, 5, 10]], aug=True, mod=14)

        matrix now has pivot in first column, first row, process is repeated with second column
        third row has invertible element in second column, so row is converted to pivot

        >>> Matrix(array=[[1, 3, 8], [0, 12, 3], [0, 1, 2]], aug=True, mod=14)

        rows are swapped and reduced

        >>> Matrix(array=[[1, 0, 2], [0, 1, 2], [0, 0, 7]], aug=True, mod=14)
        computation complete
        """

        array = self.copy()     # deep copy
        pivot_row = 0           # first pivot belongs in first row
        height = len(self)

        for j in range(len(self[0]) - int(self.augmented)):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, height):

                # if non-zero element, this row can become pivot row
                if array[i][j]:

                    # make j'th element the pivot, reducing rest of row as well
                    array[i].make_pivot(j)
                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row][:], array[i][:]
                    array.gauss(pivot_row, j)  # row reduce everything else
                    pivot_row += 1
                    break

        return array

    @property
    def is_rref(self):
        pivots = {-1}
        width = len(self[0]) - int(self.augmented)
        for i in range(len(self)):
            for j in range(width):
                e = self[i][j]
                if not e:
                    continue
                elif e == 1:
                    if j <= max(pivots):
                        return False
                    elif reduce(lambda r, c: r + 1 if c else r, self.T[j], 0) != 1:
                        return False
                    else:
                        pivots.add(j)
                        break
                else:
                    return False
        return True

    @property
    def is_consistent(self):
        if not self.augmented:
            return True
        else:
            matrix, solutions = self.separate()
            for i in range(len(solutions)):
                if not solutions[i][0] or any(matrix[i]):
                    continue
                else:
                    return False
            return True

    @property
    def rank(self):
        if self.augmented:
            return sum(1 for r in self.rref() if any(r[:-1]))
        else:
            return sum(1 for r in self.rref() if any(r))

    @property
    def null(self):
        return len(self[0]) - self.rank - int(self.augmented)

    def trace(self):
        if len(self) != len(self[0]):
            raise DimensionError(self, op=self.trace)
        else:
            return sum(self[i][i] for i in range(len(self)))

    def minor(self, index: Union[tuple, int], *args: Optional[int]):
        """
        Computes the determinant of instance matrix if no index given. If index given,
        computes the minor A1,index (referring to the resulting matrix after removing
        row 1 and column index) multiplied against the value of A[0][index] with the
        correct sign (negative if index is even [when start counting at 1] otherwise
        positive). Use A.det() for calculating determinant for efficiency, unless specific
        minors or solving for sympy.Symbol is required.

        Examples
        --------
        >>> from sympy import Symbol
        >>> x = Symbol('x')
        >>> A = rmatrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]]) - rmatrix(rows=3, identity=x)
        >>> A.minor()
        -6*x + (3 - x)*(4 - x)*(-x - 1) + 18

        >>> A = rmatrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.minor()
        6
        >>> A.det()
        6

        >>> A = rmatrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.minor(0, 1)
        27
        >>> B = rmatrix([[3, 1], [3, 4]])
        >>> B.det()
        9

        Notes
        -----
        in the first example, 'x' is a sympy.Symbol and the solution given is solvable using sympy.solve()

        in the last two examples, A.minor(1) computes the A[0, 1] * determinant of A1,2 (A with row 1 and column 2
        removed, start counting at 1) while B.det() computes the determinant of B, which is A1,2. Since the sign
        of the minor alternates, A.minor(1) returns -3 * -1 * det(A1,2) = -3 * -1 * B.det() = 27
        """

        matrix = self.copy()
        if isinstance(index, int) and len(args) == 1 and isinstance(args[0], int):
            del matrix[index]
            for row in matrix:
                del row[args[0]]
        elif isinstance(index, Sequence) and len(index) == 2 and all(isinstance(i, int) for i in index):
            del matrix[index[0]]
            for row in matrix:
                del row[index[1]]
        else:
            raise ValueError(f"{self.__class__.__name__}.{self.minor.__name__}() requires two integer indices")

        return matrix.det()

    def _cofactor(self, index: Union[tuple, int], *args: Optional[int]):
        if isinstance(index, int) and len(args) == 1 and isinstance(args[0], int):
            i, j = index, args[0]
        elif isinstance(index, tuple) and len(index) == 2 and all(isinstance(i, int) for i in index):
            i, j = index
        else:
            raise ValueError(f"{self.__class__.__name__}.{self.minor.__name__}() requires two integer indices")

        matrix = self.copy()
        del matrix[i]
        for row in matrix:
            del row[j]

        return pow(-1, i + j) * matrix.det()

    def cofactor(self):
        matrix = []
        for i in range(len(self)):
            matrix.append(rarray([]))
            for j in range(len(self[0])):
                matrix[i].append(self._cofactor(i, j))
            matrix[i].astype(self._dtype, _update=True)
        return rmatrix(matrix)

    def det(self):
        if len(self) != len(self[0]):
            raise DimensionError(self, op=self.det)
        elif len(self) == 2:
            return self[0][0] * self[1][1] - self[1][0] * self[0][1]
        # if int/Fraction, converting to np.float64 would be bad
        elif self.dtype in (int, Fraction) or any(isinstance(e, Expr) for e in self.flat):
            matrix = self.copy()
            det = 0
            for j in range(len(self[0])):
                det += matrix.minor(0, j) * pow(-1, j) * self[0][j]
            return det
        else:   # if not int/Fraction, then not worried about rounding errors so use np's for efficiency
            det = float(np_det(self._array))
            if det.is_integer():
                return int(det)
            else:
                return det

    def char_poly(self, sym=Symbol('x')):
        if len(self) != len(self[0]):
            raise DimensionError(self, op=self.char_poly)
        elif not isinstance(sym, Symbol):
            sym = Symbol(sym)

        matrix = self.copy()
        for i in range(len(matrix)):
            matrix[i][i] -= sym

        return matrix.det()

    def eigvals(self):
        def sy_to_py(s):    # converts a sympy expression to the appropriate python number type
            s = complex(s)
            if not s.imag:
                s = s.real
                if s.is_integer():
                    s = int(s)
            return s

        if len(self) == len(self[0]):
            return tuple(sy_to_py(e) for e in solve(self.char_poly()))
        else:
            raise DimensionError(self, op=self.eigvals)

    def eigvec(self, eigvals=None):
        if len(self) != len(self[0]):
            raise DimensionError(self, op=self.eigvec)
        elif eigvals is None:
            eigvals = self.eigvals()

        vectors = []
        identity = rmatrix.identity(len(self))
        for e in eigvals:
            mat = self - (identity * e)
            kern = mat.kernel()
            for v in kern:
                vectors.append(v)
        return rmatrix(vectors).transpose()

    def kernel(self):
        """Computes the basis of the kernel for the given matrix."""

        if self.augmented:
            return None

        size = len(self)  # get number of rows
        matrix = rmatrix([])
        width = len(self[0])
        blank = [0] * width
        for j in range(width):
            identity = blank[:]
            identity[j] = 1
            row = [r[j] for r in self] + identity
            matrix.append(rarray(row))

        pivot_row = 0  # first pivot belongs in first row

        for j in range(size):  # iterate only over current matrix, not attached identity matri

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(matrix)):

                # if non-zero element, this row can become pivot row
                if matrix[i][j]:

                    # make j'th element the pivot, reducing rest of row as well
                    matrix[i].make_pivot(j)
                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        matrix[i], matrix[pivot_row] = matrix[pivot_row][:], matrix[i][:]

                    matrix.gauss(pivot_row, j)  # row reduce everything else
                    pivot_row += 1

        array, kern = matrix.separate(size)  # separates original matrix from now modified identity matrix
        basis = rmatrix([])

        for i, row in enumerate(array):
            if any(row):    # all null rows in original matrix correspond to basis vector for kernel
                basis.append(kern[i])

        return basis        # basis is list of rows, transpose into standard of column vectors


class mmatrix(ABCMatrix):

    __slots__ = '_mod',

    def __init__(self, array, mod):
        if isinstance(mod, Integral):
            self._mod = int(mod)
        else:
            raise TypeError(f"modulus must be {Integral.__name__} not {type(mod)}")

        super().__init__()
        if all(isinstance(r, marray) for r in array):
            self._array = array
        elif all(isinstance(r, Sequence) and all(isinstance(e, Integral) for e in r) for r in array):
            self._array = list(marray(list(r), self._mod) for r in array)
        elif isinstance(array, Integral):
            self._array = [marray([int(array)], self._mod)]
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from 2-dimensional Sequence of Reals")

        self._dtype = int
        # set modulus of all arrays to conform with input mod
        for r in self._array:
            r.modulo(self._mod)

    @property
    def mod(self):
        return self._mod

    @staticmethod
    def identity(size):
        array = []
        row = [0] * size
        for i in range(size):
            id_row = marray(row[:], 2)
            id_row[i] = 1
            array.append(id_row)
        return mmatrix(array, 2)

    def __setitem__(self, key, value):
        if isinstance(value, marray):
            if self.mod == value.mod:
                self._array[key] = value
            else:
                raise ValueError(f"all modulus's in {self.__class__.__name__} must agree")
        elif isinstance(value, Sequence):
            self._array[key] = marray(list(value), self.mod)
        else:
            raise NotImplementedError

    def _add(a, b):
        if isinstance(b, Sequence):
            if len(a) == len(b):
                return mmatrix([x + y for x, y in zip(a, b)], a.mod)
            else:
                raise DimensionError(a, b, op=operator.add)
        elif isinstance(b, Integral):
            return mmatrix([r + b for r in a], a.mod)
        else:
            return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, Sequence):
            if len(self) == len(other):
                for i in range(len(self)):
                    self._array[i] += other[i]
                return self
            else:
                raise DimensionError(self, other, op=operator.iadd)
        elif isinstance(other, Integral):
            for r in self._array:
                r += other
            return self
        else:
            return NotImplemented

    def _sub(a, b):
        if isinstance(b, Sequence):
            if len(a) == len(b):
                return mmatrix([x - y for x, y in zip(a, b)], a.mod)
            else:
                raise DimensionError(a, b, op=operator.sub)
        elif isinstance(b, Integral):
            return mmatrix([r - b for r in a], a.mod)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Sequence):
            return self._matmul(self, other, self.mod)
        else:
            return self._mul(other)

    def __rmul__(self, other):
        if isinstance(other, Sequence):
            return self._matmul(other, self, self.mod)
        else:
            return self._mul(other)

    def _mul(a, b):
        return mmatrix([r * b for r in a], a.mod)

    @staticmethod
    def _matmul(a, b, m):
        if len(a[0]) == len(b):
            T = tuple(zip(*b))
            return mmatrix([marray([mmatrix._dot(x, y, m) for y in T], m) for x in a], m)
        else:
            raise DimensionError(a, b, op=operator.matmul)

    @staticmethod
    def _dot(a, b, m):
        return sum((x * y) % m for x, y in zip(a, b)) % m

    def __mod__(self, other):
        return mmatrix([r % other for r in self], self.mod)

    def __pos__(self):
        return mmatrix([+r for r in self], self.mod)

    def __neg__(self):
        return mmatrix([-r for r in self], self.mod)

    def append(self, value):
        if isinstance(value, marray):
            if self.mod != value.mod:
                raise ValueError(f"all modulus's in {self.__class__.__name__} must agree")
            elif len(value) == len(self[0]):
                self._array.append(value)
            else:
                raise DimensionError(self, value, op=self.append)
        elif isinstance(value, Sequence):
            try:
                self.append(marray(value, self.mod))
            except TypeError:
                raise TypeError(f"must append matrix with array")
        else:
            raise NotImplementedError

    def extend(self, values):
        if isinstance(values, mmatrix):
            if self.mod != values.mod:
                raise ValueError(f"all modulus's in {self.__class__.__name__} must agree")
            elif len(values[0]) == len(self[0]):
                self._array.extend(values.array)
            else:
                raise DimensionError(self, values, op=self.extend)
        elif isinstance(values, Sequence):
            try:
                self.extend(mmatrix(values, self.mod))
            except TypeError:
                raise TypeError(f"must extend matrix with another matrix")
        else:
            raise NotImplementedError

    def copy(self):
        return mmatrix([r.copy() for r in self], self.mod)

    def transpose(self):
        return mmatrix(list(map(lambda *a: marray([*a], self.mod), *self._array)), self.mod)

    def invert(self):
        det = self.det()
        if gcd(det, self.mod) == 1:
            return self.cofactor().transpose() * pow(det, -1, self.mod)
        else:
            raise ValueError(f"{repr(self)} does not have an inverse")

    def _compare(a, b, op):
        if isinstance(b, Integral):
            return bmatrix(list(op(r, b) for r in a))
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bmatrix(list(op(x, y) for x, y in zip(a, b)))
            else:
                raise DimensionError(a, b, op=op)
        else:
            return NotImplemented

    def insert(self, index, value):
        if isinstance(value, marray):
            if self.mod != value.mod:
                raise ValueError(f"all modulus's in {self.__class__.__name__} must agree")
            elif len(self[0]) == len(value):
                self._array.insert(index, value)
            else:
                raise DimensionError(self, value, op=self.insert)
        elif isinstance(value, Sequence):
            self.insert(index, marray(value, self.mod))
        else:
            return NotImplemented

    def __getitem__(self, item):
        if isinstance(item, slice):
            return mmatrix(self._array[item], self.mod)
        else:
            return self._array[item]

    def minor(self, index: Union[tuple, int], *args: Optional[int]):
        """
        Computes the determinant of instance matrix if no index given. If index given,
        computes the minor A1,index (referring to the resulting matrix after removing
        row 1 and column index) multiplied against the value of A[0][index] with the
        correct sign (negative if index is even [when start counting at 1] otherwise
        positive). Use A.det() for calculating determinant for efficiency, unless specific
        minors or solving for sympy.Symbol is required.

        Examples
        --------
        >>> from sympy import Symbol
        >>> x = Symbol('x')
        >>> A = rmatrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]]) - rmatrix(rows=3, identity=x)
        >>> A.minor()
        -6*x + (3 - x)*(4 - x)*(-x - 1) + 18

        >>> A = rmatrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.minor()
        6
        >>> A.det()
        6

        >>> A = rmatrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.minor(0, 1)
        27
        >>> B = rmatrix([[3, 1], [3, 4]])
        >>> B.det()
        9

        Notes
        -----
        in the first example, 'x' is a sympy.Symbol and the solution given is solvable using sympy.solve()

        in the last two examples, A.minor(1) computes the A[0, 1] * determinant of A1,2 (A with row 1 and column 2
        removed, start counting at 1) while B.det() computes the determinant of B, which is A1,2. Since the sign
        of the minor alternates, A.minor(1) returns -3 * -1 * det(A1,2) = -3 * -1 * B.det() = 27
        """

        matrix = self.copy()
        if isinstance(index, int) and len(args) == 1 and isinstance(args[0], int):
            del matrix[index]
            for row in matrix:
                del row[args[0]]
        elif isinstance(index, Sequence) and len(index) == 2 and all(isinstance(i, int) for i in index):
            del matrix[index[0]]
            for row in matrix:
                del row[index[1]]
        else:
            raise ValueError(f"{self.__class__.__name__}.{self.minor.__name__}() requires two integer indices")

        return matrix.det()

    def _cofactor(self, index: Union[tuple, int], *args: Optional[int]):
        if isinstance(index, int) and len(args) == 1 and isinstance(args[0], int):
            i, j = index, args[0]
        elif isinstance(index, tuple) and len(index) == 2 and all(isinstance(i, int) for i in index):
            i, j = index
        else:
            raise ValueError(f"{self.__class__.__name__}.{self.minor.__name__}() requires two integer indices")

        matrix = self.copy()
        del matrix[i]
        for row in matrix:
            del row[j]

        return (pow(-1, i + j) * matrix.det()) % self.mod

    def cofactor(self):
        matrix = []
        for i in range(len(self)):
            matrix.append(marray([], self.mod))
            for j in range(len(self[0])):
                matrix[i].append(self._cofactor(i, j))
        return mmatrix(matrix, self.mod)

    def det(self):
        if len(self) != len(self[0]):
            raise DimensionError(self, op=self.det)
        elif len(self) == 2:
            return (self[0][0] * self[1][1] - self[1][0] * self[0][1]) % self.mod
        else:
            matrix = self.copy()
            det = 0
            for j in range(len(self[0])):
                det += matrix.minor(0, j) * pow(-1, j) * self[0][j]
            return det % self.mod
