import operator

from numbers import *

from abc import ABCMeta, abstractmethod
from functools import reduce
from numpy import ndarray
from random import randrange
from sympy import Symbol

from numpy import invert as npinv
from numpy import array as nparray

from cryptography318.core.tools import python_number
from cryptography318.core.fraction import Fraction
from cryptography318.core.expr import Sqrt


def matrix_like(obj):
    if hasattr(obj, '__iter__') and hasattr(obj, '__len__') and hasattr(obj, '__getitem__'):
        return all(array_like(r) for r in obj)
    else:
        return False


def array_like(obj):
    if hasattr(obj, '__iter__') and hasattr(obj, '__len__') and hasattr(obj, '__getitem__'):
        return not any(hasattr(e, '__iter__') or hasattr(e, '__getitem__') for e in obj)
    else:
        return False


def has_instances(obj, cls):
    """
    Iterates over iterable object, returning True if ALL elements of iterable return
    True for isinstance(e, cls), False otherwise.
    """
    try:
        return all(isinstance(e, cls) for e in obj)
    except TypeError:
        raise TypeError(f"{obj} is not iterable")


def array(obj, mod=None):
    if isinstance(obj, ndarray):
        obj = obj.tolist()
    elif isinstance(obj, map):
        obj = list(obj)

    # checks if array is list item and non-empty and not nested
    if isinstance(obj, list) and not obj:
        arr = obj
    elif isinstance(obj, list) and not isinstance(obj[0], list):
        arr = obj
    # if nested list is flat row vector, take it
    elif isinstance(obj, list) and len(obj) == 1:
        arr = obj[0]
    elif isinstance(obj, Number) and not isinstance(obj, complex):
        arr = [python_number(obj)]
    else:
        raise TypeError("array requires input object to be non-empty list")

    if isinstance(mod, Rational):
        if mod == 2:
            return BinaryArray(arr)
        else:
            return ModArray(arr, mod)
    else:
        return Array(arr)


def matrix1(obj=None, rows=None, cols=None, rand=False, diag=None, aug=False, solution=None, mod=None):
    """Parses arguments for creation of matrix, returning an instance of Matrix object or raising exception
        if invalid arguments are given. Use this constructor instead of object constructors."""

    augmented = True if aug or solution is not None else False

    # if array given, check if list or numpy array (cannot inherit a Matrix object)
    if obj is not None:
        if isinstance(obj, map):
            obj = list(obj)
        if isinstance(obj, list):

            # empty matrix
            if not obj:
                arr = obj
            # if array given is nested list, set self.array to given, otherwise, nest it
            elif mod is not None:
                arr = list(map(lambda a: ModArray(a, mod), obj))
            elif isinstance(obj[0], list):
                arr = list(map(Array, obj))
            elif isinstance(obj[0], Array):
                arr = obj
            else:
                raise TypeError(f"{repr(obj)} is not array-like")
        elif isinstance(obj, Array):
            arr = [obj]
        else:
            raise TypeError(f"{repr(obj)} is not array-like")

    # if random, no rows or columns required, creates m x n matrix with random m, n (where not given) and
    # values random integers between -50, 50
    elif rand:
        if rows is None:
            rows = randrange(1, 10)
        if cols is None:
            cols = randrange(1, 10)
        arr = []
        for i in range(rows):
            arr.append(ModArray([], mod) if mod is not None else Array([]))
            for j in range(cols):
                arr[i].append(randrange(-50, 50))

    # if identity, matrix must be square
    elif diag is not None:
        # matrix constructed with just identity = True will return In for random n: [1, 10)
        if rows is None and cols is None:
            rows = randrange(1, 10)

        # if one of cols or rows given, make the other equal so that matrix is square
        if rows is None:
            rows = cols
        elif cols is None:
            cols = rows
        # if neither are integers, or they are not equal, raise error
        if (not isinstance(rows, int) and not isinstance(cols, int)) or rows != cols:
            raise ValueError(f"argument(s) incompatible: rows={rows}, cols={cols}")
        arr = []
        if isinstance(diag, list):
            for i in range(rows):
                arr.append(ModArray([0] * rows, mod) if mod is not None else Array([0] * rows))
                arr[i][i] = diag[i] if isinstance(diag[i], Symbol) else int(diag[i])
        else:
            for i in range(rows):
                arr.append(ModArray([0] * rows, mod) if mod is not None else Array([0] * rows))
                arr[i][i] = diag if isinstance(diag, Symbol) else int(diag)

    # if both rows and cols are ints, make empty matrix
    elif isinstance(rows, int) and isinstance(cols, int):
        arr = []
        for i in range(rows):
            arr.append(ModArray([0] * cols, mod) if mod is not None else Array([0] * cols))
    else:
        raise ValueError(f"not enough information given to construct a matrix. one of the following is required: "
                         f"valid 2-dimensional array; integer or symbol identity; bool rand; int rows "
                         f"and int cols also required if array not given")

    # if a solution is provided, attach column to end of matrix
    if solution is not None and arr:
        if not isinstance(solution, (list, Matrix, Array)):
            raise TypeError(f"{repr(solution)} is not array_like")
        if isinstance(solution[0], (list, Array)) and len(solution) == len(arr):
            for i in range(len(arr)):
                arr[i].append(solution[i][0])
        elif isinstance(solution[0], (list, Array)) and len(solution) == 1 and len(solution[0]) == len(arr):
            for i in range(len(arr)):
                arr[i].append(solution[0][i])
        elif isinstance(solution, (list, Array)) and len(solution) == len(arr):
            for i in range(len(arr)):
                arr[i].append(solution[i])
        else:
            raise TypeError(f"{repr(solution)} is incompatible with matrix of dimension {len(arr)}x{len(arr[0])}")

    if isinstance(mod, int):
        if mod == 2:
            return None
        return None
    return Matrix(arr)


def matrix(arr=None, *, rows=None, cols=None, rndm=False, diag=None, aug=False, mod=None):
    if arr is not None:
        if matrix_like(arr):
            if has_instances(arr, Array) and mod is None:
                return Matrix(arr, aug=False)
            elif has_instances(arr, ModArray) and isinstance(mod, Rational):
                # return ModMatrix
                pass
        elif arr == list():
            pass
        else:
            raise TypeError(f"{arr} is not matrix-like")


class DimensionError(Exception):
    def __init__(self, error, opt=None):
        if isinstance(error, str):
            self.message = error
        elif opt is None:
            self.message = f"incompatible dimension for attempted matrix operation: {len(error)}x{len(error[0])}"
        else:
            self.message = f"incompatible dimension(s) for matrix operations: {len(error)}x{len(error[0])} and " \
                           f"{len(opt)}x{len(opt[0])}"
        super().__init__(self.message)

    def __str__(self):
        return self.message


class ABCArray(metaclass=ABCMeta):
    __slots__ = '_array', '_dtype'

    @abstractmethod
    def __init__(self, *args, **kwargs):
        ...

    @abstractmethod
    def __repr__(self):
        ...

    @abstractmethod
    def __str__(self):
        ...

    @abstractmethod
    def __eq__(self, other):
        ...

    @abstractmethod
    def __ne__(self, other):
        ...

    @abstractmethod
    def __setitem__(self, key, value):
        ...

    @abstractmethod
    def _add(self, a, b):
        ...

    @abstractmethod
    def _sub(self, a, b):
        ...

    @abstractmethod
    def _mul(self, a, b):
        ...

    @abstractmethod
    def _floordiv(self, a, b):
        ...

    @abstractmethod
    def _truediv(self, a, b):
        ...

    @abstractmethod
    def _pow(self, power):
        ...

    @abstractmethod
    def _construct(self, arr):
        ...

    @abstractmethod
    def astype(self, dtype):
        ...

    @abstractmethod
    def append(self, item):
        ...

    @property
    def aslist(self):
        return self._array[:]

    @property
    def dtype(self):
        return self._dtype

    def copy(self):
        return self._construct(list(map(lambda e: e, self._array)))

    def __add__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._add(e, other), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._add(e1, e2), self._array, other)))
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._add(other, e), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._add(e1, e2), other, self._array)))
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._sub(e, other), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._sub(e1, e2), self._array, other)))
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._sub(other, e), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._sub(e1, e2), other, self._array)))
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._mul(e, other), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._mul(e1, e2), self._array, other)))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._mul(other, e), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._mul(e1, e2), other, self._array)))
        else:
            return NotImplemented

    def __floordiv__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._floordiv(e, other), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._floordiv(e1, e2), self._array, other)))
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda e: self._truediv(e, other), self._array)))
        elif isinstance(other, ABCArray) or array_like(other):
            return self._construct(list(map(lambda e1, e2: self._truediv(e1, e2), self._array, other)))
        else:
            return NotImplemented

    def __mod__(self, other):
        if isinstance(other, Real):
            return self._construct(list(map(lambda e: e % other, self._array)))
        else:
            return NotImplemented

    def __pow__(self, power, modulo=None):
        if isinstance(power, Real):
            if isinstance(modulo, Real):
                return self._construct(list(map(lambda e: pow(e, power, modulo), self._array)))
            else:
                return self._pow(power)
        else:
            return NotImplemented

    def _compare(self, other, op):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: int(op(r, other)), self._array)))
        elif array_like(other):
            return self._construct(list(map(lambda r1, r2: int(op(r1, r2)), self._array, other)))
        else:
            raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self._construct(self._array[item])
        return self._array[item]

    def __contains__(self, item):
        return item in self._array

    def __len__(self):
        return len(self._array)

    def __iter__(self):
        return iter(self._array)

    def index(self, item):
        return self._array.index(item)


class Array(ABCArray):

    def __init__(self, arr):
        if has_instances(arr, Number):
            types = [type(e) for e in arr]

            # most common
            if all(t is int for t in types):
                self._dtype = int
            elif all(t is int or t is Sqrt for t in types):
                self._dtype = Sqrt

            # if not all integers, this is type order resolution
            elif complex in types:
                self._dtype = complex
            elif float in types or Sqrt in types:
                self._dtype = float
            elif Fraction in types:
                self._dtype = Fraction
            else:
                # if unrecognized type, default to first
                self._dtype = types[0]

            self._array = list(map(lambda e: e if type(e) is self._dtype else self._dtype(e), arr))
        else:
            raise ValueError(f"{self.__class__.__name__} must be constructed from list of Numbers")

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array})"

    def __str__(self):
        if self._dtype in (Fraction, Sqrt):
            string = "["
            for e in self._array[:-1]:
                string += str(e) + ", "
            string += str(self._array[-1]) + "]"
            return string
        elif self._dtype is float:
            string = "["
            for e in self._array[:-1]:
                string += str(round(e, 3)) + ", "
            string += str(round(self._array[-1], 3)) + "]"
            return string
        return str(self._array)

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __setitem__(self, key, value):
        if isinstance(value, (int, float, Fraction, Sqrt)):
            self._array[key] = self._dtype(value)
        elif isinstance(value, Complex):
            if self._dtype is complex:
                self._array[key] = complex(value)
            else:
                self._array[key] = self._dtype(value.real)
        elif isinstance(value, Real):
            self._array[key] = self._dtype(float(value))
        else:
            raise TypeError(f"unable to set array key to non-Real instance: {type(value)}")

    def astype(self, dtype):
        if self._dtype is complex:
            return Array(list(map(lambda e: dtype(e.real), self._array)))
        elif dtype is Sqrt:
            return Array(list(map(lambda e: Sqrt(1, e), self._array)))
        else:
            return Array(list(map(lambda e: dtype(e), self._array)))

    def _construct(self, arr):
        return Array(arr)

    def append(self, item):
        if isinstance(item, Complex) and self._dtype is not complex:
            self._array.append(self._dtype(item.real))
        elif isinstance(item, Real):
            self._array.append(self._dtype(item))
        elif hasattr(item, '__iter__'):
            if has_instances(item, Complex) and self._dtype is not complex:
                for e in item:
                    self._array.append(self._dtype(e.real))
            elif has_instances(item, Real):
                for e in item:
                    self._array.append(e)
        else:
            raise ValueError(f"{self.__class__.__name__}.append() only accepts Real instance(s)")

    def _add(self, a, b):
        return a + b

    def _sub(self, a, b):
        return a - b

    def _mul(self, a, b):
        return a * b

    def _floordiv(self, a, b):
        return a // b

    def _truediv(self, a, b):
        return a / b

    def _pow(self, power):
        return self._construct(list(map(lambda e: pow(e, power), self._array)))


class ModArray(ABCArray):
    __slots__ = '_mod',

    def __init__(self, arr, mod):
        if has_instances(arr, Rational) and isinstance(mod, Rational):

            # all Rational instances will either be ints or can be made into Fractions
            if all(isinstance(e, int) for e in arr):
                self._dtype = int
            else:
                self._dtype = Fraction

            if mod.denominator == 1:
                mod = int(mod)
            elif not isinstance(mod, Fraction):
                mod = Fraction(mod)

            self._mod = mod
            self._array = list(map(lambda e: e % mod if type(e) is self._dtype else self._dtype(e % mod), arr))
        else:
            raise ValueError(f"{self.__class__.__name__} must be constructed from list of Rationals")

    @property
    def mod(self):
        return self._mod

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array}, {self._mod})"

    def __str__(self):
        return str(self._array)

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __setitem__(self, key, value):
        if isinstance(value, (int, Fraction)):
            self._array[key] = self._dtype(value) % self._mod
        elif isinstance(value, Rational):
            self._array[key] = self._dtype(Fraction(value)) % self._mod
        else:
            raise TypeError(f"unable to set array key to non-Rational instance: {type(value)}")

    def _construct(self, arr):
        return ModArray(arr, self._mod)

    def astype(self, dtype):
        return ModArray(list(map(lambda e: dtype(e), self._array)), self._mod)

    def append(self, item):
        if isinstance(item, Rational):
            self._array.append(item % self._mod)
        elif hasattr(item, '__iter__') and has_instances(item, Rational):
            for e in item:
                self._array.append(e % self._mod)
        else:
            raise ValueError(f"{self.__class__.__name__}.append() only accepts Rational instance(s)")

    def _add(self, a, b):
        return (a + b) % self._mod

    def _sub(self, a, b):
        return (a - b) % self._mod

    def _mul(self, a, b):
        return (a * b) % self._mod

    def _floordiv(self, a, b):
        return NotImplemented

    def _truediv(self, a, b):
        return NotImplemented

    def _pow(self, power):
        return self._construct(list(map(lambda e: e, self._array)))


class BinaryArray(ABCArray):

    def __init__(self, arr):
        if has_instances(arr, Rational):

            if all(isinstance(e, int) for e in arr):
                self._dtype = int
            else:
                self._dtype = Fraction

            self._array = list(map(lambda e: e % 2 if type(e) is self._dtype else self._dtype(e % 2), arr))
        else:
            raise ValueError(f"{self.__class__.__name__} must be constructed from list of Rationals")

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array})"

    def __str__(self):
        return str(self._array)

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __setitem__(self, key, value):
        if isinstance(value, (int, Fraction)):
            self._array[key] = self._dtype(value) % 2
        elif isinstance(value, Rational):
            self._array[key] = self._dtype(Fraction(value)) % 2
        else:
            raise TypeError(f"unable to set array key to non-Rational instance: {type(value)}")

    def _construct(self, arr):
        return BinaryArray(arr)

    def astype(self, dtype):
        if dtype in (int, Fraction):
            return BinaryArray(list(map(lambda e: dtype(e), self._array)))
        else:
            raise TypeError(f"{dtype.__name__} is not a recognized array data type")

    def append(self, item):
        if isinstance(item, Rational):
            self._array.append(item % 2)
        elif hasattr(item, '__iter__') and has_instances(item, Rational):
            for e in item:
                self._array.append(e % 2)
        else:
            raise ValueError(f"{self.__class__.__name__}.append() only accepts Rational instance(s)")

    def _add(self, a, b):
        return (a + b) % 2

    def _sub(self, a, b):
        return (a - b) % 2

    def _mul(self, a, b):
        return (a * b) % 2

    def _floordiv(self, a, b):
        return NotImplemented

    def _truediv(self, a, b):
        return NotImplemented

    def _pow(self, power):
        return self._construct(list(map(lambda e: e, self._array)))


class ABCMatrix(metaclass=ABCMeta):
    __slots__ = '_array', '_augmented', '_dtype'

    @abstractmethod
    def __init__(self, *args, **kwargs):
        ...

    @property
    @abstractmethod
    def inverse(self):
        ...

    @property
    @abstractmethod
    def T(self):
        ...

    @abstractmethod
    def __repr__(self):
        ...

    @abstractmethod
    def __str__(self):
        ...

    @abstractmethod
    def __eq__(self, other):
        ...

    @abstractmethod
    def __ne__(self, other):
        ...

    @abstractmethod
    def __setitem__(self, key, value):
        ...

    @abstractmethod
    def _construct(self, arr):
        ...

    @abstractmethod
    def astype(self, dtype):
        ...

    @abstractmethod
    def copy(self):
        ...

    @abstractmethod
    def transpose(self):
        ...

    @abstractmethod
    def append(self, item):
        ...

    @staticmethod
    @abstractmethod
    def dot(a, b):
        ...

    @property
    def aslist(self):
        return list(map(list, self._array))

    @property
    def dtype(self):
        return self._dtype

    def flatten(self):
        return [*reduce(lambda r1, r2: list(r1) + list(r2), self)] if self._array else []

    def __add__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r + other, self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            if len(self) == len(other) and len(self[0]) == len(other[0]):
                return self._construct(list(map(lambda r1, r2: r1 + r2, self._array, other)))
            else:
                raise DimensionError(self, other)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r.__radd__(other), self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            if len(self) == len(other) and len(self[0]) == len(other[0]):
                return self._construct(list(map(lambda r1, r2: r2.__radd__(r1), other, self._array)))
            else:
                raise DimensionError(self, other)
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r - other, self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            return self._construct(list(map(lambda r1, r2: r1 - r2, self._array, other)))
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r.__rsub__(other), self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            return self._construct(list(map(lambda r1, r2: r2.__rsub__(r1), other, self._array)))
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r * other, self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            if len(other[0]) == len(self):
                if isinstance(other, ABCMatrix):
                    other = other.T
                else:
                    other = list(map(lambda i: list(map(lambda r: r[i], other)), range(len(other[0]))))
                return self._construct(list(map(lambda r1, r2: self.dot(r1, r2), self._array, other)))
            else:
                raise DimensionError(self, other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r.__rmul__(other), self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            if len(other[0]) == len(self):
                return self._construct(list(map(lambda r1, r2: self.dot(r1, r2), other, self.T)))
            else:
                raise DimensionError(other, self)
        else:
            return NotImplemented

    def __floordiv__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r // other, self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            return self._construct(list(map(lambda r1, r2: r1 // r2, self._array, other)))
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: r / other, self._array)))
        elif isinstance(other, ABCMatrix) or matrix_like(other):
            return self._construct(list(map(lambda r1, r2: r1 / r2, self._array, other)))
        else:
            return NotImplemented

    def __mod__(self, other):
        if isinstance(other, Real):
            return self._construct(list(map(lambda r: r % other, self._array)))
        else:
            return NotImplemented

    def __pow__(self, power):
        if len(self) == len(self[0]):
            if isinstance(power, int) and power > 0:
                mat = self[:]
                for e in range(power - 1):
                    mat = self * mat
                return mat
            elif power == -1:
                return self.inverse
            else:
                return NotImplemented
        else:
            raise DimensionError(self)

    def _compare(self, other, op):
        if isinstance(other, Number):
            return self._construct(list(map(lambda r: op(r, other), self._array)))
        elif matrix_like(other):
            return self._construct(list(map(lambda r1, r2: op(r1, r2), self._array, other)))
        else:
            raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)

    def __getitem__(self, item):
        if isinstance(item, slice):
            try:
                arr = []
                for i in range(item.start, item.stop, item.step):
                    arr.append(self[i][:])
            except IndexError:
                raise IndexError(f"slice index out of range of matrix")
            else:
                return self._construct(arr)
        return self._array[item]

    def __contains__(self, item):
        return item in self.flatten()

    def __len__(self):
        return len(self._array)

    def __iter__(self):
        return iter(self._array)

    def index(self, item):
        i = 0
        for row in self._array:
            if item in row:
                return i, row.index(item)
            i += 1
        raise ValueError(f"{item} is not in matrix")


class Matrix(ABCMatrix):

    def __init__(self, arr, *, aug=False):
        if matrix_like(arr):
            width = len(arr[0])
            for i in range(len(arr)):
                if len(arr[i]) != width:
                    raise DimensionError(f"all rows of matrix must be same length")
                if isinstance(arr[i], Array):
                    pass
                elif isinstance(arr[i], list):
                    arr[i] = Array(arr[i])
                else:
                    raise ValueError(f"{self.__class__.__name__} must be constructed from list of Arrays")

            self._array = arr
            self._augmented = aug
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from 2-dimensional list object")

        types = [r.dtype for r in self._array]
        if all(t is int for t in types):
            self._dtype = int
        elif all(t is int or t is Sqrt for t in types):
            self._dtype = Sqrt
        elif complex in types:
            self._dtype = complex
        elif float in types or Sqrt in types:
            self._dtype = float
        elif Fraction in types:
            self._dtype = Fraction
        else:
            self._dtype = types[0]

        if self._dtype in (int, float, complex, Fraction, Sqrt):
            self._array = list(map(lambda r: r if r.dtype is self._dtype else r.astype(self._dtype), self._array))
        else:
            raise TypeError(f"invalid data type for {self.__class__.__name__}: {self._dtype}")

    @property
    def inverse(self):
        return self._construct(npinv(nparray(self.aslist)).tolist())

    @property
    def T(self):
        return self.transpose()

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array}, aug={self._augmented})"

    def __str__(self):

        types = set()
        for row in self._array:
            types.add(type(e) for e in row)
        if len(types) == 1:
            pass
        else:
            # if multiple types, complex takes priority
            if complex in types:
                self._dtype = complex
            # next if float is there and no complex, switch to float
            elif float in types:
                self._dtype = float
            # if int is there, that means a Fraction w/ sqrt num got converted to int, switch back to Fraction
            elif int in types:
                self._dtype = Fraction
            # no else case since if there are no complex, there will be no problem trying to convert back to whatever
            # self._dtype currently is
            for i in range(len(self._array)):
                self._array[i] = self._array[i].astype(self._dtype)

        if self._dtype is float:
            arr = list(map(lambda r: list(map(lambda e: round(e, 3), r)), self._array))
        else:
            arr = self.aslist

        max_len = 0
        for row in arr:
            for v in row:
                if (l := len(str(v))) > max_len:
                    max_len = l

        padding = (max_len + 1) | 1
        formatted = "["
        for i in range(l := len(arr)):
            if i == 0:
                formatted += "["
            else:
                formatted += " ["
            for j in range(len(arr[0])):
                e = str(arr[i][j])
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

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __setitem__(self, key, value):
        if isinstance(value, Array):
            self._array[key] = value if value.dtype is self._dtype else value.astype(self._dtype)
        elif isinstance(value, (list, ABCArray)):
            value = Array(value)
            self._array[key] = value if value.dtype is self._dtype else value.astype(self._dtype)
        else:
            raise TypeError(f"unable to set matrix key to non-Array instance: {type(value)}")

    def astype(self, dtype):
        if dtype in (float, int, complex, Fraction, Sqrt):
            return Matrix(list(map(lambda r: r.astype(dtype), self._array)))
        else:
            raise TypeError(f"{dtype.__name__} is not a recognized matrix data type")

    def _construct(self, arr):
        return Matrix(arr)

    def copy(self):
        return Matrix(list(map(lambda r: r[:], self._array)))

    def transpose(self):
        return self._construct(list(map(lambda i: Array(list(map(lambda r: r[i], self._array))), range(len(self[0])))))

    def append(self, item):
        if isinstance(item, Array):
            if len(item) == len(self[0]):
                self._array.append(item.astype(self._dtype))
            else:
                raise DimensionError(f"unable to append object of length {len(item)} to matrix of width {len(self[0])}")
        elif hasattr(item, '__iter__') and has_instances(item, (Array, list)):
            arr = list(map(lambda r: r[:], self._array))
            item = list(map(lambda r: r if isinstance(r, Array) else Array(r), item))
            for row in item:
                if len(item) == len(self[0]):
                    arr.append(row.astype(self._dtype))
                else:
                    raise DimensionError(
                        f"unable to append object of length {len(item)} to matrix of width {len(self[0])}"
                    )
            self._array = arr
        else:
            raise ValueError(f"{self.__class__.__name__}.append() only accepts Array instance(s) of equal length(s)")

    @staticmethod
    def dot(a, b):
        return sum(map(lambda x, y: x * y, a, b))


class ModularMatrix(ABCMatrix):
    __slots__ = '_mod',

    def __init__(self, arr, mod, *, aug=False):
        if matrix_like(arr) and isinstance(mod, Rational):
            width = len(arr[0])
            mod = int(mod) if mod.denominator == 1 else Fraction(mod)
            for i in range(len(arr)):
                if len(arr[i]) != width:
                    raise DimensionError(f"all rows of matrix must be same length")
                if isinstance(arr[i], ModArray) and arr[i].mod != mod:
                    raise ValueError(f"given modulus conflicts with row modulus")
                elif isinstance(arr[i], (Array, list)):
                    arr[i] = ModArray(arr[i], mod)
                else:
                    raise ValueError(f"{self.__class__.__name__} must be constructed from list of Arrays")

            self._augmented = aug
            self._mod = mod
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from 2-dimensional list object and "
                            f"a Rational modulus")

        dtype = int if all(r.dtype is int for r in self._array) else Fraction
        if dtype is int:
            self._array = list(map(lambda r: r if r.dtype is int else int(r), self._array))
        else:
            self._array = list(
                map(lambda r: r if r.dtype is dtype else ModArray(
                    list(map(lambda e: e if type(e) is Fraction else Fraction(e), r)), mod
                ), self._array)
            )

    @property
    def inverse(self):
        return self._construct(npinv(nparray(self.aslist)).tolist())

    @property
    def T(self):
        return self.transpose()

    def _normalize_types(self):
        dtype = int
        for row in self:
            row_types = [type(e) for e in row]
            if any(t is complex for t in row_types):
                dtype = complex
                break
            elif any(t is float for t in row_types):
                dtype = float
                break
            elif any(t is Fraction for t in row_types):
                dtype = Fraction
                break

    def __repr__(self):
        pass

    def __str__(self):
        arr = list(map(
            lambda r: list(map(
                lambda e: round(e, 3) if isinstance(e, float) else e, r
            )), self._array
        ))

        max_len = 0
        for row in arr:
            for v in row:
                if (l := len(str(v))) > max_len:
                    max_len = l

        padding = (max_len + 1) | 1
        formatted = "["
        for i in range(l := len(arr)):
            if i == 0:
                formatted += "["
            else:
                formatted += " ["
            for j in range(len(arr[0])):
                e = str(arr[i][j])
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

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __setitem__(self, key, value):
        if isinstance(value, ModArray):
            self._array[key] = value if value.dtype is self._dtype else value.astype(self._dtype)
        elif isinstance(value, (list, ABCArray)):
            value = ModArray(value, self._mod)
            self._array[key] = value if value.dtype is self._dtype else value.astype(self._dtype)
        else:
            raise TypeError(f"unable to set matrix key to non-Array instance: {type(value)}")

    def _construct(self, arr):
        return Matrix(arr)

    def astype(self, dtype):
        if dtype in (int, Fraction):
            return ModularMatrix(list(map(lambda r: r.astype(dtype), self._array)), self._mod)
        else:
            raise TypeError(f"{dtype.__name__} is not a recognized matrix data type")

    def copy(self):
        return Matrix(list(map(lambda r: r[:], self._array)))

    def transpose(self):
        return self._construct(list(map(lambda i: Array(map(lambda r: r[i], self)), range(len(self[0])))))

    def append(self, item):
        if isinstance(item, Array):
            if len(item) == len(self[0]):
                self._array.append(item)
            else:
                raise DimensionError(f"unable to append object of length {len(item)} to matrix of width {len(self[0])}")
        elif hasattr(item, '__iter__') and has_instances(item, Array):
            arr = list(map(lambda r: r[:], self._array))
            for row in item:
                if len(item) == len(self[0]):
                    arr.append(row)
                else:
                    raise DimensionError(
                        f"unable to append object of length {len(item)} to matrix of width {len(self[0])}"
                    )
            self._array = arr
        else:
            raise ValueError(f"{self.__class__.__name__}.append() only accepts Array instance(s) of equal length(s)")

    @staticmethod
    def dot(a, b):
        return sum(map(lambda x, y: x * y, a, b))
