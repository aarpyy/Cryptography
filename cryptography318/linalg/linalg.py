from cryptography318.linalg.array import Array, ArrayMod, where
from cryptography318.core.tools import dot, python_number, string_reduce, fraction, r_append
from random import randrange
from numpy import ndarray
from numpy import array as np_array
from sympy import Symbol, solve, im
from numpy.linalg import inv
from numpy.linalg import det as np_det
from math import sqrt
from functools import reduce, wraps
from itertools import combinations
from fractions import Fraction
from numbers import Number


# decorator for methods that should conserve .mod and .augmented between objects
# this decorator is used for nested methods only, since if it is used for a surface method,
# it prevents PyCharm from prompting the func w/ parentheses/parameters, so a shell function
# is used to call the actual function w/ decorator
def conserve_attributes(func):  # *CA - funcs marked with this use decorator in nested function

    @wraps(func)
    def new_func(instance, *args, **kwargs):
        if not isinstance(instance, Matrix):
            raise AttributeError(f"conserve_attributes wrapper unsupported for functions outside of Matrix class")
        result = func(instance, *args, **kwargs)
        if isinstance(result, Matrix):
            if isinstance(instance, ModularMatrix):
                result.mod = instance.mod
            result.augmented = instance.augmented
        return result

    return new_func


def array_like(obj):
    if not hasattr(obj, '__iter__'):
        return False
    for row in obj:
        if not hasattr(row, '__iter__'):
            return False
        for e in row:
            if hasattr(e, '__iter__'):
                return False
    return True


def row_vector(obj):
    """Returns True if object is a row vector, False otherwise"""

    if hasattr(obj, '__iter__') and hasattr(obj[0], '__iter__'):
        if len(obj) == 1:
            obj = obj[0]
        elif len(obj[0]) == 1:
            obj = list(transpose_obj(obj))[0]
        else:
            return False

    for r in obj:
        if isinstance(r, (list, Array, ndarray, Matrix)):
            return False
    return obj


def assert_square(obj):
    if len(obj) != len(obj[0]):
        raise AttributeError("operation unsupported for non-square matrices")


def transpose_obj(obj):
    """Transpose method that computes transpose for array-like objects without transpose
    methods."""

    return map(lambda i: list(map(lambda r: r[i], obj)), range(len(obj[0])))


def matmul(obj1, obj2, row_type=None, obj_type=None, mod=None):
    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    # attempts to transpose using built in methods, manually performs if no method exists
    transpose = getattr(obj2, 'transpose', None)
    if transpose is None:
        T = transpose_obj(obj2)
    else:
        T = transpose()
    if isinstance(other, Number) and not isinstance(other, complex):
        if row_type is ArrayMod:
            return obj_type(map(lambda r: ArrayMod(map(lambda c: dot(c, r, mod), T), mod), obj1))
        return obj_type(map(lambda r: row_type(map(lambda c: dot(c, r, mod), T)), obj1))
    return obj_type(map(lambda r: row_type(map(lambda c: dot(c, r), T)), obj1))


def is_binary_matrix(obj):
    """Returns False if non binary element in Matrix (not 0 or 1), True otherwise. Returns
    False if object is not Matrix."""

    try:
        for row in obj:
            for e in row:
                if e not in (0, 1):
                    return False
    except TypeError:
        return False
    else:
        return True


# methods reserved for binary matrices
def union(obj1, obj2, row_type=None, obj_type=None):
    """Returns the logical union of two binary matrices.

    If this operation is needed with non-matrix objects, use set.union(A, B)"""

    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    if not is_binary_matrix(obj1) and not is_binary_matrix(obj2):
        raise AttributeError(f"union between {repr(obj1)} and {repr(obj2)} is unsupported")
    if len(obj1) != len(obj2) or len(obj1[0]) != len(obj2[0]):
        raise AttributeError(f"operation unsupported for objects with dimensions {len(obj1)}x{len(obj1[0])} "
                             f"and {len(obj2)}x{len(obj2[0])}")
    if row_type is ArrayMod:
        return obj_type(map(lambda r1, r2: ArrayMod(map(lambda e1, e2: e1 ^ e2, r1, r2), 2), obj1, obj2))
    return obj_type(map(lambda r1, r2: row_type(map(lambda e1, e2: e1 ^ e2, r1, r2)), obj1, obj2))


def intersection(obj1, obj2, row_type=None, obj_type=None):
    """Returns the logical intersection of two binary matrices.

    If this operation is needed with non-Matrix objects, use set.intersection(A, B)"""

    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    if not is_binary_matrix(obj1) and not is_binary_matrix(obj2):
        raise AttributeError(f"union between {repr(obj1)} and {repr(obj2)} is unsupported")
    if len(obj1) != len(obj2) or len(obj1[0]) != len(obj2[0]):
        raise AttributeError(f"operation unsupported for objects with dimensions {len(obj1)}x{len(obj1[0])} "
                             f"and {len(obj2)}x{len(obj2[0])}")
    if row_type is ArrayMod:
        return obj_type(map(lambda r1, r2: ArrayMod(map(lambda e1, e2: e1 & e2, r1, r2), 2), obj1, obj2))
    return obj_type(map(lambda r1, r2: row_type(map(lambda e1, e2: e1 & e2, r1, r2)), obj1, obj2))


def disjunction(obj1, obj2, row_type=None, obj_type=None):
    """Returns the logical intersection of two binary matrices.

    If this operation is needed with non-Matrix objects, use set(A) - set(B)"""

    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    if not is_binary_matrix(obj1) and not is_binary_matrix(obj2):
        raise AttributeError(f"union between {repr(obj1)} and {repr(obj2)} is unsupported")
    if len(obj1) != len(obj2) or len(obj1[0]) != len(obj2[0]):
        raise AttributeError(f"operation unsupported for objects with dimensions {len(obj1)}x{len(obj1[0])} "
                             f"and {len(obj2)}x{len(obj2[0])}")
    if row_type is ArrayMod:
        return obj_type(map(lambda r1, r2: ArrayMod(map(lambda e1, e2: e1 | e2, r1, r2), 2), obj1, obj2))
    return obj_type(map(lambda r1, r2: row_type(map(lambda e1, e2: e1 | e2, r1, r2)), obj1, obj2))


def binary_inverse(obj, row_type=None, obj_type=None):
    """Returns the negative image of a binary matrix"""

    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    if not is_binary_matrix(obj):
        raise AttributeError(f"binary inverse unsupported for {repr(obj)}")

    if row_type is ArrayMod:
        return obj_type(map(lambda r: ArrayMod(map(lambda e: (e + 1) % 2, r), 2), obj))
    return obj_type(map(lambda r: row_type(map(lambda e: (e + 1) % 2, r)), obj))


def matrix(obj=None, rows=None, cols=None, rand=False, diag=None, aug=False, solution=None, mod=None):
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
                array = obj
            # if array given is nested list, set self.array to given, otherwise, nest it
            elif mod is not None:
                array = list(map(lambda a: ArrayMod(a, mod), obj))
            elif isinstance(obj[0], list):
                array = list(map(Array, obj))
            elif isinstance(obj[0], Array):
                array = obj
            else:
                raise TypeError(f"{repr(obj)} is not array-like")
        elif isinstance(obj, Array):
            array = [obj]
        else:
            raise TypeError(f"{repr(obj)} is not array-like")

    # if random, no rows or columns required, creates m x n matrix with random m, n (where not given) and
    # values random integers between -50, 50
    elif rand:
        if rows is None:
            rows = randrange(1, 10)
        if cols is None:
            cols = randrange(1, 10)
        array = []
        for i in range(rows):
            array.append(ArrayMod([], mod) if mod is not None else Array([]))
            for j in range(cols):
                array[i].append(randrange(-50, 50))

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
        array = []
        if isinstance(diag, list):
            for i in range(rows):
                array.append(ArrayMod([0] * rows, mod) if mod is not None else Array([0] * rows))
                array[i][i] = diag[i] if isinstance(diag[i], Symbol) else int(diag[i])
        else:
            for i in range(rows):
                array.append(ArrayMod([0] * rows, mod) if mod is not None else Array([0] * rows))
                array[i][i] = diag if isinstance(diag, Symbol) else int(diag)

    # if both rows and cols are ints, make empty matrix
    elif isinstance(rows, int) and isinstance(cols, int):
        array = []
        for i in range(rows):
            array.append(ArrayMod([0] * cols, mod) if mod is not None else Array([0] * cols))
    else:
        raise ValueError(f"not enough information given to construct a matrix. one of the following is required: "
                         f"valid 2-dimensional array; integer or symbol identity; bool rand; int rows "
                         f"and int cols also required if array not given")

    # if a solution is provided, attach column to end of matrix
    if solution is not None and array:
        if not isinstance(solution, (list, Matrix, Array)):
            raise TypeError(f"{repr(solution)} is not array_like")
        if isinstance(solution[0], (list, Array)) and len(solution) == len(array):
            for i in range(len(array)):
                array[i].append(solution[i][0])
        elif isinstance(solution[0], (list, Array)) and len(solution) == 1 and len(solution[0]) == len(array):
            for i in range(len(array)):
                array[i].append(solution[0][i])
        elif isinstance(solution, (list, Array)) and len(solution) == len(array):
            for i in range(len(array)):
                array[i].append(solution[i])
        else:
            raise TypeError(f"{repr(solution)} is incompatible with matrix of dimension {len(array)}x{len(array[0])}")

    if isinstance(mod, int):
        if mod == 2:
            return BinaryMatrix(array, aug=augmented)
        return ModularMatrix(array, mod, aug=augmented)
    return Matrix(array, aug=augmented)


class Matrix:
    def __init__(self, array, aug=False):
        for row in array:
            if not isinstance(row, Array):
                print(array)
                raise AttributeError(f"valid matrix must be constructed from Array object")
            for v in row:
                if isinstance(v, (list, Array, str)):
                    raise ValueError(f"matrix must be 2-dimensional array with number or symbol values. "
                                     f"given: {repr(array)}")
        self.array = array
        self.augmented = aug
        self.attr = [self.augmented]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __getitem__(self, item):
        if item == slice(None, None, None):
            return self.copy()
        return self.array[item]

    def __contains__(self, item):
        return self.array.__contains__(item)

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return iter(self.array)

    def __repr__(self):
        return f"{__class__.__name__}({self.tolist()}, aug={self.augmented})"

    def __str__(self):
        str_array = []
        max_len = 0
        for i in range(len(self)):
            str_array.append([])
            for j in range(len(self[0])):
                n = self[i][j]
                if isinstance(n, complex):
                    str_real, str_imag = string_reduce(n.real), string_reduce(n.imag) + "i"
                    if n.imag < 0:
                        str_imag = str_imag[1:]
                    if abs(n.imag) == 1:
                        str_imag = "i"
                    if n.real == 0 and n.imag == 0:
                        string = "0"
                    elif n.real == 0 and n.imag < 0:
                        string = "-" + str_imag
                    elif n.imag < 0:
                        string = str_real + " - " + str_imag
                    else:
                        string = str_real + " + " + str_imag
                    if len(string) + 1 > max_len:
                        max_len = len(string) + 1
                    str_array[i].append(string)
                elif isinstance(n, bool):
                    string = str(n)
                    if len(string) > max_len:
                        max_len = len(string)
                    str_array[i].append(string)
                else:
                    string = string_reduce(n)
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

    def __eq__(self, other):

        # if other is nested list i.e. 2-dimensional array, return False if values at each corresponding row, col
        # are not similar, True otherwise (returning True if all values within 0.1 of each other because numpy
        # rounding /dealing with long floats can result in slight variation of same matrices)
        if array_like(other) and isinstance(other[0], (list, Array, ndarray)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                return False
                # raise ValueError(f"equality between {repr(self)} and {repr(other)} is unsupported")
            for i in range(len(self)):
                for j in range(len(self[0])):
                    if abs(self[i][j] - other[i][j]) > 0.1:
                        return False
            return True

        # if not nested list, two objects will not be compared as arrays, but will return a binary matrix
        # with truth values wherever self's values are in the list
        if isinstance(other, (list, set, tuple)):
            binary_matrix = BinaryMatrix(
                list(map(lambda r: ArrayMod(map(lambda e: 1 if e in other else 0, r), 2), self))
            )

            # pivot keyword allows user to search for pivots in matrix
            if 'pivot' in other:
                # returns the union of searching for pivot
                return binary_matrix.union(self == 'pivot')
            return binary_matrix

        # if comparing to a number, return binary matrix with values = 1 wherever matrix values = number
        if isinstance(other, Number) and not isinstance(other, complex):
            return BinaryMatrix(list(map(lambda r: ArrayMod(map(lambda e: 1 if e == other else 0, r), 2), self)))

        # if user searching for pivots, return binary matrix with values = 1 wherever an entry = 1
        # that is not in a row with a preceding value != 0 or 1 (first non-zero entry in row should be pivot)
        # also checks if entry = 1 is in same column as already identified pivot
        if other == 'pivot':
            binary_matrix = matrix(rows=len(self[0]), cols=len(self), mod=2)
            adj = 1 if self.augmented else 0
            pivot_row = -1

            columns = []
            array = self.transpose()
            for i in range(len(array) - adj):
                # slices row so that if augmented, don't include the last element
                row = array[i]

                # find all instances of non-zero entries in row, if there are none, or the first one is in a column
                # with a non-zero entry, or value is not one, no more pivots can be found in matrix so return
                indices = where(row)[0]

                # if the current index is in a row with a non-zero entry preceding, it cannot be a pivot
                # checks if in column of transpose == in row of original
                index = indices[0]
                if len(indices) != 1 or index <= pivot_row or index in columns:
                    continue

                # add to list of columns every row with non-zero entry, since all pivots with preceding non-zero value
                # before are not pivots. above, the word row is referring to the column indices of the transpose of the
                # matrix, so it is actually the matrix's row, but refers to an index of a column of the transpose
                columns += indices

                # mark this row (column of transpose, row of matrix) as having a pivot
                pivot_row = index
                binary_matrix[i][index] = 1

                # if found a pivot in the last column (column of transpose, row of matrix) then there can be no more
                # pivots (since all pivots have to be below the row of the previous) so return transpose of matrix
                if index == len(array[0]) - 1:
                    return binary_matrix.transpose()
            return binary_matrix.transpose()

        # if didn't identify valid type to compare with, error thrown
        raise TypeError(f"equality between type(s) {type(self)} and {type(other)} is unsupported")

    def __ne__(self, other):

        # not equal, in all currently considered cases, is the inverse of equal
        result = self.__eq__(other)

        # if returned is a matrix, this is assumed to be a binary matrix and it is inverted
        if isinstance(result, Matrix):
            return binary_inverse(result, row_type=ArrayMod, obj_type=Matrix)

        # if matrix not returned, assumed to be bool value and is returned as not value
        return not result

    def __lt__(self, other):

        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r < other, self)))
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 < r2, self, other)))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __le__(self, other):

        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r <= other, self)))
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 <= r2, self, other)))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __gt__(self, other):

        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r > other, self)))
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 > r2, self, other)))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __ge__(self, other):

        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r >= other, self)))
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 >= r2, self, other)))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __mul__(self, other):

        # if number * matrix, return each element of matrix *= number
        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r * other, self)), *self.attr)
        # checks to see if dimensions of both objects are compatible with matrix multiplication
        if array_like(other):
            if len(self[0]) != len(other):
                raise ValueError(f"multiplication unsupported between matrix objects of dimension "
                                 f"{len(self)}x{len(self[0])} and {len(other)}x{len(other[0])}")

            return self.matmul(self, other)
        raise TypeError(f"multiplication unsupported for type(s): {type(self)} and {type(other)}")

    def __rmul__(self, other):

        # refer to __mul__ for documentation
        if isinstance(other, Number) and not isinstance(other, complex):
            return self.__mul__(other)
        if array_like(other):
            if len(other[0]) != len(self):
                raise ValueError(f"multiplication unsupported between matrix objects of dimension "
                                 f"{len(self)}x{len(self[0])} and {len(other)}x{len(other[0])}")

            return self.matmul(other, self)
        raise TypeError(f"multiplication unsupported for type(s): {type(self)} and {type(other)}")

    def __add__(self, other):

        # validates input as matrix that can be added if list-type
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 + r2, other, self)))
        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r + other, self)))
        raise TypeError(f"addition unsupported for type(s): {type(self)} and {type(other)}")

    def __radd__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            raise TypeError(f"scalar-matrix addition unsupported for scalar plus matrix")
        return self.__add__(other)

    def __sub__(self, other):
        # validates input as matrix that can be added if list-type
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 - r2, self, other)))
        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r - other, self)))
        raise TypeError(f"subtraction unsupported for type(s): {type(self)} and {type(other)}")

    def __rsub__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            raise TypeError(f"scalar-matrix subtraction unsupported for scalar minus matrix")
        # validates input as matrix that can be added if list-type
        if array_like(other):
            return self.construct(list(map(lambda r1, r2: r1 + r2, other, self)))
        raise TypeError(f"subtraction unsupported for type(s): {type(self)} and {type(other)}")

    def __floordiv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r // other, self)))
        raise TypeError(f"division unsupported for type(s): {type(self)} and {type(other)}")

    def __rfloordiv__(self, other):
        raise TypeError(f"matrix object cannot be divisor")

    def __truediv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return self.construct(list(map(lambda r: r / other, self)))
        raise TypeError(f"division unsupported for type(s): {type(self)} and {type(other)}")

    def __rtruediv__(self, other):
        raise TypeError(f"matrix object cannot be divisor")

    def __pow__(self, power, modulo=None):
        # powers currently supported are positive integers (referring to the number of times a matrix will
        # be multiplied against itself) and -1 (inverse)
        if not isinstance(power, int) or (power < 1 and power != -1):
            raise TypeError(f"exponent must be positive integer or -1")
        if len(self) != len(self[0]):
            raise ValueError(f"pow() unsupported for non-square matrices")
        if power == -1:
            return self.invert()
        product = self.copy()
        for _ in range(power - 1):
            product *= self
            if modulo is not None:
                product %= modulo
        return product

    def __mod__(self, other):
        return self.construct(list(map(lambda r: r % other, self)))

    # static methods used in matrix operations that allow subclasses to inherit more methods, using sub's statics
    @staticmethod
    def matmul(obj1, obj2):
        # attempts to transpose using built in methods, manually performs if no method exists
        transpose = getattr(obj2, 'transpose', None)
        if transpose is None:
            T = list(transpose_obj(obj2))
        else:
            T = transpose()
        return Matrix(list(map(lambda r: Array(map(lambda c: dot(c, r), T)), obj1)))

    @classmethod
    def construct(cls, *args):
        return cls(*args)

    # standard matrix methods
    def copy(self):
        """Returns exact copy of values of matrix in a new Matrix object."""

        if not self.array:
            return self.construct([], *self.attr)

        return self.construct(list(map(lambda r: r[:], self)), *self.attr)

    def __array_copy(self):
        # exact copy of array in list format
        return list(map(lambda r: r[:], self))

    def _transpose_obj(self, obj):
        return list(map(lambda i: Array(map(lambda r: r[i], obj)), range(len(obj[0]))))

    def _fraction(self):
        if self.array and all(isinstance(v, int) or (isinstance(v, float) and v.is_integer()) for v in self.flatten()):
            self.array = list(map(lambda r: Array(map(lambda e: Fraction(e), r)), self))

    def index(self, item):
        return self.array.index(item)

    def flatten(self):
        """Returns flat list containing all elements of matrix, going from left-to-right then
        top-to-bottom."""

        return [*reduce(lambda r1, r2: list(r1) + list(r2), self)] if self.array else []

    def tolist(self):
        return list(map(list, self.array))

    def append(self, row):
        if any(isinstance(e, (list, Array, ndarray)) for e in row):
            raise ValueError(f"operation unsupported for non-flat arrays. given: {repr(row)}")
        if self.array and len(row) != len(self[0]):
            raise ValueError(f"array length of {len(row)} incompatible with matrix of width {len(self[0])}")

        if isinstance(row, (ndarray, Matrix)):
            row = row.tolist()

        if isinstance(row, list):
            row = Array(row)

        if isinstance(row, Array) and not self.array:
            self.array.append(row)
        elif isinstance(row, Array) and len(row) == len(self[0]):
            self.array.append(row)
        else:
            raise TypeError(f"argument of type {type(row)} invalid for Matrix.append()")

    def invert(self):
        if len(self) != len(self[0]):
            raise AttributeError(f"inversion unsupported for non-square matrices")
        return self.construct(inv(np_array(self.tolist())).tolist(), *self.attr)  # using numpy.inv for efficiency

    def transpose(self):
        """Returns transpose of Matrix object."""

        if not self.array:
            return self.construct([], *self.attr)

        # maps indices of columns to inner map that returns the value at that index for each row
        return self.construct(list(map(lambda i: Array(map(lambda r: r[i], self)), range(len(self[0])))), *self.attr)

    def to_fraction(self):
        """Returns object with fractional instead of float values where possible."""
        return FractionMatrix(self.copy())

    def augment(self, solution):
        """Given matrix object and set of solutions, returns an augmented coefficient matrix
        with the set of solutions as the final column."""

        if self.augmented:
            raise AttributeError(f"{repr(self)} is already augmented coefficient matrix")

        self.attr = list(map(lambda e: True if isinstance(e, bool) else e, self.attr))  # only bool in self.attr is aug

        if isinstance(solution[0], list) and len(solution) == len(self):
            return self.construct(list(map(lambda r, s: r_append(r, s[0]), self, solution)), *self.attr)
        elif isinstance(solution[0], list) and len(solution[0]) == len(self):
            return self.construct(list(map(lambda r, s: r_append(r, s), self, solution[0])), *self.attr)
        elif len(solution) == len(self):
            return self.construct(list(map(lambda r, s: r_append(r, s), self, solution)), *self.attr)
        raise ValueError(f"augmentation of matrix of length {len(self)} unsupported for {repr(solution)}")

    def separate(self, index=-1):
        """Separates matrix at index, returning two matrices, first with indices [0, index) and
        second with indices [index:)."""

        self.attr = list(map(lambda e: False if isinstance(e, bool) else e, self.attr))  # only bool in self.attr is aug

        array_r = self.construct(list(map(lambda r: r[index:], self)), *self.attr)
        array_l = self.construct(list(map(lambda r: r[:index], self)), *self.attr)
        return array_l, array_r

    def remove_null_row(self, copy=False):
        """Removes rows consisting of just zeros."""

        if copy:
            return self.construct(
                reduce(lambda a, b: a if b.contains_only(0) else r_append(a, b), self, []), *self.attr
            )
        self.array = reduce(lambda a, b: a if b.contains_only(0) else r_append(a, b), self, [])

    def remove_null_column(self, copy=False):
        """Removes columns consisting of just zeros."""

        T = self.transpose()
        if self.augmented:
            T.remove_row(-1)
        if copy:
            return self.construct(
                reduce(lambda a, b: a if b.contains_only(0) else r_append(a, b), T, []), *self.attr
            ).transpose()

        self.array = list(
            self._transpose_obj(reduce(lambda a, b: a if b.contains_only(0) else r_append(a, b), T, []))
        )

    def remove_row(self, row, copy=False):

        row %= len(self)
        if copy:
            return self.construct(
                reduce(lambda a, b: a if b == row else r_append(a, self[b]), range(len(self)), []), *self.attr
            )

        self.array = reduce(
            lambda a, b: a if b == row else r_append(a, self[b]), range(len(self)), []
        )

    def remove_column(self, col, copy=False):

        T = self.transpose()
        col %= len(T)
        if copy:
            return self.construct(
                reduce(lambda a, b: a if b == col else r_append(a, T[b]), range(len(T)), []), *self.attr
            ).transpose()

        self.array = list(
            self._transpose_obj(
                reduce(lambda a, b: a if b == col else r_append(a, T[b]), range(len(T)), [])
            )
        )

    def row_reduce(self, row, col):
        """Subtracts each row besides given by given row until j'th element is zero. This function is the basis
        of Matrix.rref()"""
        for i in range(len(self)):
            if i != row:
                self[i] -= self[row] * self[i][col]

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

        array = self[:]  # deep copy

        adjust = 1 if self.augmented else 0  # if augmented don't reduce last column

        pivot_row = 0  # first pivot belongs in first row

        # if all integers, convert to fractions for more accurate rref
        # if all(isinstance(v, int) or (isinstance(v, float) and v.is_integer()) for v in self.flatten()):
        #     array = Matrix(list(map(lambda r: Array(map(lambda e: Fraction(e), r)), array)))

        for j in range(len(self[0]) - adjust):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(self)):

                # if non-zero element, this row can become pivot row
                if array[i][j] != 0:

                    # make j'th element the pivot, reducing rest of row as well
                    array[i].make_pivot(j)
                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row][:], array[i][:]
                    array.row_reduce(pivot_row, j)  # row reduce everything else
                    pivot_row += 1

                    break

        return array

    def is_rref_old(self):
        """Function that checks to see if matrix is in reduced row echelon form."""

        # adjust prevents iterating last column if augmented
        adjust = 1 if self.augmented else 0
        pivot_col, pivot_row = [], []
        for j in range(len(self[0]) - adjust):
            # iterate up from bottom row
            for i in range(len(self) - 1, -1, -1):
                e = self[i][j]
                # if non-zero entry in column with a pivot, this is not rref
                if e != 0 and j in pivot_col:
                    return False
                # if found a pivot in a row not already occupied by a pivot, then is legitimate pivot
                if e == 1 and i not in pivot_row:
                    pivot_col.append(j)
                    pivot_row.append(i)
                if e != 0 and i not in pivot_row:
                    return False
        return True

    def is_rref(self):
        columns = set()
        adj = 1 if self.augmented else 0
        for i in range(len(self)):
            for j in range(len(self[0]) - adj):
                e = self[i][j]
                if e not in (0, 1):  # if first non-zero entry in row is not 1, not in rref
                    return False
                if e == 1:
                    if columns and j <= max(columns):  # if pivot found is to left of prev pivot, not rref
                        return False
                    columns.add(j)
                    break

        for j, col in enumerate(self.transpose()):
            if j not in columns:  # if not a column w/ potential pivot, ignore
                continue
            indices = where(col != 0)[0]  # non-zero entries in the column
            if len(indices) > 1:  # should just be one, the pivot
                return False
        return True

    def is_consistent(self):
        """Function returns True if object is augmented and when in RREF has
        no non-zero solutions in rows without pivots, otherwise returns False."""

        if not self.augmented:
            raise AttributeError("This function is reserved for augmented coefficient matrices")
        matrix = self.copy()
        if not matrix.is_rref():
            matrix = matrix.rref()
        # iterate from bottom up to top row
        for i in range(len(matrix) - 1, -1, -1):
            for j in range(r := len(matrix[0])):
                e = matrix[i][j]
                # if there is a non-zero entry in this row, this row is consistent, move to next
                if e != 0 and j < r - 1:
                    break
                # if no pivot found in row and there exists non-zero entry in final column, row is inconsistent
                if e != 0 and j == r - 1:
                    return False
        return True

    def free_vars(self):
        """Function that checks if Matrix object has free variables, returning the number of free variables
        or 0 if none found."""

        pivots = self.pivot_count()
        adjust = 1 if self.augmented else 0
        return len(self[0]) - adjust - pivots

    def rank(self):
        return self.pivot_count()

    def null(self):
        return self.free_vars()

    def dimension(self):
        """Returns the dimension of the column span of the matrix. This is not necessarily the
        dimension of the domain. Use .domain_dimension to get the dimension of the domain."""

        adjust = 1 if self.augmented else 0
        return len(self[0]) - adjust

    def domain_dimension(self):
        return len(self)

    def pivot_count(self):
        matrix = self.copy()
        if not matrix.is_rref():
            matrix = matrix.rref()
        adjust = 1 if self.augmented else 0
        pivots = 0
        for i in range(len(matrix)):
            row_pivot = False
            for j in range(len(matrix[0]) - adjust):
                e = matrix[i][j]
                if e == 1:
                    row_pivot = True
            if row_pivot:
                pivots += 1
        return pivots

    def is_solvable(self):
        if not self.augmented:
            raise AttributeError("Solving matrix requires matrix to have augmented column of solutions")
        matrix = self.rref()

        # finds binary matrix representation, where each entry is 1 iff it is not 0 or a pivot
        matrix = matrix != [0, 'pivot']
        matrix = matrix.transpose()

        # removes the last row after transpose, thus removing the solutions
        matrix.remove_row(len(self[0]) - 1)

        # removes all fully-zero rows
        matrix.remove_null_row()

        # if nothing is left, this matrix has no free variables and is solvable if also consistent
        return len(matrix) == 0 and self.is_consistent()

    def solve(self):
        if not self.is_solvable():
            return False
        matrix = self.rref()
        matrix.remove_null_row()
        solutions = {}
        for i in range(len(self)):
            for j in range(c := (len(self[0]) - 1)):
                if self[i][j] == 1:
                    solutions[j] = self[i][c]
                    break
        return solutions

    def change_basis(self, basis):
        """Returns object previously in standard basis now in given basis"""

        if not isinstance(basis, Matrix):
            raise TypeError("basis must be a Matrix")
        return basis.invert() * self

    def revert_basis(self, basis):
        """Returns object previously in given basis now in standard basis"""

        if not isinstance(basis, Matrix):
            raise TypeError("basis must be a Matrix")
        return basis * self

    def column_li(self):
        """Determines if matrix consists of linearly dependent columns. Returns True if
        any column is linearly dependent, False otherwise."""

        array = self.copy()
        if not array.is_rref():
            array = array.rref()
        return array == matrix(rows=len(self), cols=len(self), diag=True)

    def to_vector(self):
        """Function that converts m x n Matrix into n columns vectors of
        length m."""

        columns = []
        for j in range(len(self[0])):
            vector = []
            for i in range(len(self)):
                vector.append([self[i][j]])
            columns.append(self.construct(vector))
        return columns

    def trace(self):
        """Calculates the tracer (diagonal sum) of a Matrix."""

        i = 0
        trace_sum = 0
        while i < len(self) and i < len(self[0]):
            trace_sum += self[i][i]
        return trace_sum

    def inner(self, other):
        """Calculates the inner product (dot product) of a vector with the instance matrix, returning
        a list of values representing the inner product of each column with the given vector."""

        if not (other := row_vector(other)) or len(self) != len(other):
            raise ValueError(f"multiplication unsupported between matrix objects of dimension "
                             f"{len(self)}x{len(self[0])} and {len(other)}x{len(other[0])}")

        # important to return python float or int instead of numpy, since numpy numbers will
        # override a matrix operator in operations where numpy object comes first
        if obj := row_vector(self):
            return python_number(dot(obj, other))
        return list(map(lambda r: python_number(dot(r, other)), obj))

    def norm(self):
        """Calculates the norm (distance from origin) of a vector"""

        if row_vector(self):
            obj = self
        elif row_vector(self[0]):
            obj = self[0]
        elif row_vector(T := self.transpose()[0]):
            obj = T
        else:
            raise ValueError(f"operation incompatible with object of dimension {len(self)}x{len(self[0])}")

        return sqrt(sum(map(lambda a: pow(a, 2), obj)))

    def orthogonal(self, others=None):
        """Determines if a list of vectors are orthogonal to each other, returning True only if
        all given vectors are orthogonal to all other given vectors.

        :param self: Matrix object, can be multi-dimensional or column vector
        :param others: list-type, containing column vectors of same length as self, can be Matrix, list, or numpy array
        """

        vectors = self.to_vector()

        if isinstance(others, Matrix):
            vectors += others.to_vector()

        # finds all two-element combinations of vectors w/o replacement
        vector_pairs = combinations(vectors, 2)

        for pair in vector_pairs:
            if Matrix.inner(pair[0], pair[1]) != 0:
                return False
        return True

    def coordinates(self, basis):
        """Returns the matrix representation of the coordinates of the given vector in the given basis."""

        # assumes instance matrix is matrix of column vectors
        matrix = self.transpose()

        if len(matrix[0]) == len(basis):
            basis = list(transpose_obj(basis))

        coord = []
        for i in range(len(matrix)):
            coord.append([])
            for j in range(len(basis)):
                coord[i].append(Matrix.inner(matrix[i], basis[j]))
        return self.construct(coord, *self.attr).transpose()

    def orthonormalize(self, steps=False):
        """Computes the orthonormal basis of a given basis using the Gram-Schmidt process. Assumes input
        basis is matrix of column vectors."""

        array = self.transpose()

        length = len(array[0])

        norm = Matrix.norm(array[0])
        e1 = array[0] / norm
        basis = [e1]
        if steps:

            # print statement to provide spacing for readability
            print()
            e1_str = str(FractionMatrix(matrix(array[0]).transpose())).split('\n')
            for j in range(length):
                if j == length // 2:
                    line = f'e1 = {fraction(1 / norm)} * '
                else:
                    line = ' ' * (8 + len(fraction(1 / norm)))
                line += e1_str[j]
                print(line)

            # print statement to provide spacing for readability
            print()
        for i, v in enumerate(array[1:]):

            # E is the sum of each already found orthonormal basis vector and its product with v
            E = sum(map(lambda e: Matrix.inner(v, e) * e, basis))

            # e_i_0 is the new basis vector that is orthogonal to all other basis vectors, but not currently
            # orthonormal, since it still must be divided by its norm
            e_i_0 = v - E
            if steps:

                # print statement to provide spacing for readability
                print()

                # list of each line of string for vector v
                v_str = str(FractionMatrix(matrix(v).transpose())).split('\n')

                # list of each dot product between v and each basis vector
                dot_str = list(map(lambda e: fraction(Matrix.inner(v, e)), basis))

                # list of each line of string for each basis vector e
                basis_str = list(map(lambda e: str(FractionMatrix(matrix(e).transpose())).split('\n'), basis))
                for j in range(length):
                    if j == length // 2:
                        line = f'e{i + 2}`= ' + v_str[j]
                        if j < length - 1:
                            line += ' '
                        line += ' -  [' + ' ' * 2
                        for k in range(d := len(dot_str)):
                            line += f'<{dot_str[k]}> * {basis_str[k][j]}' + ' '
                            if k < d - 1:
                                line += ' + '
                            else:
                                line += ' ' * 3
                    else:
                        line = ' ' * 5 + v_str[j]
                        if j < length - 1:
                            line += ' '
                        line += 4 * ' ' + '[' + ' ' * 2
                        for k in range(len(dot_str)):
                            line += (len(dot_str[k]) + 5) * ' ' + f'{basis_str[k][j]}' + ' ' * 3
                            if j < length - 1:
                                line += ' '
                    line += ']'
                    print(line)

                # print statement to provide spacing for readability
                print('\n')
                for j in range(length):
                    if j == length // 2:
                        line = f'e{i + 2}`= ' + v_str[j] + '  - '
                    else:
                        line = ' ' * 5 + v_str[j] + ' ' * 4
                    if j < length - 1:
                        line += ' '
                    sum_str = str(FractionMatrix(matrix(E).transpose())).split('\n')
                    line += sum_str[j]
                    print(line)

                # print statement to provide spacing for readability
                print('\n')
                result = str(FractionMatrix(matrix(e_i_0).transpose())).split('\n')
                for j in range(length):
                    if j == length // 2:
                        line = f'e{i + 2} = {fraction(1 / norm)} * '
                    else:
                        line = ' ' * (8 + len(fraction(1 / norm)))
                    line += result[j]
                    print(line)

                # print statement to provide spacing for readability
                print()

            norm = Matrix.norm(e_i_0)
            e_i = e_i_0 / norm
            basis.append(e_i)

        for i in range(len(basis)):
            basis[i] = list(basis[i])

        ortho_basis = matrix(basis).transpose()

        if steps:
            print('\n' + str(FractionMatrix(ortho_basis)) + '\n')

        return ortho_basis

    def eigvals(self):
        """Finds all eigenvalues for square matrix. Solves equation det(A - xI) = 0 with A
        equal to the instance matrix, I the identity, for x. The set of all solutions for x
        is analogous to the set of eigenvalues for A.

        :return: list of real or imaginary eigenvalues

        Examples
        --------
        >>> A = Matrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.eigvals()
        [1, 2, 3]

        >>> A = Matrix([[-1, -3, 1], [3, 3, 1], [3, 0, 3]])
        >>> A.eigvals()
        [0, (2.5-1.6583123951777j), (2.5+1.6583123951777j)]
        """

        if len(self) != len(self[0]):
            raise AttributeError(f"operation unsupported for non-square matrices")
        x = Symbol('x')
        return list(map(
            lambda e: int(e) if isinstance(e, float) and e.is_integer else e, map(
                lambda e: float(e) if im(e) == 0 else complex(e), solve(
                    self.char_poly(sym=x), x)
            )
        )  # inner map converts each root of solve into float if not complex
        )  # outer map converts each element of list into int if possible

    def eigvec(self, eigvals=None):
        if len(self) != len(self[0]):
            raise AttributeError(f"operation unsupported for non-square matrices")

        if eigvals is None:
            eigvals = self.eigvals()

        vectors = []
        I = matrix(rows=len(self), cols=len(self), diag=1)
        for e in eigvals:
            mat = self - (I * e)
            kern = mat.kernel().transpose()
            for v in kern:
                vectors.append(v)
        return self.construct(vectors, *self.attr).transpose()

    def minor(self, idx=None):
        """
        Computes the determinant of instance matrix if no index given. If index given,
        computes the minor A1,index (referring to the resulting matrix after removing
        row 1 and column index) multiplied against the value of A[0][index] with the
        correct sign (negative if index is even [when start counting at 1] otherwise
        positive). Use A.det() for calculating determinant for efficiency, unless specific
        minors or solving for sympy.Symbol is required.

        Examples
        --------
        >>> x = Symbol('x')
        >>> A = Matrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]]) - Matrix(rows=3, identity=x)
        >>> A.minor()
        -6*x + (3 - x)*(4 - x)*(-x - 1) + 18

        >>> A = Matrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.minor()
        6
        >>> A.det()
        6

        >>> A = Matrix([[-1, -3, 1], [3, 3, 1], [3, 0, 4]])
        >>> A.minor(1)
        27
        >>> B = Matrix([[3, 1], [3, 4]])
        >>> B.det()
        9

        Notes
        -----
        in the first example, 'x' is a sympy.Symbol and the solution given is solvable using sympy.solve()

        in the last two examples, A.minor(1) computes the A[0, 1] * determinant of A1,2 (A with row 1 and column 2
        removed, start counting at 1) while B.det() computes the determinant of B, which is A1,2. Since the sign
        of the minor alternates, A.minor(1) returns -3 * -1 * det(A1,2) = -3 * -1 * B.det() = 27
        """

        if len(self) == 2:
            return self[0][0] * self[1][1] - self[1][0] * self[0][1]

        det = 0
        if isinstance(idx, int):
            return pow(-1, idx) * self[0][idx] * self.remove_row(0, copy=True).remove_column(idx, copy=True).minor()

        sign = 1
        for j in range(len(self[0])):
            det += sign * self[0][j] * self.remove_row(0, copy=True).remove_column(j, copy=True).minor()
            sign *= -1

        if isinstance(det, float):
            return round(det) if abs(det - round(det)) < pow(10, -8) else det  # if essentially integer, return int
        return det

    def char_poly(self, sym='x'):
        """Computes the characteristic polynomial of a square matrix. Analogous
        to A.minor() for A = instance-matrix - Identity * x, for some variable x
        [typically x = sympy.Symbol('x')]."""

        assert_square(self)
        if not isinstance(sym, Symbol):
            sym = Symbol(sym)
        size = range(len(self))
        return self.construct(
            list(map(lambda i: Array(map(lambda j: self[i][j] - sym if i == j else self[i][j], size)), size))
        ).minor()

    def cross(self):
        """Given a list of n-1 vectors each of length n, returns 1 vector of length n that is
        orthogonal to all given vectors"""

        matrix = self.transpose() if len(self) == len(self[0]) + 1 else self.copy()  # if list is col vecs, transpose
        if len(self) == 1 and len(self[0]) == 2:
            return Array([self[0][1], -self[0][0]])
        return Array(map(lambda i: pow(-1, i) * matrix.remove_column(i, copy=True).minor(), range(len(matrix[0]))))

    def det(self):
        return python_number(np_det(np_array(self.tolist())))  # using numpy's det function for efficiency

    def kernel(self):
        """Computes the basis of the kernel for the given matrix."""

        if self.augmented:
            raise AttributeError(f"kernel computation uses implied set of solutions and does not support "
                                 f"augmented coefficient matrices")

        size = len(self)  # get number of rows
        array = self.copy() if self.is_rref() else self.rref()
        for j in range(l := len(self[0])):  # this loop appends identity matrix to bottom of instance matrix
            row = Array([0] * l)
            row[j] = 1
            array.append(row)

        array = array.transpose()

        pivot_row = 0  # first pivot belongs in first row

        for j in range(size):  # iterate only over current matrix, not attached identity matri

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(array)):

                # if non-zero element, this row can become pivot row
                if array[i][j] != 0:

                    # make j'th element the pivot, reducing rest of row as well
                    array[i].make_pivot(j)
                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                    array.row_reduce(pivot_row, j)  # row reduce everything else
                    pivot_row += 1

        array, kern = array.separate(size)  # separates original matrix from now modified identity matrix
        basis = Matrix([])

        for i, row in enumerate(array):
            if row.contains_only(0):  # all null rows in original matrix correspond to basis vector for kernel
                basis.append(kern[i])

        return basis.transpose()  # basis is list of rows, transpose into standard of column vectors


class ModularMatrix(Matrix):
    def __init__(self, array, mod, aug=False):
        self.mod = mod
        super().__init__(array, aug)
        self.attr = [self.mod] + self.attr

    def __repr__(self):
        return f"{__class__.__name__}({self.tolist()}, {self.mod}, aug={self.augmented})"

    @staticmethod
    def matmul(obj1, obj2):

        mod = obj1.mod if isinstance(obj1, ModularMatrix) else obj2._mod
        # attempts to transpose using built in methods, manually performs if no method exists
        transpose = getattr(obj2, 'transpose', None)
        if transpose is None:
            T = list(transpose_obj(obj2))
        else:
            T = transpose()
        return ModularMatrix(list(map(lambda r: ArrayMod(map(lambda c: dot(c, r), T), mod), obj1)), mod)

    def _transpose_obj(self, obj):
        return list(map(lambda i: ArrayMod(map(lambda r: r[i], obj), self.mod), range(len(obj[0]))))

    def _fraction(self):
        pass

    def append(self, row):
        if any(isinstance(e, (list, Array, ndarray)) for e in row):
            raise ValueError(f"operation unsupported for non-flat arrays. given: {repr(row)}")
        if self.array and len(row) != len(self[0]):
            raise ValueError(f"array length of {len(row)} incompatible with matrix of width {len(self[0])}")

        if isinstance(row, (ndarray, Matrix)):
            row = row.tolist()

        if isinstance(row, list):
            row = ArrayMod(row, self.mod)

        if isinstance(row, ArrayMod) and not self.array:
            self.array.append(row)
        elif isinstance(row, ArrayMod) and len(row) == len(self[0]) and row.mod == self.mod:
            self.array.append(row)
        else:
            raise TypeError(f"argument {repr(row)} invalid for Matrix.append()")

    def invert(self):
        raise AttributeError(f"{__class__.__name__} does not currently support inversion")

    def transpose(self):
        if not self.array:
            return self.construct([], *self.attr)

        # maps indices of columns to inner map that returns the value at that index for each row
        return self.construct(
            list(map(lambda i: ArrayMod(map(lambda r: r[i], self), self.mod), range(len(self[0])))), *self.attr
        )

    def rref(self):
        matrix = self[:]  # deep copy

        adjust = 1 if self.augmented else 0  # if augmented don't reduce last column

        pivot_row = 0  # first pivot belongs in first row
        for j in range(len(self[0]) - adjust):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(self)):

                # if non-zero element, this row can become pivot row
                if matrix[i][j] != 0 and (res := matrix[i].make_pivot(j, copy=True)):

                    # make j'th element the pivot, reducing rest of row as well
                    matrix[i] = res

                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        matrix[i], matrix[pivot_row] = matrix[pivot_row][:], matrix[i][:]

                    matrix.row_reduce(pivot_row, j)  # row reduce everything else

                    pivot_row += 1

        return matrix

    def orthonormalize(self, steps=False):
        pass

    def eigvals(self):
        pass

    def cross(self):
        pass

    def kernel(self):
        if self.augmented:
            raise AttributeError(f"kernel computation uses implied set of solutions and does not support "
                                 f"augmented coefficient matrices")

        size = len(self)  # get number of rows
        array = self.copy() if self.is_rref() else self.rref()
        for j in range(l := len(self[0])):  # this loop appends identity matrix to bottom of instance matrix
            row = ArrayMod([0] * l, self.mod)
            row[j] = 1
            array.append(row)

        array = array.transpose()
        pivot_row = 0  # first pivot belongs in first row
        for j in range(size):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(array)):

                # if non-zero element, this row can become pivot row
                if array[i][j] != 0 and (res := array[i].make_pivot(j, copy=True)):

                    # make j'th element the pivot, reducing rest of row as well
                    array[i] = res

                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                    array.row_reduce(pivot_row, j)  # row reduce everything else

                    pivot_row += 1

        mat, kern = array.separate(size)  # separates original matrix from now modified identity matrix

        if not mat.is_rref():
            print(array.transpose())
            raise ValueError(f"computation of kernel failed because matrix could not be properly row reduced. "
                             f"above is result of attempted row reduction")

        basis = ModularMatrix([], self.mod)

        for i, row in enumerate(mat):
            if row.contains_only(0):  # all null rows in original matrix correspond to basis vector for kernel
                basis.append(kern[i])

        return basis.transpose()  # basis is list of rows, transpose into standard of column vectors


class BinaryMatrix(ModularMatrix):
    def __init__(self, array, aug=False):
        super().__init__(array, 2, aug)
        self.attr = [aug]

    @staticmethod
    def matmul(obj1, obj2):
        mod = obj1.mod if isinstance(obj1, ModularMatrix) else obj2._mod
        # attempts to transpose using built in methods, manually performs if no method exists
        transpose = getattr(obj2, 'transpose', None)
        if transpose is None:
            T = list(transpose_obj(obj2))
        else:
            T = transpose()
        if isinstance(obj2, ModularMatrix):
            return ModularMatrix(list(map(lambda r: ArrayMod(map(lambda c: dot(c, r), T), mod), obj1)), obj2.mod)
        if isinstance(obj2, Matrix):
            return Matrix(list(map(lambda r: Array(map(lambda c: dot(c, r), T)), obj1)))
        return BinaryMatrix(list(map(lambda r: ArrayMod(map(lambda c: dot(c, r), T), mod), obj1)))

    def union(self, other):
        """Returns the logical union of two binary matrices.

        If this operation is needed with non-matrix objects, use set.union(A, B)"""

        if not is_binary_matrix(self) and not is_binary_matrix(other):
            raise AttributeError(f"union between {repr(self)} and {repr(other)} is unsupported")
        if len(self) != len(other) or len(self[0]) != len(other[0]):
            raise AttributeError(f"operation unsupported for objects with dimensions {len(self)}x{len(self[0])} "
                                 f"and {len(other)}x{len(other[0])}")
        return BinaryMatrix(list(map(lambda r1, r2: ArrayMod(map(lambda e1, e2: e1 ^ e2, r1, r2), 2), self, other)))

    def intersection(self, other):
        """Returns the logical intersection of two binary matrices.

        If this operation is needed with non-Matrix objects, use set.intersection(A, B)"""

        if not is_binary_matrix(self) and not is_binary_matrix(other):
            raise AttributeError(f"union between {repr(self)} and {repr(other)} is unsupported")
        if len(self) != len(other) or len(self[0]) != len(other[0]):
            raise AttributeError(f"operation unsupported for objects with dimensions {len(self)}x{len(self[0])} "
                                 f"and {len(other)}x{len(other[0])}")

        return BinaryMatrix(list(map(lambda r1, r2: ArrayMod(map(lambda e1, e2: e1 & e2, r1, r2), 2), self, other)))

    def disjunction(self, other):
        """Returns the logical intersection of two binary matrices.

        If this operation is needed with non-Matrix objects, use set(A) - set(B)"""

        if not is_binary_matrix(self) and not is_binary_matrix(other):
            raise AttributeError(f"union between {repr(self)} and {repr(other)} is unsupported")
        if len(self) != len(other) or len(self[0]) != len(other[0]):
            raise AttributeError(f"operation unsupported for objects with dimensions {len(self)}x{len(self[0])} "
                                 f"and {len(other)}x{len(other[0])}")

        return BinaryMatrix(list(map(lambda r1, r2: ArrayMod(map(lambda e1, e2: e1 | e2, r1, r2), 2), self, other)))

    def binary_inverse(self):
        """Returns the negative image of a binary matrix"""

        if not is_binary_matrix(self):
            raise AttributeError(f"binary inverse unsupported for {repr(self)}")

        return BinaryMatrix(list(map(lambda r: ArrayMod(map(lambda e: e + 1, r), 2), self)))

    def row_reduce(self, row, col):
        for i in range(len(self)):
            if self[i][col] == 1 and i != row:
                self[i] += self[row]

    def rref(self):
        matrix = self[:]  # deep copy

        adjust = 1 if self.augmented else 0  # if augmented don't reduce last column

        pivot_row = 0  # first pivot belongs in first row
        for j in range(len(self[0]) - adjust):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(self)):

                # if non-zero element, this row can become pivot row
                if matrix[i][j] == 1:

                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        matrix[i], matrix[pivot_row] = matrix[pivot_row][:], matrix[i][:]

                    matrix.row_reduce(pivot_row, j)  # row reduce everything else

                    pivot_row += 1

        return matrix

    def kernel(self):
        if self.augmented:
            raise AttributeError(f"kernel computation uses implied set of solutions and does not support "
                                 f"augmented coefficient matrices")

        size = len(self)  # get number of rows
        array = self.copy()
        for j in range(l := len(self[0])):  # this loop appends identity matrix to bottom of instance matrix
            row = ArrayMod([0] * l, 2)
            row[j] = 1
            array.append(row)

        array = array.transpose()
        pivot_row = 0  # first pivot belongs in first row
        for j in range(size):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(array)):

                # if non-zero element, this row can become pivot row
                if array[i][j] == 1:

                    if i > pivot_row:  # if pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                    array.row_reduce(pivot_row, j)  # row reduce everything else

                    pivot_row += 1

        mat, kern = array.separate(size)  # separates original matrix from now modified identity matrix

        if not mat.is_rref():
            print(array.transpose())
            raise ValueError(f"computation of kernel failed because matrix could not be properly row reduced. "
                             f"above is result of attempted row reduction")

        basis = BinaryMatrix([])

        for i, row in enumerate(mat):
            if row.contains_only(0):  # all null rows in original matrix correspond to basis vector for kernel
                basis.append(kern[i])

        return basis.transpose()  # basis is list of rows, transpose into standard of column vectors


# this class is used just for printing matrices with fractions instead of floats, no mathematical
# operations can be performed with this class currently
class FractionMatrix:
    def __init__(self, array, limit=75):

        # attempts to convert each number into a fraction before adding to matrix
        fraction_matrix = []
        for i in range(len(array)):
            fraction_matrix.append([])
            for j in range(len(array[0])):
                e = python_number(array[i][j])
                if isinstance(e, int) or (isinstance(e, float) and e.is_integer()):
                    fraction_matrix[i].append(str(int(e)))
                else:
                    fraction_matrix[i].append(fraction(e, limit))
        self.array = fraction_matrix

    def __str__(self):
        max_len = 0
        for i in range(len(self.array)):
            for j in range(len(self.array[0])):
                if len(self.array[i][j]) > max_len:
                    max_len = len(self.array[i][j])

        padding = (max_len + 3) | 1
        matrix = '['
        for i in range(c := len(self.array)):
            matrix += '['
            for j in range(len(self.array[0])):

                # having left padding at exactly this ensures all fractions have the numerators line up
                # when ignoring any negative signs
                e = str(self.array[i][j])
                if '/' in e:
                    left_pad = 1 if e[0] == '-' else 2

                # non-fractions use normal padding to center them
                else:
                    left_pad = (padding - len(e)) // 2
                right_pad = padding - len(e) - left_pad
                matrix += left_pad * ' ' + f'{e}' + ' ' * right_pad
            if i < c - 1:
                matrix += ']\n '
            else:
                return matrix + ']]'

    def __getitem__(self, item):
        return self.array[item]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return iter(self.array)
