from .array import Array, ArrayMod, where
from .tools import dot, isnumber, python_number, string_reduce, fraction, append_and_return
from random import randrange
from numpy import ndarray
from numpy import array as np_array
from sympy import Symbol, solve, im, evalf, sympify, N
from numpy.linalg import inv
from numpy.linalg import det as np_det
from math import sqrt
from functools import reduce, wraps
from itertools import combinations


def array_like(obj):
    if isinstance(obj, (list, Matrix, Array, ndarray)):
        return True
    return getattr(obj, '__getitem__', None) is not None and getattr(obj[0], '__getitem__', None) is not None


def transpose_obj(obj, row_type=None, obj_type=None):
    """Transpose method that computes transpose for array_like objects without transpose
    methods."""

    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    if row_type == ArrayMod:
        if hasattr(obj, 'mod'):
            mod = obj.mod
        elif hasattr(obj[0], 'mod'):
            mod = obj[0].mod
        else:
            raise AttributeError(f"{repr(obj)} has no associated attribute mod")
        return obj_type(map(lambda i: ArrayMod(map(lambda r: r[i], obj), mod), range(len(obj[0]))))
    return obj_type(map(lambda i: row_type(map(lambda r: r[i], obj)), range(len(obj[0]))))


def conserve_attributes(func):

    @wraps(func)
    def new_func(*args, **kwargs):
        if not isinstance(args[0], Matrix):
            raise AttributeError(f"conserve_attributes wrapper unsupported for functions outside of Matrix class")
        result = func(*args, **kwargs)
        if result is not None:
            result.mod = args[0].mod
            result.augmented = args[0].augmented
        return result

    return new_func


def matmul(obj1, obj2, row_type=None, obj_type=None, mod=None):
    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    # attempts to transpose using built in methods, manually performs if no method exists
    transpose = getattr(obj2, 'transpose', None)
    if transpose is None:
        T = []
        for j in range(len(obj2[0])):
            T.append([])
            for i in range(len(obj2)):
                T[j].append(obj2[i][j])
    else:
        T = transpose()
    if isnumber(mod):
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


class Matrix:
    def __init__(self, array=None, rows=None, cols=None, rand=False, identity=False, aug=False,
                 solution=None, mod=None, dtype: type = None):
        self.mod = mod
        self.augmented = False

        # if augmented is True or solution provided, matrix is augmented
        if aug or solution is not None:
            self.augmented = True

        # dtype allows for faster creation of matrix, assumes input is reliable and correct type
        if dtype is not None:
            if dtype not in (list, Array):
                raise TypeError(f"dtype given not recognized: {dtype}")
            if dtype == list:
                self.array = list(map(Array, array))
            else:
                self.array = array
            if self.mod is not None:
                self.array %= self.mod
            return

        # if array given, check if list or numpy array (cannot inherit a Matrix object)
        if array is not None:
            if isinstance(array, map):
                array = list(array)
            if isinstance(array, list):

                # empty matrix
                if not array:
                    self.array = array
                # if array given is nested list, set self.array to given, otherwise, nest it
                elif isinstance(array[0], list) and self.mod is not None:
                    self.array = list(map(lambda a: ArrayMod(a, self.mod), array))
                elif isinstance(array[0], list):
                    self.array = list(map(Array, array))
                elif isinstance(array[0], Array):
                    self.array = array
                else:
                    raise TypeError(f"{repr(array)} is not array_like")
            elif isinstance(array, Array):
                self.array = [array]
            else:
                raise TypeError(f"{repr(array)} is not array_like")

        # if random, no rows or columns required, creates m x n matrix with random m, n (where not given) and
        # values random integers between -50, 50
        elif rand:
            if rows is None:
                rows = randrange(1, 10)
            if cols is None:
                cols = randrange(1, 10)
            matrix = []
            for i in range(rows):
                matrix.append(ArrayMod([], self.mod) if self.mod is not None else Array([]))
                for j in range(cols):
                    matrix[i].append(randrange(-50, 50))
            self.array = matrix

        # if identity, matrix must be square
        elif identity:
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
            matrix = []
            for i in range(rows):
                matrix.append(ArrayMod([], self.mod) if self.mod is not None else Array([]))
                for j in range(cols):
                    if i == j:
                        matrix[i].append(identity if isinstance(identity, Symbol) else int(identity))
                    else:
                        matrix[i].append(0)
            self.array = matrix

        # if both rows and cols are ints, make empty matrix
        elif isinstance(rows, int) and isinstance(cols, int):
            matrix = []
            for i in range(rows):
                matrix.append(ArrayMod([0] * cols, self.mod) if self.mod is not None else Array([0] * cols))
            self.array = matrix

        # if a solution is provided, attach column to end of matrix
        if solution is not None:
            if not isinstance(solution, (list, Matrix, Array)):
                raise TypeError(f"{repr(solution)} is not array_like")
            if isinstance(solution[0], (list, Array)) and len(solution) == len(self):
                for i in range(len(self)):
                    self.array[i].append(solution[i][0])
            elif isinstance(solution[0], (list, Array)) and len(solution) == 1 and len(solution[0]) == len(self):
                for i in range(len(self)):
                    self.array[i].append(solution[0][i])
            elif isinstance(solution, (list, Array)) and len(solution) == len(self):
                for i in range(len(self)):
                    self.array[i].append(solution[i])
            else:
                raise TypeError(f"{repr(solution)} is incompatible with matrix of dimension {len(self)}x{len(self[0])}")

        if self.mod is not None:
            self.array = list(map(lambda r: r % self.mod, self.array))

        self.row_type = type(self.array[0]) if self.array else Array

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
        return f"{__class__.__name__}(array={self.array}, aug={self.augmented}, mod={self.mod})"

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
                formatted += pad_left * " " + f"{e}" + " " * pad_right
            if i == l - 1:
                formatted += "]"
            else:
                formatted += "]\n"
        return formatted + "]"

    def __eq__(self, other):

        # if other is nested list i.e. 2-dimensional array, return False if values at each corresponding row, col
        # are not similar, True otherwise (returning True if all values within 0.1 of each other because numpy
        # rounding /dealing with long floats can result in the solution to Ax = B for some x might be correct
        # solution but slightly off float values
        if array_like(other) and isinstance(other[0], (list, Array, ndarray)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError(f"equality between {repr(self)} and {repr(other)} is unsupported")
            for i in range(len(self)):
                for j in range(len(self[0])):
                    if abs(self[i][j] - other[i][j]) > 0.1:
                        return False
            return True

        # if not nested list, two objects will not be compared as arrays, but will return a binary matrix
        # with truth values wherever self's values are in the list
        if isinstance(other, (list, set, tuple)):
            binary_matrix = Matrix(map(lambda r: ArrayMod(map(lambda e: 1 if e in other else 0, r), 2), self))

            # pivot keyword allows user to search for pivots in matrix
            if 'pivot' in other:
                # returns the union of searching for pivot
                return union(binary_matrix, self == 'pivot', row_type=ArrayMod, obj_type=Matrix)
            return binary_matrix

        # if comparing to a number, return binary matrix with values = 1 wherever matrix values = number
        if isnumber(other):
            return Matrix(map(lambda r: ArrayMod(map(lambda e: 1 if e == other else 0, r), 2), self))

        # if user searching for pivots, return binary matrix with values = 1 wherever an entry = 1
        # that is not in a row with a preceding value != 0 or 1 (first non-zero entry in row should be pivot)
        # also checks if entry = 1 is in same column as already identified pivot
        if other == 'pivot':
            binary_matrix = Matrix(rows=len(self[0]), cols=len(self), mod=2)
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
                if len(indices[0]) != 1 or index <= pivot_row or index in columns:
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

        if isnumber(other):
            return Matrix(map(lambda r: r < other, self))
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 < r2, self, other))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __le__(self, other):

        if isnumber(other):
            return Matrix(map(lambda r: r <= other, self))
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 <= r2, self, other))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __gt__(self, other):

        if isnumber(other):
            return Matrix(map(lambda r: r > other, self))
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 > r2, self, other))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __ge__(self, other):

        if isnumber(other):
            return Matrix(map(lambda r: r >= other, self))
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 >= r2, self, other))
        raise TypeError(f"operation unsupported for type(s): {type(self)} and {type(other)}")

    def __mul__(self, other):

        # if number * matrix, return each element of matrix *= number
        if isnumber(other):
            return Matrix(map(lambda r: r * other, self))
        # checks to see if dimensions of both objects are compatible with matrix multiplication
        if array_like(other):
            if len(self[0]) != len(other):
                raise ValueError(f"multiplication unsupported between matrix objects of dimension "
                                 f"{len(self)}x{len(self[0])} and {len(other)}x{len(other[0])}")

            return matmul(self, other, row_type=type(self[0]), obj_type=Matrix, mod=self.mod)
        raise TypeError(f"multiplication unsupported for type(s): {type(self)} and {type(other)}")

    def __rmul__(self, other):

        # refer to __mul__ for documentation
        if isinstance(other, Matrix) or isnumber(other):
            return Matrix.__mul__(self, other)
        if array_like(other):
            if len(self[0]) != len(other):
                raise ValueError(f"multiplication unsupported between matrix objects of dimension "
                                 f"{len(self)}x{len(self[0])} and {len(other)}x{len(other[0])}")

            return matmul(other, self, row_type=type(self[0]), obj_type=Matrix, mod=self.mod)
        raise TypeError(f"multiplication unsupported for type(s): {type(self)} and {type(other)}")

    def __add__(self, other):

        # validates input as matrix that can be added if list-type
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 + r2, other, self))
        if isnumber(other):
            return Matrix(map(lambda r: r + other, self))
        raise TypeError(f"addition unsupported for type(s): {type(self)} and {type(other)}")

    def __radd__(self, other):
        if isnumber(other):
            raise TypeError(f"scalar-matrix addition unsupported for scalar + matrix")
        return self.__add__(other)

    def __sub__(self, other):
        # validates input as matrix that can be added if list-type
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 - r2, self, other))
        if isnumber(other):
            return Matrix(map(lambda r: r - other, self))
        raise TypeError(f"subtraction unsupported for type(s): {type(self)} and {type(other)}")

    def __rsub__(self, other):
        if isnumber(other):
            raise TypeError(f"scalar-matrix subtraction unsupported for scalar - matrix")
        # validates input as matrix that can be added if list-type
        if array_like(other):
            return Matrix(map(lambda r1, r2: r1 + r2, other, self))
        raise TypeError(f"subtraction unsupported for type(s): {type(self)} and {type(other)}")

    def __floordiv__(self, other):
        if isnumber(other):
            return Matrix(map(lambda r: r // other, self))
        raise TypeError(f"division unsupported for type(s): {type(self)} and {type(other)}")

    def __truediv__(self, other):
        if isnumber(other):
            return Matrix(map(lambda r: r / other, self))
        raise TypeError(f"division unsupported for type(s): {type(self)} and {type(other)}")

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
        return Matrix(map(lambda r: r % other, self))

    # pseudo-constructor class methods that return a new instance of Matrix specified by the method
    @classmethod
    def make(cls, rows, cols, aug=False, by_row=True):
        """
        Class method that takes number of rows and columns as input and returns a new instance of
        a Matrix object with the values given by user.

        :param aug : bool; determines if matrix is an augmented coefficient matrix or not.
        :param by_row : bool; determines if input from user is collected each row at a time, or each element at a time.
        """

        def get_number(n):
            if "," in n:
                c = n.split(",")
                real = float(c[0])
                imag = float(c[1])
                return complex(real, imag)
            return float(n) if "." in n else int(n)

        array = []
        for i in range(rows):
            array.append(Array([]))
            if not by_row:
                for j in range(cols):
                    n = input(f"{i + 1},{j + 1}: ")
                    array[i].append(get_number(n))
            else:
                while True:
                    try:
                        row = input(f"row {i + 1}: ")
                        values = row.split()
                        if len(values) != cols:
                            raise ValueError(f"entered {len(values)} values. this matrix is expecting {cols}")
                    except ValueError:
                        continue
                    else:
                        array[i] = Array(map(get_number, values))
                        break
        return cls(array=array, aug=aug)

    @classmethod
    def diag(cls, scalars):
        matrix = cls(rows=len(scalars), cols=len(scalars), identity=scalars)
        for i in range(len(scalars)):
            matrix[i][i] = scalars[i]
        return matrix

    @classmethod
    def make_basis(cls, in_basis=None, out_basis=None):
        """
        make_basis creates a new Matrix object which represents the change-of-basis matrix between
        the given in_basis and out_basis if both are given. in_basis determines the basis that this
        change-of-basis matrix can be applied to. out_basis determines the basis that the output of
        applying this matrix to a vector will be given in.

        ex.
            in_basis = Matrix([[1, 2], [2, -1]])
            out_basis = Matrix([[1, 1], [1, -1]])

            With just in_basis, this change of basis matrix would be the in_basis. This converts each
            input matrix from the same basis as in_basis into the standard basis for output.
            With just out_basis, this change of basis matrix would be the inverse of out_basis. Applying
            this simply returns the input matrix given in the basis of out_basis.
            With both, this change of basis matrix is the inverse of in_basis multiplied against the
            out_basis, making the matrix first convert input from in_basis into standard basis, then
            from standard basis into out_basis, except this is done in one step since both steps are
            combined into the multiplication of out_basis.invert() * in_basis.

        Parameters
        ----------
        out_basis : Matrix
            a square matrix that represents a given basis for an n-dimensional vector space, this basis will be the basis
        of matrices that this matrix will be applied to

        in_basis : Matrix
            a square matrix that represents a given basis for an n-dimensional vector space, this basis will be the basis
        of matrices that are produced by applying this matrix to another
        """
        if in_basis is None and out_basis is None:
            raise ValueError("make_basis requires at least one basis")
        if in_basis is None:
            if len(out_basis) != len(out_basis[0]):
                raise ValueError("basis must be square matrix")
            return out_basis.invert()
        if out_basis is None:
            if len(in_basis) != len(in_basis[0]):
                raise ValueError("basis must be square matrix")
            return in_basis
        if not all(length == len(in_basis) for length in [len(out_basis), len(out_basis[0]), len(in_basis[0])]):
            raise ValueError("bases must be square matrices")
        return out_basis.invert() * in_basis

    # standard matrix methods
    @conserve_attributes
    def copy(self):
        """Returns exact copy of values of matrix in a new Matrix object."""

        return Matrix(map(lambda r: r[:], self))

    def __array_copy(self):
        # exact copy of array in list format
        return list(map(lambda r: r[:], self))

    def index(self, item):
        return self.array.index(item)

    def flatten(self):
        """Returns flat list containing all elements of matrix, going from left-to-right then
        top-to-bottom."""

        return [*reduce(lambda r1, r2: list(r1) + list(r2), self)]

    def tolist(self):
        return list(map(lambda r: list(r), self))

    def append(self, row):
        if any(isinstance(e, (list, Array, ndarray)) for e in row):
            raise ValueError(f"operation unsupported for non-flat arrays. given: {repr(row)}")
        if self.array and len(row) != len(self[0]):
            raise ValueError(f"array length of {len(row)} incompatible with matrix of width {len(self[0])}")
        if isinstance(row, (ndarray, Matrix)):
            row = row.tolist()
        if isinstance(row, list) and self.row_type == ArrayMod:
            row = ArrayMod(row, self.mod)
        elif isinstance(row, list):
            row = Array(row)
        if isinstance(row, ArrayMod) and self.row_type != ArrayMod:
            raise TypeError(f"appending ArrayMod to matrix with rows of type Array is unsupported")
        elif isinstance(row, Array) and not self.array:
            self.array.append(row)
        elif isinstance(row, Array) and len(row) == len(self[0]):
            self.array.append(row)
        else:
            raise TypeError(f"argument of type {type(row)} invalid for Matrix.append()")

    def invert(self):
        if len(self) != len(self[0]):
            raise AttributeError(f"inversion unsupported for non-square matrices")
        return Matrix(inv(np_array(self.tolist())).tolist(), mod=self.mod)  # using numpy.inv for efficiency

    def transpose(self):
        """Returns transpose of Matrix object."""

        # maps indices of columns to inner map that returns the value at that index for each row
        return Matrix(map(lambda i: ArrayMod(
            map(lambda r: r[i], self), self.mod) if self.mod is not None else Array(
            map(lambda r: r[i], self)), range(len(self[0]))), mod=self.mod)

    def to_fraction(self):
        """Returns object with fractional instead of float values where possible."""
        return FractionMatrix(self.copy())

    def augment(self, solution):
        """Given matrix object and set of solutions, returns an augmented coefficient matrix
        with the set of solutions as the final column."""

        if isinstance(solution[0], list) and len(solution) == len(self):
            return Matrix(map(lambda r, s: append_and_return(r, s[0]), self, solution), aug=True, mod=self.mod)
        elif isinstance(solution[0], list) and len(solution[0]) == len(self):
            return Matrix(map(lambda r, s: append_and_return(r, s), self, solution[0]), aug=True, mod=self.mod)
        elif len(solution) == len(self):
            return Matrix(map(lambda r, s: append_and_return(r, s), self, solution), aug=True, mod=self.mod)
        raise ValueError(f"augmentation of matrix of length {len(self)} unsupported for {repr(solution)}")

    def separate(self):
        """Function used to separate an augmented coefficient matrix into
        a standard coefficient matrix and a column matrix of solutions"""

        if self.mod is not None:
            return Matrix(map(lambda r: ArrayMod(r[-1], self.mod), self), aug=False, mod=self.mod)
        return Matrix(map(lambda r: Array(r[-1]), self), aug=False, mod=self.mod)

    @conserve_attributes
    def remove_null(self, copy=False):
        """Removes all null rows and columns from matrix."""

        if copy:
            return self.remove_null_column(copy=True).remove_null_row(copy=True)

        self.remove_null_column()
        self.remove_null_row()

    @conserve_attributes
    def remove_null_row(self, copy=False):
        """Removes rows consisting of just zeros."""

        if copy:
            return Matrix(reduce(lambda a, b: a if b.contains_only(0) else append_and_return(a, b), self, []))
        self.array = reduce(lambda a, b: a if b.contains_only(0) else append_and_return(a, b), self, [])

    @conserve_attributes
    def remove_null_column(self, copy=False):
        """Removes columns consisting of just zeros."""

        T = self.transpose()
        if self.augmented:
            T.remove_row(-1)
        if copy:
            return Matrix(reduce(
                lambda a, b: a if b.contains_only(0) else append_and_return(
                    a, b), T, []), aug=self.augmented, mod=self.mod).transpose()
        self.array = transpose_obj(
            reduce(lambda a, b: a if b.contains_only(0) else append_and_return(a, b), T, []),
            row_type=self.row_type
        )

    @conserve_attributes
    def remove_row(self, row, copy=False):
        row %= len(self)
        if copy:
            return Matrix(reduce(lambda a, b: a if b == row else append_and_return(a, self[b]), range(len(self)), []))
        self.array = reduce(lambda a, b: a if b == row else append_and_return(a, self[b]), range(len(self)), [])

    @conserve_attributes
    def remove_column(self, col, copy=False):
        T = self.transpose()
        col %= len(T)
        if copy:
            return Matrix(
                reduce(lambda a, b: a if b == col else append_and_return(a, T[b]), range(len(T)), [])).transpose()
        self.array = transpose_obj(
            reduce(lambda a, b: a if b == col else append_and_return(a, T[b]), range(len(T)), []))

    def row_reduce(self, row, col):
        """Subtracts each row besides given by given row until j'th element is zero. This function is the basis
        of Matrix.rref()"""
        for i in range(len(self)):
            if i != row:
                self[i] -= self[row] * self[i][col]

    @conserve_attributes
    def rref(self):
        """Uses Guassian elimination to row reduce matrix, returning a new instance of a matrix in reduced row
        echelon form."""

        matrix = self[:]

        # if augmented don't reduce last column
        adjust = 1 if self.augmented else 0

        # first pivot belongs in first row
        pivot_row = 0
        for j in range(len(self[0]) - adjust):

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(self)):

                # if non-zero element, this row can become pivot row
                if matrix[i][j] != 0:

                    # make j'th element the pivot, reducing rest of row as well
                    matrix[i].make_pivot(j)

                    # if pivot row not already in correct position, swap
                    if i > pivot_row:
                        matrix[i], matrix[pivot_row] = matrix[pivot_row][:], matrix[i][:]

                    # row reduce everything else
                    matrix.row_reduce(pivot_row, j)
                    pivot_row += 1

        return matrix

    def is_rref(self):
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

    @conserve_attributes
    def change_basis(self, basis):
        """Returns object previously in standard basis now in given basis"""

        if not isinstance(basis, Matrix):
            raise TypeError("Basis must be a Matrix")
        return basis.invert() * self

    @conserve_attributes
    def revert_basis(self, basis):
        """Returns object previously in given basis now in standard basis"""

        if not isinstance(basis, Matrix):
            raise TypeError("Basis must be a Matrix")
        return basis * self

    def is_basis(self):
        if len(self) != len(self[0]):
            return False
        matrix = self.copy()
        if not matrix.is_rref():
            matrix = matrix.rref()
        return matrix == Matrix(rows=len(self), cols=len(self), identity=True)

    def to_vector(self):
        """Function that converts m x n Matrix into n columns vectors of
        length m."""

        columns = []
        for j in range(len(self[0])):
            vector = []
            for i in range(len(self)):
                vector.append([self[i][j]])
            columns.append(Matrix(vector))
        return columns

    def trace(self):
        """Calculates the tracer (diagonal sum) of a Matrix."""

        i = 0
        trace_sum = 0
        while i < len(self) and i < len(self[0]):
            trace_sum += self[i][i]
        return trace_sum

    def inner_prod(self, other=None):
        """Calculates the inner product (dot product) of two column vectors."""
        if other is None:
            other = self.copy()

        # converts column vectors to row vectors in form of python list, row vectors as python list
        obj1 = self if len(self) == 1 else self.transpose()
        obj2 = other if len(other) == 1 else other.transpose()

        # if nested vectors, un-nest
        if isinstance(obj1[0], list):
            obj1 = obj1[0]
        if isinstance(obj2[0], list):
            obj2 = obj2[0]

        if len(obj1) != len(obj2):
            raise ValueError("Dot product impossible with vectors of different lengths")

        # important to return python float or int instead of numpy, since numpy numbers will
        # override a matrix operator in operations where numpy object comes first
        return python_number(dot(obj1, obj2, self.mod))

    def norm(self):
        """Calculates the norm (distance from origin) of a vector"""

        obj = self if len(self) == 1 else self.transpose()

        if isinstance(obj[0], list):
            if len(obj) > 1:
                raise ValueError(f"operation incompatible with matrix of dimension {len(self)}x{len(self[0])}")
            obj = obj[0]

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
            if Matrix.inner_prod(pair[0], pair[1]) != 0:
                return False
        return True

    def coordinates(self, basis):
        """Returns the matrix representation of the coordinates of the given vector in the given basis."""

        # assumes instance matrix is matrix of column vectors
        matrix = self.transpose()

        if len(matrix[0]) == len(basis):
            basis = basis.transpose_obj()

        coord = []
        for i in range(len(matrix)):
            coord.append([])
            for j in range(len(basis)):
                coord[i].append(Matrix.inner_prod(matrix[i], basis[j]))
        return Matrix(array=coord).transpose()

    def orthonormalize(self, steps=False):
        """Computes the orthonormal basis of a given basis using the Gram-Schmidt process. Assumes input
        basis is matrix of column vectors."""

        matrix = self.transpose()

        length = len(matrix[0])

        norm = Matrix.norm(matrix[0])
        e1 = matrix[0] / norm
        basis = [e1]
        if steps:

            # print statement to provide spacing for readability
            print()
            e1_str = str(FractionMatrix(Matrix(matrix[0]).transpose())).split('\n')
            for j in range(length):
                if j == length // 2:
                    line = f'e1 = {fraction(1 / norm)} * '
                else:
                    line = ' ' * (8 + len(fraction(1 / norm)))
                line += e1_str[j]
                print(line)

            # print statement to provide spacing for readability
            print()
        for i, v in enumerate(matrix[1:]):

            # E is the sum of each already found orthonormal basis vector and its product with v
            E = sum(map(lambda e: Matrix.inner_prod(v, e) * e, basis))

            # e_i_0 is the new basis vector that is orthogonal to all other basis vectors, but not currently
            # orthonormal, since it still must be divided by its norm
            e_i_0 = v - E
            if steps:

                # print statement to provide spacing for readability
                print()

                # list of each line of string for vector v
                v_str = str(FractionMatrix(Matrix(v).transpose())).split('\n')

                # list of each dot product between v and each basis vector
                dot_str = list(map(lambda e: fraction(Matrix.inner_prod(v, e)), basis))

                # list of each line of string for each basis vector e
                basis_str = list(map(lambda e: str(FractionMatrix(Matrix(e).transpose())).split('\n'), basis))
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
                    sum_str = str(FractionMatrix(Matrix(E).transpose())).split('\n')
                    line += sum_str[j]
                    print(line)

                # print statement to provide spacing for readability
                print('\n')
                result = str(FractionMatrix(Matrix(e_i_0).transpose())).split('\n')
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

        ortho_basis = Matrix(basis).transpose()

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
                    (self - Matrix(rows=len(self), identity=x)).minor(), x)
                )
            )  # inner map converts each root of solve into float if not complex
        )  # outer map converts each element of list into int if possible

    def minor(self, index=None):
        """Computes the minor of each sub-matrix of instance matrix. Effectively computes determinant, except
        if index is given, in which case just the minor of that index is computed.
        Use matrix.det() unless solving for system of equations with sympy.Symbol

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
        if isinstance(index, int):
            return pow(-1, index) * self[0][index] * self.remove_row(
                0, copy=True).remove_column(index, copy=True).minor()
        sign = 1
        for j in range(len(self[0])):
            det += sign * self[0][j] * self.remove_row(0, copy=True).remove_column(j, copy=True).minor()
            sign *= -1
        return det

    def cross(self):
        """Given a list of n-1 vectors each of length n, returns 1 vector of length n that is
        orthogonal to all given vectors"""

        matrix = self.transpose() if len(self) == len(self[0]) + 1 else self.copy()  # if list is col vecs, transpose
        if len(self) == 1 and len(self[0]) == 2:
            return Array([self[0][1], -self[0][0]])
        return Array(map(lambda i: pow(-1, i) * matrix.remove_column(i, copy=True).minor(), range(len(matrix[0]))))

    def det(self):
        if len(self) != len(self[0]):
            raise AttributeError(f"operation unsupported for non-square matrices")

        return python_number(np_det(np_array(self.tolist())))  # using numpy's det function for efficiency


# this class is used just for printing matrices with fractions instead of floats, no mathematical
# operations can be performed with this class currently
class FractionMatrix:
    def __init__(self, array):

        # attempts to convert each number into a fraction before adding to matrix
        fraction_matrix = []
        for i in range(len(array)):
            fraction_matrix.append([])
            for j in range(len(array[0])):
                e = python_number(array[i][j])
                if isinstance(e, int) or (isinstance(e, float) and e.is_integer()):
                    fraction_matrix[i].append(str(int(e)))
                else:
                    fraction_matrix[i].append(fraction(e))
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
