import numpy
from random import randint, randrange
from warnings import warn
from math import gcd
from functools import reduce
from .tools import string_reduce, deprecated


def where(array, if_exp=None, else_exp=None):
    """Equivalent of numpy.where, accepts Matrix object as argument."""

    if isinstance(array, Matrix):
        if if_exp is None and else_exp is None:
            return numpy.where(array.array)
        return numpy.where(array.array, if_exp, else_exp)
    if isinstance(array, numpy.ndarray):
        if if_exp is None and else_exp is None:
            return numpy.where(array)
        return numpy.where(array, if_exp, else_exp)
    raise TypeError("where only accepts Matrix or numpy.ndarray objects")


def dot(obj, other, mod=None):
    """Equivalent of numpy.dot, accepts Matrix object as argument."""
    if len(obj) != len(other):
        raise ValueError(f"Unable to take product of two arrays of different length")
    if mod is not None:
        return [reduce(lambda a, b: (a + b) % mod, map(lambda x, y: (x * y) % mod, obj, other))]
    return [sum(map(lambda x, y: x * y, obj, other))]


def aslist(obj):
    """Returns object in form of Python list"""

    array = []
    if isinstance(obj[0], (list, numpy.ndarray)):
        # if normal multi-dimensional
        for i in range(len(obj)):
            array.append([])
            for j in range(len(obj[0])):
                array[i].append(obj[i][j])
        return array
    # if row vector
    for i in range(len(obj)):
        array.append(obj[i])
    return array


def isnumber(obj):
    types = (int, float, numpy.int16, numpy.int32, numpy.int64, numpy.float16,
             numpy.float32, numpy.float64)
    return isinstance(obj, types)


def is_binary_matrix(obj):
    """Returns False if non binary element in Matrix (not 0 or 1), True otherwise. Returns
    False if object is not Matrix."""

    if isinstance(obj, (Matrix, numpy.ndarray)):
        for row in obj:
            for e in row:
                if e not in [0, 1]:
                    return False
        return True
    return False


def all_elements(obj):
    """Returns all elements of list-type object in one Python list."""

    if isinstance(obj, (list, numpy.ndarray, Matrix)):
        elements = []
        for row in obj:
            if isinstance(row, list):
                for e in row:
                    elements.append(e)
            else:
                elements.append(row)
        return elements
    raise AttributeError("Cannot retrieve all elements from non-list-type object")


class Matrix:
    def __init__(self, array=None, rows=None, cols=None, rand=False, identity=False, aug=False, solution=None, mod=None):
        self.mod = mod
        self.augmented = False

        # if augmented is True or solution provided, matrix is augmented
        if aug or solution is not None:
            self.augmented = True

        # if array given, check if list or numpy array (cannot inherit a Matrix object)
        if array is not None:
            if isinstance(array, (list, numpy.ndarray)):

                # if array given is nested list, set self.array to given, otherwise, nest it
                if isinstance(array[0], (list, numpy.ndarray)):
                    self.array = numpy.array(array)
                else:
                    self.array = numpy.array([array])
            else:
                raise TypeError(f"Matrix must be type numpy.ndarray or list. Given: object of type {type(array)}")

            # if an array was given, attempt to reset type immediately
            self.reset_type()

        # if no array was given, matrix is not random or identity, and one of columns was not given matrix cannot be
        # constructed
        elif not rand and not identity and (rows is None or cols is None):
            raise ValueError("Constructor requires number of rows and columns, or valid matrix")

        # if random, no rows or columns required, creates m x n matrix with random m, n (where not given) and
        # values random integers between -50, 50
        elif rand:
            if rows is None:
                rows = randint(1, 10)
            if cols is None:
                cols = randint(2, 10)
            matrix = []
            for i in range(rows):
                matrix.append([])
                for j in range(cols):
                    matrix[i].append(randrange(-50, 50))
            self.array = numpy.array(matrix)

        # if identity, matrix must be square
        elif identity:
            if rows is not None and cols is not None and rows != cols:
                raise ValueError("Number of rows must equal number of columns in an identity matrix")

            # matrix constructed with just identity = True will return In for random n: [1, 10)
            if rows is None and cols is None:
                rows = randrange(1, 10)
                cols = rows

            # if one of cols or rows given, make the other equal so that matrix is square
            elif rows is None:
                rows = cols
            elif cols is None:
                cols = rows
            matrix = []
            for i in range(rows):
                matrix.append([])
                for j in range(cols):
                    if i == j:
                        matrix[i].append(1 if identity is True else identity)
                    else:
                        matrix[i].append(0)
            self.array = numpy.array(matrix, dtype=object)

        # if no input was given other than rows/cols, return null matrix of given size, rows/cols guaranteed
        # to not be null by this point, due to checking validity of input ~40 lines above
        else:
            matrix = []
            for i in range(rows):
                matrix.append([])
                for j in range(cols):
                    matrix[i].append(0)
            self.array = numpy.array(matrix, dtype=object)

        # if a solution is provided, attach column to end of matrix
        if solution is not None:
            if not isinstance(solution, (list, numpy.ndarray, Matrix)):
                raise TypeError(f"object {solution} is not array")
            for i in range(len(self)):
                self.array[i].append(solution[i][0])

    def __setitem__(self, key, value):
        self.array[key] = value

    def __getitem__(self, item):
        return self.array[item]

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return iter(aslist(self))

    def __repr__(self):
        return f"Matrix(array={self.array}, aug={self.augmented}, mod={self.mod})"

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
                elif isinstance(n, (bool, numpy.bool_)):
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
        if isinstance(other, (list, numpy.ndarray, Matrix)) and isinstance(other[0], (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
            for i in range(len(self)):
                for j in range(len(self[0])):
                    if abs(self[i][j] - other[i][j]) > 0.1:
                        return False
            return True

        # if not nested list, two objects will not be compared as arrays, but will return a binary matrix
        # with truth values wherever self's values are in the list
        if isinstance(other, list):
            binary_matrix = Matrix(rows=len(self), cols=len(self[0]))
            for i, row in enumerate(self):
                for j in range(len(row)):
                    e = row[j]
                    if e in other:
                        binary_matrix[i][j] = 1

            # pivot keyword allows user to search for pivots in matrix
            if 'pivot' in other:

                # returns the union of searching for pivot
                return Matrix.union(binary_matrix, self == 'pivot')
            return binary_matrix

        # if comparing to a number, return binary matrix with values = 1 wherever matrix values = number
        if isnumber(other):
            binary_matrix = []
            for i, row in enumerate(self):
                binary_matrix.append([])
                for e in row:
                    binary_matrix[i].append(1 if e == other else 0)
            return Matrix(binary_matrix)

        # if user searching for pivots, return binary matrix with values = 1 wherever an entry = 1
        # that is not in a row with a preceding value != 0 or 1 (first non-zero entry in row should be pivot)
        # also checks if entry = 1 is in same column as already identified pivot
        if other == 'pivot':
            binary_matrix = Matrix(rows=len(self), cols=len(self[0]))
            adj = 1 if self.augmented else 0

            # the reason pivot_column is kept track of but not used as a starting point for row iteration is
            # because if there are non-zero entries in front of a possible pivot, it becomes an invalid location
            pivot_column = -1
            for i in range(len(self)):
                for j in range(len(self[0]) - adj):
                    e = self[i][j]
                    if e not in [0, 1]:
                        break
                    if e == 1 and j > pivot_column:
                        binary_matrix[i][j] = 1
                        pivot_column = j
                        break
            return binary_matrix

        # if didn't identify valid type to compare with, error thrown
        raise TypeError(f"Cannot compare objects of type Matrix and type {type(other)}")

    def __ne__(self, other):

        # not equal in all currently considered cases, is the inverse of equal
        result = self.__eq__(other)

        # if returned is a matrix, this is assumed to be a binary matrix and it is inverted
        if isinstance(result, Matrix):
            return result.binary_inverse()

        # if matrix not returned, assumed to be bool value and is returned as not value
        return not result

    def __lt__(self, other):

        # uses numpy's built in array comparison methods to evaluate
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array < other)
        if isinstance(other, Matrix):
            return Matrix(self.array < other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __le__(self, other):

        # uses numpy's built in array comparison methods to evaluate
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array <= other)
        if isinstance(other, Matrix):
            return Matrix(self.array <= other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __gt__(self, other):

        # uses numpy's built in array comparison methods to evaluate
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array > other)
        if isinstance(other, Matrix):
            return Matrix(self.array > other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __ge__(self, other):

        # uses numpy's built in array comparison methods to evaluate
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array >= other)
        if isinstance(other, Matrix):
            return Matrix(self.array >= other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __mul__(self, other):

        # checks to see if dimensions of both objects are compatible with matrix multiplication
        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(self[0]) != len(other):
                raise ValueError("Number of columns in first matrix must equal number of rows in second")

        # uses numpy's matmul for all matrix multiplication
        if isinstance(other, Matrix):
            return Matrix(numpy.matmul(self.array, other.array))
        if isinstance(other, numpy.ndarray):
            return Matrix(numpy.matmul(self.array, other))
        if isinstance(other, list):
            return Matrix(numpy.matmul(self.array, numpy.array(other, dtype=object)))\

        # if number * matrix, return each element of matrix *= number
        if isnumber(other):
            matrix = self.copy()
            for i in range(len(self.array)):
                for j in range(len(self.array[0])):
                    matrix.array[i][j] *= other
            return matrix
        raise TypeError("Matrix multiplication must be done between two matrices or a matrix and a scalar")

    def __rmul__(self, other):

        # refer to __mul__ for documentation
        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(other[0]) != len(self):
                raise ValueError("Number of columns in first matrix must equal number of rows in second")
        if isinstance(other, Matrix):
            return Matrix(numpy.matmul(other.array, self.array))
        if isinstance(other, numpy.ndarray):
            return Matrix(numpy.matmul(other, self.array))
        if isinstance(other, list):
            return Matrix(numpy.matmul(numpy.array(other, dtype=object), self.array))
        if isnumber(other):
            matrix = self.copy()
            for i in range(len(self)):
                for j in range(len(self[0])):
                    matrix.array[i][j] *= other
            return matrix
        raise TypeError("Matrix multiplication must be done between two matrices or a matrix and a scalar")

    def __add__(self, other):

        # validates input as matrix that can be added if list-type
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
        if isinstance(other, (Matrix, numpy.ndarray)):
            matrix_sum = []
            for i in range(len(self)):
                matrix_sum.append([])
                for j in range(len(self[0])):
                    matrix_sum[i].append(self[i][j] + other[i][j])
            return Matrix(matrix_sum, aug=self.augmented)
        if isnumber(other):
            return Matrix(self.array + other, aug=self.augmented)
        raise TypeError("Matrix addition must be done between two matrices")

    def __radd__(self, other):
        if isnumber(other):
            raise TypeError("Cannot add matrix to number")
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
        if isinstance(other, (Matrix, numpy.ndarray)):
            matrix_sum = []
            for i in range(len(self)):
                matrix_sum.append([])
                for j in range(len(self[0])):
                    matrix_sum[i].append(self[i][j] - other[i][j])
            return Matrix(matrix_sum, aug=self.augmented)
        if isnumber(other):
            return Matrix(self.array - other, aug=self.augmented)
        raise TypeError("Matrix subtraction must be done between two matrices")

    def __rsub__(self, other):
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
        if isinstance(other, list):
            other = numpy.array(other, dtype=object)
        if isinstance(other, (Matrix, numpy.ndarray)):
            matrix_sum = []
            for i in range(len(self)):
                matrix_sum.append([])
                for j in range(len(self[0])):
                    matrix_sum[i].append(other[i][j] - self[i][j])
            return Matrix(matrix_sum, aug=self.augmented).reset_type()
        raise TypeError("Matrix subtraction must be done between two matrices")

    def __floordiv__(self, other):
        if isnumber(other):
            matrix = self.copy()
            matrix.array //= other
            return matrix
        raise TypeError("Matrix division only accepts number divisors")

    def __truediv__(self, other):
        if isnumber(other):
            matrix = self.astype(numpy.float64)
            matrix.array /= other
            return matrix
        raise TypeError("Matrix division only accepts number divisors")

    def __pow__(self, power, modulo=None):

        # powers currently supported are positive integers (referring to the number of times a matrix will
        # be multiplied against itself) and -1 (inverse)
        if not isinstance(power, int) or (power < 1 and power != -1):
            raise TypeError("Power must be real positive integer")
        if len(self) != len(self[0]):
            raise ValueError("Can only exponentiate square matrices")
        if power == -1:
            return self.invert()
        product = self.copy()
        for _ in range(power - 1):
            product = self * product
            if modulo is not None:
                product %= modulo
        product.reset_type()
        return product

    def __mod__(self, other):
        if not isnumber(other):
            raise ValueError(f"modulus must be done with number. given: {type(other)}")
        result = self.array.copy() % other
        return Matrix(result, aug=self.augmented)

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
            array.append([])
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
                            raise ValueError(f"You entered {len(values)} values, this matrix is expecting {cols}")
                    except ValueError:
                        continue
                    else:
                        array[i] = list(map(get_number, values))
                        break
        if isinstance(cls, Matrix):
            return cls(array=array, aug=cls.augmented)
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
                raise ValueError("Basis must be square matrix")
            return out_basis.invert()
        if out_basis is None:
            if len(in_basis) != len(in_basis[0]):
                raise ValueError("Basis must be square matrix")
            return in_basis
        if not all(length == len(in_basis) for length in [len(out_basis), len(out_basis[0]), len(in_basis[0])]):
            raise ValueError("Bases must be square matrix")
        return out_basis.invert() * in_basis

    # methods reserved for binary matrices
    def union(self, other):
        """Returns the logical union of two binary matrices.

        If this operation is needed with non-Matrix objects, use set.union(A, B)"""

        if not is_binary_matrix(self) and not is_binary_matrix(other):
            raise AttributeError("Union valid function only for binary matrices")

        if len(self) != len(other) or len(self[0]) != len(other[0]):
            raise AttributeError("Matrices are different dimensions. Union impossible")
        matrix = self.copy()
        for i in range(len(self)):
            for j in range(len(self[0])):
                matrix[i][j] = max((self[i][j], other[i][j]))
        return matrix

    def intersection(self, other):
        """Returns the logical intersection of two binary matrices.

        If this operation is needed with non-Matrix objects, use set.intersection(A, B)"""

        if isinstance(other, list) and not isinstance(other[0], list) and isinstance(self, list) and not isinstance(self[0], list):
            result = []
            for e in other:
                result.append(e)
            for e in self:
                result.append(e)
            return result

        if not is_binary_matrix(self) and not is_binary_matrix(other):
            raise AttributeError("Intersection valid function only for binary matrices")

        if len(self) != len(other) or len(self[0]) != len(other[0]):
            raise AttributeError("Matrices are different dimensions. Intersection impossible")
        matrix = self.copy()
        for i in range(len(self)):
            for j in range(len(self[0])):
                matrix[i][j] = 1 if self[i][j] == 1 and other[i][j] == 1 else 0
        return matrix

    def disjunction(self, other):
        """Returns the logical intersection of two binary matrices.

        If this operation is needed with non-Matrix objects, use set(A) - set(B)"""

        if not is_binary_matrix(self) and not is_binary_matrix(other):
            raise AttributeError("Disjunction valid function only for binary matrices")

        if len(self) != len(other) or len(self[0]) != len(other[0]):
            raise AttributeError("Matrices are different dimensions. Disjunction impossible")
        matrix = self.copy()
        for i in range(len(self)):
            for j in range(len(self[0])):
                matrix[i][j] = 1 if self[i][j] == 1 or other[i][j] == 1 else 0
        return matrix

    def binary_inverse(self):
        """Returns the negative image of a binary matrix"""
        obj = self.array

        if not is_binary_matrix(self):
            raise AttributeError("Binary inverse requires object to be binary matrix")

        return (self.copy() + 1) % 2

    # methods reserved for matrices with a specific modulus
    def make_pivot_mod(self, row, col=None):
        """
        Converts the first non-zero element from a row of matrix at given index into a pivot, attempting to divide
        entire row to reduce, unless division would result in float values in row, in which case a modular inverse
        is found. If column value is given, attempts to convert specific column value of given row into pivot.

        :params row: integer index of matrix corresponding to row to be operated on
        :params col: integer index of matrix corresponding to specific element in row to be converted into pivot
        """

        def mod_inv(array, e, mod):
            """Attempts to multiply entire row by modular inverse of e, if one exists. Returns True
            iff operation succeeded, else False."""
            if gcd(e, mod) == 1:
                array *= pow(int(e), -1, mod)
                array %= mod
                self[row] = array
                return True
            return False

        if self.mod is None:
            raise AttributeError("This function is reserved for matrices with a specific modulus")
        if row >= len(self):
            raise IndexError(f"Row index {row} out of bounds for Matrix with {len(self)} rows")

        self.reduce_mod()

        array = self[row].copy()
        if col is None:
            for i in range(len(array)):
                if array[i] != 0:
                    col = i
                    break

        e = array[col]
        if self.mod is None:
            raise AttributeError("This function is reserved for matrices with a specific modulus")
        res = mod_inv(array, e, self.mod)
        if res:
            return
        # if elements and mod are all 0 mod e then use integer division
        elements = all_elements(self.array) + [self.mod]
        r = gcd(*elements)
        if r > 1:
            self.mod //= r
            array //= r
            self[row] = array
            e //= r
        if e > 1:
            mod_inv(array, e, self.mod)

    def find_invertible(self):
        """Iterates through matrix over field of integers, finding values in each row that are invertible
        given their specific modulus, then returning a binary matrix, with truth values at each pivot location."""

        if self.mod is None:
            raise AttributeError("This function is reserved for matrices with a specific modulus")

        adj = 1 if self.augmented else 0

        self.reduce_mod()
        pivot_matrix = Matrix(rows=len(self), cols=len(self[0]))
        for i in range(len(self)):
            for j in range(len(self[0]) - adj):

                # if element is invertible, mark it as a possible pivot
                if gcd(self[i][j], self.mod) == 1:
                    pivot_matrix[i][j] = 1
        return pivot_matrix

    @deprecated
    def __choose_pivots_mod(self):
        """
        Function no longer has use.


        Calls find_invertible() to find all possible values that have modular inverse with row modulus,
        then decides, first through process of elimination (choosing only valid pivot in given row/column) then
        through recursion, which possible pivots would make a valid configuration of pivots in order to reduce original
        matrix into rref for solving.
        :return: set of coordinate pairs of pivots for matrix object (set of tuples)
        """

        pivot_matrix = self.find_invertible()
        adj = 1 if self.augmented else 0

        pivot_rows = []
        pivot_cols = []
        if len(set(numpy.where(pivot_matrix == 1)[1])) == len(self[0]) - adj:
            while True:

                found_pivot = False
                # finds columns with only one entry, lets that be pivot, removes column and row
                for i in range(len(pivot_matrix)):

                    # r is the list of columns with an entry == 1 from row i
                    r = set(numpy.where(pivot_matrix[i] == 1)[0]) - set(pivot_cols)

                    # if only one column to choose from, choose that column
                    if len(r) == 1:
                        p = r.pop()
                        if i in pivot_rows or p in pivot_cols:
                            continue
                        found_pivot = True
                        pivot_rows.append(i)
                        pivot_cols.append(p)

                    # print("pivot_matrix after a pass")
                    # print(pivot_matrix)
                    # print("pivots found")
                    # print(pivot_rows, "\n", pivot_cols)

                # gets transpose of binary matrix so that pivots per row can be checked easily
                pivot_t = pivot_matrix.transpose()
                for j in range(len(pivot_t)):

                    # c is the list of rows with an entry == 1 for the j'th column
                    # subtracting set of pivot rows ensures that if numpy.where found multiple columns, but
                    # some columns are invalidated by already found pivots, if statement can stil eval as True
                    c = set(numpy.where(pivot_t[j] == 1)[0]) - set(pivot_rows)

                    # if found column with one entry, use it
                    if len(c) == 1:
                        p = c.pop()
                        if j in pivot_cols or p in pivot_rows:
                            continue
                        found_pivot = True
                        pivot_rows.append(p)
                        pivot_cols.append(j)

                    # print("pivot_t after a pass")
                    # print(pivot_t)
                    # print("pivots found")
                    # print(pivot_rows, "\n", pivot_cols)

                # if found the right number of pivots, return set of tuples where each tuple is coordinate pair
                if len(pivot_cols) == len(self[0]) - adj:
                    return set(map(lambda r, c: (r, c), pivot_rows, pivot_cols))

                # if didn't find pivot in this iteration, no more pivots to find, try using recursion to choose
                if not found_pivot:
                    break

            # this is a recursive function that recurses down through rows of binary matrix, choosing pivot values,
            # it is guaranteed to return a result because the above if statement makes sure a solution is possible
            def choose_rec(piv_r, piv_c, row):
                start = 0
                while True:

                    # if found enough pivots, searching more would result in a false return since it would be searching
                    # rows that cannot have pivot, but don't need to since all are found
                    if len(piv_c) == len(pivot_matrix[0]) - adj:
                        return piv_r, piv_c

                    # if reached end of matrix, return pivot lists
                    if row == len(pivot_matrix):
                        return piv_r, piv_c

                    # make copies so actual values aren't overwritten if attempted values are wrong
                    pr, pc = piv_r[:], piv_c[:]

                    # bool for no new pivot has been found
                    found = False
                    for j in range(start, len(pivot_matrix[0]) - adj):

                        # if pivot in column already with pivot, won't work
                        if j in pc:
                            continue
                        if pivot_matrix[row][j] == 1:
                            pr.append(row)
                            pc.append(j)
                            found = True
                            break

                    # couldn't find pivot, that means this current config is unfeasible
                    if not found:
                        return False

                    result = choose_rec(pr, pc, row + 1)

                    # increment start, so that next time variable is chosen, new pivot is also chosen
                    if result is False:
                        start += 1
                        continue
                    else:
                        return result

            piv_r, piv_c = choose_rec(pivot_rows, pivot_cols, 0)
            return set(map(lambda r, c: (r, c), piv_r, piv_c))

        return False

    def reduce_mod(self):
        if self.mod is None:
            raise AttributeError("Given matrix cannot be reduced by a modulus because it does not have one")
        print(self.array, self.mod)
        self.array %= self.mod

    def rref_mod(self):
        """Uses Guassian elimination to row reduce matrix, returning a new instance of a matrix in reduced row
        echelon form. Currently nearly works, produces matrix in near-rref, only failing when an element is
        attempted to be made into a pivot but fails because it does not have a modular inverse with its modulus."""

        matrix = self.copy().astype(numpy.float64)
        matrix.reduce_mod()

        # if augmented don't reduce last column
        adj = 1 if self.augmented else 0

        pivot_row = 0
        for j in range(len(matrix[0]) - adj):
            for i in range(pivot_row, len(matrix)):

                # if non-zero element, this row can become pivot row
                if gcd(matrix[i][j], self.mod) == 1:

                    # make j'th element the pivot, reducing rest of row as well
                    matrix.make_pivot_mod(i, col=j)

                    # if pivot row not already in correct position, swap
                    if i > pivot_row:
                        matrix.swap(i, pivot_row)

                    # row reduce everything else
                    matrix.row_reduce(pivot_row, j)
                    pivot_row += 1

        matrix.reduce_mod()
        return matrix

    # standard matrix methods
    def make_pivot(self, row, col=None):
        """Converts first non-zero entry of a row into a pivot, by dividing the entire row by the entry.
        If column is provided, converts entry at index = column into a pivot. Does not return anything.

        :param row: index of row to be operated on
        :param col: column index for pivot location"""

        if row >= len(self):
            raise IndexError(f"Row index {row} out of bounds for Matrix with {len(self)} rows")
        array = self[row].copy()
        if col is None:
            for i in range(len(array)):
                if array[i] != 0:
                    col = i
                    break

        e = array[col]

        # if any elements are not 0 mod e, then use float division, otherwise use integer division
        divisible = array % e
        if divisible.any():
            array = array.astype(numpy.float64)
            self.array = self.array.astype(numpy.float64)
            array /= e
        else:
            array //= e
            array[col] = 1
        array[col] = 1
        self[row] = array
        self.reset_type()

    def append(self, row):
        def assert_len(obj):
            if len(obj) != len(self[0]):
                raise ValueError(f"Unable to append. Length should be {len(self[0])} not {len(obj)}")

        if not isinstance(row, (list, numpy.ndarray, Matrix)):
            raise TypeError("Can only append items of type list, numpy.ndarray, or Matrix to Matrix object")
        if isinstance(row, (numpy.ndarray, Matrix)):
            matrix = aslist(self)
            if isinstance(row[0], (list, numpy.ndarray)):
                for r in row:
                    assert_len(r)
                    matrix.append(r)
                self.array = numpy.array(matrix, dtype=object)
                return
            assert_len(row)
            matrix.append(list(row))
            self.array = numpy.array(matrix, dtype=object)
        if isinstance(row, list):
            new_matrix = aslist(self)
            assert_len(row)
            new_matrix.append(row)
            self.array = numpy.array(new_matrix, dtype=object)
            return

    def invert(self):
        """Inverts Matrix object using numpy.linalg.inv function."""

        array = self.array.astype(numpy.float64)
        return Matrix(numpy.linalg.inv(array), aug=self.augmented, mod=self.mod)

    def transpose(self):
        """Returns transpose of Matrix object."""

        transpose = self.array.transpose()
        return Matrix(transpose, aug=self.augmented)

    def copy(self):
        """Returns exact copy of values of matrix in a new Matrix object."""

        return Matrix(aslist(self), aug=self.augmented, mod=self.mod)

    def astype(self, data_type):
        """Returns matrix with a given data type, as specified by the numpy.ndarray dtypes. This function
        can be considered the Matrix object equivalent of numpy's astype."""

        if data_type not in (float, int, numpy.int64, numpy.float16, numpy.float32, numpy.float64):
            raise TypeError(f"Data type not recognized: {data_type}")

        result = self.array.astype(data_type)
        return Matrix(result, aug=self.augmented, mod=self.mod)

    def reset_type(self):
        """Attempts to reset matrix to smaller floats or integers if possible. If matrix has
        complex numbers they are left un-modified."""

        matrix = aslist(self.array.astype(numpy.float64))
        for i in range(len(self)):
            for j in range(len(self[0])):
                try:
                    e = matrix[i][j]
                    if "," in str(e):
                        continue
                    if e == -0:
                        matrix[i][j] = 0
                    s = str(e).split(".")
                    if len(s) == 1 or s[1] == '0':
                        matrix[i][j] = numpy.int64(matrix[i][j])
                except OverflowError:
                    continue
        self.array = numpy.array(matrix, dtype=object)

    def augment(self, solution):
        """Given matrix object and set of solutions, returns an augmented coefficient matrix
        with the set of solutions as the final column."""

        matrix_copy = aslist(self)
        if isinstance(solution[0], list):
            for i in range(len(self)):
                matrix_copy[i].append(solution[i][0])
        else:
            for i in range(len(self)):
                matrix_copy[i].append(solution[i])
        return Matrix(matrix_copy, aug=True)

    def separate(self):
        """Function used to separate an augmented coefficient matrix into
        a standard coefficient matrix and a column matrix of solutions"""

        matrix, sol = [], []
        c = len(self[0]) - 1
        for i in range(len(self)):
            matrix.append(list(self[i][:c]))
            sol.append([self[i][c]])
        return Matrix(matrix), Matrix(sol)

    def remove_null(self):
        """Removes all null rows and columns from matrix. Does not return."""

        self.remove_null_column()
        self.remove_null_row()

    def remove_null_row(self):
        """Removes rows consisting of just zeros. Does not return."""

        matrix = []
        for i in range(len(self)):

            # self[i].any() returns true if any element in the numpy array is non-zero, returning false if all
            # elements are zero means that this statement will run true if row has any non-zero element (is not null)
            if self[i].any():
                matrix.append(aslist(self[i]))
        self.array = numpy.array(matrix, dtype=object)

    def remove_null_column(self):
        """Removes columns consisting of just zeros. Does not return."""

        matrix = []
        array = self.array.transpose()
        for i in range(len(array)):

            # see above function remove_null_row for explanation of this statement
            if array[i].any():
                matrix.append(aslist(array[i]))
        self.array = numpy.array(matrix, dtype=object).transpose()

    def remove_row(self, row):
        new_matrix = []
        for i, r in enumerate(self):
            if i != row:
                new_matrix.append(r)
        self.array = numpy.array(new_matrix, dtype=object)

    def remove_column(self, col):
        new_matrix = []
        for j, c in enumerate(self.array.transpose()):
            if j != col:
                new_matrix.append(c)
        self.array = numpy.array(new_matrix, dtype=object)

    def swap(self, row1, row2):
        temp = self[row1].copy()
        self[row1] = self[row2].copy()
        self[row2] = temp

    def row_reduce(self, row, col):
        """Subtracts each row besides given by given row until j'th element is zero. This function is the basis
        of Matrix.rref()"""
        for i in range(len(self)):
            if i != row:
                self[i] -= self[row] * self[i][col]

    @deprecated
    def __rref_old(self):
        """Function that puts matrix in reduced row echelon form"""

        matrix = aslist(self.array.astype(numpy.float64))
        adjust = 1 if self.augmented else 0
        pivot_row = -1
        # iterates across columns
        for j in range((row_len := len(self[0])) - adjust):
            # iterates down columns
            pivot = False
            for i in range(col_len := len(self)):
                e = matrix[i][j]
                # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a
                # pivot
                if e != 0 and i > pivot_row:
                    # shifts entire row by factor that makes first entry = 1, therefore is pivot
                    for k in range(row_len):
                        matrix[i][k] /= e
                    # matrix[i] /= e
                    # pivot row increases from where it last was, so that new row will be one below last pivot row
                    pivot_row += 1
                    # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot
                    # row
                    if pivot_row != i:
                        temp = list(matrix[i][:])
                        matrix[i] = matrix[pivot_row][:]
                        matrix[pivot_row] = numpy.array(temp)
                        # matrix[i][:], matrix[pivot_row][:] = matrix[pivot_row][:], matrix[i][:]
                    # this row is a pivot
                    pivot = True
                    break
            if pivot:
                # iterate down through matrix, removing all non-zero entries from column with pivot
                for k in range(col_len):
                    e = matrix[k][j]
                    if k != pivot_row and e != 0:
                        for m in range(row_len):
                        # here, e represents the number of pivot row's needed to be removed to make i'th row have
                        # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                        # pivot row will make i'th row have 0 in column
                            matrix[k][m] -= matrix[pivot_row][m] * e
                        # matrix[k] -= (e * matrix[pivot_row])
        return Matrix(matrix, aug=self.augmented).reset_type()

    def rref(self):
        """Uses Guassian elimination to row reduce matrix, returning a new instance of a matrix in reduced row
        echelon form."""

        matrix = self.copy().astype(numpy.float64)

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
                    matrix.make_pivot(i, col=j)

                    # if pivot row not already in correct position, swap
                    if i > pivot_row:
                        matrix.swap(i, pivot_row)

                    # row reduce everything else
                    matrix.row_reduce(pivot_row, j)
                    pivot_row += 1

        return matrix

    def ref(self):
        """Function that puts matrix in row echelon form."""

        matrix = aslist(self.astype(numpy.float64))
        adjust = 1 if self.augmented else 0
        pivot_row = -1
        # iterates across columns
        for j in range((row_len := len(self[0])) - adjust):
            # iterates down columns
            pivot = False
            for i in range(col_len := len(self)):
                e = matrix[i][j]
                # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a
                # pivot
                if e != 0 and i > pivot_row:
                    # shifts entire row by factor that makes first entry = 1, therefore is pivot
                    # print(f"e = {e}, m[i] before shift: {m[i]}")
                    # print(shift)
                    # for k in range(row_len):
                    #     matrix[i][k] /= e
                    matrix[i] /= e
                    # print(m[i])
                    # pivot row increases from where it last was, so that new row will be one below last pivot row
                    pivot_row += 1
                    # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot
                    # row
                    if pivot_row != i:
                        matrix[i][:], matrix[pivot_row][:] = matrix[pivot_row][:], matrix[i][:]
                    # this row is a pivot
                    pivot = True
                    break
            if pivot:
                # iterate down through matrix, removing all non-zero entries from column with pivot
                for k in range(col_len):
                    e = matrix[k][j]
                    if k > pivot_row and e != 0:
                        # for m in range(row_len):
                        # here, e represents the number of pivot row's needed to be removed to make i'th row have
                        # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                        # pivot row will make i'th row have 0 in column
                            # matrix[k][m] -= matrix[pivot_row][m] * e
                        matrix[k] -= matrix[pivot_row] * e
        return Matrix(matrix, aug=self.augmented).reset_type()

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

    @deprecated
    def __is_solvable_old(self):
        """is_solvable returns True if object is (with null rows removed) a square matrix with
        no free variables, and False otherwise. If matrix is augmented, the column of solutions
        is not counted towards object being square."""

        # if not in rref, put it in that form
        if not self.is_rref():
            self.rref()
        # get copy with no null rows
        matrix = self.copy().remove_null_row()
        # any free variables in the system mean there is not a unique solution
        if matrix.free_vars() > 0:
            return False
        # if augmented, check for square excluding last column, otherwise check for square
        if self.augmented:
            if len(matrix) == len(matrix[0]) - 1:
                return matrix.is_consistent()
            return False
        if len(matrix) == len(matrix[0]):
            return matrix.is_consistent()
        return False

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

    @deprecated
    def __solve_old(self):
        """
        Solve is a function that attempts to solve a given system of linear equations.

        Function expects to be called on an augmented coefficient matrix, not necessarily
        in RREF. The function removes all null rows and attempts to solve using _solve_system.
        """
        if not self.augmented:
            raise AttributeError("This function is reserved for augmented coefficient matrices")

        if not self.is_solvable():
            return False

        # result could either be False, empty list, or list of solutions
        matrix = self.copy().remove_null_row()
        result = self.__solve_system(matrix)
        if not result:
            return result

        # sorts result by column number, adds in order from first to last column the solutions to a column
        # vector that is then returned as a Matrix
        solution = []
        for var in sorted(result):
            solution.append([result[var]])

        return Matrix(solution)

    @deprecated
    def __solve_system(self, matrix):
        """
        _solveSystem is a private helper function for public Solve that recursively
        determines solution for each variable in system of linear equations. Function
        calls itself recursively with input matrix being the given matrix minus the
        first row, until the input matrix is one row. Solving this row is then attempted,
        returning False if variable has no known solution, and a dictionary with the variable
        as the key and solution as the value otherwise. At each recursive step, the next
        variable and its solution are again added as the key and value in the dictionary.
        """

        r = len(matrix[0]) - 1

        # base case, for each non-zero entry in the row, it is added to the dictionary of solutions, if there
        # are more than one non-zero entry, the system has free variables and False is returned
        if len(matrix) == 1:
            _vars = {}
            for i in range(r):
                if matrix[0][i] != 0:
                    _vars[i] = matrix[0][i]
            if len(_vars) > 1:
                return False

            # if there is just one solution, solve for variable and return dictionary
            for x in _vars:
                _vars[x] = matrix[0][r] / _vars[x]
            return _vars

        result = self.__solve_system(matrix[1:])
        if result is False:
            return False

        # iterates through top row of matrix excluding solutions, if a value is not zero then it must be unknown or
        # in result, if it is not in result, it is truly unknown
        unknowns = []
        unknown_count = 0
        for i in range(r):
            if matrix[0][i] != 0:
                if i not in result:
                    unknown_count += 1
                unknowns.append(i)

        if unknown_count > 1:
            return False
        if unknown_count == 0:
            return result

        constant = 0
        index = None
        # iterates through all columns that have non-zero entries in this row, if they have a known
        # solution it is applied, otherwise it is considered the one true unknown and the index is marked
        for i in unknowns:
            if i in result:
                constant += result[i] * matrix[0][i]
            else:
                index = i

        # the one unknown in the row equals the row solution minus the sum of solutions of each
        # of the non-zero columns in the row
        result[index] = matrix[0][r] - constant
        return result

    def change_basis(self, basis):
        """Returns object previously in standard basis now in given basis"""

        if not isinstance(basis, Matrix):
            raise TypeError("Basis must be a Matrix")
        return basis.invert() * self

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


class LinearMap(Matrix):
    def __init__(self, array):

        # allows LinearMap to accept Matrix object as argument
        if isinstance(array, Matrix):
            super().__init__(array=array.array)
        else:
            super().__init__(array=array)

    # all binary matrix methods reserved for Matrix not LinearMap
    def union(self, other):
        pass

    def intersection(self, other):
        pass

    def disjunction(self, other):
        pass

    def binary_inverse(self):
        pass

    # all methods over field of integers reserved for Matrix not LinearMap
    def make_pivot_mod(self, row, col=None):
        pass

    def find_invertible(self):
        pass

    def reduce_mod(self):
        pass

    def rref_mod(self):
        pass

    # all methods involving augmented coefficient matrices for Matrix not LinearMap
    def is_solvable(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def is_consistent(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def solve(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def augment(self, solution):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def separate(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def is_basis(self):
        """is_basis allows a Matrix object to be a basis for its vector space. Linear maps map from one vector
        space to another (sometimes the same space but this doesn't matter here) and is not supposed to
        represent a basis of any one vector space. This method is thus ignored for linear maps."""

        pass

    # LinearMap methods not inherited or different than super class
    def dimension(self):
        """Dimension for a linear map is split into the dimension of the domain and
        dimension of the codomain, this functions serves no purpose."""

        pass

    def revert_basis(self, basis):
        """revert_basis function is incorporated in inherited change_basis function."""

        pass

    def domain_dimension(self):
        return len(self[0])

    def codomain_dimension(self):
        return len(self)

    def change_basis(self, in_basis=None, out_basis=None):
        """Function that takes as input a linear map and at least one basis and outputs
        the linear map with a change of basis.

        Linear map: an m x n matrix mapping from the standard basis vectors in m-dimensions to
        the standard basis vectors in n-dimensions.
        Notation: [S]E->E, where S, E are the linear map and standard basis vectors respectively.

        In-Basis: an m x m matrix that forms a basis in m-dimensions. This matrix, when applied
        to the linear map, gives the linear map such that the map takes as input matrices in
        this basis, and outputs them in the standard basis of n-dimensions.
        Notation: [S]B->E, where B is the In-Basis.

        Out-Basis: an n x n matrix that forms a basis in n-dimensions. This matrix, when applied
        to the linear map, gives the linear map such that the map takes as input matrices in
        the standard basis of m-dimensions and outputs a matrix in the form of the new basis.
        Notation: [S]E->C, where C is the Out-Basis"""

        if in_basis is None and out_basis is None:
            raise ValueError("Make sure to enter at least one basis")

        if not all(isinstance(e, Matrix) or e is None for e in [self, in_basis, out_basis]):
            raise ValueError("All input must be type Matrix")

        if in_basis is None:
            if len(out_basis) != len(out_basis[0]):
                raise ValueError("Basis must be square matrix")
            return LinearMap(out_basis.invert() * self)

        if out_basis is None:
            if len(in_basis) != len(in_basis[0]):
                raise ValueError("Basis must be square matrix")
            return LinearMap(self * in_basis)

        if len(in_basis) != len(in_basis[0]) or len(out_basis) != len(out_basis[0]):
            raise ValueError("Bases must be square matrix")
        # order of operations here is very important, first convert Te,e -> Tb,e converting T so that it takes input
        # in base b and output in e, then Tb,e -> Tb,b, converting T output to base b also
        return LinearMap(out_basis.invert() * (self * in_basis))

    def in_basis(self, basis):
        """This function converts a linear map into a new basis entirely, so that both input and output
        are in the given basis. This is a helper function that makes change_basis easier to work with."""

        if not all(length == len(self) for length in [len(self[0]), len(basis), len(basis[0])]):
            raise AttributeError("This function only works with square matrices for both basis and map")
        return self.change_basis(in_basis=basis, out_basis=basis)

    def map(self, other):
        """Applies linear map to matrix"""

        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(self[0]) != len(other):
                raise ValueError(f"This linear function maps from {len(self[0])} dimensions to "
                                 f"{len(self)} dimensions, given: matrix in {len(other)} dimensions")
            return self * other
        raise TypeError("Linear map can only apply to objects with type Matrix, numpy.ndarray, or list")

    def eigen_vector(self, vector):
        """Function takes as input some column vector and returns its eigenvalue if it is an
        eigenvector of the given linear map, or False if it is not an eigenvector."""

        from statistics import mean

        if len(self[0]) != len(vector):
            raise AttributeError(f"Eigenvector must have dimension {len(self[0])}. Given: vector with dimension "
                                 f"{len(vector)}")
        # maps vector with linear map
        eigen_vector = self.map(vector)
        eigen_value = []
        # iterates through column vector, checking to see if each element of the new eigen vector
        # is a scalar of the original vector, if one is zero the other must also be zero or it is not
        # an eigen vector
        for i in range(len(vector)):
            v1 = vector[i][0]
            ev1 = eigen_vector[i][0]
            if v1 == 0 or ev1 == 0:
                if v1 != ev1:
                    return False
                continue
            value = ev1 / v1
            if value not in eigen_value:
                eigen_value.append(value)

        # if all scalars were the exact same, return the eigenvalue
        if len(eigen_value) == 1:
            return eigen_value.pop()
        # if scalars were different enough, return False, but if all were within 0.1 of each other, this
        # is being considered an eigenvector, and the average of the eigenvalues is returned
        for i in range(len(eigen_value) - 1):
            e1, e2 = eigen_value[i], eigen_value[i + 1]
            if abs(e1 - e2) > 0.1:
                return False
        return numpy.round(mean(eigen_value), decimals=3)

    def is_eigen_value(self, value):
        """Function that takes in possible eigenvalue and returns True if value is eigenvalue, False otherwise.
        This only works with square maps."""

        if len(self) != len(self[0]):
            raise AttributeError("This function only works with square maps")

        # returns the determinant of M - xI where x is the potential eigenvalue
        matrix = self - Matrix(rows=len(self), cols=len(self), identity=value)
        array = matrix.array.astype(numpy.float64)
        return numpy.linalg.det(array) == 0
