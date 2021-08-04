import numpy
from random import randint
from warnings import warn
from .tools import string_reduce


def where(array, if_exp=None, else_exp=None):
    if isinstance(array, Matrix):
        if if_exp is None and else_exp is None:
            return numpy.where(array.array)
        return numpy.where(array.array, if_exp, else_exp)
    if isinstance(array, numpy.ndarray):
        if if_exp is None and else_exp is None:
            return numpy.where(array)
        return numpy.where(array, if_exp, else_exp)
    raise TypeError("where only accepts Matrix or numpy.ndarray objects")


class Matrix:
    def __init__(self, array=None, rows=None, cols=None, rand=False, identity=False, aug=False, solution=None):
        self.augmented = False
        if aug or solution is not None:
            self.augmented = True
        if array is not None:
            if isinstance(array, numpy.ndarray):
                self.array = array
            elif isinstance(array, list):
                self.array = numpy.array(array)
            else:
                raise TypeError(f"Matrix must be type numpy.ndarray or list. You gave type: {type(array)}")
        elif not rand and not identity and (rows is None or cols is None):
            raise ValueError("Constructor requires number of rows and columns, or valid matrix")
        elif rand:
            if rows is None:
                rows = randint(1, 10)
            if cols is None:
                cols = randint(1, 10)
            matrix = []
            for i in range(rows):
                matrix.append([])
                for j in range(cols):
                    matrix[i].append(randint(-50, 50))
            self.array = numpy.array(matrix)
        elif identity:
            if rows is not None and cols is not None and rows != cols:
                raise ValueError("Number of rows must equal number of columns in an identity matrix")
            if rows is None and cols is None:
                rows = randint(1, 10)
                cols = rows
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
            self.array = numpy.array(matrix)
        elif solution is not None:
            self.array = array
            for i in range(len(self)):
                self.array[i].append(solution[i][0])
        else:
            matrix = []
            for i in range(rows):
                matrix.append([])
                for j in range(cols):
                    matrix[i].append(0)
            self.array = numpy.array(matrix)

    def __setitem__(self, key, value):
        self.array[key] = value

    def __getitem__(self, item):
        return self.array[item]

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return iter(self.array)

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
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
            for i in range(len(self)):
                for j in range(len(self[0])):
                    if abs(self[i][j] - other[i][j]) > 0.5:
                        return False
            return True

        if isinstance(other, (int, float, numpy.int, numpy.float)):
            binary_matrix = []
            for i, row in enumerate(self):
                binary_matrix.append([])
                for e in row:
                    binary_matrix[i].append(True if e == other else False)
            return Matrix(binary_matrix)
        raise TypeError(f"Cannot compare objects of type Matrix and type {type(other)}")

    def __lt__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array < other)
        if isinstance(other, Matrix):
            return Matrix(self.array < other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __le__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array <= other)
        if isinstance(other, Matrix):
            return Matrix(self.array <= other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __gt__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array > other)
        if isinstance(other, Matrix):
            return Matrix(self.array > other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __ge__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float, numpy.ndarray)):
            return Matrix(self.array >= other)
        if isinstance(other, Matrix):
            return Matrix(self.array >= other.array)
        raise TypeError(f"Cannot compare between type Matrix and type {type(other)}")

    def __mul__(self, other):
        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(self[0]) != len(other):
                raise ValueError("Number of columns in first matrix must equal number of rows in second")
        if isinstance(other, Matrix):
            return Matrix(numpy.matmul(self.array, other.array)).reset_type()
        if isinstance(other, numpy.ndarray):
            return Matrix(numpy.matmul(self.array, other)).reset_type()
        if isinstance(other, list):
            return Matrix(numpy.matmul(self.array, numpy.array(other))).reset_type()
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            warn("Proper syntax assumes scalar comes before matrix in scalar-matrix multiplication", SyntaxWarning)
            matrix = self.copy()
            for i in range(len(self.array)):
                for j in range(len(self.array[0])):
                    matrix.array[i][j] *= other
            return matrix.reset_type()
        raise TypeError("Matrix multiplication must be done between two matrices or a matrix and a scalar")

    def __rmul__(self, other):
        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(other[0]) != len(self):
                raise ValueError("Number of columns in first matrix must equal number of rows in second")
        if isinstance(other, Matrix):
            return Matrix(numpy.matmul(other.array, self.array)).reset_type()
        if isinstance(other, numpy.ndarray):
            return Matrix(numpy.matmul(other, self.array)).reset_type()
        if isinstance(other, list):
            return Matrix(numpy.matmul(numpy.array(list), self.array)).reset_type()
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            matrix = self.copy()
            for i in range(len(self)):
                for j in range(len(self[0])):
                    matrix.array[i][j] *= other
            return matrix.reset_type()
        raise TypeError("Matrix multiplication must be done between two matrices or a matrix and a scalar")

    def __add__(self, other):
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
        if isinstance(other, list):
            other = numpy.array(other)
        if isinstance(other, (Matrix, numpy.ndarray)):
            matrix_sum = []
            for i in range(len(self)):
                matrix_sum.append([])
                for j in range(len(self[0])):
                    matrix_sum[i].append(self[i][j] + other[i][j])
            return Matrix(matrix_sum, aug=self.augmented).reset_type()
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            return Matrix(self.array + other, aug=self.augmented)
        raise TypeError("Matrix addition must be done between two matrices")

    def __radd__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            raise TypeError("Cannot add matrix to number")
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
        if isinstance(other, list):
            other = numpy.array(other)
        if isinstance(other, (Matrix, numpy.ndarray)):
            matrix_sum = []
            for i in range(len(self)):
                matrix_sum.append([])
                for j in range(len(self[0])):
                    matrix_sum[i].append(self[i][j] - other[i][j])
            return Matrix(matrix_sum, aug=self.augmented).reset_type()
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            return Matrix(self.array - other, aug=self.augmented)
        raise TypeError("Matrix subtraction must be done between two matrices")

    def __rsub__(self, other):
        if isinstance(other, (list, numpy.ndarray, Matrix)):
            if len(other) != len(self) or len(other[0]) != len(self[0]):
                raise ValueError("Both matrices must have the same number of rows and columns!")
        if isinstance(other, list):
            other = numpy.array(other)
        if isinstance(other, (Matrix, numpy.ndarray)):
            matrix_sum = []
            for i in range(len(self)):
                matrix_sum.append([])
                for j in range(len(self[0])):
                    matrix_sum[i].append(other[i][j] - self[i][j])
            return Matrix(matrix_sum, aug=self.augmented).reset_type()
        raise TypeError("Matrix subtraction must be done between two matrices")

    def __floordiv__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            self.array //= other
            return self
        raise TypeError("Matrix division only accepts number divisors")

    def __truediv__(self, other):
        if isinstance(other, (int, float, numpy.int, numpy.float)):
            matrix = self.astype(numpy.float64)
            matrix.array /= other
            return matrix
        raise TypeError("Matrix division only accepts number divisors")

    def __pow__(self, power, modulo=None):
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
        return product.reset_type()

    def __mod__(self, other):
        for i in range(len(self)):
            for j in range(len(self[0])):
                self[i][j] %= other
        return self.reset_type()

    @classmethod
    def make(cls, rows, cols, aug=False, by_row=True):
        """Class method that takes number of rows and columns as input and returns a new instance of
        a Matrix object with the values given by user.

        Parameters
        ----------
        aug : bool
            determines if matrix is an augmented coefficient matrix or not.
        by_row : bool
            determines if input from user is collected each row at a time, or each element at a time."""

        array = []
        if not by_row:
            for i in range(rows):
                array.append([])
                for j in range(cols):
                    n = input(f"{i + 1},{j + 1}: ")
                    if "," in n:
                        c = n.split(",")
                        real = float(c[0])
                        imag = float(c[1])
                        array[i].append(complex(real, imag))
                    else:
                        array[i].append(float(n) if "." in n else int(n))
        else:
            for i in range(rows):
                array.append([])
                while True:
                    try:
                        row = input(f"row {i + 1}: ")
                        values = row.split()
                        if len(values) != cols:
                            raise ValueError(f"You entered {len(values)} values, this matrix is expecting {cols}")
                    except ValueError:
                        continue
                    else:
                        break
                for v in values:
                    if "," in v:
                        c = v.split(",")
                        real = float(c[0])
                        imag = float(c[1])
                        array[i].append(complex(real, imag))
                    else:
                        array[i].append(float(v) if "." in v else int(v))
        return cls(array=array, aug=aug)

    @classmethod
    def diag(cls, scalars):
        matrix = cls(rows=len(scalars), cols=len(scalars), identity=True)
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

    def append(self, row):
        def assert_len(obj):
            if len(obj) != len(self[0]):
                raise ValueError(f"Unable to append. Length should be {len(self[0])} not {len(obj)}")

        if not isinstance(row, (list, numpy.ndarray, Matrix)):
            raise TypeError("Can only append items of type list, numpy.ndarray, or Matrix to Matrix object")
        if isinstance(row, (numpy.ndarray, Matrix)):
            matrix = self.array_copy()
            if isinstance(row[0], (list, numpy.ndarray)):
                for r in row:
                    assert_len(r)
                    matrix.append(r)
                self.array = numpy.array(matrix)
                return
            assert_len(row)
            matrix.append(list(row))
            self.array = numpy.array(matrix)
        if isinstance(row, list):
            new_matrix = self.array_copy()
            assert_len(row)
            new_matrix.append(row)
            self.array = numpy.array(new_matrix)
            return

    def invert(self):
        """Inverts Matrix object using numpy.linalg.inv function."""

        return Matrix(numpy.linalg.inv(self.array))

    def transpose(self):
        """Returns transpose of Matrix object."""

        matrix = []
        for j in range(len(self[0])):
            matrix.append([])
            for i in range(len(self)):
                matrix[j].append(self[i][j])
        return Matrix(matrix)

    def copy(self):
        """Returns exact copy of values of matrix in a new Matrix object."""

        return Matrix(self.array_copy(), aug=self.augmented)

    def array_copy(self):
        """Returns exact copy of values of matrix in a new list object."""

        matrix_copy = []
        for i in range(len(self)):
            matrix_copy.append([])
            for j in range(len(self[0])):
                matrix_copy[i].append(self[i][j])
        return matrix_copy

    def astype(self, data_type):
        """Returns matrix with a given data type, as specified by the numpy.ndarray dtypes:
        float16, float32, float64, or the traditional Python float or int."""

        if data_type not in [float, int, numpy.float16, numpy.float32, numpy.float64]:
            raise TypeError(f"Data type not recognized: {data_type}")

        matrix = self.array_copy()
        for i in range(len(self)):
            for j in range(len(self[0])):
                matrix[i][j] = data_type(self[i][j])
        return Matrix(matrix, aug=self.augmented)

    def reset_type(self):
        """Attempts to reset matrix to smaller floats or integers if possible. If matrix has
        complex numbers they are left un-modified."""

        matrix = self.astype(numpy.float32).array_copy()
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

        return Matrix(matrix, aug=self.augmented)

    def augment(self, solution):
        """Given matrix object and set of solutions, returns an augmented coefficient matrix
        with the set of solutions as the final column."""

        matrix_copy = self.array_copy()
        for i in range(len(self)):
            matrix_copy[i].append(solution[i][0])
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
        """Function that removes rows consisting of just zeros"""

        null_rows = []
        for i in range(c := len(self)):
            null = False
            for j in range(len(self[0])):
                if self[i][j] != 0:
                    null = False
                    break
                null = True
            if null:
                null_rows.append(i)

        clean_matrix = []
        for i in range(c):
            if i not in null_rows:
                clean_matrix.append(self[i][:])
        return Matrix(clean_matrix, aug=self.augmented)

    def remove_row(self, row):
        new_matrix = []
        for i, r in enumerate(self):
            new_matrix.append([])
            for j in range(len(r)):
                if j != row:
                    new_matrix[i].append(r[j])
        self.array = numpy.array(new_matrix)

    def rref(self):
        """Function that puts matrix in reduced row echelon form"""

        matrix = self.astype(numpy.float64).array_copy()
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

    def ref(self):
        """Function that puts matrix in row echelon form."""

        matrix = self.astype(numpy.float64).array_copy()
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
                # if there is a pivot in this row, this row is consistent, move to next
                if e == 1 and j < r - 1:
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
        """is_solvable returns True if object is (with null rows removed) a square matrix with
        no free variables, and False otherwise. If matrix is augmented, the column of solutions
        is not counted towards object being square."""

        # if not in rref, put it in that form
        if not self.is_rref():
            self.rref()
        # get copy with no null rows
        matrix = self.copy().remove_null()
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

    def solve(self):
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
        matrix = self.copy().remove_null()
        result = self._solve_system(matrix)
        if not result:
            return result

        # sorts result by column number, adds in order from first to last column the solutions to a column
        # vector that is then returned as a Matrix
        solution = []
        for var in sorted(result):
            solution.append([result[var]])

        return Matrix(solution)

    def _solve_system(self, matrix):
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

        result = self._solve_system(matrix[1:])
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
        if isinstance(array, Matrix):
            super().__init__(array=array.array)
        else:
            super().__init__(array=array)

    def is_solvable(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def is_consistent(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def solve(self):
        """Method only applies to augmented coefficient matrices. Does not apply to linear maps."""

        pass

    def _solve_system(self, matrix):
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

        matrix = self - Matrix(rows=len(self), cols=len(self), identity=value)
        return numpy.linalg.det(matrix.array) == 0
