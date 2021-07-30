import numpy
from random import randint
from warnings import warn


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
                        matrix[i].append(1)
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
        return self.array

    def __str__(self):
        """ Returns formatted string matrix """
        for i in range(len(self)):
            for j in range(len(self[0])):
                self[i][j] = numpy.round(self[i][j], decimals=3)

        max_len = 0
        for row in self.array:
            for e in row:
                if len(str(e)) > max_len:
                    max_len = len(str(e))
        padding = (max_len + 1) | 1
        formatted = "["
        for i in range(l := len(self.array)):
            if i == 0:
                formatted += "["
            else:
                formatted += " ["
            for j in range(len(self.array[0])):
                e = str(self.array[i][j])
                pad_left = (padding - len(e))//2
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
        raise TypeError(f"Cannot compare objects of type Matrix and type {type(other)}")

    def __mul__(self, other):
        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(self[0]) != len(other):
                raise ValueError("Number of columns in first matrix must equal number of rows in second")
        if isinstance(other, Matrix):
            return Matrix(numpy.matmul(self.array, other.array))
        if isinstance(other, numpy.ndarray):
            return Matrix(numpy.matmul(self.array, other))
        if isinstance(other, list):
            return Matrix(numpy.matmul(self.array, numpy.array(other)))
        if isinstance(other, int) or isinstance(other, float):
            warn("Proper syntax assumes scalar comes before matrix in scalar-matrix multiplication", SyntaxWarning)
            for i in range(len(self.array)):
                for j in range(len(self.array[0])):
                    self.array[i][j] *= other
            return self
        raise TypeError("Matrix multiplication must be done between two matrices or a matrix and a scalar")

    def __rmul__(self, other):
        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(other[0]) != len(self):
                raise ValueError("Number of columns in first matrix must equal number of rows in second")
        if isinstance(other, Matrix):
            return Matrix(numpy.matmul(other.array, self.array))
        if isinstance(other, numpy.ndarray):
            return Matrix(numpy.matmul(other, self.array))
        if isinstance(other, list):
            return Matrix(numpy.matmul(numpy.array(list), self.array))
        if isinstance(other, int) or isinstance(other, float):
            for i in range(len(self.array)):
                for j in range(len(self.array[0])):
                    self.array[i][j] *= other
            return self
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
            return Matrix(matrix_sum)
        raise TypeError("Matrix addition must be done between two matrices")

    def __radd__(self, other):
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
            return Matrix(matrix_sum)
        raise TypeError("Matrix addition must be done between two matrices")

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
            return Matrix(matrix_sum)
        raise TypeError("Matrix addition must be done between two matrices")

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
            return Matrix(matrix_sum)
        raise TypeError("Matrix addition must be done between two matrices")

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
        return product

    def __mod__(self, other):
        for i in range(len(self)):
            for j in range(len(self[0])):
                self[i][j] %= other
        return self

    def invert(self):
        """Inverts Matrix objecct using numpy.linalg.inv function."""

        return Matrix(numpy.linalg.inv(self.array))

    def transpose(self):
        """Returns transpose of Matrix object."""

        matrix = []
        for j in range(len(self[0])):
            matrix.append([])
            for i in range(len(self)):
                matrix[j].append(self[i][j])
        return Matrix(matrix)

    def make(self):
        """Function iterates through entire matrix, prompting user for input for each cell and
        subsequently changing value of cell to match input. Nothing is returned."""

        rows = len(self.array)
        cols = len(self.array[0])
        for i in range(rows):
            for j in range(cols):
                n = input(f"{i + 1},{j + 1}: ")
                if "," in n:
                    c = n.split(",")
                    real = float(c[0])
                    imag = float(c[1])
                    self.array[i][j] = complex(real, imag)
                else:
                    self.array[i][j] = float(n) if "." in n else int(n)

    def copy(self):
        """Returns exact copy of values of matrix in a new Matrix object."""

        return Matrix(self.array_copy())

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
        self.array.astype(data_type)

    def reset_type(self):
        """Attempts to reset matrix to smaller floats or integers if possible. If matrix has
        complex numbers they are left un-modified."""

        self.astype(numpy.float16)
        for i in range(len(self)):
            for j in range(len(self[0])):
                try:
                    e = self[i][j]
                    if "," in str(e):
                        continue
                    if e == -0:
                        self[i][j] = 0
                    s = str(e).split(".")
                    if len(s) == 1 or s[1] == '0':
                        self[i][j] = int(self[i][j])
                except OverflowError:
                    continue

    def augment(self, solution):
        """Given matrix object and set of solutions, returns an augmented coefficient matrix
        with the set of solutions as the final column."""

        self.augmented = True
        matrix_copy = self.array_copy()
        for i in range(len(self)):
            matrix_copy[i].append(solution[i][0])
        return Matrix(matrix_copy)

    def separate(self):
        """Function used to separate an augmented coefficient matrix into
        a standard coefficient matrix and a column matrix of solutions"""

        self.augmented = False
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
            for j in range(r := len(self[0])):
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
        return Matrix(clean_matrix)

    def rref(self):
        """Function that puts matrix in reduced row echelon form"""

        adjust = 1 if self.augmented else 0
        pivot_row = -1
        # iterates across columns
        for j in range((row_len := len(self[0])) - adjust):
            # iterates down columns
            pivot = False
            for i in range(col_len := len(self)):
                e = self[i][j]
                # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a
                # pivot
                if e != 0 and i > pivot_row:
                    # shifts entire row by factor that makes first entry = 1, therefore is pivot
                    # print(f"e = {e}, m[i] before shift: {m[i]}")
                    # print(shift)
                    for k in range(row_len):
                        self[i][k] /= e
                    # print(m[i])
                    # pivot row increases from where it last was, so that new row will be one below last pivot row
                    pivot_row += 1
                    # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot
                    # row
                    if pivot_row != i:
                        self[i][:], self[pivot_row][:] = self[pivot_row][:], self[i][:]
                    # this row is a pivot
                    pivot = True
                    break
            if pivot:
                # iterate down through matrix, removing all non-zero entries from column with pivot
                for k in range(col_len):
                    e = self[k][j]
                    if k != pivot_row and e != 0:
                        for l in range(row_len):
                            # here, e represents the number of pivot row's needed to be removed to make i'th row have
                            # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                            # pivot row will make i'th row have 0 in column
                            self[k][l] -= self[pivot_row][l] * e
        self.reset_type()
        return self

    def ref(self):
        """Function that puts matrix in row echelon form."""

        adjust = 1 if self.augmented else 0
        pivot_row = -1
        # iterates across columns
        for j in range((row_len := len(self[0])) - adjust):
            # iterates down columns
            pivot = False
            for i in range(col_len := len(self)):
                e = self[i][j]
                # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a
                # pivot
                if e != 0 and i > pivot_row:
                    # shifts entire row by factor that makes first entry = 1, therefore is pivot
                    # print(f"e = {e}, m[i] before shift: {m[i]}")
                    # print(shift)
                    for k in range(row_len):
                        self[i][k] /= e
                    # print(m[i])
                    # pivot row increases from where it last was, so that new row will be one below last pivot row
                    pivot_row += 1
                    # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot
                    # row
                    if pivot_row != i:
                        self[i][:], self[pivot_row][:] = self[pivot_row][:], self[i][:]
                    # this row is a pivot
                    pivot = True
                    break
            if pivot:
                # iterate down through matrix, removing all non-zero entries from column with pivot
                for k in range(col_len):
                    e = self[k][j]
                    if k > pivot_row and e != 0:
                        for l in range(row_len):
                            # here, e represents the number of pivot row's needed to be removed to make i'th row have
                            # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                            # pivot row will make i'th row have 0 in column
                            self[k][l] -= self[pivot_row][l] * e
        self.reset_type()
        return self

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

    def is_solvable(self):
        """is_solvable returns True if object is a square matrix with null rows removed,
        and False otherwise. If matrix is augmented, the column of solutions is not counted
        towards object being square."""

        # if not in rref, put it in that form
        if not self.is_rref():
            self.rref()
        # get copy with no null rows
        matrix = self.copy().remove_null()
        # if augmented, check for square excluding last column, otherwise check for square
        if self.augmented:
            if len(matrix) == len(matrix[0]) - 1:
                return True
            return False
        if len(matrix) == len(matrix[0]):
            return True
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
        # are more than one non-zero entry, the system has free variables and False is returend
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


class LinearMap(Matrix):
    def __init__(self, array):
        if isinstance(array, Matrix):
            super().__init__(array=array.array)
        else:
            super().__init__(array=array)

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

        inputs = [self, in_basis, out_basis]
        for e in inputs:
            if e is not None and not isinstance(e, Matrix):
                raise ValueError("All input must be type Matrix")

        if in_basis is None:
            if len(out_basis) != len(out_basis[0]):
                raise ValueError("Basis must be square matrix")
            return out_basis.invert() * self

        if out_basis is None:
            if len(in_basis) != len(in_basis[0]):
                raise ValueError("Basis must be square matrix")
            return self * in_basis

        if len(in_basis) != len(in_basis[0]) or len(out_basis) != len(out_basis[0]):
            raise ValueError("Bases must be square matrix")
        return out_basis.invert() * self * in_basis

    def map(self, other):
        """Applies linear map to matrix"""

        if isinstance(other, (Matrix, numpy.ndarray, list)):
            if len(self[0]) != len(other):
                raise ValueError(f"This linear function maps from {len(self[0])} dimensions to "
                                 f"{len(self)} dimensions, given: matrix in {len(other)} dimensions")
            return self * other
