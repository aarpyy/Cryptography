from typing import Sequence, Iterable
from functools import reduce


def bits(n):
    count = 0
    while n:
        n &= n - 1
        count += 1
    return count


# Binary-array -> array consisting of just (0, 1)
class barray:
    _masks = {}

    __slots__ = '_array', '_len',

    def __init__(self, arg, *, _len=None):
        if isinstance(arg, int) and isinstance(_len, int):
            self._array = arg
            self._len = _len
        elif isinstance(arg, Iterable):
            self._array = self.as_integer(arg)
            self._len = len(bin(self._array)) - 2
        else:
            raise TypeError(f"invalid arguments for {self.__class__.__name__}")

        # if its a non-empty array, add its length to _masks
        if self._len and self._len + 1 not in barray._masks:
            mask = 1
            # add up to len(self) + 1 so that a full mask can be retrieved by _masks[len(self)] - 1
            for i in range(1, self._len + 2):
                if i not in barray._masks:
                    barray._masks[i] = mask
                mask <<= 1

    @property
    def int(self):
        return self._array

    @staticmethod
    def mask():
        return barray._masks

    @staticmethod
    def as_integer(a):
        # converts any Integral-iterable into integer value
        value = 0
        for e in a:
            value = (value << 1) | (e & 1)
        return value

    def __str__(self):
        return "[" + ", ".join(str(e) for e in self) + "]"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.int}, _len={len(self)})"

    def __iter__(self):
        if len(self):
            mask = barray._masks[self._len]
            while mask:
                yield 1 if mask & self._array else 0
                mask >>= 1
        else:
            return

    def __add__(self, other):
        return barray(self.int ^ other.int, _len=len(self))

    def __mul__(self, other):
        return bits(self.int & other.int)

    def __eq__(self, other):
        return self.int == other.int and len(self) == len(other)

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item == slice(None, None, None):
                return +self
            elif item.step in (1, None):
                if item.start is None:
                    start = self._len
                else:
                    start = self._len - (item.start % self._len)

                if item.stop is None:
                    stop = 0
                else:
                    stop = self._len - (item.stop % self._len)

                sliced = (self._array & (barray._masks[start + 1] - 1))
                return barray(sliced >> stop, _len=start - stop)
            else:
                return barray(list(self)[item])
        else:
            return 1 if self._array & barray._masks[self._len - (item % self._len)] else 0

    def __len__(self):
        return self._len

    def __pos__(self):
        return barray(self.int, _len=len(self))

    def __lshift__(self, other):
        self._array <<= other
        self._len += 1
        if self._len + 1 not in barray._masks:
            barray._masks[self._len + 1] = barray._masks[self._len] << 1
        return self

    def __ilshift__(self, other):
        self._array <<= other
        if self._array:
            self._len += 1
            if self._len + 1 not in barray._masks:
                barray._masks[self._len + 1] = barray._masks[self._len] << 1
        return self

    def __or__(self, other):
        self._array |= other
        return self

    def __and__(self, other):
        return 1 if self._array & other else 0


class bmatrix:

    def __init__(self, array):
        self._array = array

    def __str__(self):
        if not len(self):
            return str([])
        padding = 3
        formatted = "["
        width = len(self[0])
        mask = 1 << (width - 1)
        for i in range(l := len(self)):
            if i == 0:
                formatted += "["
            else:
                formatted += " ["
            m = mask
            first_col = 1
            while m:
                e = '1' if self[i].int & m else '0'
                pad_left = (padding - len(e)) // 2
                pad_right = padding - len(e) - pad_left

                if first_col and e[0] != '-':  # sets numbers back from left [ to not squish, doesn't with negatives
                    formatted += max(pad_left, 1) * " " + f"{e}" + " " * pad_right
                else:
                    formatted += pad_left * " " + f"{e}" + " " * pad_right
                first_col = 0
                m >>= 1
            if i == l - 1:
                formatted += "]"
            else:
                formatted += "]\n"
        return formatted + "]"

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array})"

    def __iter__(self):
        return iter(self._array)

    def __setitem__(self, key, value):
        self._array[key] = value

    def __add__(self, other):
        return bmatrix(list(map(lambda a, b: a + b, self, other)))

    def __mul__(self, other):
        t = other.transpose()
        _len = len(other[0])
        # iterates through rows of first and columns of second, if dot product is odd then output is 0 else 1
        return bmatrix(list(map(lambda a: barray(reduce(lambda r, c: (r << 1) | (1 if bits(a.int & c.int) & 1 else 0),
                                                        t, 0), _len=_len), self)))

    def __rmul__(self, other):
        return barray(reduce(lambda r, c: (r << 1) | (1 if bits(other.int & c.int) & 1 else 0), self.transpose(), 0))

    def __eq__(self, other):
        return all(a.int == b.int for a, b in zip(self, other))

    def __pos__(self):
        return bmatrix([+e for e in self])

    def __getitem__(self, item):
        if isinstance(item, slice):
            return bmatrix(self._array[item])
        else:
            return self._array[item]

    def __len__(self):
        return len(self._array)

    def append(self, item):
        self._array.append(item)

    def separate(self, index=-1):
        """Separates matrix at index, returning two matrices, first with indices [0, index) and
        second with indices [index:)."""

        width = len(self[0])
        index = width - (index % width)
        size_l = width - index
        mask = barray.mask()[index + 1] - 1
        right, left = reduce(lambda r, c: ((*r[0], barray((a := self[c].int & mask), _len=index)),
                                           (*r[1], barray((self[c].int ^ a) >> index, _len=size_l))),
                             range(len(self)), ((), ()))

        return bmatrix(list(left)), bmatrix(list(right))

    def transpose(self):
        mask = barray.mask()[len(self[0])]
        array = bmatrix([])
        _len = len(self)
        while mask:
            val = 0
            for j in range(_len):
                if self._array[j].int & mask:
                    val = (val << 1) | 1
                else:
                    val <<= 1
            array.append(barray(val, _len=_len))
            mask >>= 1
        return array

    def row_reduce(self, row, col):
        for i in range(len(self)):
            if self[i].int & col and i != row:
                self[i] = self[i] + self[row]

    def rref(self):
        matrix = +self  # deep copy

        pivot_row = 0   # first pivot belongs in first row
        mask = barray.mask()[len(self[0])]
        while mask:

            # start at looking for pivot after previous pivot row
            for i in range(pivot_row, len(self)):

                # if non-zero element, this row can become pivot row
                if matrix[i].int & mask:

                    if i > pivot_row:                   # if pivot row not already in correct position, swap
                        matrix[i], matrix[pivot_row] = +matrix[pivot_row], +matrix[i]

                    matrix.row_reduce(pivot_row, mask)  # row reduce everything else
                    pivot_row += 1
            mask >>= 1

        return matrix

    def kernel(self):

        def add_rows(a, b):
            return list(x ^ y for x, y in zip(a, b))

        num_rows = len(self)        # get number of rows
        num_cols = len(self[0])     # get number of columns

        # this chunk here takes the transpose of current matrix and adds to each row one row of identity matrix
        array = []
        mask = barray.mask()[num_cols]
        index = 0
        row = [0] * num_cols
        while mask:
            extra = row[:]
            extra[index] = 1
            array.append([1 if r.int & mask else 0 for r in self] + extra)
            mask >>= 1
            index += 1

        pivot_row = 0                                   # first pivot belongs in first row
        for j in range(num_rows):
            for i in range(pivot_row, num_cols):        # start at looking for pivot after previous pivot row
                if array[i][j]:                         # if non-zero element, this row can become pivot row
                    if i > pivot_row:                   # if pivot row not already in correct position, swap
                        array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                    for k in range(num_cols):       # row reduce everything else
                        if array[k][j] and k != pivot_row:
                            array[k] = add_rows(array[k], array[pivot_row])

                    pivot_row += 1

        def separate(arr, j):
            left, right = [], []
            for i in range(len(arr)):
                left.append(arr[i][:j])
                right.append(arr[i][j:])
            return left, right

        mat, kern = separate(array, num_rows)   # separates original matrix from now modified identity matrix
        basis = []

        for i, row in enumerate(mat):
            if not any(row):                    # all null rows in original matrix correspond to basis vector for kernel
                basis.append(kern[i])

        return
