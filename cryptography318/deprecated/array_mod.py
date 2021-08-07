from cryptography318.crypto_functions import *
from cryptography318.linear_algebra import Matrix, aslist, dot
import numpy


class array_mod:
    def __init__(self, array=None, cols=None, aug=False, mod=None):
        self.augmented = aug
        self.mod = mod
        if array is None and cols is None:
            raise ValueError("Constructor must be given either valid array or number of columns")
        # inherit modulus if possible
        if isinstance(array, array_mod):
            self.mod = array.mod
            self.array = aslist(array)
        elif cols is None:
            if isinstance(array, (list, Matrix, numpy.ndarray)):
                if isinstance(array[0], (list, numpy.ndarray)) and len(array) > 1:
                    raise ValueError("Input array must be a row vector")
                self.array = aslist(array)
            else:
                raise ValueError("Input array must be a row vector")
        else:
            self.array = [0] * cols
        if self.mod is None:
            raise ValueError(f"Constructor for {__class__.__name__} requires a modulus")

    def __getitem__(self, item):
        return self.array[item]

    def __setitem__(self, key, value):
        if isinstance(value, float):
            value = int(value)
        if not isinstance(value, int):
            raise ValueError(f"Type {type(value)} unacceptable in object of type {__class__.__name__}")
        self.array[key] = value

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return iter(aslist(self))

    def __str__(self):
        return str(self.array)

    def __add__(self, other):
        self.assert_array(other)
        if isinstance(other, float):
            other = int(other)
        result = aslist(self)
        if isinstance(other, (list, numpy.ndarray, array_mod, Matrix)):
            if self.is_multi_dimensional(other):
                self.dimension_error()
            other = aslist(other)
            for i in range(len(self)):
                result[i] += other[i]
                if result[i] >= self.mod:
                    result[i] %= self.mod
            return array_mod(result, mod=self.mod)
        if isinstance(other, int):
            for i in range(len(self)):
                result[i] += other
                if result[i] >= self.mod:
                    result[i] %= self.mod
            return array_mod(result, mod=self.mod)
        self.type_error()

    def __radd__(self, other):
        if isinstance(other, (int, float)):
            raise ValueError(f"Unable to add {__class__.__name__} object to number")
        result = self.__add__(other)
        if isinstance(other, array_mod):
            result.mod = other.mod
        return result

    def __sub__(self, other):
        self.assert_array(other)
        if isinstance(other, float):
            other = int(other)
        result = aslist(self)
        if isinstance(other, (list, numpy.ndarray, array_mod, Matrix)):
            if self.is_multi_dimensional(other):
                self.dimension_error()
            other = aslist(other)
            for i in range(len(self)):
                result[i] -= other[i]
                if result[i] >= self.mod or result[i] < 0:
                    result[i] %= self.mod
            return array_mod(result, mod=self.mod)
        if isinstance(other, int):
            for i in range(len(self)):
                result[i] -= other
                if result[i] >= self.mod or result[i] < 0:
                    result[i] %= self.mod
            return array_mod(result, mod=self.mod)
        self.type_error()

    def __rsub__(self, other):
        self.assert_array(other)
        result = aslist(self)
        if isinstance(other, (list, numpy.ndarray, array_mod, Matrix)):
            if self.is_multi_dimensional(other):
                self.dimension_error()
            other = aslist(other)
            for i in range(len(self)):
                result[i] = other[i] - self[i]
                if result[i] >= self.mod or result[i] < 0:
                    result[i] %= self.mod
            return array_mod(result, mod=self.mod)
        self.type_error()

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            warn("Proper syntax assumes scalar comes before matrix in scalar-matrix multiplication", SyntaxWarning)
        # array_mod multiplication is always dot product so order does not matter thus mul == rmul
        return self.__rmul__(other)

    def __rmul__(self, other):
        result = aslist(self)
        if isinstance(other, float):
            other = int(other)
        if isinstance(other, int):
            for i in range(len(result)):
                result[i] *= other
                if result[i] >= self.mod or result[i] < 0:
                    result[i] %= self.mod
            return array_mod(result, mod=self.mod)
        if isinstance(other, (list, array_mod, numpy.ndarray, Matrix)):
            if self.is_multi_dimensional(other):
                self.dimension_error()
            other = aslist(other)
            if isinstance(other[0], list):
                temp = []
                for elem in other:
                    temp.append(elem[0])
                other = temp
            if isinstance(other, array_mod):
                return array_mod(dot(aslist(self), other, mod=other.mod), mod=other.mod)
            return array_mod(dot(aslist(self), other, self.mod), mod=self.mod)
        self.type_error()

    def __truediv__(self, other):
        return self.__floordiv__(other)

    def __floordiv__(self, other):
        if isinstance(other, float):
            other = int(other)
        r = gcd(other, self.mod)
        if r > 1:
            result = aslist(self)
            row = result[:] + [self.mod, other]
            # find gcd that can divide everything
            e = gcd(*row)
            # if no gcd, row is not divisible by number at all
            if e == 1:
                raise ValueError(f"Dividing by {other} would result in some values being converted to floats")
            # divide by gcd
            for i in range(len(self)):
                result[i] //= e
            other //= e
            # if dividing by gcd accomplished division by other, return array
            if other == 1:
                return array_mod(result, mod=self.mod//e)
            # if dividing by gcd only partly divided by other, try again with remainder of other
            return array_mod(result, mod=self.mod//e).__floordiv__(other)
        # if number is invertible mod self.mod, multiply everything by inverse
        d = pow(other, -1, self.mod)
        return self.__rmul__(d)

    def __mod__(self, other):
        result = aslist(self)
        for i in range(len(self)):
            if result[i] > other:
                result[i] %= other
        return array_mod(result, mod=self.mod)

    def make_pivot(self, col=None):
        pivot = None
        if col is None:
            pivot = col
        else:
            for i in range(len(self)):
                if self[i] != 0:
                    pivot = self[i]
                    break
        if pivot is None:
            return self.copy()
        try:
            result = self.copy() / pivot
        except ValueError:
            pass
        else:
            # if division of result by first succeeded, result will have correct modulus value as well
            return result
        return self.copy()

    def reduce(self):
        for i in range(len(self)):
            if self[i] >= self.mod or self[i] < 0:
                self[i] %= self.mod

    def copy(self):
        return array_mod(aslist(self), mod=self.mod)

    @staticmethod
    def type_error(other=None):
        if other is not None:
            raise TypeError(f"Object of type {type(other)} is incompatible for the given operation")
        raise TypeError("Object type is incompatible for the given operation")

    @staticmethod
    def dimension_error():
        raise AttributeError(f"Performing an operation with {__class__.__name__} on a "
                             f"multi-dimensional array is unsupported")

    @staticmethod
    def assert_array(other, types=(int, float)):
        standard = (list, numpy.ndarray, Matrix, array_mod)
        types += standard
        if not isinstance(other, types):
            raise ValueError(f"Type incompatible with {__class__.__name__}: {type(other)}")
        if isinstance(other, standard) and isinstance(other[0], standard) and len(other) > 1:
            array_mod.dimension_error()

    @staticmethod
    def is_multi_dimensional(other):
        standard = (list, numpy.ndarray, Matrix, array_mod)
        if isinstance(other, standard) and isinstance(other[0], standard):
            if len(other[0]) == 1:
                return False
            return len(other) > 1
        return False
