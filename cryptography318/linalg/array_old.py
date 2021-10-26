from cryptography318.core.tools import python_number
from numpy import ndarray, asarray
from math import gcd
from functools import reduce
from numbers import Number


def where(array):
    """
    Evaluates boolean values of array.

    :param array: array_like to be evaluated
    :return: tuple containing list of indices of true values, two lists if array nested

    Example
    -------

    >>> a = Array([1, 2, 3, 0, 1, 0])
    >>> where(a)
    ([0, 1, 2, 4],)

    Notes
    -----

    Limited extension of numpy,where() that provides functionality of evaluating
    array values as boolean, and returning column and row indices of truth values.
    This where, evaluates using numpy.where() but returns output as type list. Does
    not extend any further functionality from numpy.where()."""

    # where usually returns tuple of np.array of row and column indices, this converts them into python integers
    # inside a python list, so that numpy library does not extend past this function
    return tuple(map(lambda arr: list(map(lambda n: int(n), arr)), asarray(array).nonzero()))


class Array:
    def __init__(self, array):
        if isinstance(array, ndarray):
            array = array.tolist()
        elif isinstance(array, map):
            array = list(array)

        # checks if array is list item and non-empty and not nested
        if isinstance(array, list) and not array:
            self.array = array
        elif isinstance(array, list) and not isinstance(array[0], list):
            self.array = array
        # if nested list is flat row vector, take it
        elif isinstance(array, list) and len(array) == 1:
            self.array = array[0]
        elif isinstance(array, Number) and not isinstance(array, complex):
            self.array = [python_number(array)]
        else:
            raise TypeError("Array requires input object to be non-empty list")

    def __eq__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: 1 if e == other else 0, self)))
        if isinstance(other, set):
            return Array(list(map(lambda e: 1 if e in other else 0, self)))
        return not any(map(lambda e1, e2: 0 if e1 == e2 else 1, self, other))

    def __ne__(self, other):
        result = self.__eq__(other)
        if isinstance(result, Array):
            return Array(list(map(lambda e: (e + 1) % 2, result)))
        return not result

    def __lt__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: 1 if e < other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 < e2 else 0, self, other)))

    def __le__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: 1 if e <= other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 <= e2 else 0, self, other)))

    def __gt__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: 1 if e > other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 > e2 else 0, self, other)))

    def __ge__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: 1 if e >= other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 >= e2 else 0, self, other)))

    def __len__(self):
        return len(self.array)

    def __str__(self):
        return str(self.array)

    def __repr__(self):
        return f"Array({self.array})"

    def __iter__(self):
        return iter(self.array)

    def __getitem__(self, item):
        # if full slice of array, return a deep copy of array instead of list slice
        if isinstance(item, slice):
            return Array(self.array[item])
        return self.array[item]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __contains__(self, item):
        return self.array.__contains__(item)

    def __neg__(self):
        return Array(list(map(lambda e: -e, self)))

    def __add__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e + other, self)))
        elif isinstance(other, (list, Array)):
            return Array(list(map(lambda x, y: x + y, self, other)))
        raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e - other, self)))
        elif isinstance(other, (list, Array)):
            return Array(list(map(lambda x, y: x - y, self, other)))
        raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

    def __rsub__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: other - e, self)))
        elif isinstance(other, (list, Array)):
            return Array(list(map(lambda x, y: x - y, other, self)))
        raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

    def __mul__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e * other, self)))
        elif isinstance(other, (Array, list)):
            if len(self) != len(other):
                raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
            return Array(list(map(lambda x, y: x * y, self, other)))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __floordiv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e // other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x // y, self, other)))

    def __rfloordiv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e // other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x // y, other, self)))

    def __truediv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e / other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x / y, self, other)))

    def __rtruediv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return Array(list(map(lambda e: e / other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x / y, other, self)))

    def __pow__(self, power, modulo=None):
        return Array(list(map(lambda e: pow(e, power, modulo), self)))

    def __mod__(self, other):
        return Array(list(map(lambda e: e % other, self)))

    def __abs__(self):
        return Array(list(map(abs, self)))

    def __bool__(self):
        return bool(self.array)

    def append(self, item):
        self.array.append(item)

    def copy(self):
        return Array(self.array[:])

    def mod(self, mod):
        return ArrayMod(self.array[:], mod)

    def index(self, item):
        return self.array.index(item)

    def to_ndarray(self):
        return ndarray(self.array)

    def make_pivot(self, index=None, copy=False):
        if index is None:
            index = where(self)[0][0]
        if copy:
            return self.__truediv__(self[index])
        self.array = list(map(lambda n: n / self[index], self.array))

    def shift_elements(self, shift, copy=False):
        """Performs a logical shift of elements in array.

        Examples
        --------
        >>> a = Array([1, 2, 3, 4])
        >>> a.shift_elements(shift=-1)
        [2, 3, 4, 1]

        >>> a.shift_elements(shift=2)
        [4, 1, 2, 3]

        >>> a = Array([1, 2, 3, 4])
        >>> b = a.shift_elements(shift=-1, copy=True)
        >>> b
        [2, 3, 4, 1]

        >>> b.shift_elements(3)
        [3, 4, 1, 2]

        :param shift: integer which determines how many indices to shift by
        :param copy: boolean dictates if instance is shifted or if new instance is created, shifted, then returned
        :return: shifted array if copy is True, otherwise nothing is returned
        """

        if copy:
            array = self.copy()
            array.shift_elements(shift, copy=False)
            return array
        shift = -shift % len(self)
        if shift > 0:
            self.array = self.array[shift:] + self.array[:shift]

    def contains_only(self, elements=None):
        if elements is None:
            elements = (0, 1)
        if getattr(elements, '__iter__', None) is not None:
            for e in self.array:
                if e not in elements:
                    return False
        else:
            for e in self.array:
                if e != elements:
                    return False
        return True


class ArrayMod(Array):
    def __init__(self, array, mod):
        super().__init__(array)
        self.mod = mod
        self.array = list(map(lambda e: e % mod, self.array))

    def __getitem__(self, item):
        # if full slice of array, return a deep copy of array instead of list slice
        if isinstance(item, slice):
            return ArrayMod(self.array[item], self.mod)
        return self.array[item]

    def __eq__(self, other):
        if isinstance(other, ArrayMod) and other.mod != self.mod:
            return False
        if isinstance(other, str) and other == 'inv':
            return list(map(lambda e: 1 if gcd(e, self.mod) == 1 else 0, self))
        if isinstance(other, Number) and not isinstance(other, complex) and other < 0:
            other %= self.mod
        return super().__eq__(other)

    def __ne__(self, other):
        if isinstance(other, ArrayMod) and other.mod != self.mod:
            return True
        if isinstance(other, str) and other == 'inv':
            return list(map(lambda e: 0 if gcd(e, self.mod) == 1 else 1, self))
        if isinstance(other, Number) and not isinstance(other, complex) and other < 0:
            other %= self.mod
        return super().__ne__(other)

    def __lt__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex) and other < 0:
            other %= self.mod
        return super().__lt__(other)

    def __le__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex) and other < 0:
            other %= self.mod
        return super().__le__(other)

    def __gt__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex) and other < 0:
            other %= self.mod
        return super().__gt__(other)

    def __ge__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex) and other < 0:
            other %= self.mod
        return super().__ge__(other)

    def __repr__(self):
        return f"ArrayMod({self.array}, mod={self.mod})"

    def __neg__(self):
        return ArrayMod(list(map(lambda e: -e % self.mod, self)), self.mod)

    def __add__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return ArrayMod(list(map(lambda e: (e + other) % self.mod, self)), self.mod)
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")
        return ArrayMod(list(map(lambda x, y: (x + y) % self.mod, self, other)), self.mod)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return ArrayMod(list(map(lambda e: (e - other) % self.mod, self)), self.mod)
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

        return ArrayMod(list(map(lambda x, y: (x - y) % self.mod, self, other)), self.mod)

    def __rsub__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return ArrayMod(list(map(lambda e: (other - e) % self.mod, self)), self.mod)
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")
        return ArrayMod(list(map(lambda x, y: (y - x) % self.mod, self, other)), self.mod)

    def __mul__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            return ArrayMod(list(map(lambda e: (e * other) % self.mod, self)), self.mod)
        elif len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return ArrayMod(list(map(lambda x, y: (x * y) % self.mod, self, other)), self.mod)

    def __floordiv__(self, other: int):
        """Attempts to multiply array by modular inverse of given number or divide by gcd, raising an
        error if no inverse exists or no greatest common divisor greater than 1 exists. If a divisor
        partially divides the array, and the remainder after partial division is invertible given the modulus,
        then both steps will be performed, first dividing the entire array by the gcd of all elements and divisor,
        then multiplying array by modular inverse of remaining divisor after initial division (see examples
        for more details on how this works/its purpose).

        Examples
        --------

        [1]

        >>> a = ArrayMod([2, 4, 6, 8, 10], 14)
        >>> repr(a // 6)
        ArrayMod([5, 10, 1, 6, 11], mod=14)

        [2]

        >>> a = ArrayMod([2, 4, 6, 8, 10], 14)
        >>> a // 4
        ValueError: division cannot be performed because 2 is not invertible mod 14 or it does not divide the array
        >>> repr(a // 2)
        ArrayMod([1, 2, 3, 4, 5], mod=14)

        Notes
        -----
        The full trace of dividing ArrayMod([2, 4, 6, 8, 10], 14) by 6 is as follows: first, 2 is found to be
        the gcd of the divisor and all elements of the array. Both the divisor (6) and the array are divided
        by this gcd, and then floor divison is attempted once more, with the goal of multiplying the array
        by the modular inverse of the remaining divisor. In this case, this is possible, since the resulting array
        after floor divison by 2 is ArrayMod([1, 2, 3, 4, 5], 14) and the remaining divisor (3) is invertible
        mod 14. The array is then mutliplied by the modular inverse, and final array is as shown in example 1.
        This is built like this for the specific reason of being able to convert an array value to 1, even if
        it is not directly invertible with the given modulus. In the final array, after dividing
        ArrayMod([2, 4, 6, 8, 10], 14) by 6, the 3rd element is 1, since division by 6 achieves the goal of
        converting all elements of value 6 into 1, either through direct multiplication by inverse, or through
        floor division followed by multiplication. This is integral in providing support for the function
        ArrayMod.make_pivot(), as it allows for a pivot to be constructed even if it requires more than one
        standard operator.
        """

        if isinstance(other, Number) and not isinstance(other, complex):
            if gcd(other, self.mod) == 1:
                inverse = pow(other, -1, self.mod)
                return ArrayMod(list(map(lambda e: (e * inverse) % self.mod, self)), self.mod)

            # gcd of other and elements of row is checked, if row can be divided it is done, and remainder of
            # divisor is then used to try and divide fully again, this allows for division by a number not
            # invertible given the modulus but
            if (d := gcd(other, *self.array)) > 1:
                return ArrayMod(list(map(lambda e: e // d, self)), self.mod) // (other // d)
            raise ValueError(f"division cannot be performed because {other} is not invertible mod {self.mod} or "
                             f"it does not divide the array")
        raise ValueError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

    def __rfloordiv__(self, other):
        raise ValueError(f"unsupported operation for object of type: {__class__.__name__}")

    def __truediv__(self, other):
        if isinstance(other, Number) and not isinstance(other, complex):
            if (m := gcd(*(*self.array, self.mod, other))) == 1:
                raise ValueError(f"{repr(self)} is not divisible by {other}")
            array = self.copy()
            array.mod = self.mod // m
            return array // other

        raise ValueError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

    def __rtruediv__(self, other):
        raise ValueError(f"unsupported operation for object of type: {__class__.__name__}")

    def __pow__(self, power, modulo=None):
        mod = self.mod if modulo is None else modulo
        return ArrayMod(list(map(lambda e: pow(e, power, mod), self)), self.mod)

    def __mod__(self, other=None):
        """Computes each element of the array mod an integer. If no input is given, will compute mod
        the mod attribute of the instance. If an input is given, the array is converted to a new
        instance of ArrayMod with mod attribute equal to input integer."""

        mod = self.mod if other is None else other
        return ArrayMod(list(map(lambda e: e % mod, self)), mod)

    def __abs__(self):
        return ArrayMod(list(map(lambda e: e % self.mod, self)), self.mod)  # abs redundant, in case of negatives: mod

    def copy(self):
        return ArrayMod(self.array[:], self.mod)

    def make_pivot(self, index=None, copy=False):
        """Attempts to divide instance array by integer given using floor division. If complete
        division fails, instead of throwing an error, False is returned. If division successful,
        result of division is returned."""

        if index is None:

            # iterates over list backwards, if current value is invertible return it, otherwise
            # return prev value, initial value as False so if value never gets changed, can be evaluated
            # as bool to determine if any elements in array are invertible
            g = gcd(*self.array)
            pivot = reduce(lambda i, c: c if (
                    gcd(c, self.mod) == 1 or (g > 1 and gcd(c//g, self.mod) == 1)
            ) else i, self.array[::-1], False)
            if not pivot:
                raise ValueError(f"no elements of {repr(self)} divide the array")
        else:
            pivot = self[index]

        try:
            if copy:
                return self // pivot
            array = self // pivot
            self.array = array.array
        except ValueError:
            return False
