from .tools import isnumber, python_number
from numpy import ndarray, asarray
from math import gcd


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
        elif isnumber(array):
            self.array = [python_number(array)]
        else:
            raise TypeError("Array requires input object to be non-empty list")

    def __eq__(self, other):
        if isnumber(other):
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
        if isnumber(other):
            return Array(list(map(lambda e: 1 if e < other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 < e2 else 0, self, other)))

    def __le__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: 1 if e <= other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 <= e2 else 0, self, other)))

    def __gt__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: 1 if e > other else 0, self)))
        return Array(list(map(lambda e1, e2: 1 if e1 > e2 else 0, self, other)))

    def __ge__(self, other):
        if isnumber(other):
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
        if item == slice(None, None, None):
            return Array(self.array[:])
        return self.array[item]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __contains__(self, item):
        return self.array.__contains__(item)

    def __neg__(self):
        return Array(list(map(lambda e: -e, self)))

    def __add__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: e + other, self)))
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

        return Array(list(map(lambda x, y: x + y, self, other)))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: e - other, self)))
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

        return Array(list(map(lambda x, y: x - y, self, other)))

    def __rsub__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: other - e, self)))
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

        return Array(list(map(lambda x, y: x - y, other, self)))

    def __mul__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: e * other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x * y, self, other)))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __floordiv__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: e // other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x // y, self, other)))

    def __rfloordiv__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: e // other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x // y, other, self)))

    def __truediv__(self, other):
        if isnumber(other):
            return Array(list(map(lambda e: e / other, self)))
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return Array(list(map(lambda x, y: x / y, self, other)))

    def __rtruediv__(self, other):
        if isnumber(other):
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

    def append(self, item):
        self.array.append(item)

    def copy(self):
        return Array(self.array[:])

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

    def __eq__(self, other):
        if isinstance(other, ArrayMod) and other.mod != self.mod:
            return False
        if isinstance(other, str) and other == 'inv':
            return list(map(lambda e: 1 if gcd(e, self.mod) == 1 else 0, self))
        if isnumber(other) and other < 0:
            other %= self.mod
        return super().__eq__(other)

    def __ne__(self, other):
        if isinstance(other, ArrayMod) and other.mod != self.mod:
            return True
        if isinstance(other, str) and other == 'inv':
            return list(map(lambda e: 0 if gcd(e, self.mod) == 1 else 1, self))
        if isnumber(other) and other < 0:
            other %= self.mod
        return super().__ne__(other)

    def __lt__(self, other):
        if isnumber(other) and other < 0:
            other %= self.mod
        return super().__lt__(other)

    def __le__(self, other):
        if isnumber(other) and other < 0:
            other %= self.mod
        return super().__le__(other)

    def __gt__(self, other):
        if isnumber(other) and other < 0:
            other %= self.mod
        return super().__gt__(other)

    def __ge__(self, other):
        if isnumber(other) and other < 0:
            other %= self.mod
        return super().__ge__(other)

    def __repr__(self):
        return f"ArrayMod({self.array}, mod={self.mod})"

    def __neg__(self):
        return ArrayMod(list(map(lambda e: -e % self.mod, self)), self.mod)

    def __add__(self, other):
        if isnumber(other):
            return ArrayMod(list(map(lambda e: (e + other) % self.mod, self)), self.mod)
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")
        return ArrayMod(list(map(lambda x, y: (x + y) % self.mod, self, other)), self.mod)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isnumber(other):
            return ArrayMod(list(map(lambda e: (e - other) % self.mod, self)), self.mod)
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

        return ArrayMod(list(map(lambda x, y: (x - y) % self.mod, self, other)), self.mod)

    def __rsub__(self, other):
        if isnumber(other):
            return ArrayMod(list(map(lambda e: (other - e) % self.mod, self)), self.mod)
        if not isinstance(other, (list, Array)):
            raise TypeError(f"unsupported operation for type(s): {type(self)} and {type(other)}")
        return ArrayMod(list(map(lambda x, y: (y - x) % self.mod, self, other)), self.mod)

    def __mul__(self, other):
        if isnumber(other):
            return ArrayMod(list(map(lambda e: (e * other) % self.mod, self)), self.mod)
        if len(self) != len(other):
            raise AttributeError(f"unsupported operation for objects of length(s): {len(self)} and {len(other)}")
        return ArrayMod(list(map(lambda x, y: (x * y) % self.mod, self, other)), self.mod)

    def __floordiv__(self, other):
        """Attempts to multiply array by modular inverse of given number, raising an error if no inverse exists."""

        if isnumber(other):
            if gcd(other, self.mod) == 1:
                inverse = pow(other, -1, self.mod)
                return ArrayMod(list(map(lambda e: (e * inverse) % self.mod, self)), self.mod)
            raise ValueError(f"division cannot be performed because {other} is not invertible mod {self.mod}")
        raise ValueError(f"unsupported operation for type(s): {type(self)} and {type(other)}")

    def __rfloordiv__(self, other):
        raise ValueError(f"unsupported operation for object of type: {__class__.__name__}")

    def __truediv__(self, other):
        if isnumber(other):
            if gcd(other, self.mod) == 1:
                inverse = pow(other, -1, self.mod)
                return ArrayMod(list(map(lambda e: (e * inverse) % self.mod, self)), self.mod)
            elif (g := gcd(*(other, self.mod, *self))) > 1:
                div_array = list(map(lambda e: e // g, self))
                return ArrayMod(div_array, self.mod // g).__floordiv__(other // g)
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
        pass

    def copy(self):
        return ArrayMod(self.array[:], self.mod)

    def make_pivot(self, index=None, copy=False):
        if index is None:
            index = where(self)[0][0]

        pivot = self[index]
        if gcd(pivot, self.mod) == 1:
            if copy:
                return self.__mul__(inv := pow(pivot, -1, self.mod))
            self.array = list(map(lambda n: n * inv, self.array))
        raise ValueError("array is not invertible at the given index")

