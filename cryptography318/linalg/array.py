from .arrayabc import ABCArray, DimensionError, get_dtype
from typing import Sequence, Iterable
from numbers import Integral, Real
from math import gcd, inf
from sympy.core.expr import Expr
from functools import reduce
import operator


# Binary-array -> array consisting of just (0, 1)
class bitarray(ABCArray):

    _masks = {}

    __slots__ = '_len',

    def __init__(self, arg, *, _len=None):
        super().__init__()
        if isinstance(arg, int) and isinstance(_len, int):
            self._array = arg
            self._len = _len
        elif isinstance(arg, Sequence):
            self._array = self.as_integer(arg)
            self._len = len(arg)
        else:
            raise TypeError(f"invalid arguments for binary array")

        if self._len and self._len not in bitarray._masks:
            mask = 1
            for i in range(1, self._len):
                if i in bitarray._masks:
                    continue
                bitarray._masks[i] = mask
                mask <<= 1

        # dtype always integers since either 0 or 1
        self._dtype = int

    @property
    def int(self):
        return self._array

    def invert(self):
        return bitarray(self.int ^ ((bitarray._masks[len(self)] << 1) - 1), _len=len(self))

    def copy(self):
        return +self

    @property
    def bits(self):
        return sum(self)

    @property
    def array(self):
        return list(self)

    @staticmethod
    def as_integer(a):
        # converts any Integral-iterable into integer value
        return int(''.join(str(e & 1) for e in a), 2) if a else 0  # empty list returns as 0

    def __str__(self):
        return "[" + ", ".join(str(e) for e in self) + "]"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.int}, _len={len(self)})"

    def __iter__(self):
        if len(self):
            mask = 1 << (len(self) - 1)
            while mask:
                yield 1 if mask & self._array else 0
                mask >>= 1
        else:
            return

    def __setitem__(self, key, value):
        if key < -len(self) or key > len(self) - 1:
            raise IndexError(f"{self.__class__.__name__} index out of range")

        value &= 1
        # first shift mask so that its 1 is at beginning of list, then shift by key to get to index
        mask = bitarray._masks[len(self) - (key % len(self))]

        # get the bit at the index
        _key_bit = mask & self._array
        # if the bit at index: key == value, then no need to change anything
        if bool(_key_bit) == bool(value):
            pass
        else:
            # _mask is 1 at index: key, 0's everywhere else, so XOR makes all 1's in self stay 1's
            # if index: key is 1, then XOR will perform 1 ^ 1 which sets to 0 which is our value (we know since
            # if it was 1 first if statement would have run)
            # if index: key is 0 then XOR will perform 1 ^ 0 which sets to 1 which is our value
            self._array ^= mask

    def __add__(self, other):
        if isinstance(other, Integral):
            if other & 1:
                return self.invert()
            else:
                return +self
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                if isinstance(other, bitarray):
                    return bitarray(self._array ^ other.int, _len=len(self))
                else:
                    return bitarray(self._array ^ self.as_integer(other), _len=len(self))
            else:
                raise ValueError(f"unsupported operation for arrays of length(s): {len(self)} and {len(other)}")
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Integral):
            if other & 1:
                return self.invert()
            else:
                return +self
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                return bitarray(self._array ^ self.as_integer(other), _len=len(self))
            else:
                raise ValueError(f"unsupported operation for arrays of length(s): {len(self)} and {len(other)}")
        else:
            return NotImplemented

    # -1 == 1 and -0 == 0 mod 2 so subtraction is the same as addition
    __sub__, __rsub__ = __add__, __radd__

    def __mul__(self, other):
        if isinstance(other, bitarray):
            if len(self) == len(other):
                return bitarray(self.int & other.int, _len=len(self))
            else:
                raise ValueError(f"unsupported operation for arrays of length(s): {len(self)} and {len(other)}")
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                # if a is 1 and b & 1 is 1, then product is 1, otherwise its 0
                return bitarray(list(map(lambda a, b: 1 if a and b & 1 else 0, self, other)))
            else:
                raise ValueError(f"unsupported operation for arrays of length(s): {len(self)} and {len(other)}")
        elif isinstance(other, Integral):
            # any value % 2 == 0 will result in full 0 array,
            if not other & 1:
                return bitarray(0, _len=len(self))
            # otherwise it will just return same array
            else:
                return +self
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Sequence):
            if len(other) == len(self):
                # if a is 1 and b & 1 is 1, then product is 1, otherwise its 0
                return bitarray(list(map(lambda a, b: 1 if a & 1 and b else 0, other, self)))
            else:
                raise ValueError(f"unsupported operation for arrays of length(s): {len(self)} and {len(other)}")
        elif isinstance(other, Integral):
            # any value % 2 == 0 will result in full 0 array,
            if not other & 1:
                return bitarray(0, _len=len(self))
            # otherwise it will just return same array
            else:
                return +self
        else:
            return NotImplemented

    def dot(self, other):

        def bits(n):
            count = 0
            while n:
                n &= n - 1
                count += 1
            return count

        if isinstance(other, bitarray):
            return bits(self.int & other.int)
        elif isinstance(other, Sequence):
            return bits(self.int & self.as_integer(other))
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Sequence):
            if len(self) == len(other):
                if isinstance(other, bitarray):
                    # if any are not the same value, XOR will return positive integer value, so not that is False
                    return not bool(self._array ^ other.int)
                else:
                    return all(e1 == e2 for e1, e2 in zip(self, other))
            else:
                return False
        elif isinstance(other, int):
            if other == 1:
                return +self
            elif other == 0:
                return self.invert()
            else:
                return bitarray(0, _len=len(self))
        else:
            return NotImplemented

    def _compare(a, b, op):
        if isinstance(b, Integral):
            b &= 1
            return bitarray([op(e, b) for e in a])
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bitarray([op(x, y) for x, y in zip(a, b)])
            else:
                raise ValueError(f"unsupported operation for arrays of length(s): {len(a)} and {len(b)}")
        else:
            return NotImplemented

    def __pos__(self):
        return bitarray(self._array, _len=len(self))

    def __neg__(self):
        # -1 == 1 and -0 == 0 so no changes from __pos__()
        return bitarray(self._array, _len=len(self))

    def append(self, value):
        self._array = (self._array << 1) | (value & 1)
        self._len += 1
        if self._len not in bitarray._masks:
            bitarray._masks[self._len] = bitarray._masks[self._len - 1] << 1

    def reverse(self):
        val = self._array
        self._array = 0
        for _ in range(len(self)):
            self._array <<= 1
            self._array |= val & 1
            val >>= 1

    def extend(self, values):
        if isinstance(values, Sequence):
            _len = len(values)
            if isinstance(values, bitarray):
                self._array = (self._array << _len) | values.int
            else:
                self._array = (self._array << _len) | self.as_integer(values)
            if (k := self._len + _len) not in bitarray._masks:
                mask = bitarray._masks[len(self)]
                for i in range(self._len + 1, k):
                    mask <<= 1
                    bitarray._masks[i] = mask
            self._len = k
        else:
            return NotImplemented

    def pop(self, index=None):
        if not len(self):
            raise IndexError(f"pop from empty {self.__class__.__name__}")
        elif index is None:
            mask = bitarray._masks[len(self)]
            k = self._array & mask
            self._array &= mask - 1
        else:
            k = None
            mask = bitarray._masks[len(self)]
            val = 0
            for i in range(len(self)):
                v = self._array & mask
                mask >>= 1
                if index == i:
                    k = v
                else:
                    val <<= 1
                    val |= 1 if v else 0
            if k is None:
                raise IndexError(f"{self.__class__.__name__} index out of range")
            self._array = val
        self._len -= 1
        # if no index given, k will be some large value if the popped value is 1 (ex. [1, 0, 0].pop(), k = 100 = 4)
        return 1 if k else 0

    def remove(self, value):
        found = 0
        val = 0
        if not len(self):
            raise ValueError(f"unable to remove from empty {self.__class__.__name__}")
        elif value == 1:
            mask = 1 << (len(self) - 1)
            while mask:
                v = self._array & mask
                mask >>= 1
                if not found and v:
                    found = 1
                else:
                    val <<= 1
                    val |= 1 if v else 0
        elif not value:
            mask = 1 << (len(self) - 1)
            while mask:
                v = self._array & mask
                mask >>= 1
                if not found and not v:
                    found = 1
                else:
                    val <<= 1
                    val |= 1 if v else 0
        if not found:
            raise ValueError(f"{value} is not in {self.__class__.__name__}")
        self._array = val
        self._len -= 1

    def insert(self, index, value):
        if not len(self):
            self._len = 1
            self._array = value & 1
        elif index == 0:
            # only need to do anything if value inserted is 1
            if value & 1:
                self._array |= 1 << len(self)
            self._len += 1
        elif 0 < index < len(self):
            mask = (1 << (len(self) - index - 1)) - 1
            lower = self._array & mask
            upper = self._array ^ lower

            if value & 1:
                self._array = (upper << 1) | (1 << (len(self) - index)) | lower
            else:
                self._array = (upper << 1) | lower

            self._len += 1
        elif index == len(self):
            self.append(value)
        else:
            raise IndexError(f"{self.__class__.__name__} index out of range")

    def __getitem__(self, item):
        if isinstance(item, slice):
            # use builtin list slicing since a lot easier to deal w/ and understand
            return bitarray(list(self)[item])
        elif isinstance(item, int):
            if -len(self) <= item < len(self):
                return 1 if self._array & (1 << (len(self) - item - 1)) else 0
            else:
                raise IndexError(f"{self.__class__.__name__} index out of range")
        else:
            raise TypeError(f"{self.__class__.__name__} indices must be integers or slices, not {type(item)}")

    def __delitem__(self, item):
        self.remove(item)

    def __len__(self):
        return self._len


# Binary-array -> array consisting of just (0, 1)
class binarray(ABCArray):

    def __init__(self, array):
        super().__init__()
        if isinstance(array, list):
            self._array = array
        elif isinstance(array, binarray):
            self._array = array.array
        elif isinstance(array, Integral):
            self._array = [int(array)]
        elif isinstance(array, Iterable):
            self._array = list(array)
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from list or array"
                            f" instance")

        self._dtype = int
        self._array = list(e & 1 for e in self._array)

    def __str__(self):
        return str(self._array)

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array})"

    @property
    def complement(self):
        return binarray([e ^ 1 for e in self._array])

    def copy(self):
        return binarray(self._array[:])

    def __iter__(self):
        return iter(self._array)

    def __setitem__(self, key, value):
        if isinstance(value, int):
            self._array[key] = value & 1
        else:
            raise ValueError(f"all {self.__class__.__name__} items must be integers not {type(value)}")

    def __add__(self, other):
        if isinstance(other, Integral):
            other = int(other) & 1
            if other:
                return self.complement
            else:
                return self.copy()
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                if isinstance(other, binarray):
                    return binarray([x ^ y for x, y in zip(self, other)])
                else:
                    return binarray([(x + y) & 1 for x, y in zip(self, other)])
            else:
                raise DimensionError(self, other, op=operator.add)
        else:
            return NotImplemented

    # sub/add are same with binary array
    __radd__, __sub__, __rsub__ = __add__, __add__, __add__

    def __mul__(self, other):
        if isinstance(other, Integral):
            other = int(other) & 1
            if other:
                return self.copy()
            else:
                return binarray([0] * len(self))
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return binarray([x * y for x, y in zip(self, other)])
            else:
                raise DimensionError(self, other, op=operator.mul)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def dot(self, other):
        return sum(y for x, y in zip(self, other) if x)

    def __eq__(self, other):
        if isinstance(other, Integral):
            return binarray([1 if e == other else 0 for e in self])
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return binarray([1 if x == y else 0 for x, y in zip(self, other)])
            else:
                return False
        else:
            return NotImplemented

    def _compare(a, b, op):
        if isinstance(b, Integral):
            return binarray([1 if op(e, b) else 0 for e in a])
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return binarray([1 if op(x, y) else 0 for x, y in zip(a, b)])
            else:
                raise DimensionError(a, b, op=op)
        else:
            return NotImplemented

    def __pos__(self):
        return binarray([+e for e in self])

    def __neg__(self):
        return binarray([-e for e in self])

    def append(self, value):
        if isinstance(value, Integral):
            self._array.append(int(value) & 1)
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.append.__name__}() is not implemented for "
                                      f"non-Integral values")

    def reverse(self):
        self._array.reverse()

    def extend(self, values):
        if all(isinstance(v, Integral) for v in values):
            self._array.extend(int(e) & 1 for e in values)
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.extend.__name__}() is not implemented for "
                                      f"non-Integral values")

    def pop(self, index=None):
        return self._array.pop()

    def remove(self, value):
        if isinstance(value, Integral):
            self._array.remove(int(value))
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.remove.__name__}() is not implemented for "
                                      f"non-Integral values")

    def insert(self, index, value):
        if isinstance(value, Integral):
            self._array.insert(index, int(value) & 1)
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.insert.__name__}() is not implemented for "
                                      f"non-Integral values")

    def __getitem__(self, item):
        if isinstance(item, slice):
            # if slice, return smaller marray, otherwise just return the value
            return binarray(self._array[item])
        else:
            return self._array[item]

    def __delitem__(self, item):
        del self._array[item]

    def __len__(self):
        return len(self._array)


# Real-array -> array consisting of Real instances
class rarray(ABCArray):

    def __init__(self, array):
        super().__init__()
        if isinstance(array, Iterable):
            self._array = list(array)
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from Iterable")

        _types = set()
        for v in self._array:
            if isinstance(v, Expr):
                continue
            elif isinstance(v, float) and v.is_integer():
                _types.add(int)
            else:
                _types.add(type(v))

        # array could be non-empty but _types empty if consisting only of sympy Expr's
        self._dtype = get_dtype(*_types) if self._array and _types else int
        self.astype(self._dtype, _update=True)

    def astype(self, dtype, *, _update=False):
        if _update:
            self._array = [e if isinstance(e, Expr) else dtype(e) for e in self]
            self._dtype = dtype
        else:
            return rarray([e if isinstance(e, Expr) else dtype(e) for e in self])

    def copy(self):
        return rarray(self._array[:])

    def __str__(self):
        return str(self._array)

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array})"

    def __iter__(self):
        return iter(self._array)

    def __setitem__(self, key, value):
        self._array[key] = value
        if not isinstance(value, (self._dtype, Expr)):
            self._dtype = get_dtype(self._dtype, type(value))
            self.astype(self._dtype, _update=True)

    def __add__(self, other):
        if isinstance(other, Real):
            other = float(other)
            return rarray(list(map(lambda a: a + other, self._array)))
        elif isinstance(other, Expr):
            return rarray(list(map(lambda a: a + other, self._array)))
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                return rarray(list(map(lambda a, b: a + b, self._array, other)))
            else:
                raise DimensionError(self, other, op=operator.add)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Real):
            other = float(other)
            return rarray(list(map(lambda a: a - other, self._array)))
        elif isinstance(other, Expr):
            return rarray(list(map(lambda a: a - other, self._array)))
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                return rarray(list(map(lambda a, b: a - b, self._array, other)))
            else:
                raise DimensionError(self, other, op=operator.sub)
        else:
            return NotImplemented

    def __rsub__(self, other):
        return other + -self

    def __mul__(self, other):
        if isinstance(other, Real):
            other = float(other)
            return rarray(list(map(lambda a: a * other, self._array)))
        elif isinstance(other, Expr):
            return rarray(list(map(lambda a: a * other, self._array)))
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                return rarray(list(map(lambda a, b: a * b, self._array, other)))
            else:
                raise DimensionError(self, other, op=operator.mul)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def dot(self, other):
        if len(self) == len(other):
            return sum(x * y for x, y in zip(self, other))
        else:
            raise DimensionError(self, other, op=self.dot)

    def __truediv__(self, other):
        if isinstance(other, Real):
            other = float(other)
            if other.is_integer():
                other = int(other)
            return rarray([e / other for e in self])
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return rarray([x / y for x, y in zip(self, other)])
            else:
                raise DimensionError(self, other, op=operator.truediv)
        else:
            return NotImplemented

    def __floordiv__(self, other):
        if isinstance(other, Real):
            return rarray(list(e // float(other) for e in self))
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return rarray([x // float(y) for x, y in zip(self, other)])
            else:
                raise DimensionError(self, other, op=operator.floordiv)
        else:
            return NotImplemented

    def __mod__(self, other):
        if isinstance(other, Real):
            other = float(other)
            return rarray([e % other for e in self])
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return rarray([x % float(y) for x, y in zip(self, other)])
            else:
                raise DimensionError(self, other, op=operator.mod)
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Sequence):
            if len(self) == len(other):
                return all(a == b for a, b in zip(self._array, other))
            else:
                return False
        elif isinstance(other, Real):
            return bitarray([1 if e == other else 0 for e in self._array])
        else:
            return NotImplemented

    def _compare(a, b, op):
        if isinstance(b, Real):
            return bitarray([1 if op(e, b) else 0 for e in a])
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bitarray([1 if op(x, y) else 0 for x, y in zip(a, b)])
            else:
                raise DimensionError(a, b, op=op)
        else:
            return NotImplemented

    def __abs__(self):
        return rarray(list(map(lambda a: abs(a), self._array)))

    def __pos__(self):
        return rarray(list(e for e in self._array))

    def __neg__(self):
        return rarray(list(-e for e in self._array))

    def insert(self, index, value):
        if isinstance(value, (Real, Expr)):
            self._array.insert(index, value)
            if not isinstance(value, (self._dtype, Expr)):
                self._dtype = get_dtype(self._dtype, type(value))
                self.astype(self._dtype, _update=True)
        else:
            raise NotImplementedError

    def __getitem__(self, item):
        return self._array[item]

    def __delitem__(self, item):
        del self._array[item]

    def __len__(self):
        return len(self._array)

    def append(self, value):
        if isinstance(value, Real):
            value = float(value)
            if value.is_integer():
                value = int(value)
            if not self._array:
                self._array.append(value)
                self._dtype = type(value)
            else:
                self._array.append(self._dtype(value))
        elif isinstance(value, Expr):
            self._array.append(value)
        else:
            raise NotImplementedError

    def reverse(self):
        self._array.reverse()

    def extend(self, values):
        if isinstance(values, Sequence):
            _types = set()
            for v in values:
                if isinstance(v, (Expr, self._dtype)):
                    continue
                elif isinstance(v, float) and v.is_integer():
                    _types.add(int)
                else:
                    _types.add(type(v))

            # _types only has non-expr/non-dtype types, so if they are all expr/same as dtype then this doesn't run
            if _types:
                self._dtype = get_dtype(self._dtype, *_types)
                self.astype(self._dtype, _update=True)
            self._array.extend(values)
        else:
            return NotImplemented

    def pop(self, index=None):
        return self._array.pop(index)

    def remove(self, value):
        del self[value]

    def __iadd__(self, other):
        if isinstance(other, Real):
            other = self._dtype(other)
            for i in range(len(self)):
                self._array[i] += other
            return self
        elif isinstance(other, Expr):
            for i in range(len(self)):
                self._array[i] += other
            return self
        elif isinstance(other, Sequence):
            return self + other
        else:
            return NotImplemented

    def make_pivot(self, index=None):
        if index is None:
            index = reduce(lambda r, c: c if self[c] and c < r else r, range(len(self)), inf)

        value = self[index]
        self._array = [e / value for e in self._array]


# Modulus-array -> array with a specific modulus consisting of integers
class marray(ABCArray):

    __slots__ = '_mod',

    def __init__(self, array, modulus):
        super().__init__()
        if isinstance(array, list):
            self._array = array
        elif isinstance(array, ABCArray):
            self._array = array.array
        elif isinstance(array, Integral):
            self._array = [int(array)]
        elif isinstance(array, Iterable):
            self._array = list(array)
        else:
            raise TypeError(f"{self.__class__.__name__} must be constructed from list or array"
                            f" instance")

        self._dtype = int
        if not all(isinstance(e, (Integral, Expr)) for e in self) or not isinstance(modulus, int):
            raise TypeError(f"{self.__class__.__name__} must be constructed from list of integers and integer modulus")
        else:
            self._array = [e % modulus if isinstance(e, Expr) else int(e) % modulus for e in self]
        self._mod = modulus

    @property
    def mod(self):
        return self._mod

    def modulo(self, m):
        self._mod = m

    def __str__(self):
        return str(self._array)

    def __repr__(self):
        return f"{self.__class__.__name__}({self._array}, {self._mod})"

    def copy(self):
        return rarray(self._array[:])

    def __iter__(self):
        return iter(self._array)

    def __setitem__(self, key, value):
        if isinstance(value, Integral):
            self._array[key] = int(value) % self.mod
        elif isinstance(value, Expr):
            self._array[key] = value % self.mod
        else:
            raise ValueError(f"all {self.__class__.__name__} items must be integers not {type(value)}")

    def __add__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            return marray([e + other for e in self], self._mod)
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return marray([x + y for x, y in zip(self, other)], self._mod)
            else:
                raise DimensionError(self, other, op=operator.add)
        else:
            return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            return marray([e - other for e in self], self._mod)
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return marray([x - y for x, y in zip(self, other)], self._mod)
            else:
                raise DimensionError(self, other, op=operator.sub)
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            return marray([other - e for e in self], self._mod)
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return marray([x - y for x, y in zip(other, self)], self._mod)
            else:
                raise DimensionError(self, other, op=operator.sub)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            return marray([e * other for e in self], self._mod)
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return marray([x * y for x, y in zip(self, other)], self._mod)
            else:
                raise DimensionError(self, other, op=operator.mul)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def dot(self, other):
        if len(self) == len(other):
            return sum((x * y) % self._mod for x, y in zip(self, other)) % self._mod
        else:
            raise DimensionError(self, other, self.dot)

    def __truediv__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            g = gcd(other, self._mod)
            if g == 1:
                return self * pow(other, -1, self._mod)
            elif (k := gcd(*self._array, g)) > 1:
                return marray([e // k for e in self], self._mod // k) / (other // k)
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                if all(gcd(e, self.mod) == 1 for e in other):
                    return marray([x * pow(y, -1, self._mod) for x, y in zip(self, other)], self._mod)
                else:
                    raise ValueError(f"{self} is not divisible by {other}")
            else:
                raise DimensionError(self, other, op=operator.truediv)

        else:
            return NotImplemented

    def __floordiv__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            if gcd(other, self._mod) == 1:
                return self * pow(other, -1, self._mod)
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                if all(gcd(e, self.mod) == 1 for e in other):
                    return marray([x * pow(y, -1, self._mod) for x, y in zip(self, other)], self._mod)
                else:
                    raise ValueError(f"{self} is not divisible by {other}")
            else:
                raise DimensionError(self, other, op=operator.floordiv)
        else:
            return NotImplemented

    def __mod__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            return marray([e % other for e in self], self.mod)
        elif isinstance(other, Sequence):
            if not all(isinstance(e, Integral) for e in other):
                raise ValueError(f"all modulus must be Integrals")
            elif len(self) == len(other):
                return marray([x % int(y) for x, y in zip(self, other)], self.mod)
            else:
                raise DimensionError(self, other, op=operator.mod)
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Integral):
            return bitarray([1 if e == other else 0 for e in self])
        elif isinstance(other, Sequence):
            if len(self) == len(other):
                return all(x == y for x, y in zip(self, other))
            else:
                return False
        else:
            return NotImplemented

    def _compare(a, b, op):
        if isinstance(b, Integral):
            return bitarray([1 if op(e, b) else 0 for e in a])
        elif isinstance(b, Sequence):
            if len(a) == len(b):
                return bitarray([1 if op(x, y) else 0 for x, y in zip(a, b)])
            else:
                raise DimensionError(a, b, op=op)
        else:
            return NotImplemented

    def __pos__(self):
        return marray([+e for e in self], self._mod)

    def __neg__(self):
        return marray([-e for e in self], self._mod)

    def append(self, value):
        if isinstance(value, Integral):
            self._array.append(int(value) % self.mod)
        elif isinstance(value, Expr):
            self._array.append(value % self.mod)
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.append.__name__}() is not implemented for "
                                      f"non-Integral values")

    def reverse(self):
        self._array.reverse()

    def extend(self, values):
        if all(isinstance(v, (Integral, Expr)) for v in values):
            self._array.extend(e % self.mod if isinstance(e, Expr) else int(e) % self.mod for e in values)
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.extend.__name__}() is not implemented for "
                                      f"non-Integral values")

    def pop(self, index=None):
        return self._array.pop()

    def remove(self, value):
        if isinstance(value, Integral):
            self._array.remove(int(value))
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.remove.__name__}() is not implemented for "
                                      f"non-Integral values")

    def insert(self, index, value):
        if isinstance(value, Integral):
            self._array.insert(index, int(value) % self.mod)
        elif isinstance(value, Expr):
            self._array.insert(index, value % self.mod)
        else:
            raise NotImplementedError(f"{self.__class__.__name__}.{self.insert.__name__}() is not implemented for "
                                      f"non-Integral values")

    def __getitem__(self, item):
        if isinstance(item, slice):
            # if slice, return smaller marray, otherwise just return the value
            return marray(self._array[item], self._mod)
        else:
            return self._array[item]

    def __delitem__(self, item):
        del self._array[item]

    def __len__(self):
        return len(self._array)

    def __iadd__(self, other):
        if isinstance(other, Integral):
            other = int(other)
            for i, e in enumerate(self):
                self._array[i] = (e + other) % self.mod
            return self
        elif isinstance(other, Sequence):
            return self + other
        else:
            return NotImplemented
