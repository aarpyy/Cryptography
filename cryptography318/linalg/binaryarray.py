from abc import abstractmethod
from collections import Iterable, Sequence
from typing import Iterator, overload
from functools import reduce


class BinaryArray(Sequence, Iterable):

    _full_mask = {}

    __slots__ = '_value', '_len'

    def __new__(cls, arg, *, _len=None):
        self = super().__new__(cls)
        if isinstance(arg, int) and isinstance(_len, int):
            self._value = arg
            self._len = _len if isinstance(_len, int) else len(bin(arg)[2:])
        elif isinstance(arg, list):
            self._value = int(''.join(str(e) for e in arg), 2)
            self._len = len(arg)
        else:
            raise TypeError(f"invalid arguments for binary array")

        if self._len not in BinaryArray._full_mask:
            BinaryArray._full_mask[self._len] = pow(2, self._len) - 1
        return self

    @property
    def int(self):
        return self._value

    @property
    def invert(self):
        return BinaryArray(self.int ^ BinaryArray._full_mask[self._len], _len=self._len)

    @property
    def bits(self):
        return reduce(lambda r, c: (r + 1) if c == '1' else r, self, 0)

    def __str__(self):
        return "[" + " ".join(str(e) for e in self) + "]"

    def __len__(self):
        return self._len

    def __iter__(self) -> Iterator:
        _arr = bin(self.int)[2:]
        _arr = '0' * (self._len - len(_arr)) + _arr
        return (int(e) for e in _arr)

    def __getitem__(self, item):
        if isinstance(item, slice):
            # default values: start default is '0' which in our case is 2**self._len - 1
            if item.start is None:
                start = self._len
            else:
                start = -item.start % self._len

            if item.stop is None:
                stop = 0
            else:
                stop = -item.stop % self._len

            if start not in BinaryArray._full_mask:
                BinaryArray._full_mask[start] = pow(2, start) - 1
            if stop not in BinaryArray._full_mask:
                BinaryArray._full_mask[stop] = pow(2, stop) - 1

            if item.step in (1, None):
                _slice = self.int & BinaryArray._full_mask[start]
                return BinaryArray(_slice >> stop, _len=start - stop)
            else:
                raise TypeError(f"{self.__class__.__name__} does not support slicing with step other than 1")
        elif isinstance(item, int):
            # invert index since index of 0 means & w/ 2**(len-1)
            item = -(item + 1) % self._len
            if item in BinaryArray._full_mask:
                # _full_mask[item] + 1 = pow(2, item), so if bool(result) is true, then there is a 1 at index=item
                return int(bool(self.int & (BinaryArray._full_mask[item] + 1)))
            else:
                k = pow(2, item)
                BinaryArray._full_mask[item] = k - 1
                return int(bool(self.int & k))
        else:
            raise TypeError(f"index must be int or slice")

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return all(e1 == e2 for e1, e2 in zip(self, other))
        else:
            return False

    def __add__(self, other):
        if isinstance(other, BinaryArray):
            if len(other) == len(self):
                return BinaryArray(self.int ^ other.int, _len=self._len)
            else:
                raise ValueError(f"different length arrays")
        elif isinstance(other, int):
            if other % 2:
                return self.invert
            else:
                return +self
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                return BinaryArray(self.int ^ int(''.join(str(e) for e in other), 2), _len=self._len)
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (BinaryArray, Sequence, int)):
            return self + other
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, BinaryArray):
            if len(other) == len(self):
                return BinaryArray(self.int & other.int, _len=self._len)
            else:
                raise ValueError(f"different length arrays")
        elif isinstance(other, Sequence):
            if len(other) == len(self):
                return BinaryArray(self.int & int(''.join(str(e) for e in other), 2), _len=self._len)

    def __pos__(self):
        return BinaryArray(self.int, _len=self._len)

    def __neg__(self):
        return self.invert
