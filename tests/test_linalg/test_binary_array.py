from cryptography318.linalg.binaryarray import *
from typing import Sequence


def test_constructor():
    b = BinaryArray(4, _len=5)
    assert b == [0, 0, 1, 0, 0]

    b = BinaryArray([1, 0, 0, 1, 0])
    assert b == [1, 0, 0, 1, 0]
    assert len(b) == 5

    print(f"b.next(): {b.__next__()}")


if __name__ == '__main__':
    a = BinaryArray([1, 0, 0, 1, 0, 1])
    b = BinaryArray([0, 1, 0, 1, 1, 1])
    c = BinaryArray([1, 0, 1, 0, 0, 0])
    d = ([1, 0, 0, 1, 0, 1])
    print(b[2:4])
