import tempfile
from pathlib import Path

import pytest

from cryptography318 import primesieve


def test_search():
    # 10 is between 7 and 11 which are the 4th and 5th primes
    assert primesieve.search(10) == (4, 5)
    assert primesieve.search(23) == (9, 9)

    with pytest.raises(ValueError):
        primesieve.search(0)


def test_extend():
    primesieve.extend(100)
    assert primesieve.data[-1] == 97

    primesieve.extend(10)
    assert primesieve.data[-1] == 97

    primesieve.extend(1000)
    assert primesieve.data[-1] == 997


def test_load():
    # Get temp path
    temp_path = Path(tempfile.gettempdir())
    temp_file = temp_path / "primesieve.txt"
    with open(temp_file, "w") as outfile:
        outfile.write("2\n3\n5\n7\n11\n13\n17\n19\n23\n29\n31\n37\n41\n43\n47\n53\n59\n61\n67\n71\n73\n79\n83\n89\n")

    primesieve.load(temp_file, overwrite=True)
    assert primesieve.data[-1] == 89

    # Now try appending
    with open(temp_file, "w") as outfile:
        outfile.write("89\n97\n101\n103\n107\n109\n113\n127\n131\n137\n139\n149\n151\n157\n163\n167\n173\n179")

    primesieve.load(temp_file)
    assert primesieve.data[-1] == 179

    # Now try partial
    with open(temp_file, "w") as outfile:
        outfile.write("101\n103\n107\n109\n113\n127\n131\n137\n139\n149\n151\n157\n163\n167\n173\n179")

    primesieve.load(temp_file)
    assert primesieve.data[-1] == 179

    # Now try discontinuous which should throw error
    with open(temp_file, "w") as outfile:
        outfile.write("211\n223\n227\n229\n233\n239\n241\n251\n257\n263\n269\n271\n277\n281\n283\n293\n307\n311\n313")

    with pytest.raises(ValueError):
        primesieve.load(temp_file)


def test_primerange():
    assert [*primesieve.primerange(10)] == [2, 3, 5, 7]
    assert [*primesieve.primerange(10, 20)] == [11, 13, 17, 19]


def test_primepi():
    assert primesieve.primepi(10) == 4
    assert primesieve.primepi(100) == 25
    assert primesieve.primepi(0) == 0


def test_getitem():
    # Assert retrieving nth prime
    assert primesieve[1] == 2
    assert primesieve[2] == 3

    assert primesieve[:10] == [2, 3, 5, 7, 11, 13, 17, 19, 23]
    assert primesieve[2:10] == [3, 5, 7, 11, 13, 17, 19, 23]


def test_contains():
    assert 2 in primesieve
    assert 3 in primesieve
    assert 4 not in primesieve
    assert -1 not in primesieve


def test_iter():
    # Can't just assert equal since primesieve is a generator that infinitely yields primes
    assert [p for p, _ in zip(primesieve, range(10))] == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
