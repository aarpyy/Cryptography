import json
from math import ceil, log, log10, log2, prod


class Details:
    __slots__ = ("_details",)

    def __init__(self, **kwargs):
        self._details = kwargs

    def __str__(self):
        return json.dumps(self._details)

    def add_details(self, name, value):
        self._details[name] = value

    def clear_details(self):
        for name in self._details:
            self._details[name] = None


def binary_search(a, key, *, start=0, end=None, exist=True):
    """
    Given a sorted sequence, find key in the sequence if it exists,
    or the index at which key would be inserted to in order to keep
    sort order. This is an attempted replication of Java's binarySearch.

    If exist is set to False, then if the value is found in the list,
    None will be returned. This configuration is useful when we want
    to force a unique list and only insert a value if it does not exist.

    :param a: sorted sequence
    :param key: value to search
    :param start: index to start search at
    :param end: index to stop search at
    :param exist: if value is allowed to exist in sequence
    :return: index value should exist at
    """
    low = start
    high = end or (len(a) - 1)
    while low <= high:
        mid = (low + high) // 2
        if a[mid] < key:
            low = mid + 1
        elif a[mid] > key:
            high = mid - 1
        elif exist:
            return mid
        else:
            return None

    return high + 1


def smooth_factor(n, factors):
    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while not (d := divmod(n, f))[1]:
            n = d[0]
            exp[i] += 1

    if abs(n) == 1:
        return exp

    return None


def eval_power(exp, primes):
    """Calculates the value of a list of powers of primes. If only p is given, assumes list of primes to be
    from 2 to the largest prime <= p. If list of exponents does not match the powers of the continuous ascending
    list of primes, this will compute incorrectly."""

    # Raises each prime to the corresponding power in list exp, then reduces that list with multiplication
    return prod(pow(p, e) for p, e in zip(primes, exp, strict=True))


def from_base(lst, base):
    multiplier = 1
    acc = 0
    for e in lst:
        acc += e * multiplier
        multiplier *= base
    return acc


def extended_gcd(*args):
    """
    Computes the Euclidean algorithm for finding the greatest common divisor for multiple integer values
    :param args: integer values
    :return:
    """

    # actual extended gcd function
    def ext_gcd(a, b):
        if a == 0:
            return b, 0, 1

        _g, x, y = ext_gcd(b % a, a)
        return _g, y - (b // a) * x, x

    # if just two arguments, return normal extended gcd
    if len(args) == 2:
        return ext_gcd(args[0], args[1])

    g, u, v = ext_gcd(args[0], args[1])  # gcd of first two values; u, v s.t. values sum to gcd
    values = [u, v]
    for e in args[2:]:
        g, u1, v1 = ext_gcd(g, e)
        # u1 s.t. u1 * gcd prev two values + v1 * curr val = gcd, so u1 applied to prev values
        values = [e * u1 for e in values]
        values.append(v1)
    return g, *values


def n_digits(n, radix=10):
    """
    Determines the minimum number of digits required to represent n in the given radix, defaults to 10
    :param n: Integer value
    :param radix: Radix of n
    :return: Number of digits required to represent n
    """

    match radix:
        case 10:
            return ceil(log10(n))
        case 2:
            return ceil(log2(n))
        case _:
            return ceil(log(n) / log(radix))
