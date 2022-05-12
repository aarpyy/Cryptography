from math import prod
from typing import Sequence


def smooth_factor(n, factors):
    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while not (d := divmod(n, f))[1]:
            n = d[0]
            exp[i] += 1

    if abs(n) == 1:
        return exp
    else:
        return None


def eval_power(exp, primes):
    """Calculates the value of a list of powers of primes. If only p is given, assumes list of primes to be
    from 2 to largest prime <= p. If list of exponents does not match the powers of the continuous ascending
    list of primes, this will compute incorrectly."""

    # Raises each prime to the corresponding power in list exp, then reduces that list with multiplication
    return prod(map(lambda p, e: pow(p, e), primes, exp))


def from_base(lst, base):
    multiplier = 1
    acc = 0
    for e in lst:
        acc += e * multiplier
        multiplier *= base
    return acc


def extended_gcd(*args):
    # actual extended gcd function
    def ext_gcd(a, b):
        if a == 0:
            return b, 0, 1

        _g, x, y = ext_gcd(b % a, a)
        return _g, y - (b // a) * x, x

    # if just two arguments, return normal extended gcd
    if len(args) == 2:
        return ext_gcd(args[0], args[1])

    g, u, v = ext_gcd(args[0], args[1])  # gcd of first two args; u, v s.t. args sum to gcd
    values = [u, v]
    for e in args[2:]:
        g, u1, v1 = ext_gcd(g, e)
        # u1 s.t. u1 * gcd prev two values + v1 * curr val = gcd, so u1 applied to prev values
        values = [e * u1 for e in values]
        values.append(v1)
    return g, *values


def n_digits(n):
    return len(str(int(n)))


def _where(iterable, cond):
    for i, e in enumerate(iterable):
        if cond(e):
            yield i


def where(iterable, cond=None):
    """
    Returns index(es) where condition function returns true. Does not search recursively
    if given a non-flat list.

    Note
    ----
    This has the same name as numpy.where() and is intended as a similar function but is not
    the same in what it does.


    :param iterable: iterable to search
    :param cond: condition function
    :return: tuple of index(es)
    """

    if cond is None:
        def cond(x):
            return bool(x)

    return tuple(_where(iterable, cond))


def shape(o):
    """
    Returns the shape of nested sequences as a tuple.

    :param o: potentially not flat sequence
    :return: shape of all non-jagged portions of sequence
    """

    # Non-sequences don't have shape, return empty tuple
    if isinstance(o, (int, float)):
        return tuple()

    # Get shapes of all sequences
    s = [shape(e) for e in o if isinstance(e, Sequence)]

    # If not all the elements were sequences, we are done going deeper just return length
    if len(s) != len(o):
        return len(o),

    # Get length of most recent call to shape, if they are all the same then return
    x = set(a[0] for a in s)
    if len(x) == 1:
        return len(o), *s[0]
    else:
        return len(o),
