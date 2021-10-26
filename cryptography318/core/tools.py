# The purpose of this file is to easily provide tools for use in this package that are used in multiple files

from warnings import warn, simplefilter
from functools import wraps, reduce
from math import sqrt, gcd
from numpy import int16, int32, int64, float16, float32, float64
from fractions import Fraction as PyFraction
from numbers import *
from typing import SupportsRound


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""

    @wraps(func)
    def new_func(*args, **kwargs):
        simplefilter('always', DeprecationWarning)  # turn off filter
        warn_message = f"Call to deprecated function {func.__name__}."
        if func.__name__ == "ModularInverse":
            warn_message = warn_message[:-1] + f", instead use pow({args[0]}, -1, {args[1]})."
        if func.__name__ == "GCD":
            warn_message = warn_message[:-1] + f", instead use math.gcd({args[0]}, {args[1]})."
        warn(warn_message,
             category=DeprecationWarning,
             stacklevel=2)
        simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)

    return new_func


def string_reduce(n):
    """Reduces numbers in string form to shortest possible number for printing as matrix. Numbers
    are reduces such that floats with no floating point value are printed as integers and 0's with
    a negative prefix lose the prefix (ex. 2.0 would return as 2, but 2.01 would return as 2.01)"""
    if not isinstance(n, SupportsRound) or isinstance(n, (complex, PyFraction)):
        return str(n)

    n = python_number(round(n, 3))

    if isinstance(n, int):
        return str(n)
    elif n.is_integer():
        return str(int(n))
    return str(n)


def isnumber(*args):
    return NotImplemented


def r_append(obj: list, item) -> list:
    """Appends item to object and returns object. Helper function for lambda's
    involving appending."""

    obj.append(item)
    return obj


def join_dict(a: dict, b: dict, *args: dict, _null=0) -> dict:
    """Joins multiple dictionaries in a way that sums values of shared keys. Assumes all values
    are Numbers that support + method.

    :param a: first dictionary
    :param b: second dictionary
    :param args: any further dictionaries
    :param _null: default value for retrieving dictionary value
    """

    def join(dict1, dict2):
        for key in dict2:
            dict1[key] = dict1.get(key, _null) + dict2[key]
        return dict1

    return reduce(lambda x, y: join(x, y), args, join(a, b))


def replace_all(string: str, values: str, replace: str = '') -> str:
    """Replaces all instances of items from string values with string replace, default is to remove all items
    from values."""

    return reduce(lambda r, c: r.replace(c, replace), values, string)


def read_mm_int(fname='mathematica_numbers.txt'):
    """Reads file containing integer in syntax of mathematica's integer notation. Returns actual value."""

    with open(fname, "r") as f:
        lines = f.readlines()
        f.close()

    integers = {}
    var = None
    for i, line in enumerate(lines[2:]):
        if '=' in line:
            var_info = line.split('=')
            var = var_info[0][:-1]
            num = var_info[1]
            integers[var] = replace_all(num, '\\\n;')
        else:
            integers[var] += replace_all(line, '\\\n;')

    for var in integers:
        integers[var] = int(integers[var])

    return integers


def python_number(number):
    """Returns Python version of given number, instead of numpy's version. Useful in ensuring correct
    operations are performed with matrices (numpy.int64 * Matrix -> numpy.ndarray, not Matrix). Converts
    strings to Python floats."""

    if isinstance(number, str):
        number = float(number)

    elif hasattr(number, 'evalf'):
        number = number.evalf()

    if isinstance(number, int):
        return number
    elif isinstance(number, Rational):
        if number.denominator == 1:
            return number.numerator
        else:
            number = float(number)
    elif isinstance(number, Real):
        number = float(number)
    elif isinstance(number, Complex):
        return complex(number)

    if isinstance(number, float) and number.is_integer():
        return int(number)
    else:
        return number


def number_to_integer(number):
    """Converts non-integer number into integer if it is integer-valued, otherwise
    returns the original value."""
    if isinstance(number, int):
        return number
    if isinstance(number, Complex) and number.imag == 0:
        number = number.real
    if isinstance(number, float) and number.is_integer():
        return int(number)
    if isinstance(number, Rational) and number.denominator == 1:
        return int(number.numerator)
    return number


def fraction(n, limit=75):
    """Converts given number into a fraction if denominator < limit. Fractions are
    represented by strings."""

    n = python_number(n)

    if isinstance(n, int) or (isinstance(n, float) and n.is_integer()):
        return str(n)

    for i in range(2, limit):
        x = round(n * i, 5)
        y = round(n * sqrt(i), 5)

        if isinstance(x, float) and x.is_integer():
            x = int(x)
        if isinstance(y, float) and y.is_integer():
            y = int(y)

        if isinstance(x, int):
            if x == 0:
                return '0'
            return f'{x}/{i}'
        if isinstance(y, int):
            if y == 0:
                return '0'
            return f'{y}/sqrt({i})'

    return string_reduce(n)


def lcm(*args):

    def _lcm(a, b):
        g = gcd(a, b)
        x = a // g
        return x * b

    return reduce(lambda i, c: _lcm(i, c), args)


def evaluate(obj):
    if isinstance(obj, str):
        return eval(obj, {'sqrt': sqrt})

    array = []
    for i in range(len(obj)):
        array.append([])
        for j in range(len(obj[0])):
            array[i].append(eval(obj[i][j], {'sqrt': sqrt}))
    return array


def dot(obj, other, mod=None):
    """Equivalent of numpy.dot, accepts Matrix object as argument."""

    if len(obj) != len(other):
        raise ValueError(f"unable to take product of two arrays of different length")
    if isinstance(mod, int):
        return reduce(lambda a, b: (a + b) % mod, map(lambda x, y: (x * y) % mod, obj, other))
    return sum(map(lambda x, y: x * y, obj, other))
