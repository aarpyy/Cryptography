# The purpose of this file is to easily provide tools for use in this package that are used in multiple files

from warnings import warn, simplefilter
from functools import wraps, reduce
from math import sqrt
import numpy
from sympy import Symbol, evalf, re, im


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

    n = python_number(round(n, 3))

    if isinstance(n, int):
        return str(n)
    elif n.is_integer():
        return str(int(n))
    return str(n)


def append_and_return(obj, item):
    obj.append(item)
    return obj


def join_dict(*args: dict) -> dict:
    """Joins multiple dictionaries in a way that sums values of shared keys. Assumes all values
    support + method."""

    def join(dict1, dict2):
        for key in dict2:
            if key in dict1:
                dict1[key] += dict2[key]
            else:
                dict1[key] = dict2[key]
        return dict1

    def update(dict1, dict2):
        dict1.update(dict2)
        return dict1

    return reduce(lambda a, b: update(a, b) if not any(k in b for k in a) else join(a, b), args)


def replace_all(string, values, replace=''):
    """Replaces all instances of items from string values with string replace, default is to remove all items
    from values."""

    for v in values:
        string = string.replace(v, replace)
    return string


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

    if isinstance(number, (numpy.float16, numpy.float32, numpy.float64)):
        number = float(number)
    elif isinstance(number, (numpy.int16, numpy.int32, numpy.int64)):
        number = int(number)
    if isinstance(number, float) and number.is_integer():
        return int(number)
    return number


def isnumber(obj):
    types = (int, float, numpy.int16, numpy.int32, numpy.int64, numpy.float16,
             numpy.float32, numpy.float64)
    return isinstance(obj, types)


def fraction(n, limit=75):
    """Converts given number into a fraction if denominator < limit. Fractions are
    represented by strings."""

    n = python_number(n)

    for i in range(2, limit):
        x = round(n * i, 5)
        y = round(n * sqrt(i), 5)
        if isinstance(x, float) and x.is_integer():
            x = int(x)
        if isinstance(y, float) and y.is_integer():
            y = int(y)
        if isinstance(x, int):
            return f'{x}/{i}'
        if isinstance(y, int):
            return f'{y}/sqrt({i})'
    return string_reduce(n)


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
        raise ValueError(f"Unable to take product of two arrays of different length")
    if mod is not None:
        return reduce(lambda a, b: (a + b) % mod, map(lambda x, y: (x * y) % mod, obj, other))
    return sum(map(lambda x, y: x * y, obj, other))
