# The purpose of this file is to easily provide tools for use in this package that are not
# related to the mathematical content

from warnings import warn, simplefilter
from functools import wraps, reduce
from numpy import round


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
    s = str(n).split(".")
    if len(s) == 1 or s[1] != '0':
        return str(round(n, decimals=3))
    return s[0]


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


def read_mm_int(fname='mathematica_numbers.txt'):
    """Reads file containing integer in syntax of mathematica's integer notation. Returns actual value."""

    with open(fname, "r") as f:
        lines = f.readlines()
        f.close()

    integers = {}
    var = None
    for i, line in enumerate(lines):
        if i == len(lines) - 1:
            integers[var] += line[:-1]
        elif '=' in line:
            var_info = line.split('=')
            var = var_info[0][:-1]
            num = var_info[1]
            integers[var] = num[:-2]
            continue
        else:
            integers[var] += line[:-2]

    for var in integers:
        integers[var] = int(integers[var])

    return integers
