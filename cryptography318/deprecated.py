"""The purpose of this file is to easily provide warnings to old code that should not be used
but is left in the project."""

from warnings import warn, simplefilter
from functools import wraps


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
