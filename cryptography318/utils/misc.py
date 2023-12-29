import operator


def as_int(n):
    try:
        if isinstance(n, bool):
            raise TypeError
        return operator.index(n)
    except TypeError:
        raise TypeError(f"{n} is not an integer")


def extended_gcd(*args):
    """
    Computes the Euclidean algorithm for finding the greatest common divisor for multiple integer values
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
