from math import prod


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
