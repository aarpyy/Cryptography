def smooth_factor(n, factors):
    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while n % f == 0:
            n //= f
            exp[i] += 1

    if abs(n) == 1:
        return exp
    else:
        return None


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
