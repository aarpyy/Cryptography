from math import ceil, log10, log2, log


def smooth_factor(n, factors):
    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while (d := divmod(n, f))[1] == 0:
            n = d[0]
            exp[i] += 1

    if abs(n) == 1:
        return exp

    return None


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


def choose_f(digits):
    """
    Choose size of factor base. Hard coded values for sizes of n.

    :param digits: Number of digits in n
    :return: Size of factor base, f
    """
    if digits < 38:
        return 4200
    elif digits < 40:
        return 5600
    elif digits < 42:
        return 7000
    elif digits < 44:
        return 8400
    elif digits < 48:
        return 13000
    elif digits < 52:
        return 16000
    elif digits < 56:
        return 29000
    elif digits < 60:
        return 60000
    elif digits < 70:
        return 100000
    elif digits < 80:
        return 350000
    else:
        return 900000


def choose_m(digits):
    if digits < 45:
        return 30000
    elif digits < 52:
        return 65536
    elif digits < 88:
        return 196608
    else:
        return 589824
