from math import isqrt


def integer_nth_root(n, k):
    return round(pow(n, 1 / k))


def is_square(n):
    """
    Replacement for Sympy's is_square function that follows almost
    explicitly the routine outlined in the link. The major difference
    is that due to the speed of Python's integer, we don't need
    to pre-mod our value to increase the speed of future mods, for each
    test we can mod n directly.

    References
    ----------
    https://mersenneforum.org/showpost.php?p=110896

    :param n: integer
    :return: if a * a == n for some integer a
    """
    if n < 0:
        return False
    elif n in (0, 1):
        return True

    m = n & 127  # n % 128
    if (m * 0x8bc40d7d) & (m * 0xa1e2f5d1) & 0x14020a:
        return False

    m = n % 63
    if (m * 0x3d491df7) & (m * 0xc824a9f9) & 0x10f14008:
        return False

    m = n % 25
    if (m * 0x1929fc1b) & (m * 0x4c9ea3b2) & 0x51001005:
        return False

    m = 0xd10d829a * (n % 31)
    if m & (m + 0x672a5354) & 0x21025115:
        return False

    m = n % 23
    if (m * 0x7bd28629) & (m * 0xe7180889) & 0xf8300:
        return False

    m = n % 19
    if (m * 0x1b8bead3) & (m * 0x4d75a124) & 0x4280082b:
        return False

    m = n % 17
    if (m * 0x6736f323) & (m * 0x9b1d499) & 0xc0000300:
        return False

    m = n % 11
    if (m * 0xabf1a3a7) & (m * 0x2612bf93) & 0x45854000:
        return False

    m = isqrt(n)
    return m * m == n