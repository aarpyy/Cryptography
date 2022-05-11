from math import log, isqrt, floor, sqrt, gcd
from typing import Union

from .curve import *
from cryptography318.prime.prime import isprime, primesieve


def safe_weierstrass(p):
    """
    Computes a non-singular Weierstrass elliptic curve safe for elliptic
    curve cryptography, and a point (x, y) that lies on the curve.

    :param p: composite integer
    :return: Weierstrass curve and point
    """

    while 1:
        a = rndm.randrange(p)
        x = rndm.randrange(p)
        y = rndm.randrange(p)

        # generate b based on point
        b = (pow(y, 2, p) - pow(x, 3, p) - a * x) % p

        # if singular curve, try again
        if 4 * pow(a, 3) + 27 * pow(b, 2) == 0:
            continue
        curve = Weierstrass(a, b, p)
        return curve, WeierstrassPoint(x, y, curve)


def safe_montgomery(p):
    """
    Implementation of Suyama's parameterization of Montgomery curves.

    Curve E over FF(q) for some prime q constructed s.t. E has order = 12 * k for some
    positive integer k. This helps in finding some q factor of N s.t. q is a product of small
    primes, in the computational process of Lenstra's ECM.

    References
    ----------
    Refer to section 2.3 p.4 of Gaj.pdf

    :param p: composite integer
    :return: Montgomery curve and point, or int if non-trivial factor found
    """
    sigma = rndm.randrange(6, p - 1)
    sig_sq = sigma * sigma
    u = (sig_sq - 5) % p
    v = (4 * sigma) % p
    u_sub_v = v - u
    u_3 = u * u * u
    try:
        a_sig = ((u_sub_v * u_sub_v * u_sub_v) * (3 * u + v) * pow(4 * u_3 * v, -1, p) - 2) % p
    except ValueError:

        # If composite m for FF(m) then this curve is likely intended for Lenstra's ECM and gcd is factor
        return gcd(4 * u_3 * v, p)

    x, z = u_3 % p, (v * v * v) % p
    E = Montgomery(a_sig, p)
    P = MontgomeryPoint(x, z, E)
    return E, P


def z_1_montgomery(p):
    """
    Implementation of Suyama's parameterization of Montgomery curves as above except that the
    coordinates x, z are optimized such that x = x / z, z = z / z, resulting in z = 1 allowing
    for the use of add_z1 and ladder_z1 methods which saves a multiplication.

    :param p: composite integer
    :return: Montgomery curve and point, or int if non-trivial factor found
    """
    sigma = rndm.randrange(6, p - 1)
    sig_sq = sigma * sigma
    u = (sig_sq - 5) % p
    v = (4 * sigma) % p
    u_sub_v = v - u
    u_3 = u * u * u
    try:
        # Value error occurs if a modular inverse does not exist
        u_3_inv = pow(4 * u_3 * v, -1, p)
    except ValueError:

        # If composite p then lack of inverse means gcd is greater than 1
        return gcd(4 * u_3 * v, p)

    a_sig = ((u_sub_v * u_sub_v * u_sub_v) * (3 * u + v) * u_3_inv - 2) % p

    # here, z would be v^3, but instead x = x / z and z = z / z so x = u^3 / v^3 and z = 1
    # this allows for one less multiplication to be used in add and ladder methods
    u_div_v = u_3 * u * 4 * u_3_inv  # 4u^4 / 4u^3v = u / v
    x = u_div_v * u_div_v * u_div_v % p  # u^3 / v^3
    E = Montgomery(a_sig, p)
    P = MontgomeryPoint(x, 1, E)
    return E, P


def ecm_mont_basic(N, *, B1=None, B2=None, retry=50):
    """
    A more basic version of the below elliptic factorization over a Montgomery form elliptic curve.
    This method opts for a direct multiplication approach for phase 2. Phase 1 is computed
    identically to find Q = kP. If kP missed the order of E(FF(q)), for some q factor of N, by a small
    integer, then we compute some multiple of kP to find the factor. In this phase 2, the integer
    multiplied against kP is the z coord of the product of Q with each prime p, B1 < p < B2. This results
    in a sometimes quicker calculation if N is small enough, but in general the below function should
    be used for elliptic curve factorization.
    """
    if B1 is None:
        from math import e
        B1 = (int(pow(e, isqrt(int((log(N) * log(log(N))) / 2)))) | 1) + 1  # | 1 + 1 ensures even

    if B2 is None:
        B2 = B1 * 100

    primesieve.extend(B2)

    # with initial point P, get k to find initial Q = kP
    k = 1
    for p in primesieve.range(B1 + 1):
        k *= floor(log(B1, p))

    trials = 0
    while trials <= retry:
        trials += 1

        # get montgomery curve for phase 1
        res = safe_montgomery(N)
        if isinstance(res, int):
            return res

        E, P = res
        Q = P.ladder(k)

        # if found factor already, return, otherwise continue to phase 2
        q = gcd(Q.z, N)
        if 1 < q < N:
            return q

        d = 1

        for p in primesieve.range(B1, B2):
            d = (d * Q.ladder(p).z) % N

        q = gcd(d, N)
        if 1 < q < N:
            return q

    return None


def ecm_mont(N, *, B1=None, B2=None, retry=50):
    """
    Performs Lenstra's elliptic curve factorization method with elliptic curve E
    in Montgomery form over FF(N) with Suyama's parameterization.

    References
    ----------
    [1] Gaj K. et al. Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware
    https://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf

    :param N: number to be factored
    :param B1: limit for phase 1; start for phase 2
    :param B2: limit for phase 2
    :param retry: number of trials to be run
    :return: one non-trivial factor of N or None
    """
    if B1 is None:
        from math import e
        q = isqrt(N)
        B1 = int(pow(e, sqrt(log(q) * log(log(q)) * 0.5)))
    elif B1 & 1 == 1:
        B1 += 1

    if B2 is None:
        B2 = B1 * 100
    elif B2 & 1 == 1:
        B2 += 1

    primesieve.extend(B2)

    # with initial point P, get k to find initial Q = kP
    k = 1
    for p in primesieve.range(2, B1):
        k *= floor(log(B1, p))

    # Type hints for more helpful indexing by PyCharm
    D: int
    S: list[Union[int, MontgomeryPoint]]
    E: Montgomery
    P: MontgomeryPoint
    Q: MontgomeryPoint
    Q_0: MontgomeryPoint
    R: MontgomeryPoint
    D: MontgomeryPoint
    Q_D: MontgomeryPoint

    D = isqrt(B2)
    S = [0] * (D + 1)

    for _ in range(retry):

        # Get montgomery curve for phase 1
        res = safe_montgomery(N)

        # If in finding curve a factor was found, return it
        if isinstance(res, int):
            return res

        E, P = res
        Q_0 = P.ladder(k)

        # If found factor already, return, otherwise continue to phase 2
        q = gcd(Q_0.z, N)
        if 1 < q < N:
            return q

        # Begin Phase 2

        # Initialize full table S (algorithm 5) where S = [Q_0, 2 * Q_0, 3 * Q_0, ..., D/2 * Q_0]
        Q = Q_0
        Q_2 = Q_0.double()
        S[0] = Q
        S[1] = Q.add(Q_2, Q_0)
        for i in range(2, D + 1):
            S[i] = S[i - 1].add(Q_2, S[i - 2])

        d = 1

        R = Q_0.ladder(B1)          # B * Q
        T = Q_0.ladder(B1 - D)      # (B - D) * Q
        Q_D = S[D]
        for i in range(1, D + 1):
            d = (d * (R.x * S[i].z - R.z * S[i].x))

            # Swap T, R to keep track of difference between points for montgomery addition
            T, R = R, R.add(Q_D, T)     # R = (B + kD)Q + DQ, T = (B + (k-1) * D)Q = diff

        q = gcd(d, N)
        if 1 < q < N:
            return q

    return None


def ecm_weierstrass(N, *, B=None, retry=50):
    """
    Lenstra's Elliptic Curve Factorization method over a short Weierstrass curve
    with parameters a, b s.t. 4a^3 + 27b^2 != 0
    """

    if B is None:
        from math import e
        L = pow(e, isqrt(int(log(N) * log(log(N)))))
        B = int(pow(L, 1 / pow(2, .5)))

    trials = 0
    while trials <= retry:
        trials += 1

        E, P = safe_weierstrass(N)

        try:
            for j in range(1, B):
                P *= j
        except ValueError as exc:

            # This should always be a valid number, but in case it isn't ignore this one and go next curve
            digit = str(exc).strip()
            if digit.isdigit():
                # If an error was thrown, a number was inverted that is not invertible, take the gcd
                k = int(digit)
                q = gcd(k, N)
                if 1 < q < N:
                    return q

    return None


def lenstra_ecm(N, *, B=None, retry=50):
    """
    General purpose function to compute Lenstra's elliptic curve factorization method
    on a composite integer N. If number is a 32-bit integer, use a Weierstrass curve
    to find a factor. If this succeeds, return the factor. If this fails, or number is
    greater than a 32-bit integer, use a Montgomery curve to find a factor.
    """
    if isprime(N):
        return N
    if N < 0x100000000:
        f = ecm_weierstrass(N, B=B, retry=retry)
        if f:
            return f
    return ecm_mont(N, B1=B, retry=retry)
