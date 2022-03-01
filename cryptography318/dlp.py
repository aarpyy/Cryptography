import random
from math import isqrt, gcd, sqrt
from functools import reduce

from prime import primesieve
from utils import smooth_factor, from_base


def baby_step_giant_step(g, h, p, order=None):
    """Function attempts to solve DLP using classic baby-step-giant-step algorithm."""
    if order is None:
        order = p - 1

    n = isqrt(order) + 1

    # find lists A = g^i B = h * g^-jn s.t. A[i] == B[j] for some indices i, j, this collision allows us to solve
    A = list(map(lambda e: pow(g, e, p), range(n)))
    B = list(map(lambda e: (h * pow(g, -e * n, p)) % p, range(n)))

    # convert to set for intersection calculation
    a, b = set(A), set(B)
    U = a.intersection(b)

    # if empty set, no collisions found
    if not U:
        return None

    # otherwise, find first indices of match and use to solve
    match = U.pop()
    i, j = A.index(match), B.index(match)
    return i + (j * n)


def pollard_rho_dlp(g, h, p, order=None):
    """Uses Pollard's Rho algorithm for logarithms to solve given discrete log problem. Function will run
    indefinitely until a solution is found.

    :param g: integer base
    :param h: integer solution to g^x for some x
    :param p: integer prime modulus
    :param order: integer order of g (smallest integer s.t. pow(g, order, p) == 1
    :return: solution to g^x = h for integer x"""

    xstate = (1, 0, 0)
    ystate = (1, 0, 0)

    if order is None:
        order = p - 1

    while True:
        xstate = calculate_state(xstate, g, h, p)
        ystate = calculate_state(calculate_state(ystate, g, h, p), g, h, p)

        if xstate[0] == ystate[0]:
            try:

                # try to return result right away, fails if beta value is not invertible mod q
                return (xstate[1] - ystate[1]) * pow(ystate[2] - xstate[2], -1, order)
            except ValueError:

                # try to reduce entire equation by gcd, then try to invert again
                alpha = xstate[1] - ystate[1]
                beta = ystate[2] - xstate[2]
                e = gcd(alpha, beta, order)
                if e > 1:
                    alpha //= e
                    beta //= e
                    mod = order // e

                    # if reduced until invertible, find solution
                    if gcd(beta, mod) == 1:
                        log_g_h = alpha * pow(beta, -1, mod)

                        # current solution is mod q // e, but real solution could be mod any increment of
                        # q up until its original value, find real solution by checking
                        for i in range(e - 1):
                            if pow(g, log_g_h, p) == h % p:
                                return log_g_h
                            log_g_h += mod
                continue


def index_calculus_dlp(g, h, p):
    """Function attempts to solve DLP through index calculus algorithm, computing a series of smaller dlp's
    used to solve the larger."""

    raise NotImplementedError("ICM currently in development.")

    from math import e, log

    B = int(pow(e, sqrt((log(p) * log(log(p))) / 2)))

    print(f"B: {B}")

    primesieve.extend(B)
    primes = primesieve[:B]

    print(f"primesieve extended")

    # currently brute forces solutions to log_g_x with x for each prime <= B
    logs = []
    for n in primes:
        lg = pollard_rho_dlp(g, n, p)
        print(f"log {n} mod {p} = {lg}")
        logs.append(lg)

    k = 0
    while 1:
        k += 1
        x = (h * pow(g, -k, p)) % p
        if (exponents := smooth_factor(x, primes)) is not None:
            # reduce sums list returned by map, which returns list of products of exponents and logs
            return reduce(lambda a, b: a + b, list(map(lambda i, l: i * l, exponents, logs))) + k


def calculate_state(state, g, h, p):
    """This function is a helper method for Pollard's Rho algorithm for solving DLP's."""

    x, alpha, beta = state[0], state[1], state[2]

    if 0 <= x < p // 3:
        x *= g
        if x >= p:
            x %= p
        alpha += 1
        if alpha >= p - 1:
            alpha %= (p - 1)
        return x, alpha, beta
    elif p // 3 <= x < 2 * p // 3:
        x = pow(x, 2, p)
        alpha *= 2
        beta *= 2
        if x >= p:
            x %= p
        if alpha >= p - 1:
            alpha %= (p - 1)
        if beta >= p - 1:
            beta %= (p - 1)
        return x, alpha, beta
    elif 2 * p // 3 <= x < p:
        x *= h
        if x >= p:
            x %= p
        beta += 1
        if beta >= p - 1:
            beta %= (p - 1)
        return x, alpha, beta


def pohlig_hellman(g, h, p, q, exp, prog=False):
    """Function attempts to solve DLP mod p where order of g, q, is some prime raised to a power."""

    X = []
    if prog:
        print("Starting process for X0")

    r = pow(g, pow(q, exp - 1), p)
    X0 = baby_step_giant_step(r, pow(h, pow(q, exp - 1), p), p, q)
    X.append(X0)

    if prog:
        print(f"Found X0 = {X0}\n")

    for i in range(1, exp):
        if prog:
            print(f"Starting process for X{i}")
        exp_term = from_base(X[::-1], q)
        h_term = pow(h * pow(pow(g, exp_term, p), -1, p), pow(q, exp - i - 1), p)
        Xi = baby_step_giant_step(r, h_term, p, q)
        X.append(Xi)
        if prog:
            print(f"Found X{i} = {Xi}\n")

    return from_base(X[::-1], q)


if __name__ == "__main__":
    import prime
    p = prime.randprime(pow(10, 6), pow(10, 7))
    g = random.randrange(2, p)
    x = random.randrange(2, p - 1)
    h = pow(g, x, p)
    y = index_calculus_dlp(g, h, p)
    print(f"g^y: {pow(g, y, p)}, g^x: {h}")
