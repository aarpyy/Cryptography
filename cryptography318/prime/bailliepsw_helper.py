"""
This file is directly from https://github.com/smllmn/baillie-psw excluding jacobi. It is
included in a project otherwise consisting of entirely original work because my focus
was largely on factoring, but since the Miller-Rabin primality test is only deterministic
up to 3317044064679887385961981 and I needed to confirm that integers larger than that limit
were composite before factoring, I included this project as a further primality test.
"""


def D_chooser(candidate):
    """Choose a D value suitable for the Baillie-PSW test - from internet"""

    D = 5
    while jacobi(D, candidate) != -1:
        D += 2 if D > 0 else -2
        D *= -1
    return D


def U_V_subscript(k, n, U, V, P, Q, D):
    """Helper function that returns suitable integer values
    for Lucas pseudoprime test - from internet"""

    k, n, U, V, P, Q, D = map(int, (k, n, U, V, P, Q, D))
    digits = list(map(int, str(bin(k))[2:]))
    subscript = 1
    for digit in digits[1:]:
        U, V = U*V % n, (pow(V, 2, n) - 2*pow(Q, subscript, n)) % n
        subscript <<= 1
        if digit == 1:
            if not (P*U + V) & 1:
                if not (D*U + P*V) & 1:
                    U, V = (P*U + V) >> 1, (D*U + P*V) >> 1
                else:
                    U, V = (P*U + V) >> 1, (D*U + P*V + n) >> 1
            elif not (D*U + P*V) & 1:
                U, V = (P*U + V + n) >> 1, (D*U + P*V) >> 1
            else:
                U, V = (P*U + V + n) >> 1, (D*U + P*V + n) >> 1
            subscript += 1
            U, V = U % n, V % n
    return U, V


def LucasPseudoPrime(n, D, P, Q):
    """Perform the Lucas probable prime test - from internet
    Test performed with D, P, Q, s.t. P = 1,
    Q = (1 - D) / 4, and D is an integer that satisfies Jacobi(D / n) = -1"""

    U, V = U_V_subscript(n+1, n, 1, P, P, Q, D)

    if U != 0:
        return False

    d = n + 1
    s = 0
    while not d & 1:
        d >>= 1
        s += 1

    U, V = U_V_subscript(n+1, n, 1, P, P, Q, D)

    if not U:
        return True

    for r in range(s):
        U, V = (U*V) % n, (pow(V, 2, n) - 2*pow(Q, d*(2**r), n)) % n
        if not V:
            return True

    return False


def jacobi(a, n):
    """
    Function that returns Jacobi symbol, assuming n is odd - original
    --
    If n is odd prime: the Jacobi symbol for (a / n) = 0 iff
    a reduces to 0 mod n, = 1 iff a is a quadratic residue
    s.t. there exists some c where c^2 = a mod n,
    = -1 iff there does not exist some c where c^2 = a mod n
    """

    # checks to see if n is odd
    if not n & 1:
        raise ValueError(f"modulus 'n' must be odd")

    # if a is 0 or 1, automatically return a
    if a in (0, 1):
        return a

    # if a is not already mod n, reduce it mod n
    elif a != a % n:
        return jacobi(a % n, n)

    # if a is even, it can be reduced until it is odd, multiplying
    # the Legendre symbol for (a / n) for each iteration of a being even
    elif not a & 1:

        # checks if n = +/- 1 mod 8, multiplying by the result
        # of the Legendre symbol (a / n) = 1
        # otherwise, multiply result by Legendre symbol (a / n) = -1
        return jacobi(a >> 1, n) if n % 8 in (1, 7) else -jacobi(a >> 1, n)

    else:

        # if a is odd, return Jacobi of (n / a)
        # if a and n both reduce to 3 mod 4, return the negative
        # otherwise, return (n / a)
        return -jacobi(n % a, a) if a % 4 == 3 == n % 4 else jacobi(n % a, a)
