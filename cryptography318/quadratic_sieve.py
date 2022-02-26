from math import e, sqrt, log, isqrt, ceil, prod, gcd, log2
from prime import primesieve, randprime, quadratic_residue, sqrt_mod
from utils import smooth_factor
import linalg

required_relations_ratio = 1.05
relations_found = 0
sieve_index = 0
sieve_array = []
sqrt_n = 0
trial_error = 20


def exp_value(exp, primes):
    """Calculates the value of a list of powers of primes. If only p is given, assumes list of primes to be
    from 2 to largest prime <= p. If list of exponents does not match the powers of the continuous ascending
    list of primes, this will compute incorrectly."""

    # raises each prime to the corresponding power in list exp, then reduces that list with multiplication
    return prod(map(lambda p, e: pow(p, e), primes, exp))


def sieve(fb, s1, s2, lp):

    global sieve_index, sieve_array

    size = len(fb)
    for i in range(10):
        for j in range(size):
            p = fb[j]
            i1 = s1[j] + sieve_index * p
            if p == 2:
                if i1 >= len(sieve_array):
                    sieve_array += [0] * ((i1 + 1) - len(sieve_array))
                sieve_array[i1] += lp[j]
            else:
                i2 = s2[j] + sieve_index * p
                imax = max(i1, i2)
                if imax >= len(sieve_array):
                    sieve_array += [0] * ((imax + 1) - len(sieve_array))
                sieve_array[i1] += lp[j]
                sieve_array[i2] += lp[j]
        sieve_index += 1


def trial_division(primes, b, relations_x, relations_div):
    global sqrt_n, trial_error, relations_found

    for x in range(len(sieve_array)):
        if not x or sieve_array[x] >= log2(2 * x * sqrt_n) - trial_error:
            b_x = b + x
            a = b_x * b_x - sqrt_n
            if (factors := smooth_factor(a, primes)) is not None:
                relations_x.append(a)
                relations_div.append(factors)
                relations_found += 1

    return relations_x, relations_div


def solve_matrix(bases, exp, primes, n):
    mod2 = []
    for i in range(len(exp[0])):
        mod2.append([])
        for j in range(len(exp)):
            mod2[i].append(exp[j][i] % 2)

    kernel = linalg.binary_kernel(mod2)

    for vector in kernel:  # iterate over basis of kernel
        e = map(lambda x: x // 2, linalg.vecmatmul(vector, exp))
        a = 1
        for j, k in zip(vector, bases):
            if j:
                a *= k

        b = exp_value(e, primes)
        p, q = gcd(a + b, n), gcd(a - b, n)

        print(f"p: {p}, q: {q}")

        if 1 < p < n:
            return p
        if 1 < q < n:
            return q

    return None


def quadratic_sieve(n, F=None):
    if F is None:
        F = pow(e, sqrt(log(n) * log(log(n)) / 2))

    primesieve.extend(F)
    primes = primesieve[:F]

    try:
        b = ceil(sqrt(n))
    except OverflowError:
        b = isqrt(n) + 1

    global sqrt_n
    sqrt_n = b

    t_sqrt = []
    log_p = []
    soln1 = []
    soln2 = []
    factor_base = []
    for p in primes:
        if quadratic_residue(n, p):
            factor_base.append(p)

            t = sqrt_mod(n % p, p)
            t_sqrt.append(t)
            log_p.append(round(log(p, 2)))
            soln1.append((t - b) % p)
            soln2.append((-t - b) % p)

    required_relations = int(len(factor_base) * required_relations_ratio)

    print(f"required relations: {required_relations}")

    smooth_factors = []
    smooth_input = []

    while 1:
        while 1:
            sieve(factor_base, soln1, soln2, log_p)
            smooth_input, smooth_factors = trial_division(primes, b, smooth_input, smooth_factors)
            if relations_found >= required_relations:
                break
        factor = solve_matrix(smooth_input, smooth_factors, primes, n)
        if factor is not None:
            return factor


if __name__ == "__main__":
    lower = pow(10, 6)
    upper = pow(10, 7)
    A = randprime(lower, upper)
    B = randprime(lower, upper)
    print(f"{A} * {B} = {(N := A * B)}")
    factor = quadratic_sieve(N, 1000)
    if factor is not None:
        print(f"factor of n: {factor}; n / factor = {N // factor}")
