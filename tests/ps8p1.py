from cryptography318.crypto_functions import *
import numpy
import random
from functools import reduce


# given
g = 5
h = 676468030
p = 717994967
B = 20
# primes <= B
primes = PrimesLT(B)


def find_log_g():
    smooth_nums = []
    # to start, we know the solution for log_g_5 = 1
    exp_known = [0] * 8
    exp_known[2] = 1
    while True:
        i = random.randrange(2, p-1)
        c = pow(g, i, p)
        if BSmoothQ(c, factors=primes):
            factors = factor_base_exp(c, primes)
            unknown = 0
            new_known = exp_known[:]
            for j in range(len(factors)):
                # if new random i produces smooth number with only one unknown log_g_x, add it to matrix
                # so that can solve for that x
                if factors[j] != 0 and exp_known[j] == 0:
                    new_known[j] = 1
                    unknown += 1
            # if more than one unknown in new smooth number, don't add it
            if unknown > 1:
                continue
            exp_known = new_known[:]
            factors.append(i)
            smooth_nums.append(factors)

        if len(smooth_nums) == 10:
            break

    # this is the system I solved to find each log written in homework
    system = Matrix(smooth_nums, aug=True).astype(numpy.float64)
    print(system)


find_log_g()


# wrote out dictionary with keys as primes and values as solutions log_g so that I can use for next step
logs = {2: 386603800, 3: 584671346, 5: 1, 7: 138593831, 11: 183202078, 13: 163170104, 17: 697282118, 19: 8417856}

"""Here, I am finding some k s.t. h * g^-k is 20-smooth"""
k = 0
while True:
    k += 1
    x = (h * pow(g, -k, p)) % p
    if BSmoothQ(x, B):
        exponents = factor_base_exp(x, primes)
        print(k)
        """
        map is taking in a list of exponents (the powers of each prime < 20 that equal our 20-smooth number) and a list of log values (log_g_x for x = each prime < 20) and multiplying the log value by the exponent
        
        reduce is summing the list returned by map, which when added to k returns result of log_g_h
        
        since this print-to-pdf doesn't include terminal, my result for g^x = h is x = 9683732 mod p-1
        """
        print(reduce(lambda a, b: a + b, list(map(lambda i, n: i * logs[n], exponents, logs))) + k)
        break
