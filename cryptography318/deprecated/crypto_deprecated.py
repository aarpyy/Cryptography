
def combinations_cumulative(found, source):
    """Function takes in set of combinations of k-choices without replacement and returns a new
    set of combinations of k + 1 choices without replacement.

    :param found: list of already found combinations, cannot be empty
    :param source: list of source for combinations"""

    new_combination = []
    for choice in found:
        for e in source:
            if e not in choice:
                temp = choice[:]
                temp.append(e)
                add = True
                for c in new_combination:
                    different = False
                    i = 0
                    while i < len(temp) and not different:
                        t = temp[i]
                        i += 1
                        if t not in c:
                            different = True
                            continue
                    if not different:
                        add = False
                        break
                if add:
                    new_combination.append(temp)

    return new_combination


def DSA(D, S1, S2, g, p, q, A):
    S2_inv = pow(S2, -1, q)
    V1 = (D * S2_inv) % q
    V2 = (S1 * S2_inv) % q
    return ((pow(g, V1, p) * pow(A, V2, p)) % p) % q == S1




def _factor_with_known(p, q, n):
    """Helper function for all integer factoring functions, which further factors integer given known factors.

    :param p: integer that divides n, not necessarily prime
    :param q: same as p
    :param n: integer to be factored
    :return: dictionary. keys: all primes factors of n, values: powers of prime factors"""

    # if inputs are non-trivial, try to factor more, if trivial, ignore
    fact_p = factor(p) if p not in [1, n] else {}
    fact_q = factor(q) if q not in [1, n] else {}

    # if failed to factor p or q further, add them to dictionary as non-prime factors
    if fact_q is None:
        fact_q = {q: 1}
    if fact_p is None:
        fact_p = {p: 1}

    factors_known = join_dict(fact_q, fact_p)

    factors = {}
    for f in factors_known:
        factors[f] = 0

    for f in factors:
        while n % f == 0:
            n //= f
            factors[f] += 1

    if n == 1:
        return factors
    if isprime(n):
        return join_dict(factors, {n: 1})

    more_factors = factor(n)
    return join_dict(factors, {n: 1}) if more_factors is None else join_dict(more_factors, factors)

