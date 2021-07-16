import random, multiprocessing, time
from crypto_functions import *


def MRHelper(candidate, n, k=40):
    if MillerRabinPrimality(n, k) and candidate % n == 0:
        print("{} is a factor of {}".format(n, candidate))


def ConfirmPrime1(n):
    """Uses infinitely deterministic primality test, checking for factors of all
    primes <= square root of candidate"""
    if KnownPrime(n) is not None:
        return KnownPrime(n)

    processes = []
    upper = (math.isqrt(n) + 1) | 1
    first = (upper // 3) | 1
    second = (2 * first) | 1
    for i in range(3, first, 2):
        p = multiprocessing.Process(target=lambda x: MillerRabinPrimality(x, 5), args=(i,))
        processes.append(p)
        p.start()
    for i in range(first, second, 2):
        p = multiprocessing.Process(target=lambda x: MillerRabinPrimality(x, 5), args=(i,))
        processes.append(p)
        p.start()
    for i in range(second, upper, 2):
        p = multiprocessing.Process(target=lambda x: MillerRabinPrimality(x, 5), args=(i,))
        processes.append(p)
        p.start()

    for process in processes:
        process.join()


if __name__ == '__main__':
    start1 = time.time()
    for _ in range(10):
        p = RandomPrime(pow(2, 40), certainty=5)
        ConfirmPrime1(p)
    print(f"Multiprocessing took {time.time()-start1:.2f}s")
    start2 = time.time()
    for _ in range(10):
        p = RandomPrime(pow(2, 40), certainty=5)
        ConfirmPrime(p)
    print(f"Single-threaded processing took {time.time()-start2:.2f}s")

