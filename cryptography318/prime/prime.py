from collections import UserList
from itertools import islice
from math import isqrt, prod
from pathlib import Path
from random import choice, randrange

from cryptography318.prime.bailliepsw_helper import D_chooser, LucasPseudoPrime
from cryptography318.utils.utils import binary_search, Details


isprime_details = Details(**{
    "Miller-Rabin": None,
    "Baillie-PSW": None,
    "Direct": None,
})


def get_details():
    global isprime_details
    return str(isprime_details)


class Sieve(UserList[int]):
    """
    Unbound list of primes starting at 2. Object is iterable, index-able, print-able, searchable,
    and extendable. Unless another specific use is required, import primesieve object as a global
    variable for use. Intended to help with primality tests and factoring integers.
    """

    def __init__(self, data=None):
        if data is None:
            data = [2, 3, 5, 7, 11, 13]
        super().__init__(data)

    def search(self, value, *args):
        """
        Search for one or more values in the list. This is essentially
        list.index() except when our list does not contain the value,
        we will return the indices of the bordering primes.

        :param value:
        :param args:
        :return:
        """

        # If we are just searching for the one value, much easier
        if not args:
            try:
                # Try to just find it in the list
                return self.data.index(value)
            except ValueError:
                # We either got a value error because the data is within the bounds of the list
                # but not a prime OR it is outside the bounds of the list and possibly prime.
                # If outside the bounds, we can't do anything so raise a ValueError
                if value > self.data[-1] or value < self.data[0]:
                    raise ValueError(f"{{{value}}} is not contained within the primesieve")

                # Otherwise we know the value was within the bounds of the list but composite
                # so lets find and return the indices of the bordering primes
                i = binary_search(self.data, value, exist=False)
                if i is None:
                    raise ValueError(f"{{{value}}} is not contained within the primesieve")
                return i - 1, i

        # sorted() gets it back into a list, but now we've made sure to remove duplicates and add value
        values = sorted({value, *args})

        # If the largest value is larger than our upper bound or the smallest value is smaller than our
        # lower bound then we cannot return information for all values so raise an error
        if values[-1] > self.data[-1] or values[0] < self.data[0]:
            raise ValueError(f"{set(values)} are not contained within the primesieve")

        indices = []
        composite = []
        comp_indices = []

        last_index = 0

        # Try to get indices using index, add all that fail to our list, as well as the last prime
        # we saw before the fail
        for j, v in enumerate(values):
            try:
                last_index = self.data.index(v)
                indices.append(last_index)
            except ValueError:
                # Keep track of where we should insert this into indices
                comp_indices.append(j)
                # Add the composite value and the index of the last prime before we searched for this one
                composite.append((v, last_index))

        # If all values being searched for were primes, just return the values
        if len(composite) == 0:
            return tuple(indices)

        for k, (v, j) in enumerate(composite):
            i = binary_search(self.data, v, start=j, exist=False)
            if i is None:
                raise ValueError(f"{{{value}}} is not contained within the primesieve")
            indices.insert(comp_indices[k], (i - 1, i))

        return tuple(indices)

    def extend(self, n):
        """
        Extends the list to include all primes up to ``n`` inclusive.
        If ``n`` is composite, list extends to include the smallest
        prime larger than ``n``.

        :param n: upper bound
        :return: None
        """
        if n <= self.tail:
            return

        p = self.tail
        while True:
            p = next_prime(p)
            self.data.append(p)
            if p > n:
                return

    def load(self, file=None, overwrite=False):
        """
        Loads primes from text file into sieve.
        :param file: Reads from primes.txt local to project
        :param overwrite: If provided list of primes should replace existing list
        """

        if file is None:
            fp = Path(__file__).parent.absolute().joinpath("primes.txt")
        else:
            fp = Path(file)

        if not fp.is_file():
            raise FileNotFoundError(f"Unable to locate file {fp}")

        with open(fp, "r") as infile:

            # If we are overwriting and just accepting this as the list
            if overwrite:
                self.data = []
                for line in infile:
                    self.data.append(int(line.strip()))
                return

            first_prime = int(infile.readline().strip())

            # If we are starting from first prime, read all and then check to see if it read a higher prime
            if first_prime == 2:
                data = [2]
                for line in infile:
                    data.append(int(line.strip()))
                if data[-1] > self.data[-1]:
                    self.data = data
                return

            # If first prime is somewhere in list, find where it is and start adding from there
            if first_prime < self.data[-1]:
                # If this throws error then first_prime is not real prime
                start_index = self.data.index(first_prime)
                data = [first_prime]
                for line in infile:
                    data.append(int(line.strip()))

                if data[-1] > self.data[-1]:
                    self.data = self.data[:start_index] + data
                return

            # If first prime is our last prime, simple, just start appending to our list
            if first_prime == self.data[-1]:
                for line in infile:
                    self.data.append(int(line.strip()))
                return

            # If first prime read is next sequentially, we can accept it
            if next_prime(self.data[-1]) == first_prime:
                self.data.append(first_prime)
                for line in infile:
                    self.data.append(int(line.strip()))
                return

            # Otherwise we have discontinuous list and have to throw error
            raise ValueError(f"List of primes starting with {first_prime} is not continuous with current primesieve")

    def range(self, a, b=None):
        """
        Returns all primes within the given range. If a lower bound is not given
        the range starts at 2.

        :param a: first bound
        :param b: upper bound
        :return: all primes in range [a, b)
        """
        if b is None:
            b = a
            a = 2

        # Make sure that we've got a valid sized range
        if b <= a:
            return None

        # Extend sieve up to upper bound inclusive
        self.extend(b)

        # Get indices of upper and lower bounds
        i, j = self.search(a, b)
        if isinstance(i, tuple):
            i = i[1]

        # If it's a tuple, lets take the lower since we need to add 1 later since if it's not a tuple
        # we want to make sure we include that prime
        if isinstance(j, tuple):
            j = j[0]

        return islice(self.data, i, j + 1)

    @property
    def list(self):
        return self.data

    @property
    def tail(self):
        return self.data[-1]


primesieve = Sieve()


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
    elif n & 1:
        return False

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


def miller_rabin(n, k=40, *, details=False):
    """
    MRPrimality test reduces n - 1 to a power of 2 and an odd number, then
    tests if random `a` is a witness of n's composite-ness, testing with
    k random a's
    """

    d = n - 1
    r = 0
    while not d & 1:
        r += 1
        d >>= 1

    for _ in range(k):
        if not _mr_test(d, n):
            if details:
                global isprime_details
                isprime_details.add_details("Miller-Rabin", f"{d} is a witness to {n}'s composite-ness")
            return False

    return True


def _mr_test(d, n):
    """Helper function for miller_rabin which uses previously found d to
    check if random `a` is a witness to n's composite-ness"""

    a = randrange(2, n - 1)
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True

    # doubles d every time until d returns to original n-1 value
    while d != n - 1:
        x = pow(x, 2, n)
        d <<= 1

        if x == 1:
            return False
        elif x == n - 1:
            return True
    return False


def _miller_rabin_base_a(a, n):
    """Miller Rabin test with specific base of a"""

    if a >= n:
        a %= n

    if not a:
        return True

    q = n - 1
    k = 0
    while not q & 1:
        q >>= 1
        k += 1

    a = pow(a, q, n)
    if a == 1 or a == n - 1:
        return True
    for _ in range(k):
        # If we found any a^2 = -1 mod n then we know `a` is not a witness to n's compositeness
        if a == -1 or a == n - 1:
            return True

        # If we found an a^2 = 1 mod n where a != +/- 1 then a is definitely composite
        elif a == 1:
            return False
        a = pow(a, 2, n)

    return False


def miller_rabin_bases(bases, n, *, details=False):
    """Helper function that allows for a list of witnesses to be tested
    using MillerRabin_base_a function"""
    global isprime_details

    for a in bases:
        if not _miller_rabin_base_a(a, n):
            if details:
                isprime_details.add_details("Miller-Rabin", f"{a} is witness to {n}'s compositeness")
            return False
    if details:
        isprime_details.add_details("Miller-Rabin", f"{n} is probably prime")
    return True


def baillie_psw(n, mr=True):
    """
    Perform the Baillie-PSW probabilistic primality test on candidate.

    :param n: prime candidate
    :param mr: if Miller-Rabin test base 2 should be used
    :return:
    """

    # Check divisibility by a short list of primes less than 50
    if (res := known_prime(n)) is not None:
        return res

    # Now perform the Miller-Rabin primality test base 2
    if mr and not _miller_rabin_base_a(2, n):
        return False

    # Checks if number has square root
    if is_square(n):
        return False

    # Finally perform the Lucas primality test
    D = D_chooser(n)
    return LucasPseudoPrime(n, D, 1, (1 - D) // 4)


def known_prime(n):
    """Helper function, confirming prime candidate is not easily known"""

    known_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
                    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
                    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443,
                    449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
                    587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
                    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
                    853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983,
                    991, 997]

    for p in known_primes:
        if n == p:
            return True
        elif n % p == 0:
            return False
    return None


def isprime(n, *, details=False):
    """
    IsPrime function returns False iff the prime-candidate is composite, and True
    if the prime-candidate is probably prime.

    Uses deterministic variants of the Miller-Rabin Primality test, which, through
    the use of specific bases and ranges, can deterministically return True iff
    candidate is prime for n < 3317044064679887385961981. For all larger n,
    there is no  known set of bases that makes the MR test deterministic. Thus, a
    SPRP-test consisting of a Strong Lucas Pseudo-prime test and a Miller-Rabin
    test with 20 random bases `a`, s.t. 1 < a < n is used to determine if candidate is
    probably prime.
    """

    global isprime_details
    if details:
        isprime_details.clear_details()

    if n < 2:
        if details:
            isprime_details.add_details("Direct", f"{n} < 2")
        return False

    elif n < 10:
        if details:
            isprime_details.add_details("Direct", f"{n} < 10")
        return bool([0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0][n])

    # check for odds
    elif not n & 1:
        if details:
            isprime_details.add_details("Direct", f"{n} is even")
        return False

    # check for all other instances n != 6k +/- 1
    elif not n % 3:
        if details:
            isprime_details.add_details("Direct", f"{n} is divisible by 3")
        return False

    # This step is pretty useless unless primesieve is being used for something else or is
    # being purposefully generated, since it is constructed only with first 6 primes
    global primesieve
    if n in primesieve:
        if details:
            isprime_details.add_details("Direct", f"{n} is known prime")
        return True
    elif n < 2047:
        return miller_rabin_bases([2], n, details=details)
    elif n < 1373653:
        return miller_rabin_bases([2, 3], n, details=details)
    elif n < 9080191:
        return miller_rabin_bases([31, 73], n, details=details)
    elif n < 1050535501:
        return miller_rabin_bases([336781006125, 9639812373923155], n, details=details)
    elif n < 3215031751:
        return miller_rabin_bases([2, 3, 5, 7], n, details=details)
    elif n < 4759123141:
        return miller_rabin_bases([2, 7, 61], n, details=details)
    elif n < 1122004669633:
        return miller_rabin_bases([2, 13, 23, 1662803], n, details=details)
    elif n < 55245642489451:
        return miller_rabin_bases([2, 141889084524735, 1199124725622454117, 11096072698276303650], n, details=details)
    elif n < 7999252175582851:
        return miller_rabin_bases([2, 4130806001517, 149795463772692060, 186635894390467037, 3967304179347715805], n, details=details)
    elif n < 18446744073709551616:
        return miller_rabin_bases([2, 325, 9375, 28178, 450775, 9780504, 1795265022], n, details=details)
    elif n < 318665857834031151167461:
        return miller_rabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37], n, details=details)
    elif n < 3317044064679887385961981:
        return miller_rabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41], n, details=details)
    else:
        res = miller_rabin(n, k=40, details=details)
        if not res:
            return False

        # If miller rabin didn't say it was composite, then we should note that we used Baillie-PSW
        if details:
            isprime_details.add_details("Baillie-PSW", True)
        return baillie_psw(n, mr=False)


def randprime(a: int, b: int = None):
    """Uses combination of Miller-Rabin and Baillie-PSW primality tests to generate random prime

    Note
    ----
    If no lower bound is specified, 2 will never be generated, since bounds of sampling are set [3, b)

    :param a: integer starting point of range for random prime
    :param b: integer stopping point of range for random prime (exclusive)
    """

    # determines if user entered a lower and upper limit or just an upper
    if b is None:
        b, a = a, 3

    base_2 = a == 2
    a = a | 1  # If base_2, `a` isn't used, if not base 2 then even `a` is not prime so adjust to odd affects nothing

    global primesieve
    if b <= primesieve.tail:
        start, stop = primesieve.search(a), primesieve.search(b)
        if isinstance(start, tuple):
            start = start[1]
        if isinstance(stop, tuple):
            stop = stop[0]
        return choice(primesieve[start:stop])

    # If base_2 is True then it uses 2 as a base and increments by 1 (default) for generating random int
    # If base != 2, generates random int starting at lower limit, incrementing by 2
    if base_2:
        def rprime():
            return randrange(2, b)
    else:
        def rprime():
            return randrange(a, b, 2)

    while True:
        prime = rprime()
        if isprime(prime):
            return prime


def confirm_prime(n):
    """Uses deterministic AKS (Agrawal-Kayal-Saxena) primality test which
    returns True if-and-only-if n is prime"""

    if n < 2:
        return False
    elif n == 2:
        return True

    global primesieve
    if n in primesieve:
        return True

    # generates the n-th row of Pascal's triangle, if any of the coefficients != 0 mod n, n is not prime
    for k in range(1, (n + 1) // 2):
        res = 1
        if k > (n - k):
            k = n - k
        for i in range(k):
            res *= n - i
            res //= i + 1
        if res % n != 0:
            return False
    return True


def next_prime(n):
    """Returns first probable prime after number given"""

    if n < 2:
        return 2
    elif n < 11:
        return [2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11][n]
    elif n < primesieve.tail:
        i = primesieve.search(n)
        if isinstance(i, tuple):
            i = i[1]
        else:
            i += 1
        return primesieve[i]

    # Ensures that n starts at the nearest 6k + 1

    # `a` is the closest 6k + 1 to n
    a = n - (n % 6) + 1
    if a <= n:

        # If a <= n, try 6k + 5, only return if that's greater than n
        a += 4
        if a > n and isprime(a):
            return a

        # Otherwise, get `a` to 6k + 1, if this is prime, it's guaranteed > n so return
        a += 2
        if isprime(a):
            return a

        # Otherwise, since 6k + 1 above n didn't work, set n to 6k + 5 for below loop
        n = a + 4
    elif isprime(a):

        # If a > n and is prime, just return (this case only runs when n % 6 == 0 and n + 1 is prime)
        return a
    else:

        # Otherwise, start off n at a + 4 which is (6k + 1) + 4
        n = a + 4

    assert n % 6 == 5

    # Iterate up through each 6k +/- 1
    while True:
        if isprime(n):
            return n
        n += 2
        if isprime(n):
            return n
        n += 4


def prev_prime(n):
    """Returns first prime before number given"""

    if n < 3:
        raise ValueError(f"No primes exist < {n}")

    if n < 11:
        return [0, 0, 0, 0, 3, 3, 5, 5, 7, 7, 7][n]

    if n <= primesieve.tail:
        i = primesieve.search(n)
        if isinstance(i, tuple):
            i = i[0]
        else:
            i -= 1
        return primesieve[i]

    # ensures that n starts at the nearest 6k - 1 below
    r = n % 6
    if not r:
        n -= 1
    elif r == 1:
        n -= 2
    elif r <= 4:
        n -= r
        if isprime(n + 1):
            return n + 1
        n -= 1
    elif r == 5:
        n -= 4
        if isprime(n):
            return n
        n -= 2

    while True:
        if isprime(n):
            return n
        n -= 4
        if isprime(n):
            return n
        n -= 2

        # If we are sure that n passed was larger than 2, and we somehow end up below 2, just return 2
        if n < 2:
            return 2


def prime_range(a, b=None):
    """Constructs list of a <= primes < b"""
    if b is None:
        b = a
        a = 1

    if b < 2:
        return []

    global primesieve
    if b < primesieve.tail:
        primesieve.extend(b)

    if a < 2:
        a = 2
    indices = primesieve.search(a, b)
    start, stop = indices[0], indices[1]
    if isinstance(start, tuple):
        start = start[1]
    if isinstance(stop, tuple):
        stop = stop[1]
    return primesieve[start:stop]


def sqrt_mod(a, p):
    """Finds a solution for x to equation x^2 = a (mod p). If a solution is returned, a second
    solution s2 will also exist where s2 = -x (mod p)."""

    if not a:
        return 0
    elif not quadratic_residue(a, p):
        return None

    mod8 = p % 8
    if mod8 == 1:
        q = p - 1
        s = 0
        while not q & 1:
            q >>= 1
            s += 1

        z = randrange(2, p)
        while not quadratic_non_residue(z, p):
            z = randrange(2, p)

        m = s
        c = pow(z, q, p)
        t = pow(a, q, p)
        r = pow(a, (q + 1) // 2, p)

        while True:
            if t == 0:
                return 0
            if t == 1:
                return r

            i = 0
            x = t
            while x != 1:
                x = pow(x, 2, p)
                i += 1

            b = pow(c, pow(2, m - i - 1), p)
            c = pow(b, 2, p)
            m = i

            t = (t * c) % p
            r = (r * b) % p
    elif mod8 == 5:
        a2 = a + a
        v = pow(a2, (p - 5) // 8, p)
        i = (a2 * v * v) % p
        return (a * v * (i - 1)) % p
    else:
        return pow(a, (p + 1) // 4, p)


def lift_sqrt(root, n, modulus, q=None):
    """
    Given integer root that is the modular square root of ``n
    % modulus`` compute and return the square root of ``n % modulus * q``.
    That is, if q is not given, this function computes the square root
    of ``n % modulus ** 2``.

    Note
    ----
    If q is not given, modulus must be prime. If q is given it must be
    prime and modulus must be a prime power of q.

    References
    ----------
    This algorithm is Hensel's Lemma.

    :param root: square root of n mod modulus
    :param n: integer in FF(modulus)
    :param modulus: modulus of field
    :param q: prime
    :return: square root of n mod modulus * q
    """
    if q is None:
        q = modulus

    s = ((n - root * root) // modulus) * pow(root + root, -1, q)
    return (root + s * modulus) % (modulus * q)


def quadratic_residue(a, p):
    """Returns True if n is a quadratic residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == 1


def quadratic_non_residue(a, p):
    """Returns True if n is a quadratic non-residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == p - 1


def chinese_remainder(values, moduli):
    # Initializes lists of moduli, mod = product of all moduli
    mod = prod(moduli)

    # Maps list of moduli and their inverses to x and y respectively
    x, y = [], []
    for m in moduli:
        mi = mod // m
        x.append(mi)
        y.append(pow(mi, -1, m))

    # Accumulates product of number and moduli and their inverses
    acc = 0
    for i in range(len(values)):
        acc = (acc + values[i] * x[i] * y[i]) % mod

    return acc
