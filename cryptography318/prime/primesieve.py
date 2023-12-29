from bisect import bisect
from collections import UserList
from itertools import count

from utils.misc import as_int


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

    def search(self, n):
        """
        Search for one or more values in the list. This is essentially
        list.index() except when our list does not contain the value,
        we will return the indices of the bordering primes.

        :param n:
        :return:
        """

        n = as_int(n)
        if n < 2:
            raise ValueError("Value must be greater than 1")

        # We either got a value error because the data is within the bounds of the list
        # but not a prime OR it is outside the bounds of the list and possibly prime.
        # If outside the bounds, we can't do anything so raise a ValueError
        if n > self.data[-1]:
            self.extend(n)

        # Otherwise we know the value was within the bounds of the list but composite
        # so lets find and return the indices of the bordering primes
        i = bisect(self.data, n)

        # If the value we are searching for exists, return that index exactly
        if self.data[i - 1] == n:
            return i, i
        return i, i + 1

    def extend(self, n):
        """
        Extends the list to include all primes up to ``n`` inclusive.
        If ``n`` is composite, list extends to include the smallest
        prime larger than ``n``.

        :param n: upper bound
        :return: None
        """
        n = int(n)
        if n <= self.data[-1]:
            return

        # Snippet from sympy.ntheory.generate.extend
        maxbase = int(n ** 0.5) + 1
        self.extend(maxbase)

        # Create a new sieve starting from sqrt(n)
        begin = self.data[-1] + 1
        newsieve = [*range(begin, n + 1)]

        # Now eliminate all multiples of primes in [2, sqrt(n)]
        for p in self.primerange(maxbase):
            # Start counting at a multiple of p, offsetting
            # the index to account for the new sieve's base index
            startindex = (-begin) % p
            for i in range(startindex, len(newsieve), p):
                newsieve[i] = 0

        # Merge the sieves
        self.data.extend(x for x in newsieve if x)

    def extend_to(self, i):
        """ Extends the list to include at least the first i primes """
        while len(self.data) < i:
            self.extend(self.data[-1] * 2)

    def load(self, fp, overwrite=False):
        """
        Loads primes from text file into sieve.
        :param fp: Reads from primes.txt local to project
        :param overwrite: If provided list of primes should replace existing list
        """

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

            # Otherwise we have discontinuous list and have to throw error
            raise ValueError(f"List of primes starting with {first_prime} is not continuous with current primesieve")

    def primerange(self, a, b=None):
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
        _, i = self.search(a)
        maxi = len(self.data) + 1
        while i < maxi:
            p = self.data[i - 1]
            if p >= b:
                break
            yield p
            i += 1

    def primepi(self, n):
        """
        Returns the number of primes less than or equal to n.
        """
        try:
            n = as_int(n)
        except ValueError:
            return 0
        if n < 2:
            return 0
        self.extend(n)
        a, b = self.search(n)
        return a

    def __contains__(self, n):
        """ Check if n is prime """
        try:
            n = as_int(n)
        except ValueError:
            return False
        if n < 2:
            return False
        a, b = self.search(n)
        return a == b

    def __iter__(self):
        """ Iterate over all primes """

        # Use __getitem__ to arbitrarily extend sieve
        for i in count(1):
            yield self[i]

    def __getitem__(self, n):
        """ Return the n-th prime number """
        if isinstance(n, slice):
            self.extend_to(n.stop)
            start = n.start or 1
            return self.data[start - 1: n.stop - 1: n.step]

        n = as_int(n)
        if n < 1:
            raise ValueError("Value must be greater than 0")
        if n >= len(self.data):
            self.extend_to(n)
        return self.data[n - 1]

    @property
    def list(self):
        return self.data

    @property
    def tail(self):
        return self.data[-1]


primesieve = Sieve()
