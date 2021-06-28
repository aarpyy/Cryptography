import random, time, statistics, math, numpy


def StringToNum(s, base=128):
    result, i = 0, 0
    for c in s[::-1]:
        result += (base ** i) * ord(c)
        i += 1
    return result


def NumToString(n, base=128):
    string = ''
    while n > 0:
        index = 0
        while True:
            if pow(base, index + 1) >= n:
                break
            index += 1

        coeff = pow(base, index)
        letter = n // coeff
        n -= coeff * letter
        string += chr(letter)
    return string


# from internet
def ExtendedGCD(a, b):
    # Base Case
    if a == 0:
        return b, 0, 1

    gcd, x1, y1 = ExtendedGCD(b % a, a)

    # Update x and y using results of recursive
    # call
    x = y1 - (b // a) * x1
    y = x1

    return gcd, x, y


def GCD(a, b):
    return ExtendedGCD(a, b)[0]


def ModularInverse(x, m):
    k = ExtendedGCD(x, m)
    if k[0] == 1:
        return k[1]
    return None


# from internet
def IsPrime(num):
    # 2 is the only even prime, checks for 2 first
    if num == 2:
        return True
    # num & 1 returns intersection of binary 1 and binary num
    # if num is odd, it will always intersect with 1 and return True
    # not num & 1 filters all evens to return False, otherwise check below
    if not num & 1:
        return False
    # checks if fermat's little theroem works, will sometimes produce psuedo-primes
    return pow(2, num - 1, num) == 1


def ConfirmPrime(n):
    if n == 2:
        return True
    elif not n & 1:
        return False
    elif n % 3 == 0 or n % 5 == 0:
        return False
    else:
        sq = math.floor(math.sqrt(n))
        primes = []
        for num in range(3, sq):
            if IsPrime(num):
                primes.append(num)

        for num in primes:
            if n % num == 0:
                return False
        return True


def PercentChar(s):
    # function that returns the percent of letter characters in string s
    percent = 0
    for c in s:
        if 65 <= ord(c) <= 90 or 97 <= ord(c) <= 122:
            percent += 1
    return percent / len(s)


def RandomPrime(*args):
    # determines if user entered a lower range for prime
    base_2 = False
    if len(args) == 2:
        # if user entered lower base =/= 2, adjusts to nearest odd number
        lower = args[0]
        # if user entered lower base == 2, sets base_2 to True
        if lower == 2:
            base_2 = True
        elif lower % 2 == 0:
            lower += 1
        upper = args[1]
    # default lower base = 3, since most times function used for large primes, 2 not desired
    else:
        lower = 3
        upper = args[0]

    # if base_2, uses 2 as a base and increments by 1 (default) for generating random int
    if base_2:
        while True:
            prime = random.randrange(2, upper)
            if IsPrime(prime):
                return prime

    # if base =/= 2, generates random int starting at lower limit, incrementing by 2
    while True:
        prime = random.randrange(lower, upper, 2)
        if IsPrime(prime):
            return prime


def MakeMatrix(rows, cols):
    # gets input from user for rows of matrix to be made
    # splits into imaginary and real parts when necessary
    mat = []
    for i in range(rows):
        mat.append([])
        for j in range(cols):
            n = input(f"{i + 1},{j + 1}: ")
            if "," in n:
                c = n.split(",")
                real = int(c[0])
                imag = int(c[1])
                z = complex(real, imag)
                mat[i].append(z)
            else:
                mat[i].append(int(n))
    return numpy.array(mat)


def MultiplyMatrix():
    print("Enter dimensions for Matrix A: ")
    rowsA = int(input("Rows: "))
    colsA = int(input("Columns: "))

    A = MakeMatrix(rowsA, colsA)

    print("Enter dimensions for Matrix B: ")
    print(f"Rows: {colsA}")
    colsB = int(input("Columns: "))

    B = MakeMatrix(colsA, colsB)

    M = numpy.matmul(A, B)
    return M


def InvertMatrix():
    size = int(input("Enter matrix dimensions nxn: "))
    A = MakeMatrix(size, size)

    M = numpy.linalg.inv(A)
    return M


def MakeChineseRemainder():
    # gets lists of solutions and moduli from user for chinese remainder
    nums, mods = [], []
    equations = int(input("How many equations: "))
    for i in range(equations):
        print("x = ", end="")
        x = int(input())
        print("mod ", end="")
        m = int(input())
        print()
        nums.append(x)
        mods.append(m)
    return nums, mods


def ChineseRemainder(nums, mods):
    # initializes lists of moduli, M = product of all moduli
    M = 1
    alt_mods, mods_inverse = [], []

    for m in mods:
        M *= m
        mi = M // m
        alt_mods.append(mi)
        mods_inverse.append(ModularInverse(mi, m))

    # accumulates product of mi, mi-inverse, and original number
    acc = 0
    for i in range(len(nums)):
        acc = (acc + nums[i] * alt_mods[i] * mods_inverse[i]) % M

    return acc
