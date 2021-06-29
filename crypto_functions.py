import random, time, statistics, math, numpy


def StringToNum(s, base=128):
    result, i = 0, 0
    for c in s[::-1]:
        result += (base ** i) * ord(c)
        i += 1
    return result


def NumToString(n, base=128):
    string = ''
    chars = toBase(n, base)
    for c in chars:
        string += chr(c)
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
    return pow(x, -1, m)


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
    for m in mods:
        M *= m

    # maps list of moduli and their inverses to x and y respectively
    x, y = [], []
    for m in mods:
        mi = M // m
        x.append(mi)
        y.append(pow(mi, -1, m))

    # accumulates product of numbers and moduli and their inverses
    acc = 0
    for i in range(len(nums)):
        acc = (acc + nums[i] * x[i] * y[i]) % M

    return acc


def BSGS(g, h, p, prog=False, N=None):
    if N is None:
        N = p-1
    n = math.ceil(math.sqrt(N))

    if prog:
        increment = n // 25
        print("Starting Part 1: Creating Lists")
        A, B = {}, {}
        count = 0
        for i in range(n):
            if count >= increment:
                count = 0
                print("-", end="")
            count += 1
            j = i*n
            A[pow(g, i, p)] = i
            B[(h * pow(g, -1*j, p)) % p] = j

        print("\nDone With Part 1. Starting Part 2: Finding Matches")

        count = 0
        for e in A:
            if count >= increment:
                count = 0
                print("-", end="")
            count += 1
            if e in B:
                print()
                return A[e] + B[e]
        print()

    else:
        A, B = {}, {}
        for i in range(n):
            j = i * n
            A[pow(g, i, p)] = i
            B[(h * pow(g, -1 * j, p)) % p] = j

        for e in A:
            if e in B:
                return A[e] + B[e]
    return None


def toBase(n, base):
    # finds highest exponent of base
    exp = 0
    while True:
        if n < pow(base, exp + 1):
            break
        exp += 1

    # itereates starting at highest exponent down to 0, creates list of 'base' base
    list_coeff = [0 for _ in range(exp + 1)]
    ind = 0
    for i in range(exp, -1, -1):
        k = n // pow(base, i)
        list_coeff[ind] = k
        n -= k * pow(base, i)
        ind += 1

    return list_coeff

def fromBase(lst, base):
    acc = 0
    l = len(lst)
    for i in range(l):
        acc += lst[i] * pow(base, l-i-1)
    return acc


def PohligHellman(g, h, p, q, exp, prog=False):
    X = []

    r = pow(g, pow(q, exp-1), p)
    X0 = BSGS(r, pow(h, pow(q, exp-1), p), p, prog, q)
    X.append(X0)

    if prog:
        print(f"Found X0 = {X0}\n"
              f"Starting process for X1\n")

    for i in range(1, exp):
        exp_term = fromBase(X[::-1], q)
        h_term = pow(h * pow(pow(g, exp_term, p), -1, p), pow(q, exp-i-1), p)
        Xi = BSGS(r, h_term, p, prog, q)
        X.append(Xi)
        if prog:
            print(f"Found X{i} = {Xi}\n"
                  f"Starting process for X{i+1}\n")

    return fromBase(X[::-1], q)




