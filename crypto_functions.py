import random, time, statistics
from matmult import *


def StringToNum(s):
    n = ''
    for i in range(len(s) - 1, -1, -1):
        n += s[i]

    result = 0
    for i in range(len(n)):
        result += (128 ** i) * ord(n[i])
    return result


def NumToString(n, base=128):
    s = ''
    while True:
        index = 0
        while True:
            if pow(base, index) >= n:
                index -= 1
                break
            index += 1

        k = pow(base, index)
        m = int(n // k)
        n -= m * k
        s += chr(m)
        if n == 0:
            return s


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


def IsPrime(num):
    if num == 2:
        return True
    if not num & 1:
        return False
    return pow(2, num - 1, num) == 1


def PercentChar(s):
    percent = 0
    for c in s:
        if 65 <= ord(c) <= 90 or 97 <= ord(c) <= 122:
            percent += 1
    return percent / len(s)


def RandomPrime(*args):
    base_2 = False
    if len(args) == 2:
        lower = args[0]
        if lower == 2:
            base_2 = True
            lower += 1
        elif lower % 2 == 0:
            lower += 1
        upper = args[1]
    else:
        lower = 3
        upper = args[0]

    if base_2:
        while True:
            prime = random.randrange(2, upper)
            if IsPrime(prime):
                return prime
    else:
        while True:
            prime = random.randrange(lower, upper, 2)
            if IsPrime(prime):
                return prime


def MakeMatrix(rows, cols):
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


def MakeChineseRemainder():
    nums = []
    mods = []
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
    M = 1
    mods_inverse = []
    alt_mods = []

    for m in mods:
        M *= m
        mi = M // m
        alt_mods.append(mi)
        mods_inverse.append(ModularInverse(mi, m))

    acc = 0
    for i in range(len(nums)):
        acc = (acc + nums[i] * alt_mods[i] * mods_inverse[i]) % M

    return acc
