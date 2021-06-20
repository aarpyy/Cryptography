import random


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


def Inverse(x, m):
    return ExtendedGCD(x, m)[1]


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
