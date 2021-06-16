def StringToNum(s):
    n = ''
    for i in range(len(s) - 1, -1, -1):
        n += s[i]

    result = 0
    for i in range(len(n)):
        result += (128 ** i) * ord(n[i])
    return result


def NumToString(n):
    s = ''
    while True:
        index = 0
        while True:
            if (128 ** index) < n:
                index += 1
            else:
                index -= 1
                break

        k = 128 ** index
        m = n // k
        char = m * k
        n -= char
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