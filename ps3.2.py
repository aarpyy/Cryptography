def StringToNum(s):
    n = ''
    for i in range(len(s)-1, -1, -1):
        n += s[i]

    result = 0
    for i in range(len(n)):
        result += (128**i) * ord(n[i])
    return result


def PowerMod(x, y, m):
    result = x
    for _ in range(y-1):
        result = (result * x) % m
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


def Inverse(x, m):
    for i in range(m):
        if (x * i) % m == 1:
            return i


p = 7954763
g = 4
a = 1234092
A = PowerMod(g, a, p)

k = 6234043
ma = StringToNum("Hi!")
c1a = PowerMod(g, k, p)
c2a = (ma * PowerMod(A, k, p)) % p

print(f"\nProblem 2. a)\n"
      f"-------------\n"
      f"Alice publishes:\n"
      f"A: {A}\n\n"
      f"Bobs sends:\n"
      f"c1: {c1a}\n"
      f"c2: {c2a}\n")

c1b = 7724479
c2b = 3951884

c1bi = Inverse(c1b, p)
mb = (PowerMod(c1bi, a, p) * c2b) % p

print(f"\nProblem 2. b)\n"
      f"-------------\n"
      f"Bobs' message:\n"
      f"m: {NumToString(mb)}\n")