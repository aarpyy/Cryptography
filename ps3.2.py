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

def NumToString(n, p):
    s = ''
    while n > 0:
        index = 0
        while True:
            if n >= (128 ** index):
                index += 1
            else:
                index -= 1
                break
        k = 128 ** index
        print(n, k, index)
        char = (n % k)
        n -= k * char
        s += chr(char)
    return s

def Inverse(x, m):
    for i in range(m):
        if (x * i) % m == 1:
            return i

p = 7954763
g = 4
a = 1234092
# A = PowerMod(g, a, p)
#
# k = 6234043
# m = StringToNum("Hi!")
# c1 = PowerMod(g, k, p)
# c2 = (m * PowerMod(A, k, p)) % p

A = 1575862
# c1 = 5138015
# c2 = 3284434

c1 = 7724479
c2 = 3951884

print(f"A: {A}\n"
      f"c1: {c1}\n"
      f"c2: {c2}")


c1i = Inverse(c1, p)
m = (PowerMod(c1i, a, p) * c2) % p
print(NumToString(1193121, p))
