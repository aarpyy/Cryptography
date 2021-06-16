import time, math

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


def gcdExtended(a, b):
    # Base Case
    if a == 0:
        return b, 0, 1

    gcd, x1, y1 = gcdExtended(b % a, a)

    # Update x and y using results of recursive
    # call
    x = y1 - (b // a) * x1
    y = x1

    return gcd, x, y

a=49373499544578370101774047474458596291858628914848986697049154431838471330308629870034771349828765608506214150646075883127724969423387080923462618892131720056973983178666529419981335809893468246677965122600792388743123353
p=279760100626694627279092765868653174038469183666115163297748251165513105527406250598899180452349653311619118936729470268005858893110819678598861045418161109304147233015344698372446823092219881400249824039335535334245324399
h=215925731920882621780955531531474284542499806848480889459342752927318993381101610129486236944309350492375503279282310460788289345110865810956201729670178612655011476341155587241111062935153949275411989371869694010076887076

ai = gcdExtended(a, p-1)[1]
m = pow(h, ai, p)
print(NumToString(m))