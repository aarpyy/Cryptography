from crypto_functions import *

alt_mods, mods_inverse = [], []

mods = [RandomPrime(1000) for _ in range(5)]
M = 1
for m in mods:
    M *= m

start1 = time.time()
count = 0
while count < 50:
    alt_mods, mods_inverse = [], []
    for m in mods:
        mi = M // m
        alt_mods.append(mi)
        mods_inverse.append(pow(mi, -1, m))
    count += 1

print(f"For loop took: {time.time() - start1:.10f}s")

start2 = time.time()
count = 0
while count < 50:
    x = list(map(lambda a: M // a, mods))
    y = list(map(lambda a, b: pow(a, -1, b), x, mods))
    count += 1
print(f"Map took: {time.time() - start2:.10f}s")
print(x, y)

print(alt_mods, mods_inverse)

print(ChineseRemainder([4, 14, 19], [63, 142, 97]))

