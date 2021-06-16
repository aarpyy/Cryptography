import time
from statistics import mean

def PowerMod(x, y, m):
    result = x
    for _ in range(y-1):
        result = (result * x) % m
    return result


p = 7954763
a = 1234092
b = 6234043
g = 4

start = time.time()

A = PowerMod(g, a, p)
B = PowerMod(g, b, p)
k = PowerMod(A, b, p)

print(f"A: {A}\n"
      f"B: {B}\n"
      f"key: {k}")
print(f"PowerMod took {time.time() - start}s")