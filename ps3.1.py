import time
from statistics import mean


p = 7954763
a = 1234092
b = 6234043
g = 4

start = time.time()

A = pow(g, a, p)
B = pow(g, b, p)
k = pow(A, b, p)

print(f"A: {A}\n"
      f"B: {B}\n"
      f"key: {k}")
print(f"This took {time.time() - start:.2f}s")