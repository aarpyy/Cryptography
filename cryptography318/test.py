from prime import *
import time


start1 = time.time()
for _ in range(500):
    ConfirmPrime2(random.randrange(pow(2, 10)))

print(f"2 took: {time.time()-start1:.2f}s")

start2 = time.time()
for _ in range(500):
    ConfirmPrime(random.randrange(pow(2, 10)))

print(f"1 took: {time.time()-start2:.2f}s")

start1 = time.time()
for _ in range(500):
    ConfirmPrime2(random.randrange(pow(2, 10)))

print(f"2 took: {time.time()-start1:.2f}s")
