from crypto_functions import *

# start1 = time.time()
# p = RandomPrime(pow(2, 1000))
# print(p)
# print(f"K = 40 took: {time.time() - start1:.2f}s")


# print(MillerRabinPrimality(53))
# print(Jacobi(19, 53))
p = RandomPrime(pow(2, 1000))
print(p)
print(BailliePSW_Primality(p))

