from crypto_functions import *
import time

def Compare(func1, func2, iter):
    inputs = []
    for _ in range(iter):
        inputs.append(RandomPrime(pow(2, 1000)))
    start1 = time.time()
    f1 = list(map(func1, inputs))
    print(f"{func1} took: {time.time()-start1:.2f}s")

    start2 = time.time()
    f2 = list(map(func2, inputs))
    print(f"{func2} took: {time.time()-start2:.2f}s")

    for ind in range(len(f1)):
        i, j = f1[ind], f2[ind]
        if not i or not j:
            print(f'Psuedo-prime produced: {inputs[ind]}')
            print("Function 1 failed") if not j else print("Function 2 failed")
        ind += 1




Compare(MillerRabinPrimality, lambda x: MillerRabinPrimality(x, 100), 500)
