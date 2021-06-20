from crypto_functions import *

p = 7954763
g = 4
a = 1234092
A = pow(g, a, p)

k = 6234043
ma = StringToNum("Hi!")
c1a = pow(g, k, p)
c2a = (ma * pow(A, k, p)) % p

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
mb = (pow(c1bi, a, p) * c2b) % p

print(f"\nProblem 2. b)\n"
      f"-------------\n"
      f"Bobs' message:\n"
      f"m: {NumToString(mb)}\n")