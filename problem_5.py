cipher = "SAHVPBJWU**XTDMY*EOZIFQ*G*"
for i in range(26):
    for j in (range(26 - i)):
        print(cipher[j], end="")
    print()
    for k in (range(i, 26)):
        print(cipher[k], end = "")
    print()
