s = "XPMHKIDHKPNGRHXGPVFIMHUVJWXPZVMHVROFWFWTUPPDGJCKH"
alph = [chr(a+65) for a in range(26)]
dict = {}
count = 0
for c in alph:
    dict[c] = count
    count += 1

def encrypt(c):
    return chr((((19 * (ord(c) - 65)) + 9) % 26) + 65)


def bad_encrypt(c):
    return chr((((6 * (ord(c) - 65)) + 9) % 26) + 65)


# pairs = {}
# for c in alph:
#     n = bad_encrypt(c)
#     if n in pairs:
#         pairs[n].append(c)
#     else:
#         pairs[n] = [c]
#
# print(pairs)
#
# def func(c):
#     return ((6 * (ord(c) - 65)) + 9) % 26
#
# print(func('N'))
# print(ord('N') - 65)
# print(ord('A') - 65)

key = {}
for c in alph:
    key[encrypt(c)] = c

for k in key:
    print(f"{k} --> {key[k]}")

for c in s:
    print(key[c], end="")
print()


k = int(input("Enter: "))