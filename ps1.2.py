s = "UZQSOVUOHXMOPVGPOZPEVSGZWSZOPFPESXUDBMETSXAIZ" \
    "VUEPHZHMDZSHZOWSFPAPPDTSVPQUZWYMXUZUHSX" \
    "EPYEPOPDZSZUFPOMBZWPFUPZHMDJUDTMOHMQ"
a = "itwasdisclosedyesterdaythatseveralinformalbut" \
    "directcontactshavebeenmadewithpolitical" \
    "representativesofthevietconginmoscow"

def confirm(cipher_text):
    # length is number of chars in keyword, r is number of chars in last row
    # k is min number of chars in each col, rows is max number of chars in each col
    solved = False
    length = 4
    while not solved:
        if length == 14:
            break
        r = 26 % length
        k = 26 // length
        rows = k + min(r, 1)

        # adds letters to columns to form cipher-block
        columns = []
        l, h, b = 0, 0, 0
        for i in range(r):
            columns.append([])
            for j in range(l, rows + l):
                columns[i].append(cipher_text[j])
                h = j
            l = h + 1
            b = i

        for m in range(b + 1, length - r + b + 1):
            columns.append([])
            for j in range(l, k + l):
                columns[m].append(cipher_text[j])
                h = j
            l = h + 1

        row = []
        for i in range(length):
            row.append(columns[i][1])

        solved = True
        for i in range(len(row) - 1):
            one, two = ord(row[i]), ord(row[i+1])
            if 65 <= one <= 90 and 65 <= two < 90 and one > two:
                solved = False

        length += 1
        if solved:
            # iterates through columns, printing out cipher-block
            for i in range(len(columns[0])):
                for j in range(len(columns)):
                    if i < len(columns[j]):
                        print(columns[j][i], end=" ")
                print()



# creates dictionary matching each letter in plaintext to its corresponding letter in cipher-text
d = {}
for i in range(len(a)):
    if a[i].upper() not in d:
        d[a[i].upper()] = s[i]

# manual matches made by looking at missing letters in cipher-block
d["J"] = "C"
d["Q"] = "N"
d["K"] = "L"
d["X"] = "K"
d["Z"] = "R"

# creates list of uppercase alphabet
alph = [chr(x + 65) for x in range(26)]

# prints alphabet for comparison of plaintext to cipher-text
for c in alph:
    print(c, end=" ")
print()

# creates list for comparison of cipher-text to plaintext, prints out
cipher_text = []
for c in alph:
    if c in d:
        cipher_text.append(d[c])
        print(d[c], end=" ")
    else:
        cipher_text.append("_")
        print("_", end=" ")
print("\n\n")




confirm(cipher_text)
