SAMPLE = "LPTBT DPKIKDLTIB KB NST QVFPL ITKJVWU FMTBB YNIQ K DVOPTI LPKL V" \
         "B LN BKU LPTU DNSXTU K QTKSVSF GML LPTS YINQ HPKL VB ZSNHS NY DKOLKVS " \
         "ZVJJ V DNMWJ SNL BMOONBT PVQ DKOKGWT NY DNSBLIMDLVSF KSU NY LPT QNIT " \
         "KGBLIMBT DIUOLNFIKOPB V QKJT MO QU QVSJ KL NSDT LPKL LPVB HKB NY K BVQ" \
         "OWT BOTDVTB BMDP PNHTXTI KB HNMWJ KOOTKI LN LPT DIMJT VSLTWWTDL NY L" \
         "PT BKVWNI KGBNWMLTWU VSBNWMGWT HVLPNML LPT ZTU"
ORIGINAL = SAMPLE.replace(" ", "")
KEY = {'A': None, 'B': 'S', 'C': None, 'D': 'C', 'E': None, 'F': 'G', 'G': 'B', 'H': 'W', 'I': 'R', 'J': 'D', 'K': 'A',
       'L': 'T', 'M': 'U', 'N': 'O', 'O': 'P', 'P': 'H', 'Q': 'M', 'R': None, 'S': 'N', 'T': 'E', 'U': 'Y', 'V': 'I',
       'W': 'L', 'X': 'V', 'Y': 'F', 'Z': 'K'}


def print_row(row, t):
    for char in row:
        if char == ' ':
            print(" ", end='')
        elif t[char] is None:
            print("_", end='')
        else:
            print(t[char], end='')


def sort(lst):
    for i in range(len(lst)):
        l = i
        m = (0, 0)
        for j in range(i, len(lst)):
            if m[1] < lst[j][1]:
                m = lst[j]
                l = j
        lst[i], lst[l] = lst[l], lst[i]
    return lst


# makes dictionary of each letter and frequency in string 's'
singles = {}
for c in ORIGINAL:
    if c in singles:
        singles[c] += 1
    else:
        singles[c] = 1

# makes dictionary of each set of bigrams and their frequency
doubles = {}
for i in range(len(ORIGINAL) - 1):
    b = ORIGINAL[i] + ORIGINAL[i + 1]
    if b in doubles:
        doubles[b] += 1
    else:
        doubles[b] = 1

# makes dictionary of each set of trigrams and their frequency
triples = {}
for i in range(len(ORIGINAL) - 2):
    b = ORIGINAL[i] + ORIGINAL[i + 1] + ORIGINAL[i + 2]
    if b in triples:
        triples[b] += 1
    else:
        triples[b] = 1

single_freq = []
for k in singles:
    single_freq.append((k, singles[k]))

double_freq = []
for k in doubles:
    if doubles[k] > 3:
        double_freq.append((k, doubles[k]))

triple_freq = []
for k in triples:
    if triples[k] > 2:
        triple_freq.append((k, triples[k]))

print(single_freq := sort(single_freq), "\n")
print(double_freq := sort(double_freq), "\n")
print(triple_freq := sort(triple_freq), "\n")

count = 0
to_edit = []
for i in range(len(SAMPLE)):
    if count == 50:
        print(SAMPLE[i])
        to_edit.append(SAMPLE[i])

        print_row(to_edit, KEY)
        print('\n')

        to_edit = []
        count = 0
    else:
        print(SAMPLE[i], end='')
        to_edit.append(SAMPLE[i])
    count += 1

print()
print_row(to_edit, KEY)
print("\n")
x = [chr(a + 65) for a in range(26)]
for c in KEY:
    if KEY[c] in x:
        x.remove(KEY[c])
print(f"Missing letters: {x}")
