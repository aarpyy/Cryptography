SAMPLE = "LPTBT DPKIKDLTIB KB NST QVFPL ITKJVWU FMTBB YNIQ K DVOPTI LPKL V" \
         "B LN BKU LPTU DNSXTU K QTKSVSF GML LPTS YINQ HPKL VB ZSNHS NY DKOLKVS " \
         "ZVJJ V DNMWJ SNL BMOONBT PVQ DKOKGWT NY DNSBLIMDLVSF KSU NY LPT QNIT " \
         "KGBLIMBT DIUOLNFIKOPB V QKJT MO QU QVSJ KL NSDT LPKL LPVB HKB NY K BVQ" \
         "OWT BOTDVTB BMDP PNHTXTI KB HNMWJ KOOTKI LN LPT DIMJT VSLTWWTDL NY L" \
         "PT BKVWNI KGBNWMLTWU VSBNWMGWT HVLPNML LPT ZTU"
ORIGINAL = SAMPLE.replace(" ", "")

key = {'A': None, 'B': 'S', 'C': None, 'D': 'C', 'E': None, 'F': 'G', 'G': 'B', 'H': 'W',
       'I': 'R', 'J': 'D', 'K': 'A', 'L': 'T', 'M': 'U', 'N': 'O', 'O': 'P', 'P': 'H',
       'Q': 'M', 'R': None, 'S': 'N', 'T': 'E', 'U': 'Y', 'V': 'I', 'W': 'L', 'X': 'V',
       'Y': 'F', 'Z': 'K'}

dictionary = {}
for i in range(26):
    dictionary[chr(i+65)] = 0

print(dictionary)
for x in ORIGINAL:
    dictionary[x] += 1

count = 0

print(dictionary)
to_add = []
for c in ORIGINAL:
    count+=1
    print(c, end ="")
    to_add.append(c)
    if count == 50:
        print()
        for a in to_add:
            if key[a] is not None:
                print(key[a], end="")
            else:
                print("_", end="")
        to_add = []
        count = 0
        print()
print()
for a in to_add:
    if key[a] is not None:
        print(key[a], end="")
    else:
        print("_", end="")

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

bigrams_dict = {}
for x in range(len(ORIGINAL) - 1):
    s = ORIGINAL[x] + ORIGINAL[x+1]
    if s in bigrams_dict:
        bigrams_dict[s] += 1
    else:
        bigrams_dict[s] = 1

bigrams = []
for c in bigrams_dict:
    if bigrams_dict[c] > 4:
        bigrams.append((c, bigrams_dict[c]))
print(sort(bigrams))

trigrams_dict = {}
for y in range(len(ORIGINAL) - 2):
    s = ORIGINAL[y] + ORIGINAL[y + 1] + ORIGINAL[y +2]
    if s in trigrams_dict:
        trigrams_dict[s] += 1
    else:
        trigrams_dict[s] = 1

trigrams = []
for d in trigrams_dict:
    if trigrams_dict[d] > 1:
        trigrams.append((d, trigrams_dict[d]))
print(sort(trigrams))












