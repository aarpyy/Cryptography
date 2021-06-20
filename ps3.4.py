from crypto_functions import *
from statistics import mean

c1 = 4516304024640477153195106845293721117799090654728271883418515743052367794616879147604327199514612019313711520115245892276791096207746360519457897707418484363985930012464303009581430152261880976655550127023091136429974827930905265087972147266160039932299733788576778804954788110000273140780041426445719302674660
c2 = 1114121861501427347654329789397047458702244739466649313965505455544143262051012746257276977369330404278120401764504419446588908438048320056902456792504440989837342289544534632131663016421557598857151249648977316475848388831137081255986402796295140318635044811756402911956259639712229718031668723318603801882210451438162
gcd = GCD(c1, c2)

def Sort(lst):
    for i in range(len(lst)):
        k = i
        max = lst[i]
        for j in range(i, len(lst)):
            if lst[j][1] > max[1]:
                max = lst[j]
                k = j
        lst[k], lst[i] = lst[i], lst[k]
    return lst


def MaximizeChar():
    lst = []
    for d in range(1, 10):
        k = gcd // d
        m1 = NumToString(c1 // k)
        m2 = NumToString(c2 // k)
        p1 = PercentChar(m1)
        p2 = PercentChar(m2)
        pavg = mean((p1, p2))
        lst.append((d, pavg))

    lst = Sort(lst)
    d = lst[0][0]
    k = gcd // d
    m1 = NumToString(c1 // k)
    m2 = NumToString(c2 // k)
    return d, m1, m2, lst[0]





# k = gcd//3
# m1 = NumToString(c1//k)
# m2 = NumToString(c2//k)
# print(m1)
# print(m2)

print(MaximizeChar())
