def from_base(lst, base):
    multiplier = 1
    acc = 0
    for e in lst:
        acc += e * multiplier
        multiplier *= base
    return acc
