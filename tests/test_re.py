import re

s = re.compile(r"""
    \A\s*
    (?P<coeff>\d*)
    (?:\.(?P<coeff_dec>\d*))?
    \s*?
    \*?
    \s*?
    (?:sqrt\(
        (?P<val>\d+))?
    (?:\.(?P<decimal>\d*))?
    \)
""", re.VERBOSE | re.IGNORECASE)


_RATIONAL_FORMAT = re.compile(r"""
    \A\s*                      # optional whitespace at the start, then
    (?P<sign>[-+]?)            # an optional sign, then
    (?=\d|\.\d)                # lookahead for digit or .digit
    (?P<num>\d*)               # numerator (possibly empty)
    (?:                        # followed by
       (?:/(?P<denom>\d+))?    # an optional denominator
    |                          # or
       (?:\.(?P<decimal>\d*))? # an optional fractional part
       (?:E(?P<exp>[-+]?\d+))? # and optional exponent
    )
    \s*\Z                      # and optional whitespace to finish
""", re.VERBOSE | re.IGNORECASE)


def print_re(regex):
    groups = ('sign', 'num', ('denom', 'decimal', 'exp'))

    to_print = ""
    for g in groups:
        if isinstance(g, str):
            v = regex.group(g)
            if v:
                to_print += f'{g}: {v}; '
        else:
            d = regex.group(g[0])
            if d is not None:
                to_print += f'{g[0]}: {d}; '
            else:
                for e in g[1:]:
                    if (d := regex.group(e)) is not None:
                        to_print += f'{e}: {d}; '
    print(to_print)


# x = _RATIONAL_FORMAT.match("1.2E5")
x = s.match("2 * sqrt(1.1)")
print(x)
print(x, x.group('val'), x.group('decimal'), x.groups(), x.group(), x.group() == 'sqrt(12)')
vals = tuple(map(lambda g: '0' if g is None else g, x.groups()))
coeff = '.'.join(vals[:2])
value = '.'.join(vals[2:])
print(coeff, value)
# print(re.search('(?<=sqrt)\d+', 'sqrt12'))
