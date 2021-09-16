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


def test_fraction():
    _FRACTION_WITH_SR_FORMAT = re.compile(r"""
        \A\s*                      # optional whitespace at the start, then
        (?P<sign>[-+]?)            # an optional sign, then             # lookahead for digit or .digit
        (?P<value>\d*)               # numerator (possibly empty)
        (?:sqrt\((?P<radicand>\d+))?
        \)
        (?:/(?P<denom>\d+))?    
    """, re.VERBOSE | re.IGNORECASE)
    test = "sqrt(2)"
    print(_FRACTION_WITH_SR_FORMAT.match(test))


if __name__ == '__main__':
    test_fraction()
