from cryptography318.crypto_functions import *
from cryptography318.prime import randprime, isprime
from cryptography318.elliptic import *
from cryptography318.factor import *
import pytest
from timeit import timeit
from sympy.ntheory.ecm import Point as syPoint
from sympy import factorint
from sympy.ntheory.factor_ import _factorint_small


def test_point():
    for _ in range(50):
        p = randprime(pow(10, 4), pow(10, 5))
        E = EllipticCurve.safe_curve(p)
        P = E.point(x=randrange(2, p))
        assert E.on_curve(P)


def test_safe_curve_and_point():
    for _ in range(50):
        E, P = EllipticCurve.safe_curve_and_point(randprime(pow(10, 3), pow(10, 4)))
        assert isinstance(E, EllipticCurve) and isinstance(P, Elliptic) and E.on_curve(P)


def test_elliptic_add():
    for _ in range(500):
        p = randprime(pow(10, 2), pow(10, 3))
        E, P = EllipticCurve.safe_curve_and_point(p)
        Q = P + P
        R = Q - P
        assert R == P and R + P == Q


def test_elliptic_mult():
    for _ in range(50):
        p = randprime(pow(10, 2), pow(10, 3))
        E, P = EllipticCurve.safe_curve_and_point(p)
        R = P * 2
        order = 2
        while R != P:
            order += 1
            R += P

        assert P * order == P

        # if you can factor order, try multiplying part way by order, then completely to see if result is original
        if not isprime(order):
            fact = factor(order)
            factors = []
            for f in fact:
                factors.append(pow(f, fact[f]))

            x = factors[0]
            y = prod(factors[1:])

            Q = P * x
            R = Q * y

            assert R == P


def test_elliptic_bsgs():
    E = EllipticCurve(14, 19, 3623)
    P = E.point(6, 730)
    Q = P * 947
    assert 947 == elliptic_bsgs(P, Q)


def test_string():
    E = EllipticCurve(0, 1, 180039337167793501897175493739004038762395261681316852988742367480124359904807289130328068843833873066421702742460640101217376600809338199468636198555098044628998299482972323944983045884234913875508525458735536722067910377840070662914565334732578373261984514856548773711583945448134428877109037377605127434409832936443635849865584735606971650245470719854253303747830625905599682675415954316175111606488355254401541327809795822857218931391498811975792908442781893479296187087742281713894913694970210301919109047958435358258637147)
    s = "Some say the world will end in fire,\n" \
        "Some say in ice.\n" \
        "From what I've tasted of desire,\n" \
        "I hold with those who favor fire. \n" \
        "But if it had to perish twice\n" \
        "I think I know enough of hate\n" \
        "To say that for destruction ice\n" \
        "Is also great\n" \
        "And would suffice."
    e = string_to_elliptic(E, s)
    assert elliptic_to_string(e) == s

    for _ in range(50):
        p = randprime(pow(10, 3), pow(10, 4))
        E = EllipticCurve.safe_curve(p)

        # get some form of string non-reliant on Elliptic curves
        n = randrange(2, p)
        s = int_to_string(n)

        p = string_to_elliptic(E, s)
        assert elliptic_to_string(p) == s


def test_all():
    test_safe_curve_and_point()
    test_elliptic_add()
    test_elliptic_mult()
    test_string()
    test_point()


def test_efficiency():
    p = 851820817
    curve, point = Weierstrass.safe_curve_and_point(p)
    E = EllipticCurve(curve.a, curve.b, p)
    P = E.point(point.x, point.y)
    time_1 = timeit(lambda: P * 2, number=100000)

    P = syPoint(11, 16, 7, 851820817)

    time_2 = timeit(lambda: P.double(), number=100000)

    time_3 = timeit(lambda: point * 2, number=100000)

    print(f"time mine: {time_1:.2f}s")
    print(f"time sy: {time_2:.2f}s")
    print(f"time new point: {time_3:.2f}s")
    print(f"sy = 100%; mine = {(time_1 / time_2) * 100:.2f}%; new point = {(time_3 / time_2) * 100:.2f}%")


def test_weierstrass():
    e = Weierstrass(2, 3, 13)
    print(e.__slots__, e._p)


def test_compare_mont():
    E = Montgomery(7, 29)
    p1 = MontgomeryPoint(11, 16, E)
    p2 = MontgomeryPoint(13, 10, E)

    syp1 = syPoint(11, 16, 7, 29)
    syp2 = syPoint(13, 10, 7, 29)

    it = 1000000

    time_add = timeit(lambda: p2.add(p1, p1), number=it)
    time_syadd = timeit(lambda: syp2.add(syp1, syp1), number=it)

    time_double = timeit(lambda: p1.double(), number=it)
    time_sydouble = timeit(lambda: syp1.double(), number=it)

    time_ladder = timeit(lambda: p1.ladder(3), number=it)
    time_syladder = timeit(lambda: syp1.mont_ladder(3), number=it)

    print(f"time to add: {(time_add / time_syadd) * 100:.2f}%")
    print(f"time to double: {(time_double / time_sydouble) * 100:.2f}%")
    print(f"time to multiply: {(time_ladder / time_syladder) * 100:.2f}%")


def test_compare_all(test_accurate=False):
    if test_accurate:
        for _ in range(10000):
            p = randprime(pow(2, 10), pow(2, 14))
            N = p + 4
            if isprime(N):
                continue
            km, kw = ecm_mont(N), ecm_weierstrass(N)
            if km is None or kw is None:
                print(N, km, kw)
                continue
            assert not N % km and not N % kw

    time_mont, time_mont2 = 0, 1
    power = 32
    count = 0
    while 1:
        N = 7
        while isprime(N):
            N = randprime(pow(2, power), pow(2, power + 2)) + 4
        print(N)

        time_mont = timeit(lambda: ecm_mont(N), number=100)
        print(f"time 1: {time_mont}")
        time_mont2 = timeit(lambda: ecm_mont_large(N), number=100)
        print(f"time 2: {time_mont2}")

        print(f"mont: {(time_mont / time_mont2) * 100:.2f}%")
        if time_mont > time_mont2:
            if count > 10:
                break
            count += 1
        else:
            power += 1

    N = randprime(pow(2, 15), pow(2, 16)) + 4
    if not isprime(N):
        print(N)
        time_factor = timeit(lambda: factor_small({}, N, 2 ** 15), number=10000)
        time_sy = timeit(lambda: _factorint_small({}, N, 2 ** 15, 300), number=10000)
        print(f"time comparison: {(time_factor / time_sy) * 100:.2f}%")


def inc_mul(p, n):
    b = p
    m = 1
    while not n % b:
        b *= p
        m += 1
    return m - 1


if __name__ == '__main__':
    # test_compare_all()
    time_m = timeit(lambda: multiplicity(3, 54), number=pow(10, 7))
    time_i = timeit(lambda: inc_mul(3, 54), number=pow(10, 7))
    print(f"efficiency: {(time_m / time_i) * 100:2f}%")
