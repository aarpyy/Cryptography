from cryptography318.crypto_functions import *
from cryptography318.prime import RandomPrime, IsPrime
import pytest


def test_safe_curve_and_point():
    for _ in range(50):
        E, P = EllipticCurve.safe_curve_and_point(RandomPrime(pow(10, 3), pow(10, 4)))
        assert isinstance(E, EllipticCurve) and isinstance(P, Elliptic)


def test_elliptic_add():
    for _ in range(500):
        p = RandomPrime(pow(10, 2), pow(10, 3))
        E, P = EllipticCurve.safe_curve_and_point(p)
        Q = P + P
        R = Q - P
        assert R == P and R + P == Q


def test_elliptic_mult():
    for _ in range(50):
        p = RandomPrime(pow(10, 2), pow(10, 3))
        E, P = EllipticCurve.safe_curve_and_point(p)
        R = P * 2
        order = 2
        while R != P:
            order += 1
            R += P

        assert P * order == P

        # if you can factor order, try multiplying part way by order, then completely to see if result is original
        if not IsPrime(order):
            fact = FactorInt(order)
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


@pytest.mark.skip
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
    e = StringToElliptic(E, s)
    assert EllipticToString(e) == s


if __name__ == '__main__':
    test_safe_curve_and_point()
    test_elliptic_add()
    test_elliptic_mult()
