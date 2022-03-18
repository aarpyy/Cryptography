from cryptography318 import siqs, randprime, factor, pollard_rho_factor, pollard_p1

from pathlib import Path


test_root = Path(__file__).parent.absolute()


def test_siqs():
    a = randprime(pow(10, 17), pow(10, 18))
    b = randprime(pow(10, 17), pow(10, 18))
    file = test_root.parent.joinpath("cryptography318").joinpath("primes.txt")
    if file.is_file():
        c = siqs(a * b, fp=str(file), loud=False)
    else:
        c = siqs(a * b, loud=False)
    assert c == a or c == b


def test_rho():
    a = randprime(pow(10, 12), pow(10, 13))
    b = randprime(pow(10, 12), pow(10, 13))
    c = pollard_rho_factor(a * b)
    assert c == a or c == b


def test_p1():
    a = randprime(pow(10, 12), pow(10, 13))
    b = randprime(pow(10, 12), pow(10, 13))
    c = randprime(pow(10, 2), pow(10, 3))
    d = pollard_p1(a * b * c)
    assert d in (a, b, c)


def test_factor():
    a = randprime(pow(10, 12), pow(10, 13))
    b = randprime(pow(10, 12), pow(10, 13))
    c = randprime(pow(10, 2), pow(10, 3))
    factors = factor(a * b * b * c * c)
    assert factors == {a: 1, b: 2, c: 2}


def test_all():
    test_siqs()
    test_rho()
    test_p1()
    test_factor()


if __name__ == "__main__":
    test_all()
