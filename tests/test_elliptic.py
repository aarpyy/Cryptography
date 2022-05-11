from cryptography318.factor import ecm_mont
from cryptography318.factor.elliptic import ecm_weierstrass, ecm_mont_basic
from cryptography318.prime import randprime, isprime
from timeit import timeit


def test_ecm_mont():
    a = randprime(pow(10, 25), pow(10, 26))
    b = randprime(pow(10, 25), pow(10, 26))
    c = randprime(pow(10, 4), pow(10, 5))
    d = ecm_mont(a * b, B1=960, B2=5700)
    assert d in (a, b, c)


def test_ecm_weierstrass():
    a = randprime(pow(10, 7), pow(10, 8))
    b = randprime(pow(10, 7), pow(10, 8))
    f = ecm_weierstrass(a * b)
    if f is not None:
        assert f in (a, b)


def test_compare_ecm():
    count = 0
    a = 20
    while 1:
        N = randprime(pow(2, a - 1), pow(2, a)) + 4
        if isprime(N):
            continue

        time_mont = timeit(lambda: ecm_mont_basic(N), number=50)
        time_we = timeit(lambda: ecm_mont(N), number=50)
        print(time_mont / time_we)
        if time_we < time_mont:
            if count == 10:
                print(N, a)
                break
            count += 1
        else:
            count = 0
        a += 1
        print(count, a)


if __name__ == "__main__":
    test_ecm_weierstrass()
