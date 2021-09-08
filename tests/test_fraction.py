from cryptography318.fraction import *


def test_constructor():
    # for i in range(1, 100):
    #     x = Fraction(1, sqrt(i))
    #     try:
    #         assert abs(x - 1 / sqrt(i)) < pow(10, -9)
    #     except AssertionError:
    #         print(i, x, x.eval(), 1/sqrt(i))
    x = Fraction(6, 25)
    y = Fraction(8, 15)
    print(x)
    print(y)
    print(x + y)
    print(x + Fraction(4, 25))


if __name__ == '__main__':
    test_constructor()
