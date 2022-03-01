from cryptography318.crypto_functions import from_base


def test_from_base():
    frm_bse = lambda lst, base: sum(map(lambda i, n: n * pow(base, i), range(len(lst)), lst))

    array = [1, 1, 1]
    assert frm_bse(array, 2) == from_base(array, 2)
