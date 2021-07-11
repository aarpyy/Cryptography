from crypto_functions import *
# from twilio.rest import Client
#
# account_sid = 'AC236a9c8d47fcdacbfda0511b52170843'
# auth_token = 'c1d15c40890e12452819a6e4217ee64b'
# client = Client(account_sid, auth_token)
# number = "+15039577675"


while True:
    q = RandomPrime(pow(2, 500))
    p = (2 * q) + 1
    if MillerRabinPrimality(p, 80) and MillerRabinPrimality(pow(q, -1, p), 80):
        print(q, p)
        break


