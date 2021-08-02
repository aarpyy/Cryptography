//
// Created by Andrew Carpenter on 8/2/21.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


unsigned long long power_mod(unsigned long long b, long double e, unsigned long long m) {
    long double x = (long double) b;
    long double result = powl(x, e);
    if (result < m) return result;
    return result % m
}


int miller_rabin_base_a(unsigned long long n, unsigned long long a) {
    if (a >= n) a %= n;
    if (a == 0) return 1;

    long double q = n - 1;
    int k = 0;

    while (q % 2 == 0) {
        q >>= 1;
        k++;
    }

    a = power_mod(a, q, n);
    if (a == 1 || a == n - 1) return 1;

    for (int i = 0; i < k; i++) {
        if (a == n - 1) return 1;
        else if (a == 1) return 0;
        a = power_mod(a, 2, n);
    }
    return 0;
}


int miller_rabin_bases(unsigned long long n, unsigned long long bases [], int size) {
    for (int i = 0; i < size; i++) {
        if (miller_rabin_base_a(n, bases[i]) == 0) return 0;
    }
    return 1;
}


int is_prime(unsigned long long n) {
    if (n < 2047)
        return miller_rabin_bases(n, {2}, 1)
    if (n < 1373653)
        return MillerRabin_bases(n, {2, 3}, 2)
    if (n < 9080191)
        return MillerRabin_bases(n, {31, 73}, 2)
    if (n < 1050535501)
        return MillerRabin_bases(n, {336781006125, 9639812373923155}, 2)
    if (n < 3215031751)
        return MillerRabin_bases(n, {2, 3, 5, 7}, 4)
    if (n < 4759123141)
        return MillerRabin_bases(n, {2, 7, 61}, 3)
    if (n < 1122004669633)
        return MillerRabin_bases(n, {2, 13, 23, 1662803}, 4)
    if (n < 55245642489451)
        return MillerRabin_bases(n, {2, 141889084524735, 1199124725622454117, 11096072698276303650}, 4)
    if (n < 7999252175582851)
        return MillerRabin_bases(n, {2, 4130806001517, 149795463772692060, 186635894390467037, 3967304179347715805}, 5)
    if (n < 18446744073709551616)
        return MillerRabin_bases(n, {2, 325, 9375, 28178, 450775, 9780504, 1795265022}, 7)
    if (n < 318665857834031151167461)
        return MillerRabin_bases(n, {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}, 12)
    if (n < 3317044064679887385961981)
        return MillerRabin_bases(n, {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41}, 13)
    return -1;
}
