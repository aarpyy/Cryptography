//
// Created by Andrew Carpenter on 8/2/21.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>


unsigned long long power_mod(unsigned long long b, int e, unsigned long long m) {
    unsigned long long result = b;
    for (int i = 1; i < e; i++) {
        result *= b;
        if (result > m) result %= m;
    }
    return result;
}


int miller_rabin_base_a(unsigned long long n, unsigned long long a) {
    if (a >= n) a %= n;
    if (a == 0) return 1;

    unsigned long long q = n - 1;
    int k = 0;

    while (q % 2 == 0) {
        q /= 2;
        k++;
    }
    a = power_mod(a, (int)q, n);
    if (a == 1 || a == n - 1) return 1;

    for (int i = 0; i < k; i++) {
        if (a == n - 1) {
            return 1;
        }
        else if (a == 1) {
            return 0;
        }
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
    // numbers >= 1122004669633 will return -1 automatically here
    if (n < 2047) {
        unsigned long long base[] = {2};
        return miller_rabin_bases(n, base, 1);
    }
    if (n < 1373653) {
        unsigned long long base[] = {2, 3};
        return miller_rabin_bases(n, base, 2);
    }
    if (n < 9080191) {
        unsigned long long base[] = {31, 73};
        return miller_rabin_bases(n, base, 2);
    }
    if (n < 1050535501) {
        unsigned long long base[] = {336781006125, 9639812373923155};
        return miller_rabin_bases(n, base, 2);
    }
    if (n < 3215031751) {
        unsigned long long base[] = {2, 3, 5, 7};
        return miller_rabin_bases(n, base, 4);
    }
    if (n < 4759123141) {
        unsigned long long base[] = {2, 7, 61};
        return miller_rabin_bases(n, base, 3);
    }
    if (n < 1122004669633) {
        unsigned long long base[] = {2, 13, 23, 1662803};
        return miller_rabin_bases(n, base, 4);
    }
    return -1;
}


void test_power_mod() {
    unsigned long long base = 5;
    unsigned long long mods [] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 10; i++) {
        if (base == mods[i]) continue;
        unsigned long long result = power_mod(base, (int) (mods[i] - 1), mods[i]);
        assert (result == 1);
    }
}


int main() {
    time_t start = time(NULL);
    int result = is_prime(539049392);
    if (result) printf("n is prime\n");
    else printf("n is not prime\n");
    printf("This took %ds\n", (int)(time(NULL) - start));
    return 0;
}
