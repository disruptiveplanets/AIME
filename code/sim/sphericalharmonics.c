#include "sphericalharmonics.h"


int factorial(int i) {
    if (i == 0) {
        return 1;
    }
    return i * factorial(i-1);
}
int choose(int a, int b) {
    return factorial(a) / factorial(b) / factorial(a - b);
}
double assoc_lagrange(uint l, int m, double x) {
    double sum = 0;
    for(int k = m; k <= l; k++) {
        sum += factorial(k) / factorial(k - m) * pow(x, k-m) *
        choose(l, k) * choose((l + k - 1) / 2, l);
    }
    return sum * pow(-1, m) * pow(2, l) * pow(1 - x * x, m/2.0);
}

double real_spherical_harmonic(uint l, int m, double theta, double phi) {
    double pre = pow(-1, m) * sqrt((2 * l + 1) / (2 * PI) *
                      factorial(l - abs(m)) / factorial(l + abs(m))) *
        assoc_lagrange(l, m, cos(theta));
    if (m < 0) {
        return pre * sin(abs(m) * phi);
    }
    if (m > 0) {
        return pre * cos(abs(m) * phi);
    }
    return pre / sqrt(2);
}
