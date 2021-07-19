#include "wignerd.hpp"

// z-y-z Euler angle convention

// https://en.wikipedia.org/wiki/Wigner_D-matrix
uint fact(uint i) {
    if (i == 0) return 1;
    return i * fact(i-1);
}

DMatGen::DMatGen(double alpha, double beta, double gamma)
    : alpha(0, alpha), gamma(0, gamma) {
    cs = cos(beta/2);
    sn = sin(beta/2);
}

cdouble DMatGen::operator()(uint l, int mp, int m) {
    auto d_pos = small_d_state.find({l, mp, m});
    if (d_pos == small_d_state.end()) {
        small_d_state.insert({{l, mp, m}, wigner_small_d(l, mp, m)});
        d_pos = small_d_state.find({l, mp, m});
    }
    return std::exp(-(double)mp * alpha) * std::exp(-(double)m * gamma)
        * d_pos->second;
}

double DMatGen::wigner_small_d(uint j, int mp, int m) {
    double pre = sqrt(fact(j+mp) * fact(j-mp) * fact(j+m) * fact(j-m));
    double sum = 0;
    for(uint s = max(0, m-mp); s <= min(j+m, j-mp); s++) {
        sum += (sign(mp-m+s) * pow(cs, 2 * j + m - mp - 2 * s)
        * pow(sn, mp - m + 2 * s)) /
        (fact(j + m - s) * fact(s) * fact(mp - m + s) * fact(j - mp - s));
    }
    return pre * sum;
}
