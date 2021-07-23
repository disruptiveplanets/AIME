#include "wignerd.hpp"

// z-y-z Euler angle convention

// https://en.wikipedia.org/wiki/Wigner_D-matrix
uint fact(uint i) {
    if (i == 0) return 1;
    return i * fact(i-1);
}

uint choose(uint a, uint b) {
    uint num = 1;
    for (uint i = 0; i < b; i++) {
        num *= a - i;
    }
    return num / fact(b);
}

double gen_choose(double a, uint b) {
    double num = 1;
    for (uint i = 0; i < b; i++) {
        num *= a - i;
    }
    return num / fact(b);
}

DMatGen::DMatGen(double alpha, double beta, double gamma)
    : alpha(0, alpha), gamma(0, gamma) {
    sn = sin(beta/2);
    cs = cos(beta/2);
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
        sum += (parity(mp-m+s) * pow(cs, 2 * j + m - mp - 2 * s)
        * pow(sn, mp - m + 2 * s)) /
        (fact(j + m - s) * fact(s) * fact(mp - m + s) * fact(j - mp - s));
    }
    return pre * sum;
}

cdouble slm_c(uint l, int m, double r, double costheta, double phi) {
    return parity(m) * fact(l - m) / pow(r, l+1) * ylm_c(l, m, costheta, phi);
}

cdouble ylm_c(uint l, int m, double costheta, double phi) {
    double sum = 0;
    for (uint k = abs(m); k <= l; k++) {
        sum += fact(k) / (double)fact(k - abs(m)) * pow(costheta, k - abs(m))
            * choose(l, k) * gen_choose((l + k - 1) / 2.0, l);
    }
    cdouble out = pow(2, l) * pow(1 - costheta * costheta, abs(m) / 2.0) * sum;

    if (m >= 0) {
        return out * std::complex(cos(m * phi), sin(m * phi));
    }
    else {
        return (double)parity(m) * fact(l + m) / fact(l - m) * out
            * std::complex(cos(m * phi), sin(m * phi));
    }
}
