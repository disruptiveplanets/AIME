#include "wignerd.hpp"

// z-y-z Euler angle convention

// https://en.wikipedia.org/wiki/Wigner_D-matrix



int fact(int i) {
    if (i == 0) return 1;
    return i * fact(i-1);
}


int choose(int a, int b) {
    int num = 1;
    for (int i = 0; i < b; i++) {
        num *= a - i;
    }
    return num / fact(b);
}

double gen_choose(double a, int b) {
    double num = 1;
    for (int i = 0; i < b; i++) {
        num *= a - i;
    }
    return num / fact(b);
}

DMatGen::DMatGen(double alpha, double beta, double gamma)
    : alpha(0, alpha), gamma(0, gamma) {
    sn = sin(beta/2);
    cs = cos(beta/2);
}

cdouble DMatGen::operator()(int l, int mp, int m) const {
    return my_exp(-(double)mp * alpha) * my_exp(-(double)m * gamma)
        * wigner_small_d(l, mp, m);
}

double DMatGen::wigner_small_d(int j, int mp, int m) const {
    double pre = sqrt(fact(j+mp) * fact(j-mp) * fact(j+m) * fact(j-m));
    double sum = 0;
    for(int s = max_me(0, m-mp); s <= min_me(j+m, j-mp); s++) {
        sum += ((double)parity(mp-m+s) * pow(cs, 2 * j + m - mp - 2 * s)
        * pow(sn, mp - m + 2 * s)) /
        (fact(j + m - s) * fact(s) * fact(mp - m + s) * fact(j - mp - s));
    }
    return pre * sum;
}

cdouble slm_c(int l, int m, double r, double costheta, double phi) {
    // About 0.33 s
    return (double)parity(m) * fact(l - m) / pow(r, l+1) * ylm_c(l, m, costheta, phi);
}

cdouble ylm_c(int l, int m, double costheta, double phi) {
    double sum = 0;
    //for (int k = abs(m); k <= l; k++) {
    double power = 1;
    for (int u = 0; u <= l - abs(m); u++) {
        sum += power / fact(l - u - abs(m)) * gen_choose((l + u + abs(m) - 1) / 2.0, l);
        power *= costheta / (u + 1);
    }
    cdouble out = pow(2, l) * pow(1 - costheta * costheta, abs(m) / 2.0) * sum * fact(l);

    if (m >= 0) {
        return out * cdouble(cos(m * phi), sin(m * phi));
    }
    else {
        return (double)parity(m) * fact(l + m) / fact(l - m) * out
            * cdouble(cos(m * phi), sin(m * phi));
    }
}

std::ostream& operator<<(std::ostream& os, const cdouble& c) {
    if (c.i >= 0) {
        os << c.r << " + " << c.i << " i";
    }
    else {
        os << c.r << " - " << abs(c.i) << " i";
    }
    return os;
}
cdouble operator*(double d, cdouble c) { return cdouble(d * c.r, d * c.i); }
cdouble operator/(double d, cdouble c) { return cdouble(d * c.r / c.norm2(), -d * c.i / c.norm2()); }
cdouble my_exp(cdouble c) {return exp(c.r) * cdouble(cos(c.i), sin(c.i)); }
cdouble my_sqrt(double d) {
    if (d < 0) {return cdouble (0, sqrt(-d));}
    return cdouble(sqrt(d), 0);
}
cdouble my_pow(cdouble c, double p) {
    return my_exp(p * cdouble(0.5 * log(c.norm2()), atan2(c.i, c.r)));
}
cdouble operator+(double d, cdouble c) {
    return cdouble(d + c.r, c.i);
}
cdouble operator-(double d, cdouble c) {
    return cdouble(d - c.r, -c.i);
}
