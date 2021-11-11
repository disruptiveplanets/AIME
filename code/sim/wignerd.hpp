#pragma once

#include <math.h>
#include <map>
#include <tuple>
#include <iostream>

#define max_me(a, b) ((a) < (b) ? (b) : (a))
#define min_me(a, b) ((a) > (b) ? (b) : (a))
#define parity(a) ((a) % 2 == 0 ? (1) : (-1))

using uint = unsigned int;



struct cdouble {
    cdouble(double r, double i) : r(r), i(i) {}
    cdouble(double r) : r(r), i(0) {}
    cdouble() {}

    double norm2() const {return r * r + i * i;}
    cdouble conj() const {return cdouble(r, -i);}
    cdouble operator+(double d) const { return cdouble(r + d, i);}
    cdouble operator-(double d) const { return cdouble(r - d, i);}
    cdouble operator+(cdouble c) const { return cdouble(r + c.r, i + c.i);}
    cdouble operator-(cdouble c) const { return cdouble(r - c.r, i - c.i);}
    cdouble operator*(double d) const { return cdouble(r * d, i * d);}
    cdouble operator/(double d) const { return cdouble(r / d, i / d);}
    void operator+=(cdouble c) { r += c.r; i += c.i; }
    void operator-=(cdouble c) { r -= c.r; i -= c.i;}
    void operator*=(double d) { r *= d; i *= d;}
    void operator/=(double d) { r /= d; i /= d;}

    cdouble operator*(cdouble c) const { return cdouble(r * c.r - i * c.i, i * c.r + r * c.i);}
    void operator*=(cdouble c) { double newr = r * c.r - i * c.i; double newi = i * c.r + r * c.i; r = newr; i = newi;}

    double r;
    double i;
};

uint fact(uint i);
uint choose(uint a, uint b);
double gen_choose(double a, uint b);

class DMatGen {
public:
    DMatGen(double alpha, double beta, double gamma);
    cdouble operator()(uint l, int m, int mp) const;
    cdouble db(uint l, int m, int mp) const;

private:
    cdouble premul;
    double cs;// cos(beta / 2)
    double sn;// sin(beta / 2)
    cdouble alpha;
    cdouble gamma;
};

cdouble slm_c(uint l, int m, double r, double costheta, double phi);
cdouble ylm_c(uint l, int m, double costheta, double phi);

std::ostream& operator<<(std::ostream& os, const cdouble& c);
cdouble operator*(double d, cdouble c);
cdouble operator/(double d, cdouble c);
cdouble operator+(double d, cdouble c);
cdouble operator-(double d, cdouble c);
cdouble my_exp(cdouble c);
cdouble my_sqrt(double d);
cdouble my_pow(cdouble c, double p);
