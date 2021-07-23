#pragma once

#include <complex>
#include <math.h>
#include <map>
#include <tuple>
#include <iostream>

#define max_me(a, b) ((a) < (b) ? (b) : (a))
#define min_me(a, b) ((a) > (b) ? (b) : (a))
#define parity(a) ((a) % 2 == 0 ? (1) : (-1))

using uint = unsigned int;
using cdouble = std::complex<double>;

uint fact(uint i);
uint choose(uint a, uint b);
double gen_choose(double a, uint b);

class DMatGen {
public:
    DMatGen(double alpha, double beta, double gamma);
    cdouble operator()(uint l, int m, int mp);

    double wigner_small_d(uint j, int mp, int m);

private:
    std::map<std::tuple<int, int, int>, double> small_d_state;
    cdouble premul;
    double cs;// cos(beta / 2)
    double sn;// sin(beta / 2)
    cdouble alpha;
    cdouble gamma;
};

cdouble slm_c(uint l, int m, double r, double costheta, double phi);
cdouble ylm_c(uint l, int m, double costheta, double phi);
