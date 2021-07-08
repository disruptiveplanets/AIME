#pragma once

#include <complex>
#include <math.h>
#include <map>
#include <tuple>
#include <iostream>

#define max(a, b) ((a) < (b) ? (b) : (a))
#define min(a, b) ((a) > (b) ? (b) : (a))
#define sign(a) ((a) % 2 == 0 ? (1) : (-1))

using uint = unsigned int;
using cdouble = std::complex<double>;

uint fact(uint i);

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
