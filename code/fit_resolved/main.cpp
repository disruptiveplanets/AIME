#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <complex>
#include <vector>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"
#include "../sim/triangle.hpp"

using cdouble = std::complex<double>;

std::vector<double> simulate(double cadence, const std::vector<double> jlms_raw,
    const std::vector<double> mlms_raw, double spin, double impact_parameter,
    double speed) {

    int maxjl = sqrt(jlms_raw.size()) - 1;
    int maxml = sqrt(mlms_raw.size()+4) - 1;

    std::vector<cdouble> jlms;
    auto jlm_raw_iter = jlms_raw.begin();
    for (int l = 0; l <= maxjl; l++){
        if (l != 0) {
            for (int m = 1; m <= l; m++) {
                jlms.push_back({*jlm_raw_iter++, *jlm_raw_iter++});
            }
        }
        jlms.push_back({*jlm_raw_iter++, 0});
    }

    std::vector<cdouble> mlms;
    auto mlm_raw_iter = mlms_raw.begin();
    for (int l = 2; l <= maxml; l++){
        for (int m = 1; m <= l; m++) {
            mlms.push_back({*mlm_raw_iter++, *mlm_raw_iter++});
        }
        mlms.push_back({*mlm_raw_iter++, 0});
    }

    Asteroid asteroid(jlms, mlms, spin, impact_parameter, speed);

    // Run asteroid
    std::vector<double> resolved_data;
    asteroid.simulate(cadence, resolved_data);


    return resolved_data;// Unfortunately, copying the data is necessary.
}

PYBIND11_MODULE(asteroids, m) {
    m.doc() = "Asteroid simulation wrapper to be used by the python module"
        "emcee"; // optional module docstring

    m.def("simulate", &simulate, "Simulates a flyby");
}
