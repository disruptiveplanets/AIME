#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "../sim/backend.hpp"
#include "../../../code/sim/algebra.hpp"

std::vector<double> simulate(double cadence, const std::vector<double> parameters,  double radius, double initial_precess,
    double distance_ratio_cut=0, double enforce_drc=false, double velocity_mul=1) {

    std::vector<cdouble> halfklms;
    auto klm_raw_iter = parameters.begin();
    const double initial_roll = *klm_raw_iter++;
    halfklms.push_back({1, 0}); // K00
    halfklms.push_back({0, 0}); // K11
    halfklms.push_back({0, 0}); // K10
    halfklms.push_back({*klm_raw_iter++, 0}); // K22
    halfklms.push_back({0, 0}); // K20
    halfklms.push_back({*klm_raw_iter++, 0}); // K20
    for (int l = 3; l <= ASTEROIDS_MAX_K; l++){
        for (int m = 1; m <= l; m++) {
            halfklms.push_back({*klm_raw_iter++, *klm_raw_iter++});
        }
        halfklms.push_back({*klm_raw_iter++, 0});
    }


    // Make full jlms
    cdouble klms[(ASTEROIDS_MAX_K+1) * (ASTEROIDS_MAX_K+1)];
    uint k = 0;
    for (uint l = 0; l <= ASTEROIDS_MAX_K; l++) {
        for (int m = -l; m <= (int)l; m++) {
            assert(k < (ASTEROIDS_MAX_K+1) * (ASTEROIDS_MAX_K+1));
            if (m < 0) {
                klms[k] = halfklms[l * (l + 1) / 2 + l - abs(m)].conj()
                    * (double)parity(m);
            }
            else {
                klms[k] = halfklms[l * (l + 1) / 2 + l - abs(m)];
            }
            k++;
        }
    }

    Asteroid asteroid(klms, radius, initial_roll, initial_precess, distance_ratio_cut, enforce_drc, velocity_mul);

    // Run asteroid
    std::vector<double> resolved_data;
    asteroid.simulate(cadence, resolved_data);


    return resolved_data;// Unfortunately, copying the data is necessary.
}

PYBIND11_MODULE(NAME, m) {
    m.doc() = "Asteroid simulation wrapper to be used by the python module"
        "emcee"; // optional module docstring

    m.def("simulate", &simulate, "Simulates a flyby");
}
