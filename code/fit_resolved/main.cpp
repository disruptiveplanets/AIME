#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"

std::vector<double> simulate(double cadence, const std::vector<double> jlms_raw,
    const std::vector<double> parameters,  double radius, double spinx, double spiny, double spinz,
    double perigee, double speed, double central_mu, double central_radius,
    double distance_ratio_cut=0, double enforce_drc=false,
    double velocity_mul=1) {

    std::vector<cdouble> halfjlms;
    auto jlm_raw_iter = jlms_raw.begin();
    for (int l = 0; l <= ASTEROIDS_MAX_J; l++){
        if (l != 0) {
            for (int m = 1; m <= l; m++) {
                halfjlms.push_back({*jlm_raw_iter++, *jlm_raw_iter++});
            }
        }
        halfjlms.push_back({*jlm_raw_iter++, 0});
    }

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
    cdouble jlms[(ASTEROIDS_MAX_J+1) * (ASTEROIDS_MAX_J+1)];
    cdouble klms[(ASTEROIDS_MAX_K+1) * (ASTEROIDS_MAX_K+1)];
    uint j = 0;
    for (uint l = 0; l <= ASTEROIDS_MAX_J; l++) {
        for (int m = -l; m <= (int)l; m++) {
            assert(j < (ASTEROIDS_MAX_J+1) * (ASTEROIDS_MAX_J+1));
            if (m < 0) {
                jlms[j] = halfjlms[l * (l + 1) / 2 + l - abs(m)].conj()
                    * (double)parity(m);
            }
            else {
                jlms[j] = halfjlms[l * (l + 1) / 2 + l - abs(m)];
            }
            j++;
        }
    }
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


    Asteroid asteroid(jlms, klms, radius, Vector3({spinx, spiny, spinz}),
        initial_roll, perigee, speed, central_mu, central_radius,
        distance_ratio_cut, enforce_drc, velocity_mul);

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
