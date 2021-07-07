#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"
#include "../sim/triangle.hpp"

std::vector<double> simulate(double cadence, int L, int n, int m,
    const std::vector<double> clms,
    const std::vector<double> densities, double spin, double impact_parameter,
    double speed, double central_mass) {


    Asteroid asteroid(L, n, m, clms, densities, spin,
        impact_parameter, speed, central_mass);

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
