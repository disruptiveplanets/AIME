#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include "algebra.hpp"
#include "triangle.hpp"

#define DELTA_T_MIN 1
#define INTEGRAL_LIMIT_FRAC 1.0e-5 
    // Torque at closest approach divided by torque at start of sim.
#define G 6.67408e-11
#define _DEBUG

using uint = unsigned int;

class Asteroid {
public:
    Asteroid(std::string filename);

    int simulate(std::ofstream&& resolved, std::ofstream&& unresolved);

private:
    void make_chunks();
    void calculate_moi();
    void calculate_mass();
    Vector3 get_com();
    void recenter();
    void set_pos(double impact_parameter);

    Vector3 get_torque();
    void update_position(double dt);
    void update_orientation(double dt);

private:
    uint L; // Max degree of harmonics
    uint n; // Number of cube subdivisions on outermost shell
    uint m; // Number of shells
    double mass;
    Matrix3 moi;
    std::vector<Chunk> chunks;
    std::vector<Triangle> triangles;
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    double closest_approach;
};