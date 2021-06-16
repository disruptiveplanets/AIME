#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include "algebra.hpp"
#include "triangle.hpp"

#define ONE_SECOND_TORQUE 1e-8
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
    double edge_dist; // Limit of the integration region
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    double mean_density;// Used to generate rough parameters of the asteroid
    double radius;// Used to generate rough parameters of the asteroid
    double closest_approach;
};
