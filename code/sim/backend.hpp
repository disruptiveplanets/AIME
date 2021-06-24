#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include "algebra.hpp"
#include "triangle.hpp"

#define DT_SLOPE 5e-6// 5e-6
#define ONE_SECOND_TORQUE 1e-8// 1e-8
#define INTEGRAL_LIMIT_FRAC 1.0e-5// 5e-6
    // Torque at closest approach divided by torque at start of sim.
#define max(a, b) ((a) > (b) ? (a) : (b))
#define G 6.67408e-11
//#define _DEBUG

using uint = unsigned int;

class Asteroid {
public:
    Asteroid(std::string filename);

    int simulate(double cadence, std::ofstream&& resolved, std::ofstream&& unresolved);
    void draw(std::string filename, Vector3 axis) const;

private:
    void make_chunks();
    void calculate_moi();
    void calculate_mass();
    Vector3 get_com() const;
    void recenter();
    void set_pos(double impact_parameter);

    Vector3 get_torque();
    void update_position(double dt);
    void update_orientation(double dt);
    Vector3 get_rot_ang_mom();
    Matrix3 global_to_inertial() const;
    Matrix3 inertial_to_global() const;

private:
    uint L; // Max degree of harmonics
    uint n; // Number of cube subdivisions on outermost shell
    uint m; // Number of shells
    double mass;
    Matrix3 moi;
    Matrix3 moiInverse;
    std::vector<Chunk> chunks;
    std::vector<Triangle> triangles;
    double edge_dist; // Limit of the integration region
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    Quaternion orientation;
    double time;

    // Shape features
    double mean_density;// Used to generate rough parameters of the asteroid
    double radius;// Used to generate rough parameters of the asteroid

    // Orbital factors
    double energy;
    Vector3 ang_mom;
    double closest_approach;
    double excess_vel;
    double impact_parameter;

    #ifdef _DEBUG
    double max_quat_mag;
    #endif
};
