#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <deque>

#include "algebra.hpp"
#include "triangle.hpp"

#define MAX_DT 2.0
#define MIN_DT 1.0e-0
#define DT_POWER 1// 5M

/*#define MAX_DT 3.0
#define MIN_DT 1.0e-0
#define DT_POWER 3*/// 5M
// Sort this out










#define INTEGRAL_LIMIT_FRAC 1.0e-5// 5e-6
#define NUM_THREADS 4
    // Torque at closest approach divided by torque at start of sim.
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define G 6.67408e-11
//#define _DEBUG

using uint = unsigned int;

class Asteroid {
public:
    Asteroid(int L, int n, int m, const std::vector<double>& clms,
        const std::vector<double>& densities, double spin,
        double impact_parameter, double speed, double central_mass);

    int simulate(double cadence, std::vector<double>& resolved_data);
    void draw(std::string filename, int axis) const;

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
    std::deque<Chunk> chunks;
    std::vector<Triangle> triangles;
    double edge_dist; // Limit of the integration region
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    Quaternion orientation;
    double time;
    std::array<double, 3> moi_evals;
    std::array<Vector3, 3> moi_evecs;

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


    // Variables defined globally for computational efficiency
    Vector3 accel;
    Vector3 Omega;
    Matrix3 inv_mat;
    Matrix3 moiGlobal;
    Matrix3 moiGlobalInverse;;
    Vector3 torque;
    Vector3 omegaDot;
    Quaternion d_quat;
};
