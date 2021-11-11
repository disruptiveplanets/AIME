#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <deque>
#include <assert.h>
#include <stdexcept>
#include <chrono>

#include "algebra.hpp"
#include "wignerd.hpp"

#define MAX_DT 2.0
#define MIN_DT 1.0e-0
#define DT_POWER 1// 5M
#define GM 3.986004418e14
#define RADIUS 6371000
#define POSITION_DT 1.0

#define INTEGRAL_LIMIT_FRAC 1.0e-3// 5e-6
#define NUM_THREADS 4
    // Torque at closest approach divided by torque at start of sim.
#define _DEBUG

using uint = unsigned int;

class Asteroid {
public:
    Asteroid(const cdouble* jlms, const cdouble* klms, const double asteroid_radius,
        Vector3 spin, double initial_roll, double perigee, double speed, double distance_ratio_cut);

    int simulate(double cadence, std::vector<double>& resolved_data);

private:
    cdouble jlm(uint l, int m) const;
    cdouble klm(uint l, int m) const;
    cdouble klmc(uint l, int m) const;
    void calculate_moi(double initial_roll);
    void calculate_poses();
    bool extract_pos(double time, Vector3& position, Vector3& velocity);
    Vector3 extract_spin(Vector3 angles, Vector3 momenta);
    void set_pos(double speed);
    Vector3 get_torque();
    void get_derivatives(Vector3 position, const Vector3 angles, const Vector3 momenta, Vector3& dangles, Vector3& dmomenta);
    Matrix3 local_to_global(Vector3 angles);

private:
    const cdouble* jlms;
    const cdouble* klms;
    cdouble klmcs[(ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1)];
    double asteroid_radius;
    double distance_ratio_cut;

    Matrix3 moi, inv_moi;
    double edge_dist; // Limit of the integration region
    double mu;
    double time;
    std::array<double, 3> moi_evals;
    std::array<Vector3, 3> moi_evecs;

    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;

    // Orbital factors
    double pericenter_pos;
    double excess_vel;
    Vector3 initial_spin;
    double expire_time;
    double initial_roll;

    #ifdef _DEBUG
    double max_quat_mag;
    #endif

    int max_time_index;
};
