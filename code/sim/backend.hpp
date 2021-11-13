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

#define MIN_DT 1.0e-0
#define MAX_DT 1000
#define DT_POWER 1// 5M
#define POSITION_DT 1.0

#define INTEGRAL_LIMIT_FRAC 1.0e-3// 5e-6
#define NUM_THREADS 4
    // Torque at closest approach divided by torque at start of sim.
//#define TEXT_DEBUG

class Asteroid {
public:
    Asteroid(const cdouble* jlms, const cdouble* klms, const double asteroid_radius,
        Vector3 spin, double initial_roll, double perigee, double speed,
        double central_mu, double central_radius, double distance_ratio_cut);

    int simulate(double cadence, std::vector<double>& resolved_data);

private:
    cdouble jlm(int l, int m) const;
    cdouble klm(int l, int m) const;
    cdouble klmc(int l, int m) const;
    void calculate_moi(double initial_roll);
    Vector3 get_torque();
    void get_derivatives(Vector3 position, Vector3 spin, Quaternion quat, Vector3& dspin, Quaternion& dquat);
    void calculate_poses();
    bool extract_pos(double time, Vector3& position, Vector3& velocity);
    Vector3 extract_spin(Vector3 angles, Vector3 momenta);

private:
    const cdouble* jlms;
    const cdouble* klms;
    cdouble klmcs[(ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1)];
    double asteroid_radius;
    double distance_ratio_cut;

    Vector3 moi, inv_moi;
    Vector3 torque;
    double avg_moi;
    double mu;
    double central_radius;
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;

    // Orbital factors
    double pericenter_pos;
    double excess_vel;
    Vector3 initial_spin;
    double expire_time;
    double initial_roll;

    //double max_quat_mag;
    /// Throws a segfault when this is not present.
};
