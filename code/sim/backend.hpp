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
#include <complex>

#include "algebra.hpp"
#include "wignerd.hpp"




#define ASTEROIDS_MAX_J 0
#define ASTEROIDS_MAX_K 3




#define MAX_DT 2.0
#define MIN_DT 1.0e-0
#define DT_POWER 1// 5M
#define GM 3.986004418e14
#define RADIUS 6371000

#define INTEGRAL_LIMIT_FRAC 1.0e-3// 5e-6
#define NUM_THREADS 4
    // Torque at closest approach divided by torque at start of sim.
//#define _DEBUG

using uint = unsigned int;
using cdouble = std::complex<double>;

class Asteroid {
public:
    Asteroid(const cdouble* jlms, const cdouble* klms, const double asteroid_radius,
        Vector3 spin, double initial_roll, double impact_parameter, double speed, double distance_ratio_cut);

    int simulate(double cadence, std::vector<double>& resolved_data);

private:
    cdouble jlm(uint l, int m) const;
    cdouble klm(uint l, int m) const;
    void calculate_moi(double initial_roll);
    void set_pos(double impact_parameter);
    Vector3 get_torque();
    void update_position(double dt);
    void update_orientation(double dt);

private:
    const cdouble* jlms;
    const cdouble* klms;
    double asteroid_radius;
    double distance_ratio_cut;

    Vector3 moi, inv_moi;
    double edge_dist; // Limit of the integration region
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    Quaternion orientation;// Global to local
    double time;
    std::array<double, 3> moi_evals;
    std::array<Vector3, 3> moi_evecs;

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
    Vector3 torque;
    Vector3 omegaDot;
    Quaternion d_quat;
    std::array<double, 3> angles;
    Vector3 rot_pos;
    cdouble x_torque, y_torque, z_torque;
    double rot_pos_r, rot_pos_ct, rot_pos_p;
    cdouble nowjlm;
    cdouble prelpmp;
};
