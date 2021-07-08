#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <deque>
#include <assert.h>
#include <complex>

#include "algebra.hpp"
#include "triangle.hpp"
#include "wignerd.hpp"

#define MAX_DT 2.0
#define MIN_DT 1.0e-0
#define DT_POWER 1// 5M

#define INTEGRAL_LIMIT_FRAC 1.0e-5// 5e-6
#define NUM_THREADS 4
    // Torque at closest approach divided by torque at start of sim.
#define G 6.67408e-11
#define _DEBUG

using uint = unsigned int;
using cdouble = std::complex<double>;

class Asteroid {
public:
    Asteroid(const std::vector<cdouble> jlms, const std::vector<cdouble> mlms,
        double spin, double impact_parameter, double speed);

    int simulate(double cadence, std::vector<double>& resolved_data);

private:
    cdouble mlm(uint l, int m);
    cdouble nowmlm(uint l, int m);
    void set_nowmlm(uint l, int m, cdouble val);
    cdouble jlm(uint l, int m);
    void calculate_moi();
    void set_pos(double impact_parameter);
    void update_mlms();

    Vector3 get_torque();
    void update_position(double dt);
    void update_orientation(double dt);
    Vector3 get_rot_ang_mom();
    Matrix3 global_to_inertial() const;
    Matrix3 inertial_to_global() const;

private:
    uint maxml;
    uint maxjl;

    const std::vector<cdouble> mlms;
    std::vector<cdouble> nowmlms;
    const std::vector<cdouble> jlms;

    Matrix3 moi;
    Matrix3 moiInverse;
    double edge_dist; // Limit of the integration region
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    Quaternion orientation;
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
    Matrix3 inv_mat;
    Matrix3 moiGlobal;
    Matrix3 moiGlobalInverse;;
    Vector3 torque;
    Vector3 omegaDot;
    Quaternion d_quat;
};
