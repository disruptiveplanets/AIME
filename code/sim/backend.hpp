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

#define MIN_DT 10.0 // 11_420
#define MAX_DT 20.0
#define POSITION_DT 1.0

#define INTEGRAL_LIMIT_FRAC 1.0e-3// 5e-6

//#define TEXT_DEBUG
#define SEGFAULT_DEBUG
//#define FIRST_ORDER

#define SPACE_MPP (1)
#define SPACE_MP (SPACE_MPP * (2 * (ASTEROIDS_MAX_K + 1) + 1))
#define SPACE_LP (SPACE_MP * (2 * (ASTEROIDS_MAX_K + 1) + 1))
#define SPACE_M (SPACE_LP * (ASTEROIDS_MAX_K + 1))
#define SPACE_L (SPACE_M * (2 * (ASTEROIDS_MAX_J + 1) + 1))
#define COEFF_LENGTH (SPACE_L * (ASTEROIDS_MAX_J + 1))

class Asteroid {
public:
    Asteroid(const cdouble* jlms, const cdouble* klms, const double asteroid_radius,
        Vector3 spin, double initial_roll, double perigee, double speed,
        double central_mu, double central_radius, double distance_ratio_cut, bool enforce_drc);

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
    void compute_coefficient(int l, int m, int lp, int mp, int mpp, 
        cdouble& coeff_x, cdouble& coeff_y, cdouble& coeff_z) const;
    void make_coefficients();
    void get_coefficients(int l, int m, int lp, int mp, int mpp, cdouble mul,
        cdouble& torque_x, cdouble& torque_y, cdouble& torque_z) const;

private:
    const cdouble* jlms;
    const cdouble* klms;
    cdouble klmcs[(ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1)];
    double asteroid_radius;
    double distance_ratio_cut;
    bool enforce_drc;

    cdouble coeffs_x[COEFF_LENGTH];
    cdouble coeffs_y[COEFF_LENGTH];
    cdouble coeffs_z[COEFF_LENGTH];

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
};
