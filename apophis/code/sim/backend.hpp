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

#include "../../../code/sim/algebra.hpp"
#include "../../../code/sim/wignerd.hpp"

#define MIN_DT 10.0 // 11_420
#define MAX_DT 20.0
#define POSITION_DT 1.0

#define INTEGRAL_LIMIT_FRAC_START 1.0e-3 // 1.0e-12
#define INTEGRAL_LIMIT_FRAC_END 1.0e-5

//#define TEXT_DEBUG
#define SEGFAULT_DEBUG
//#define FIRST_ORDER

#define SPACE_MPP (1)
#define SPACE_MP (SPACE_MPP * (2 * (ASTEROIDS_MAX_K + 1) + 1))
#define SPACE_LP (SPACE_MP * (2 * (ASTEROIDS_MAX_K + 1) + 1))
#define SPACE_M (SPACE_LP * (ASTEROIDS_MAX_K + 1))
#define SPACE_L (SPACE_M * (2 * (ASTEROIDS_MAX_J + 1) + 1))
#define COEFF_LENGTH (SPACE_L * (ASTEROIDS_MAX_J + 1))

// Before modifying these, remember to rerun the angle program
#define P_PSI (264.178 * 3600)
#define P_PHI (27.38547 * 3600)
// The rest are just Apophis / coordinate numbers
#define L_LON (275 * PI / 180)
#define L_LAT (-85 * PI / 180)
#define NORTH_LON (-PI / 2)
#define NORTH_LAT (66.5607083333 * PI / 180)
#define ARG_PERI (32.56762207107416 * PI / 180)
#define INCLINATION (162.8561543134850 * PI / 180)
#define LON_ASC_NODE (152.0019490904037 * PI / 180)
#define OBLATENESS 1.082627e-3
#define PERIGEE 3.801148649129103e07 // In meters
#define ECCENTRICITY 4.253901873212790
#define MU 3.986004418e14 // SI
#define CENTRAL_RADIUS 6371071.027 // SI


class Asteroid {
public:
    Asteroid(const cdouble* klms, const double asteroid_radius, double initial_roll, double initial_precess,
        double distance_ratio_cut, bool enforce_drc, double velocity_mul);

    int simulate(double cadence, std::vector<double>& resolved_data);

private:
    cdouble jlm(int l, int m) const;
    cdouble klm(int l, int m) const;
    cdouble klmc(int l, int m) const;
    void calculate_moi();
    void calculate_jlms();
    Vector3 get_torque();
    void get_derivatives(Vector3 position, Vector3 spin, Quaternion quat, Vector3& dspin, Quaternion& dquat);
    void calculate_poses();
    bool extract_pos(double time, Vector3& position, Vector3& velocity);
    Vector3 extract_spin(Vector3 angles, Vector3 momenta);
    void initialize_rotation(Quaternion& quat, Vector3& spin) const;
    void compute_coefficient(int l, int m, int lp, int mp, int mpp, 
        cdouble& coeff_x, cdouble& coeff_y, cdouble& coeff_z) const;
    void make_coefficients();
    void get_coefficients(int l, int m, int lp, int mp, int mpp, cdouble mul,
        cdouble& torque_x, cdouble& torque_y, cdouble& torque_z) const;

private:
    cdouble* jlms;
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
    double central_radius;
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;

    // Orientation factors
    double l_lon;
    double l_lat;
    double initial_roll;
    double initial_precess;

    // Orbital factors
    double excess_vel;
    double expire_time;
    double start_time;
    double velocity_mul;

    //double max_quat_mag;
};


double get_tilt(double roll);
Quaternion get_ecliptic_quat();
