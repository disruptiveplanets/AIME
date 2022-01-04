#include "backend.hpp"

#define my_min(a, b) ((a) < (b) ? (a) : (b))

const int max_k = (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1);
const int max_j = (ASTEROIDS_MAX_J + 1) * (ASTEROIDS_MAX_J + 1);

Asteroid::Asteroid(const cdouble* jlms, const cdouble* klms, double asteroid_radius,
    Vector3 spin, double initial_roll, double perigee, double speed, double central_mu,
    double central_radius, double distance_ratio_cut) :
    jlms(jlms), klms(klms), asteroid_radius(asteroid_radius), distance_ratio_cut(distance_ratio_cut),
    mu(central_mu), central_radius(central_radius), pericenter_pos(perigee), excess_vel(speed), initial_spin(spin),
    initial_roll(initial_roll) {

    for (int i = 0; i < (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1); i++) {
        klmcs[i] = klms[i].conj();
    }

    calculate_moi(initial_roll);
    calculate_poses();

    /*#ifdef TEXT_DEBUG
    std::cout<< "Klms: " << std::endl;
    for (int l = 0; l <= ASTEROIDS_MAX_K; l++){
        for (int m = -l; m<= l; m++) {
            std::cout << l << ' ' << m << ": " << klm(l, m) << std::endl;
        }
    }
    std::cout<< "Jlms: " << std::endl;
    for (int l = 0; l <= ASTEROIDS_MAX_J; l++){
        for (int m = -l; m<= l; m++) {
            std::cout << l << ' ' << m << ": " << jlm(l, m) << std::endl;
        }
    }

    std::cout << "MOI (local): " << moi << std::endl;
    std::cout << std::endl;
    #endif*/

    make_coefficients();
}

cdouble Asteroid::jlm(int l, int m) const {
    #ifdef SEGFAULT_DEBUG
    if (l * l + l + m < 0 || l * l + l + m > max_j) {
        std::cout << "Jlm array exceeded" << std::endl;
        throw std::runtime_error("Jlm array exceeded.");
    }
    #endif
    return jlms[l * l + l + m];
}

cdouble Asteroid::klm(int l, int m) const {
    if (abs(m) > l) {return 0;}
    #ifdef SEGFAULT_DEBUG
    if (l * l + l + m < 0 || l * l + l + m > max_k) {
        std::cout << "Klm array exceeded" << std::endl;
        throw std::runtime_error("Klm array exceeded.");
    }
    #endif
    return klms[l * l + l + m];
}

cdouble Asteroid::klmc(int l, int m) const {
    if (abs(m) > l) {return 0;}
    #ifdef SEGFAULT_DEBUG
    if (l * l + l + m < 0 || l * l + l + m > max_k) {
        std::cout << "Klmc array exceeded" << std::endl;
        throw std::runtime_error("Klm array exceeded.");
    }
    #endif
    return klmcs[l * l + l + m];
}

void Asteroid::compute_coefficient(int l, int m, int lp, int mp, int mpp, 
    cdouble& coeff_x, cdouble& coeff_y, cdouble& coeff_z) const {

    cdouble prelm = pow(central_radius, l) * jlm(l, m) * mu / 2.0;
    cdouble prelpmp = prelm * parity(lp) * pow(asteroid_radius, lp - 2);
    cdouble mppfactor = prelpmp * sqrt(fact(lp - mpp) * fact(lp+mpp) / (double)(fact(lp-mp) * fact(lp+mp)));

    coeff_x = mppfactor * (cdouble(0, lp - mpp + 1) * klm(lp, mpp - 1)
        + cdouble(0, lp + mpp + 1) * klm(lp, mpp + 1));

    coeff_y = mppfactor * (-(lp - mpp + 1) * klm(lp, mpp - 1)
        + (lp + mpp + 1) * klm(lp, mpp + 1));

    coeff_z = mppfactor * cdouble(0, 2) * mpp * klm(lp, mpp);
}

void Asteroid::make_coefficients() {
    for (int l = 0; l <= ASTEROIDS_MAX_J; l++) {
        for (int m = -l; m <= l; m++) {
            for (int lp = 2; lp <= ASTEROIDS_MAX_K; lp++) {
                for (int mp = -lp; mp <= lp; mp++) {
                    for (int mpp = -lp; mpp <= (int)lp; mpp++) {
                        const int index = l * SPACE_L + (m + ASTEROIDS_MAX_J) * SPACE_M + lp * SPACE_LP +
                            (mp + ASTEROIDS_MAX_K) * SPACE_MP + (mpp + ASTEROIDS_MAX_K) * SPACE_MPP;
                        compute_coefficient(l, m, lp, mp, mpp, 
                            coeffs_x[index],
                            coeffs_y[index],
                            coeffs_z[index]);
                    }
                }
            }
        }
    }
}

void Asteroid::get_coefficients(int l, int m, int lp, int mp, int mpp, cdouble mul,
        cdouble& torque_x, cdouble& torque_y, cdouble& torque_z) const {
    const int index = l * SPACE_L + (m + ASTEROIDS_MAX_J) * SPACE_M + lp * SPACE_LP +
        (mp + ASTEROIDS_MAX_K) * SPACE_MP + (mpp + ASTEROIDS_MAX_K) * SPACE_MPP;
    #ifdef SEGFAULT_DEBUG
    if (index >= COEFF_LENGTH) {
        std::cout << "Coeff array exceeded" << std::endl;
        throw std::runtime_error("Coeff array exceeded.");
    }
    #endif
    torque_x += mul * coeffs_x[index];
    torque_y += mul * coeffs_y[index];
    torque_z += mul * coeffs_z[index];
}

void Asteroid::calculate_poses() {
    #ifdef TEXT_DEBUG
    auto start = std::chrono::high_resolution_clock::now();
    #endif

    double pericenter_vel = sqrt(excess_vel * excess_vel + 2 * mu / pericenter_pos);
    const double cutoff_dist_squared = pericenter_pos * pericenter_pos
        * pow(INTEGRAL_LIMIT_FRAC, -2/3.0);// Go a little over for safety
    #ifdef TEXT_DEBUG
    std::cout << "cutoff distance: " << sqrt(cutoff_dist_squared) << std::endl;
    std::cout << "periapsis: " << pericenter_pos << std::endl;
    #endif
    Vector3 velocity = Vector3({0, pericenter_vel, 0});
    Vector3 position = Vector3({pericenter_pos, 0, 0});

    double time;
    for (time = 0; position.mag2() <= cutoff_dist_squared; time += POSITION_DT) {
        Vector3 accel = -mu / position.mag2() * position / position.mag();
        velocity += POSITION_DT * accel;
        position += POSITION_DT * velocity;
        positions.push_back(position);
        velocities.push_back(velocity);
    }
    
    // Add one extra position data point for the sake of interpolation
    Vector3 accel = -mu / position.mag2() * position / position.mag();
    velocity += POSITION_DT * accel;
    position += POSITION_DT * velocity;
    positions.push_back(position);
    velocities.push_back(velocity);


    expire_time = time;

    #ifdef TEXT_DEBUG
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Position time " << duration.count()/ 1000000.0 << std::endl;
    std::cout << "Expire time (h): " << expire_time / 3600 << std::endl;
    #endif
}

bool Asteroid::extract_pos(double time, Vector3& position, Vector3& velocity) {
    bool flip = false;
    if (time < 0) {
        flip = true;
        time *= -1;
    }
    int time_index = time / POSITION_DT;
    if (time > expire_time) {
        return false;
    }
    double mixer = (time - time_index * POSITION_DT) / POSITION_DT;

    #ifdef SEGFAULT_DEBUG
    if (time_index + 1 >= (int)positions.size()) {
        std::cout << "Position array length exceeded: " << time_index << ' ' << positions.size() << "; " << time << ' ' << expire_time << std::endl;
    }
    #endif

    position = positions[time_index + 1] * mixer + positions[time_index] * (1 - mixer);
    velocity = velocities[time_index + 1] * mixer + velocities[time_index] * (1 - mixer);
    if (flip) {
        position = Vector3({position[0], -position[1], 0});
        velocity = Vector3({-velocity[0], velocity[1], 0});
    }
    return true;
}

void Asteroid::calculate_moi(double initial_roll) {
    double Ixx = (2/3.0 * klm(2,0) - 2.0 * klm(2, -2) - 2.0 * klm(2, 2)).r + 2/3.0;
    double Iyy = (2/3.0 * klm(2,0) + 2.0 * klm(2, -2) + 2.0 * klm(2, 2)).r + 2/3.0;
    double Izz = (- 4 / 3.0 * klm(2, 0)).r + 2/3.0;
    // These mois are really moi per radius^2 per M.

    if (Izz <= 0 || Iyy <= 0 || Ixx <= 0) {
        throw std::runtime_error("Moment of inertia values were negative."
        );
    }

    if (abs(Izz) <= abs(Ixx) || abs(Izz) <= abs(Iyy)) {
        throw std::runtime_error("Moment of inertia was not maximized along z.");
    }

    if (Izz >= Ixx + Iyy || Iyy >= Ixx + Izz || Ixx >= Izz + Iyy) {
        throw std::runtime_error("Triangle inequality violated");
    }

    moi = Vector3({Ixx, Iyy, Izz});
    avg_moi = (abs(moi[0]) + abs(moi[1]) + abs(moi[2])) / 3;
    inv_moi = Vector3({1.0 / Ixx, 1.0 / Iyy, 1.0 / Izz});
}

void Asteroid::get_derivatives(Vector3 position, Vector3 spin, Quaternion quat, Vector3& dspin, Quaternion& dquat) {
    
    #ifndef FIRST_ORDER
    Vector3 angles = quat.euler_angles();
    DMatGen dgen = DMatGen(angles[0], angles[1], angles[2]);
    double pos_r = position.mag();
    double pos_ct = position[2] / position.mag();
    double pos_p = atan2(position[1], position[0]);
    cdouble x_torque = 0;
    cdouble y_torque = 0;
    cdouble z_torque = 0;
    // Up until now, takes about 0.1 s, L=2

    /*for (int l = 0; l <= ASTEROIDS_MAX_J; l++) {
        for (int m = -l; m <= l; m++) {
            cdouble prelm = pow(central_radius, l) * jlm(l, m);
            for (int lp = 2; lp <= ASTEROIDS_MAX_K; lp++) {
                for (int mp = -lp; mp <= lp; mp++) {
                    cdouble prelpmp = prelm * parity(lp) * pow(asteroid_radius, lp - 2)
                        * slm_c(l+lp, m+mp, pos_r, pos_ct, pos_p).conj();

                    for (int mpp = -lp; mpp <= (int)lp; mpp++) {
                        cdouble mppfactor = dgen(lp, mp, mpp).conj() * prelpmp * sqrt(fact(lp - mpp) * fact(lp+mpp) / (double)(fact(lp-mp) * fact(lp+mp)));

                        x_torque += mppfactor * (cdouble(0, lp - mpp + 1) * klm(lp, mpp - 1)
                            + cdouble(0, lp + mpp + 1) * klm(lp, mpp + 1));

                        y_torque += mppfactor * (-(lp - mpp + 1) * klm(lp, mpp - 1)
                            + (lp + mpp + 1) * klm(lp, mpp + 1));

                        z_torque += mppfactor * cdouble(0, 2) * mpp * klm(lp, mpp);
                    }
                }
            }
        }
    }
    torque = Vector3({x_torque.r, y_torque.r, z_torque.r}) * mu / 2.0;*/

    for (int l = 0; l <= ASTEROIDS_MAX_J; l++) {
        for (int m = -l; m <= l; m++) {
            for (int lp = 2; lp <= ASTEROIDS_MAX_K; lp++) {
                for (int mp = -lp; mp <= lp; mp++) {
                    cdouble prelpmp = slm_c(l+lp, m+mp, pos_r, pos_ct, pos_p).conj();
                    for (int mpp = -lp; mpp <= (int)lp; mpp++) {
                        cdouble mppfactor = dgen(lp, mp, mpp).conj() * prelpmp;
                        get_coefficients(l, m, lp, mp, mpp, mppfactor, x_torque, y_torque, z_torque);
                    }
                }
            }
        }
    }
    torque = Vector3({x_torque.r, y_torque.r, z_torque.r});
    
    #else

    double beta = acos(-position[2] / position.mag());
    double gamma = atan2(-position[1], position[0]);
    Matrix3 r = quat.matrix().transpose();
    Quaternion q = Quaternion(cos(beta / 2), 0, sin(beta / 2), 0)
        * Quaternion(cos(gamma / 2), 0, 0, sin(gamma / 2));
    Matrix3 b = q.matrix();
    //std::cout << b.transpose() * position / position.mag() << std::endl;
    Matrix3 diag_moi = Matrix3({moi[0], 0, 0, 0, moi[1], 0, 0, 0, moi[2]});

    //std::cout << r.transpose() * diag_moi * r << std::endl;

    Matrix3 moi_z_at_planet = b.transpose() * r * diag_moi * r.transpose() * b;
    Vector3 torque2_z_at_planet = 3 * mu / pow(position.mag(), 3) * Vector3({-moi_z_at_planet(1, 2), moi_z_at_planet(0, 2), 0});
    Vector3 torque = r.transpose() * b * torque2_z_at_planet;

    if (torque.is_nan()) {
        std::cout << "Torque was nan" << std::endl;
        throw std::runtime_error("Torque was nan");
    }
    // This torque is really torque per radius^2 per M.
    #endif

    dspin = Vector3({
        (torque[0] + (moi[1] - moi[2]) * spin[1] * spin[2]) * inv_moi[0],
        (torque[1] + (moi[2] - moi[0]) * spin[2] * spin[0]) * inv_moi[1],
        (torque[2] + (moi[0] - moi[1]) * spin[0] * spin[1]) * inv_moi[2],
    });


    dquat = 0.5 * Quaternion(0, spin[0], spin[1], spin[2]) * quat;
}


int Asteroid::simulate(double cadence, std::vector<double>& resolved_data) {
    int frames = 0;
    int cadence_index = -expire_time / cadence-1;
    double dt = MIN_DT;
    Quaternion dquat1, dquat2, dquat3, dquat4;
    Vector3 dspin1, dspin2, dspin3, dspin4, position, velocity;


    double alpha = initial_roll;
    double beta = acos(initial_spin[2] / initial_spin.mag());
    double gamma = atan2(initial_spin[1], -initial_spin[0]);

    Quaternion quat = Quaternion(cos(alpha / 2), 0, 0, sin(alpha / 2))
        * Quaternion(cos(beta / 2), 0, sin(beta / 2), 0)
        * Quaternion(cos(gamma / 2), 0, 0, sin(gamma / 2));
    Vector3 spin = Vector3::z() * initial_spin.mag();

    const double expire_distance = distance_ratio_cut * distance_ratio_cut * pericenter_pos * pericenter_pos;

    double time;
    for (time = MIN_DT-expire_time; time+MIN_DT < expire_time; time += dt) {
        if (!extract_pos(time, position, velocity)) {
            throw std::runtime_error("Simulation ran out of bounds");
        }
        if (distance_ratio_cut > 0 && position.mag2() > expire_distance &&
            Vector3::dot(position, velocity) > 0) {
            break;
        }

        get_derivatives(position, spin, quat, dspin1, dquat1);
        const double scale = (pow(position.mag2() / (pericenter_pos * pericenter_pos), 3.0 / 2) - 1.0) * INTEGRAL_LIMIT_FRAC;
        dt = MIN_DT + (MAX_DT - MIN_DT) * scale;

        extract_pos(time+dt/2, position, velocity);
        get_derivatives(position, spin + dt / 2 * dspin1, quat + dt / 2 * dquat1, dspin2, dquat2);
        get_derivatives(position, spin + dt / 2 * dspin2, quat + dt / 2 * dquat2, dspin3, dquat3);
        extract_pos(time+dt, position, velocity);
        get_derivatives(position, spin + dt * dspin3, quat + dt * dquat3, dspin4, dquat4);

        spin += dt / 6 * (dspin1 + 2 * dspin2 + 2 * dspin3 + dspin4);
        quat += dt / 6 * (dquat1 + 2 * dquat2 + 2 * dquat3 + dquat4);

        quat /= quat.mag();

        while (int(time / cadence) > cadence_index) {
            Vector3 global_spin = quat.rotate(spin);
            resolved_data.push_back(global_spin[0]);
            resolved_data.push_back(global_spin[1]);
            resolved_data.push_back(global_spin[2]);
            cadence_index++; //cadence_index = int(time / cadence)
        }
        frames++;
    }

    #ifdef TEXT_DEBUG
    std::cout << "Simulation took " << time << " seconds." << std::endl;
    #endif

    return frames;
}
