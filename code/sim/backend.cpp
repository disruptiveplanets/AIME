#include "backend.hpp"

const int max_k = (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1);
const int max_j = (ASTEROIDS_MAX_J + 1) * (ASTEROIDS_MAX_J + 1);

Asteroid::Asteroid(const cdouble* jlms, const cdouble* klms, double asteroid_radius,
    Vector3 spin, double initial_roll, double perigee, double speed,
    double distance_ratio_cut) :
    jlms(jlms), klms(klms), asteroid_radius(asteroid_radius), distance_ratio_cut(distance_ratio_cut),
    velocity(Vector3({0, 0, speed})), spin(spin), perigee(perigee) {

    klmcs = new cdouble[(ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1)];
    for (int i = 0; i < (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1); i++) {
        klmcs[i] = klms[i].conj();
    }

    calculate_moi(initial_roll);
    set_pos(speed);

    #ifdef _DEBUG
    std::cout<< "Klms: ";
    for (int i = 0; i < (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1); i++) {
        std::cout << klms[i] << ' ';
    }
    std::cout << std::endl;

    std::cout << "Jlms: ";
    for (int i = 0; i < (ASTEROIDS_MAX_J + 1) * (ASTEROIDS_MAX_K + 1); i++) {
        std::cout << jlms[i] << ' ';
    }
    std::cout << std::endl;

    std::cout << "PARAMETERS:" << std::endl;
    std::cout << "MOI (local): " << moi << std::endl;
    std::cout << "Start position: " << position << std::endl;
    std::cout << "Start velocity: " << velocity << std::endl;
    std::cout << "Start spin: " << spin << std::endl;
    std::cout << std::endl;
    #endif
}

cdouble Asteroid::jlm(uint l, int m) const {
    if (l * l + l + m < 0 || l * l + l + m > max_j) {
        throw std::runtime_error("Jlm array exceeded.");
    }
    return jlms[l * l + l + m];
}

cdouble Asteroid::klm(uint l, int m) const {
    if (abs(m) > l) {return 0;}
    if (l * l + l + m < 0 || l * l + l + m > max_k) {
        throw std::runtime_error("Klm array exceeded.");
    }
    return klms[l * l + l + m];
}

cdouble Asteroid::klmc(uint l, int m) const {
    if (abs(m) > l) {return 0;}
    if (l * l + l + m < 0 || l * l + l + m > max_k) {
        throw std::runtime_error("Klm array exceeded.");
    }
    return klmcs[l * l + l + m];
}

void Asteroid::set_pos(double speed) {
    // Asteroid enters at +x and is traveling towards -x, with offset in +y
    // direction.
    edge_dist = perigee * pow(INTEGRAL_LIMIT_FRAC, -1/3.0);
    double semi_major_axis = GM / speed / speed;
    double eccentricity = perigee / semi_major_axis + 1;
    impact_parameter = semi_major_axis * sqrt(eccentricity * eccentricity - 1);
    ang_mom = Vector3::x() * impact_parameter * speed;
    double semi_latus_rectum = ang_mom.mag2() / GM;
    double nu = -acos((semi_latus_rectum / edge_dist - 1) / eccentricity);
    position0 = Vector3({0, cos(nu), sin(nu)}) * edge_dist;
    position0 *= 1 - EPSILON;
        // So that it reads as just inside the allowed region
    velocity0 = sqrt(GM / semi_latus_rectum) * Vector3({0, -sin(nu), eccentricity + cos(nu)});
    energy = 0.5 * speed * speed;
}

void Asteroid::calculate_moi(double initial_roll) {
    double Ixx = (2/3.0 * klm(2,0) - 2.0 * klm(2, -2) - 2.0 * klm(2, 2)).r + 2/3.0;
    double Iyy = (2/3.0 * klm(2,0) + 2.0 * klm(2, -2) + 2.0 * klm(2, 2)).r + 2/3.0;
    double Izz = (- 4 / 3.0 * klm(2, 0)).r + 2/3.0;
    // These mois are really moi per radius^2 per M.

    if (Izz < 0 || Iyy < 0 || Ixx < 0) {
        throw std::runtime_error(
            "Moment of inertia values were negative."
        );
    }

    if (abs(Izz) < abs(Ixx) || abs(Izz) < abs(Iyy)) {
        throw std::runtime_error(
            "Moment of inertia was not maximized along z.");
    }

    if (Izz > Ixx + Iyy || Iyy > Ixx + Izz || Ixx > Izz + Iyy) {
        std::runtime_error("Triangle inequality violated");
    }

    moi = Vector3({Ixx, Iyy, Izz});
    inv_moi = Vector3({1.0 / Ixx, 1.0 / Iyy, 1.0 / Izz});

    Vector3 nspin = spin / spin.mag();

    double A = sqrt((1 - nspin[2]) / 2);
    double B = sqrt((nspin[0] * nspin[0] + nspin[1] * nspin[1])
        / (2 * (1 - nspin[2])));
    double phi = atan2(nspin[0], nspin[1]) - initial_roll;
    orientation0 = Quaternion(B*sin(phi), A * sin(initial_roll),
        A * cos(initial_roll), B * cos(phi));

    orientation0 = orientation0.inverse();
    orientation0 /= orientation0.mag();
    if (orientation0.is_nan()) {
        orientation0 = Quaternion(cos(initial_roll), 0, 0, sin(initial_roll));
    }
    spin0 = Vector3::z() * spin.mag();
}

int Asteroid::simulate(double cadence, std::vector<double>& resolved_data) {
    double time = 0;
    int frames = 0;
    int cadence_index = -1;
    double dt = 10;
    Vector3 spin = spin0;
    Quaternion orientation = orientation0;
    Vector3 d_spin1, d_orientation1;
    Vector3 d_spin2, d_orientation2;
    Vector3 d_spin3, d_orientation3;
    Vector3 d_spin4, d_orientation4;

    for (double time = 0;; time += time_step) {
        get_diffs(time, spin, orientation, d_spin1, d_orientation1);
        get_diffs(time + dt / 2, spin + dt * d_spin1 / 2, orientation + dt * d_orientation1 / 2,
            d_spin2, d_orientation2);
        get_diffs(time + dt / 2, spin + dt * d_spin2 / 2, orientation + dt * d_orientation2 / 2,
            d_spin3, d_orientation3);
        get_diffs(time + dt, spin + dt * d_spin3, orientation + dt * d_orientation3,
            d_spin4, d_orientation4);
        spin += dt / 6.0 * (d_spin1 + 2 * d_spin2 + 2 * d_spin3 + d_spin4)
        orientation += dt / 6.0 * (d_orientation1 + 2 * d_orientation2 + 2 * d_orientation3 + d_orientation4)

        if (int(time / cadence) > cadence_index) {
            Vector3 global_spin = orientation.matrix().transpose() * spin;
            resolved_data.push_back(global_spin[0]);
            resolved_data.push_back(global_spin[1]);
            resolved_data.push_back(global_spin[2]);
            cadence_index = int(time / cadence);
            if (distance_ratio_cut >= 0 && Vector3::dot(position, velocity) > 0
                && position.mag() / perigee >= distance_ratio_cut) {
                break;
            }
            //std::cout << spin[0] <<' '<< spin[1]<<' ' << spin[2] << ' ' << dt << std::endl;
        }
    }

    #ifdef _DEBUG
    std::cout << "Simulation took " << time << " seconds." << std::endl;
    std::cout << "Maximum dquaternion magnitude (want 0) " << max_quat_mag
        << std::endl;
    #endif

    return frames;
}

void generate_positions() {
    #ifdef _DEBUG
    auto start = high_resolution_clock::now();
    #endif
    Vector3 position = position0;
    Vector3 velocity = velocity0;
    Vector3 accel;
    for (int i = 0; position.mag() < edge_dist; i++) {
        accel = -GM / pow(position.mag(), 3) * position;
        velocity += accel * POSITION_DT;
        position += velocity * POSITION_DT;
        positions.push_back(position);
        velocities.push_back(velocity);
        max_time_index = i;
    }
    #ifdef _DEBUG
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Position time" << std::endl;
    #endif
}

bool extract_position(double time, Vector& position, Vector& velocity) {
    int time_index = time / POSITION_DT;
    if time_index >= max_time_index {
        return false;
    }
    double mixer = (time - time_index * POSITION_DT) / POSITION_DT
    position = positions[time_index] * mixer + positions[time_index] * (1 - mixer)
    velocity = velocity[time_index] * mixer + velocity[time_index] * (1 - mixer)
    return true;
}

Vector3 Asteroid::get_diff_orientation(double time, Vector3 spin, Quaternion orientation,
    Vector3& d_spin, Quaternion& d_orientation) {
    if not extract_position(time, position, velocity) {
        // Hit simulation end
        return;
    }
    const std::vector<double, 3> angles = orientation.euler_angles();
    const DMatGen dgen = DMatGen(angles[0], angles[1], angles[2]);
    const Vector3 rot_pos = orientation.matrix() * -position;
    const double rot_pos_r = position.mag();
    const double rot_pos_ct = rot_pos[2] / rot_pos.mag();
    const double rot_pos_p = atan2(rot_pos[1], rot_pos[0]);
    double x_torque = 0;
    double y_torque = 0;
    double z_torque = 0;

    for (uint l = 0; l <= ASTEROIDS_MAX_J; l++) {
        for (int m = -l; m <= (int)l; m++) {
            nowjlm = 0;
            for (int mpp = -l; mpp <= (int)l; mpp++) {
                nowjlm += sqrt(fact(l-mpp) * fact(l+mpp))
                    * dgen(l,m,mpp).conj() * jlm(l,mpp);
            }
            nowjlm *= (double)parity(l) / sqrt(fact(l-m) * fact(l+m)) * pow(RADIUS, l);
            for (uint lp = 2; lp <= ASTEROIDS_MAX_K; lp++) {
                for (int mp = -lp; mp <= (int)lp; mp++) {
                    prelpmp = nowjlm * slm_c(l + lp, m + mp, rot_pos_r, rot_pos_ct, rot_pos_p)
                        * pow(asteroid_radius, lp - 2);
                    x_torque += prelpmp * ((double)(lp - mp + 1) * klmc(lp, mp-1)
                        + (double)(lp + mp + 1) * klmc(lp, mp+1));
                    y_torque += prelpmp * ((double)(lp - mp + 1) * klmc(lp, mp-1)
                        - (double)(lp + mp + 1) * klmc(lp, mp+1));
                    z_torque += prelpmp * 2.0 * (double)mp * klmc(lp, mp);
                }
            }
        }
    }
    // This torque is really torque per radius^2 per M.

    d_spin = Vector3({
        (x_torque.i * -0.5 * GM + (moi[1] - moi[2]) * spin[1] * spin[2]) * inv_moi[0],
        (y_torque.r * -0.5 * GM + (moi[2] - moi[0]) * spin[2] * spin[0]) * inv_moi[1],
        (z_torque.i * -0.5 * GM + (moi[0] - moi[1]) * spin[0] * spin[1]) * inv_moi[2],
    });


    /// HERE"S WHERE THE PROBLEM COMES IN: THE ANGULAR CHANGE IS A SECOND DERIVATIVE
    spin += dt * d_spin;
    d_quat = 0.5 * Quaternion(0, spin[0], spin[1], spin[2]) * orientation;

    return d_orientation
}
