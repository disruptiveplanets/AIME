#include "backend.hpp"

Asteroid::Asteroid(const std::vector<cdouble> jlms,
    const std::vector<cdouble> klms, Vector3 spin, double initial_roll,
    double impact_parameter, double speed) :
    jlms(jlms), klms(klms), velocity(Vector3({0, 0, speed})),
    spin(spin), mu(jlms[0].real() * G) {

    maxkl = sqrt(klms.size()) - 1;
    maxjl = sqrt(jlms.size()) - 1;
    assert((maxkl + 1) * (maxkl + 1) == klms.size());
    assert((maxjl + 1) * (maxjl + 1) == jlms.size());

    calculate_moi(initial_roll);
    set_pos(impact_parameter);

    //for(auto j : jlms) {std::cout << j << ' ';} std::cout << std::endl;
    //for(auto j : klms) {std::cout << j << ' ';} std::cout << std::endl;

    #ifdef _DEBUG
    std::cout << "PARAMETERS:" << std::endl;
    std::cout << "MOI (local): " << moi << std::endl;
    std::cout << "Start position: " << position << std::endl;
    std::cout << "Start velocity: " << velocity << std::endl;
    std::cout << "Start spin: " << spin << std::endl;
    std::cout << "Central mu: " << mu << std::endl;
    std::cout << std::endl;
    #endif
}

cdouble Asteroid::jlm(uint l, int m) const {
    return jlms[l * l + l + m];
}

cdouble Asteroid::klm(uint l, int m) const {
    if (abs(m) > l) {return 0;}
    return klms[l * l + l + m];
}

void Asteroid::set_pos(double b) {
    // Asteroid enters at +x and is traveling towards -x, with offset in +y
    // direction.
    edge_dist = b * pow(INTEGRAL_LIMIT_FRAC, -1/3.0);
    position = Vector3({0, b, -sqrt(edge_dist * edge_dist - b * b)});
    position *= 1 - EPSILON;
        // So that it reads as just inside the allowed region

    energy = 0.5 * velocity.mag2() - mu / position.mag();
    ang_mom = Vector3::cross(position, velocity);
    double velocity_periapsis = mu / ang_mom.mag() +
        sqrt(mu * mu / ang_mom.mag2() + 2 * energy);
    closest_approach = ang_mom.mag() / velocity_periapsis;
    excess_vel = sqrt(2 * energy);
    impact_parameter = ang_mom.mag() / excess_vel;
}

void Asteroid::calculate_moi(double initial_roll) {
    double sph = 2/5.0 * klm(0, 0).real();
    double Ixx = (2/3.0 * klm(2,0) - 2.0 * klm(2, -2) - 2.0 * klm(2, 2)).real()
        + sph;
    double Iyy = (2/3.0 * klm(2,0) + 2.0 * klm(2, -2) + 2.0 * klm(2, 2)).real()
        + sph;
    double Izz = (- 4 / 3.0 * klm(2, 0)).real()
        + sph;

    assert(abs(Izz) > abs(Ixx) && abs(Izz) > abs(Iyy));

    moi = Vector3({Ixx, Iyy, Izz});

    Vector3 nspin = spin / spin.mag();

    double A = sqrt((1 - nspin[2]) / 2);
    double B = sqrt((nspin[0] * nspin[0] + nspin[1] * nspin[1])
        / (2 * (1 - nspin[2])));
    double phi = atan2(nspin[0], nspin[1]) - initial_roll;
    orientation = Quaternion(B*sin(phi), A * sin(initial_roll),
        A * cos(initial_roll), B * cos(phi));

    orientation = orientation.inverse();
    orientation /= orientation.mag();
    if (orientation.is_nan()) {
        orientation = Quaternion(cos(initial_roll), 0, 0, sin(initial_roll));
    }
    spin = Vector3::z() * spin.mag();
}

int Asteroid::simulate(double cadence, std::vector<double>& resolved_data) {
    time = 0;
    int frames = 0;
    int cadence_index = -1;
    double dt;
    double close_power = pow(closest_approach, DT_POWER);
    double denom = 1.0 / (pow(edge_dist, DT_POWER) - close_power);

    for (;position.mag() < edge_dist; frames++) {
        dt = MIN_DT + (MAX_DT - MIN_DT) *
            (pow(position.mag(), DT_POWER) - close_power) * denom;
        update_orientation(dt);
        update_position(dt);
        time += dt;

        /*if (frames == 200) {
            int i = 0;
            std::cout << time << std::endl;
            for(i = 0; i < 1e15; i++);
            std::cout << i << std::endl;
        }*/

        if (int(time / cadence) > cadence_index) {
            Vector3 global_spin = orientation.matrix().transpose() * spin;
            resolved_data.push_back(global_spin[0]);
            resolved_data.push_back(global_spin[1]);
            resolved_data.push_back(global_spin[2]);
            cadence_index = int(time / cadence);
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

Vector3 Asteroid::get_torque() {
    angles = orientation.euler_angles();
    DMatGen dgen = DMatGen(angles[0], angles[1], angles[2]);
    rot_pos = orientation.matrix() * -position;
    rot_pos_r = rot_pos.mag();
    rot_pos_ct = rot_pos[2] / rot_pos.mag();
    rot_pos_p = atan2(rot_pos[1], rot_pos[0]);
    x_torque = 0;
    y_torque = 0;
    z_torque = 0;
    for (uint l = 0; l <= maxjl; l++) {
        for (int m = -l; m <= (int)l; m++) {
            nowjlm = 0;
            for (int mpp = -l; mpp <= (int)l; mpp++) {
                nowjlm += sqrt(fact(l-mpp) * fact(l+mpp))
                    * std::conj(dgen(l,m,mpp)) * jlm(l,mpp);
            }
            nowjlm *= (double)parity(l) / sqrt(fact(l-m) * fact(l+m));
            for (uint lp = 2; lp <= maxkl; lp++) {
                for (int mp = -lp; mp <= (int)lp; mp++) {
                    prelpmp = nowjlm * slm_c(l + lp, m + mp,
                        rot_pos_r, rot_pos_ct, rot_pos_p);
                    x_torque += prelpmp * ((double)(lp - mp + 1)
                        * std::conj(klm(lp, mp-1))
                        + (double)(lp + mp + 1)
                        * std::conj(klm(lp, mp+1)));
                    y_torque += prelpmp * ((double)(lp - mp + 1)
                        * std::conj(klm(lp, mp-1))
                        - (double)(lp + mp + 1)
                        * std::conj(klm(lp, mp+1)));
                    z_torque += prelpmp * 2.0 * (double)mp
                        * std::conj(klm(lp, mp));
                }
            }
        }
    }
    //std::cout << x_torque << ' ' << y_torque << ' ' << z_torque << std::endl;
    return Vector3({-std::move(x_torque.imag()),
        std::move(y_torque.real()),
        std::move(-z_torque.imag())})
        * -0.5 * G;
}

void Asteroid::update_position(double dt) {
    accel = -mu / pow(position.mag(), 3) * position;
    velocity += accel * dt;
    position += velocity * dt;
}

void Asteroid::update_orientation(double dt) {
    torque = get_torque();

    omegaDot = Vector3({
        (torque[0] + (moi[1] - moi[2]) * spin[1] * spin[2]) / moi[0],
        (torque[1] + (moi[2] - moi[0]) * spin[2] * spin[0]) / moi[1],
        (torque[2] + (moi[0] - moi[1]) * spin[0] * spin[1]) / moi[2],
    });

    spin += dt * omegaDot;

    d_quat = 0.5 * Quaternion(0, spin[0], spin[1], spin[2]) * orientation;
    orientation += d_quat * dt;

    #ifdef _DEBUG
    max_quat_mag = max(max_quat_mag, d_quat.mag());
    #endif

    orientation /= orientation.mag();
}
