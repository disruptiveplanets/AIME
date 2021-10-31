#include "backend.hpp"

const int max_k = (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1);
const int max_j = (ASTEROIDS_MAX_J + 1) * (ASTEROIDS_MAX_J + 1);

Asteroid::Asteroid(const cdouble* jlms, const cdouble* klms, double asteroid_radius,
    Vector3 spin, double initial_roll, double perigee, double speed,
    double distance_ratio_cut) :
    jlms(jlms), klms(klms), asteroid_radius(asteroid_radius), distance_ratio_cut(distance_ratio_cut),
    velocity(Vector3({0, 0, speed})), spin(spin), perigee(perigee) {

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
    position = Vector3({0, cos(nu), sin(nu)}) * edge_dist;
    position *= 1 - EPSILON;
        // So that it reads as just inside the allowed region
    velocity = sqrt(GM / semi_latus_rectum) * Vector3({0, -sin(nu), eccentricity + cos(nu)});
    energy = 0.5 * speed * speed;
}

void Asteroid::calculate_moi(double initial_roll) {
    double Ixx = (2/3.0 * klm(2,0) - 2.0 * klm(2, -2) - 2.0 * klm(2, 2)).real() + 2/3.0;
    double Iyy = (2/3.0 * klm(2,0) + 2.0 * klm(2, -2) + 2.0 * klm(2, 2)).real() + 2/3.0;
    double Izz = (- 4 / 3.0 * klm(2, 0)).real() + 2/3.0;
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
    //double close_power = pow(perigee, DT_POWER);
    //double denom = 1.0 / (pow(edge_dist, DT_POWER) - close_power);

    for (;position.mag() < edge_dist; frames++) {
        dt = 0.5;
        //MIN_DT + (MAX_DT - MIN_DT) * (pow(position.mag(), DT_POWER) - close_power) * denom;
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

Vector3 Asteroid::get_torque() {
    angles = orientation.euler_angles();
    DMatGen dgen = DMatGen(angles[0], angles[1], angles[2]);
    rot_pos = orientation.matrix() * -position;
    rot_pos_r = position.mag();
    rot_pos_ct = rot_pos[2] / rot_pos.mag();
    rot_pos_p = atan2(rot_pos[1], rot_pos[0]);
    x_torque = 0;
    y_torque = 0;
    z_torque = 0;
    for (uint l = 0; l <= ASTEROIDS_MAX_J; l++) {
        for (int m = -l; m <= (int)l; m++) {
            nowjlm = 0;
            for (int mpp = -l; mpp <= (int)l; mpp++) {
                nowjlm += sqrt(fact(l-mpp) * fact(l+mpp))
                    * std::conj(dgen(l,m,mpp)) * jlm(l,mpp);
            }
            nowjlm *= (double)parity(l) / sqrt(fact(l-m) * fact(l+m)) * pow(RADIUS, l);
            for (uint lp = 2; lp <= ASTEROIDS_MAX_K; lp++) {
                for (int mp = -lp; mp <= (int)lp; mp++) {
                    prelpmp = nowjlm * slm_c(l + lp, m + mp, rot_pos_r, rot_pos_ct, rot_pos_p)
                        * pow(asteroid_radius, lp - 2);
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
        -std::move(z_torque.imag())})
        * -0.5 * GM;
    // This torque is really torque per radius^2 per M.
}

void Asteroid::update_position(double dt) {
    accel = -GM / pow(position.mag(), 3) * position;
    velocity += accel * dt;
    position += velocity * dt;
}

void Asteroid::update_orientation(double dt) {
    torque = get_torque();

    omegaDot = Vector3({
        (torque[0] + (moi[1] - moi[2]) * spin[1] * spin[2]) * inv_moi[0],
        (torque[1] + (moi[2] - moi[0]) * spin[2] * spin[0]) * inv_moi[1],
        (torque[2] + (moi[0] - moi[1]) * spin[0] * spin[1]) * inv_moi[2],
    });

    /*if (isnan(omegaDot[0])) {
        std::cout << spin << std::endl;
        std::cout << torque << std::endl;
        std::cout << orientation << std::endl;
        std::cin >> time;
    }*/

    spin += dt * omegaDot;

    d_quat = 0.5 * Quaternion(0, spin[0], spin[1], spin[2]) * orientation;
    orientation += d_quat * dt;

    /*if (orientation.mag() == 0) {
        std::cout << spin << std::endl;
        std::cout << torque << std::endl;
        std::cout << orientation << std::endl;
        std::cin >> time;
    }*/

    #ifdef _DEBUG
    max_quat_mag = max_me(max_quat_mag, d_quat.mag());
    #endif

    orientation /= orientation.mag();
}
