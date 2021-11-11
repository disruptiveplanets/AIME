#include "backend.hpp"

#define sign(a) (((a) % 2 == 0) ? 1 : -1)

const int max_k = (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1);
const int max_j = (ASTEROIDS_MAX_J + 1) * (ASTEROIDS_MAX_J + 1);

Asteroid::Asteroid(const cdouble* jlms, const cdouble* klms, double asteroid_radius,
    Vector3 spin, double initial_roll, double perigee, double speed,
    double distance_ratio_cut) :
    jlms(jlms), klms(klms), asteroid_radius(asteroid_radius), distance_ratio_cut(distance_ratio_cut),
    pericenter_pos(perigee), excess_vel(speed), initial_spin(spin), initial_roll(initial_roll) {

    for (int i = 0; i < (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1); i++) {
        klmcs[i] = klms[i].conj();
    }

    calculate_moi(initial_roll);
    calculate_poses();


    #ifdef _DEBUG
    std::cout<< "Klms: ";
    for (int i = 0; i < (ASTEROIDS_MAX_K + 1) * (ASTEROIDS_MAX_K + 1); i++) {
        std::cout << klms[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << "Jlms: ";
    for (int i = 0; i < (ASTEROIDS_MAX_J + 1) * (ASTEROIDS_MAX_J + 1); i++) {
        std::cout << jlms[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << "MOI (local): " << moi << std::endl;
    std::cout << std::endl;
    #endif
}

cdouble Asteroid::jlm(uint l, int m) const {
    #ifdef _DEBUG
    if (l * l + l + m < 0 || l * l + l + m > max_j) {
        throw std::runtime_error("Jlm array exceeded.");
    }
    #endif
    return jlms[l * l + l + m];
}

cdouble Asteroid::klm(uint l, int m) const {
    //if (abs(m) > l) {return 0;}
    #ifdef _DEBUG
    if (l * l + l + m < 0 || l * l + l + m > max_k) {
        throw std::runtime_error("Klm array exceeded.");
    }
    #endif
    return klms[l * l + l + m];
}

cdouble Asteroid::klmc(uint l, int m) const {
    //if (abs(m) > l) {return 0;}
    #ifdef _DEBUG
    if (l * l + l + m < 0 || l * l + l + m > max_k) {
        throw std::runtime_error("Klm array exceeded.");
    }
    #endif
    return klmcs[l * l + l + m];
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

    moi = Matrix3({Ixx, 0, 0, 0, Iyy, 0, 0, 0, Izz});
    inv_moi = Matrix3({1.0 / Ixx, 0, 0, 0, 1.0 / Iyy, 0, 0, 0, 1.0 / Izz});
}

void Asteroid::calculate_poses() {
    #ifdef _DEBUG
    auto start = std::chrono::high_resolution_clock::now();
    #endif

    double pericenter_vel = sqrt(excess_vel * excess_vel + 2 * GM / pericenter_pos);
    const double cutoff_dist_squared = pericenter_pos * pericenter_pos * pow(INTEGRAL_LIMIT_FRAC, -2/3.0);
    Vector3 velocity = Vector3({0, pericenter_vel, 0});
    Vector3 position = Vector3({pericenter_pos, 0, 0});

    double time;
    for (time = 0; position.mag2() < cutoff_dist_squared; time += POSITION_DT) {
        Vector3 accel = -GM / position.mag2() * position / position.mag();
        velocity += POSITION_DT * accel;
        position += POSITION_DT * velocity;
        positions.push_back(position);
        velocities.push_back(velocity);
    }
    expire_time = time;

    #ifdef _DEBUG
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Position time " << duration.count()/ 1000000.0 << std::endl;
    std::cout << "Expire time (h): " << expire_time / 3600 << std::endl;
    #endif
}

Vector3 Asteroid::extract_spin(Vector3 angles, Vector3 momenta) {
    // Local
    double ca = cos(angles[0]);
    double sa = sin(angles[0]);
    double cb = cos(angles[1]);
    double sb = sin(angles[1]);
    Matrix3 a_inv = 1 / sb * Matrix3({
        -ca * cb, -sa * cb, sb,
        -sa * sb, ca * sb, 0,
        ca, sa, 0});
    return inv_moi * a_inv.transpose() * momenta;
    /// Check to make sure the multiplication is done in the correct order.
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
    position = positions[time_index] * mixer + positions[time_index] * (1 - mixer);
    velocity = velocities[time_index] * mixer + velocities[time_index] * (1 - mixer);
    if (flip) {
        position = Vector3({position[0], -position[1], 0});
        position = Vector3({-position[0], position[1], 0});
    }
    return true;
}
Matrix3 Asteroid::local_to_global(Vector3 angles) {
    const double ca = cos(angles[0]);
    const double sa = sin(angles[0]);
    const double cb = cos(angles[1]);
    const double sb = sin(angles[1]);
    const double cc = cos(angles[2]);
    const double sc = sin(angles[2]);

    return Matrix3({cc, -sc, 0, sc, cc, 0, 0, 0, 1})
        * Matrix3({cb, 0, sb, 0, 1, 0, -sb, 0, cb})
        * Matrix3({ca, -sa, 0, sa, ca, 0, 0, 0, 1});
}

void Asteroid::get_derivatives(Vector3 position, const Vector3 angles, const Vector3 momenta, Vector3& dangles, Vector3& dmomenta) {
    const double ca = cos(angles[0]);
    const double sa = sin(angles[0]);
    const double cb = cos(angles[1]);
    const double sb = sin(angles[1]);
    const double pos_costheta = position[2] / position.mag();
    const double pos_phi = atan2(position[1], position[0]);
    DMatGen wignerd = DMatGen(angles[0], angles[1], angles[2]);
    const Matrix3 a_inv = 1 / sb * Matrix3({
        -ca * cb, -sa * cb, sb,
        -sa * sb, ca * sb, 0,
        ca, sa, 0});

    const Matrix3 a_inv_x= 1 / sb * Matrix3({
        sa * cb, -ca * cb, 0,
        -ca * sb, -sa * sb, 0,
        -sa, ca, 0});
    const Matrix3 a_inv_y = 1 / sb / sb * Matrix3({
        ca, sa, 0,
        0, 0, 0,
        -ca * cb, -sa * cb, 0});


    const Vector3 inv_transpose_p = inv_moi * a_inv.transpose() * momenta;
    dangles = a_inv * inv_transpose_p;
    dmomenta = Vector3({Vector3::dot(momenta, a_inv_x * inv_transpose_p),
        Vector3::dot(momenta, a_inv_y * inv_transpose_p), 0});

    cdouble dpx;
    cdouble dpy;
    cdouble dpz;

    for (int l = 0; l < ASTEROIDS_MAX_J; l++) {
        for (int m = l; m <= l; m++) {
            for (int lp = 0; lp < ASTEROIDS_MAX_K; lp++) {
                for (int mp = -lp; mp <= lp; mp++) {
                    cdouble slm_now = slm_c(l+lp, m+mp, position.mag(), pos_costheta, pos_phi);
                    for (int mpp = -lp; mpp <= lp; mpp++) {
                        cdouble premul = -GM * position.mag() * sign(lp + mp + mpp) * slm_now * pow(RADIUS, l)
                        * pow(asteroid_radius, lp) * fact(lp+mpp) / fact(lp+mp);
                        cdouble wd = wignerd(l, mp, mpp);
                        dpx += premul * cdouble(0, m) * angles[0] * wd;
                        dpy += premul * wignerd.db(l, mp, mpp);
                        dpz += premul * cdouble(0, mp) * angles[2] * wd;
                    }
                }
            }
        }
    }
    //std::cout << momenta << ' ' << dmomenta << std::endl;
    if (dangles.is_nan() || dmomenta.is_nan()) {
        throw std::runtime_error("Nan encountered");
    }
    //std::cout << dpx << ' ' << dpy << ' ' << dpz << ' ' << std::endl;
    /// Update momentum
}

int Asteroid::simulate(const double cadence, std::vector<double>& resolved_data) {
    double beta = acos(initial_spin[2] / initial_spin.mag());
    double gamma = atan2(initial_spin[1], initial_spin[0]);
    Vector3 angles = Vector3({initial_roll, beta, gamma});

    Vector3 momenta = 2 / 3.0 * (-2 * klm(2, 0).r + 1) * initial_spin.mag() * Vector3({1, 0, cos(angles[1])});


    momenta += Vector3::y() * momenta[2] / 100;


    Vector3 dangles, dmomenta;
    Vector3 position, velocity;
    double dt = 1;
    int cadence_index = -expire_time / cadence - 1;
    int frames = 0;
    for (double time = dt-expire_time; time < expire_time; time += dt) {
        if (!extract_pos(time, position, velocity)) {
            throw std::runtime_error("Simulation ran out of bounds");
        }
        if (distance_ratio_cut > 0 &&
            position.mag2() > distance_ratio_cut * pericenter_pos &&
            Vector3::dot(position, velocity) > 0) {
            break;
        }

        // Update
        get_derivatives(position, angles, momenta, dangles, dmomenta);
        angles += dangles * dt;
        momenta += dmomenta * dt;

        if (int(time / cadence) > cadence_index) {

            Vector3 global_spin = local_to_global(angles) * extract_spin(angles, momenta);
            resolved_data.push_back(global_spin[0]);
            resolved_data.push_back(global_spin[1]);
            resolved_data.push_back(global_spin[2]);
            cadence_index = int(time / cadence);
            //std::cout << spin[0] <<' '<< spin[1]<<' ' << spin[2] << ' ' << dt << std::endl;
        }
        frames ++;
    }

    return frames;
}
