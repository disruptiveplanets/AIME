#include "backend.hpp"

Asteroid::Asteroid(const std::vector<cdouble> jlms,
    const std::vector<cdouble> mlms, double spin,
    double impact_parameter, double speed) :
    mlms(mlms), jlms(jlms), velocity(Vector3({0, 0, speed})),
    spin(Vector3({spin, 0, 0})), mu(jlms[0].real() * G) {

    maxml = -3/2.0 + sqrt(1 + 8 * (mlms.size()+3)) / 2;
    maxjl = -3/2.0 + sqrt(1 + 8 * jlms.size()) / 2;
    assert((maxml + 2) * (maxml + 1) / 2 == mlms.size()+3);
    assert((maxjl + 2) * (maxjl + 1) / 2 == jlms.size());

    nowmlms = mlms;// copy

    calculate_moi();
    set_pos(impact_parameter);

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

cdouble Asteroid::mlm(uint l, int m) {
    if (m < 0) {
        return std::conj(mlm(l, -m)) * (double)sign(m);
    }
    return mlms[l * (l + 1) / 2 + l - m - 3]; // l=0 and l=1 are not included
}

cdouble Asteroid::nowmlm(uint l, int m) {
    if (m < 0) {
        return std::conj(nowmlm(l, -m)) * (double)sign(m);
    }
    return nowmlms[l * (l + 1) / 2 + l -  m - 3]; // l=0 and l=1 are not included
}

void Asteroid::set_nowmlm(uint l, int m, cdouble val) {
    if (m < 0) {
        return set_nowmlm(l, -m, std::conj(val) * (double)sign(m));
    }
    nowmlms[l * (l + 1) / 2 + l - m - 3] = val;
}

cdouble Asteroid::jlm(uint l, int m) {
    if (m < 0) {
        return std::conj(jlm(l, -m)) * (double)sign(m);
    }
    return jlms[l * (l + 1) / 2 + l - m];
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

    Vector3 Omega = ang_mom / position.mag2();
    spin -= Omega;
}

void Asteroid::calculate_moi() {
    double Ixx = (-(2.0 * mlm(2, 0) - 1.0) / 3.0
        + 2.0 * mlm(2, -2) + 2.0 * mlm(2, 2)).real();
    double Iyy = (-(2.0 * mlm(2, 0) - 1.0) / 3.0
        - 2.0 * mlm(2, -2) - 2.0 * mlm(2, 2)).real();
    double Izz = (4.0* mlm(2, 0) / 3.0 + 1/3.0).real();
    double Ixz = (-mlm(2, 1) + mlm(2, -1)).real();
    double Iyz = -(mlm(2, 1) + mlm(2, -1)).imag();
    double Ixy = -2.0 * (mlm(2, 2) + mlm(2, -2)).imag();
    moi = Matrix3({ Ixx, Ixy, Ixz,
                    Ixy, Iyy, Iyz,
                    Ixz, Iyz, Izz,});

    moi_evals = moi.get_evals();
    moi_evecs = moi.get_symmetric_evecs(moi_evals);
    moiInverse = Matrix3::symmetric_invert(moi_evals, moi_evecs);

    moi = Matrix3::symmetric_reconstruct(moi_evals, moi_evecs);
        // Do if there are issues with moi * moiInverse != identity.

    // Sort the evals and evecs
    double temp_val;
    Vector3 temp_vec;

    if (moi_evals[2] > moi_evals[1]) {
        temp_val = moi_evals[1];
        temp_vec = moi_evecs[1];
        moi_evals[1] = moi_evals[2];
        moi_evecs[1] = moi_evecs[2];
        moi_evals[2] = temp_val;
        moi_evecs[2] = temp_vec;
    }
    if (moi_evals[1] > moi_evals[0]) {
        temp_val = moi_evals[1];
        temp_vec = moi_evecs[1];
        moi_evals[1] = moi_evals[0];
        moi_evecs[1] = moi_evecs[0];
        moi_evals[0] = temp_val;
        moi_evecs[0] = temp_vec;
    }
    if (moi_evals[2] > moi_evals[1]) {
        temp_val = moi_evals[1];
        temp_vec = moi_evecs[1];
        moi_evals[1] = moi_evals[2];
        moi_evecs[1] = moi_evecs[2];
        moi_evals[2] = temp_val;
        moi_evecs[2] = temp_vec;
    }

    // Set the spin to the principle axis with the greatest moi
    spin = moi_evecs[0] / moi_evecs[0].mag() * spin.mag();

    orientation = Quaternion::identity();
    //moi = Matrix3::identity();
    //moiInverse = Matrix3::identity();
}

Matrix3 Asteroid::global_to_inertial() const {
    double theta = atan2(position[1], -position[2]);
    return Matrix3::rotation_x(theta);
}

Matrix3 Asteroid::inertial_to_global() const {
    double theta = atan2(position[1], -position[2]);
    return Matrix3::rotation_x(-theta);
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

        if (int(time / cadence) > cadence_index) {
            resolved_data.push_back(spin[0]);
            resolved_data.push_back(spin[1]);
            resolved_data.push_back(spin[2]);
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

Vector3 Asteroid::get_rot_ang_mom() {
    Vector3 Omega = ang_mom / position.mag2();
    Matrix3 inv_mat = orientation.inverse().matrix();
    Matrix3 moiGlobal = orientation.matrix() * moi * inv_mat;
    return global_to_inertial() * (moiGlobal * (spin + Omega));
}

Vector3 Asteroid::get_torque() {
    double tx = 0;
    double ty = 0;
    double tz = 0;
    cdouble back;
    cdouble now;
    double pre;
    for (int l = 0; l <= maxjl; l++) {
        for (int m = -l; m <= l; m++) {
            for (int lp = max(abs(m), 2); lp <= maxml; lp++) {
                pre = sign(l + m) * fact(l + lp)
                    / pow(position.mag(), l + lp + 1);
                back = std::conj(jlm(l, m)) * nowmlm(lp, m-1);
                now = std::conj(jlm(l, m)) * nowmlm(lp, m);
                tx += pre * (lp - m + 1) * back.real();
                ty += pre * (lp - m + 1) * back.imag();
                tz += pre * 2 * m * now.imag();
            }
        }
    }

    inv_mat = orientation.inverse().matrix();
    moiGlobal = orientation.matrix() * moi * inv_mat;
    std::cout << 3 * G * jlm(0, 0).real() / pow(position.mag(), 3)
        * Vector3({-moiGlobal(2, 1), moiGlobal(2, 0), 0}) << std::endl;

    std::cout << G * Vector3({tx, ty, tz}) << std::endl;

    //return G * Vector3({tx, ty, tz});
    return 3 * G * jlm(0, 0).real() / pow(position.mag(), 3)
        * Vector3({-moiGlobal(2, 1), moiGlobal(2, 0), 0});
}

void Asteroid::update_mlms() {
    std::array<double, 3> angles = orientation.euler_angles();
    DMatGen d_generator(angles[0], angles[1], angles[2]);

    for (int l = 2; l <= maxml; l++) {
        for (int m = -l; m <= l; m++) {
            cdouble newmlm = 0;
            for (int mp = -l; mp <= -l; mp++) {
                newmlm += sqrt(fact(l - mp) * fact(l + mp)) * sign(m + mp)
                    * d_generator(l, m, mp) * mlm(l, mp);
            }
            set_nowmlm(l, m, newmlm / sqrt(fact(l - m) * fact(l + m)));
        }
    }
}

void Asteroid::update_position(double dt) {
    accel = -mu / pow(position.mag(), 3) * position;
    velocity += accel * dt;
    position += velocity * dt;
}

void Asteroid::update_orientation(double dt) {
    update_mlms();
    Omega = ang_mom / position.mag2();

    inv_mat = orientation.inverse().matrix();
    moiGlobal = orientation.matrix() * moi * inv_mat;
    moiGlobalInverse = orientation.matrix() * moiInverse * inv_mat;

    Vector3 torque = get_torque();

    omegaDot = moiGlobalInverse * (torque - Vector3::cross(
        Omega + spin, moiGlobal * (Omega + spin)))
        + 2 * Omega * Vector3::dot(position, velocity) / position.mag2()
        - Vector3::cross(Omega, spin);

    spin += dt * omegaDot;

    d_quat = 0.5 * Quaternion(0, spin[0], spin[1], spin[2])
        * orientation;
    orientation += d_quat * dt;

    #ifdef _DEBUG
    max_quat_mag = max(max_quat_mag, d_quat.mag());
    #endif

    orientation /= orientation.mag();
}
