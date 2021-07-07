#include "backend.hpp"

Asteroid::Asteroid(int L, int n, int m, const std::vector<double>& clms,
    const std::vector<double>& densities, double spin, double impact_parameter,
    double speed, double central_mass) :
    L(L), n(n), m(m), velocity(Vector3({0, 0, speed})),
    spin(Vector3({spin, 0, 0})), mu(G * central_mass),
    radius(clms[0]) {

    int num_chunks = n * n * (m + 1) * (2 * m + 1) / m;
    make_chunks();

    // Density
    triangles.reserve(12 * num_chunks);
    mean_density = 0;
    for(int i = 0; i < densities.size(); i++) {
        if (i >= chunks.size()) {
            std::cout << "Too few chunks." << std::endl;
        }
        mean_density += densities[i];
        chunks[i].shape(densities[i], L, clms, triangles);
    }
    mean_density /= chunks.size();

    calculate_mass();
    recenter();
    calculate_moi();
    set_pos(impact_parameter);

    #ifdef _DEBUG
    std::cout << "PARAMETERS:" << std::endl;
    std::cout << "L=" << L << ", m="<<m<<", n=" << n << std::endl;
    std::cout << "Mass: " << mass << std::endl;
    std::cout << "MOI (local): " << moi << std::endl;
    std::cout << "Start position: " << position << std::endl;
    std::cout << "Start velocity: " << velocity << std::endl;
    std::cout << "Start spin: " << spin << std::endl;
    std::cout << "Central mu: " << mu << std::endl;
    std::cout << std::endl;
    #endif
}

void Asteroid::draw(std::string filename, int axis_num) const {
    // Write data about how to draw an asteroid to a text file (.ast) to be read
    // by a python script and displayed with matplotlib.

    // Project everything along axis. Axis must be normalized
    Vector3 axis = moi_evecs[axis_num];
    Vector3 up = Vector3::z();
    if (Vector3::cross(up, axis).mag() < EPSILON) {
        up = Vector3::x();
    }
    up = up - Vector3::dot(up, axis) * axis;
    up /= up.mag();

    std::vector<double> depths;
    std::vector<std::string> lines;
    depths.reserve(triangles.size());
    lines.reserve(triangles.size());

    for (Triangle const& triangle : triangles) {
        if (!triangle.is_edge()){
            continue;
        }
        std::array<Vector3, 3> corners = triangle.get_corners();
        Vector3 center = (corners[0] + corners[1] + corners[2]) / 3;
        Vector3 norm = Vector3::cross(corners[1] - corners[0],
            corners[2] - corners[0]);
        norm /= norm.mag();
        if (Vector3::dot(norm, axis) < 0) {
            continue;
        }
        double depth = Vector3::dot(center, axis);

        std::string line = "";
        for(Vector3 const& c : corners) {
            Vector3 offset = c - Vector3::dot(c, axis) * axis;
            if (offset.mag() == 0){
                line += "0 0 ";
                continue;
            }
            double theta = acos(Vector3::dot(offset, up) / offset.mag());
            if (Vector3::dot(Vector3::cross(offset, up), axis) < 0) {
                theta = 2 * PI - theta;
            }
            double r = offset.mag();
            line += std::to_string(r * cos(theta)) + " "
                    + std::to_string(r * sin(theta)) + " ";
        }
        line += std::to_string(triangle.get_density()) + " "
                + std::to_string(Vector3::dot(norm, axis));

        // Insert the line and depth in order from low depth to high
        if (lines.size() == 0) {
            depths.push_back(depth);
            lines.push_back(line);
        }
        int i;
        for (i = 0; i < depths.size(); i++){
            if (depths[i] > depth) {
                depths.insert(depths.begin() + i, depth);
                lines.insert(lines.begin() + i, line);
                break;
            }
        }
        if (i == depths.size()){
            depths.push_back(depth);
            lines.push_back(line);
        }
    }

    std::ofstream output;
    output.open(filename);
    for (std::string& line : lines) {
        output << line << '\n';
    }
    output.close();
}

void Asteroid::make_chunks() {
    for (int f = 0; f < 6; f++) {
        for (int s = 0; s < m; s++){
            const int length = n * (s + 1) / m;
            const double alpha_above = double(s + 1) / m;
            const double alpha_below = double(s) / m;
            for (int a = 0; a < length; a++){
                for (int b = 0; b < length; b++) {
                    chunks.push_back(
                        Chunk(alpha_above, alpha_below, length, f, a, b));
                }
            }
        }
    }
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
    double Ixx = 0;
    double Iyy = 0;
    double Izz = 0;
    double Ixz = 0;
    double Ixy = 0;
    double Iyz = 0;
    for (Triangle& t : triangles) {
        Ixx += t.get_Isame(0);
        Iyy += t.get_Isame(1);
        Izz += t.get_Isame(2);
        Ixy += t.get_Idiff(0, 1);
        Ixz += t.get_Idiff(0, 2);
        Iyz += t.get_Idiff(1, 2);
    }
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

void Asteroid::calculate_mass() {
    mass = 0;
    for (Triangle const& t : triangles) {
        mass += t.get_mass();
    }
}

void Asteroid::recenter() {
    Vector3 com = get_com();
    for (Triangle& t : triangles) {
        t.recenter(-com);
    }

    /*#ifdef _DEBUG
    std::cout << "New COM: " << get_com() << " should be zero." << std::endl;
    #endif*/
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

            if (isnan(spin[0])) {
                std::cout << mass << std::endl;
                std::cout << moi << std::endl;
                std::cout << moiInverse << std::endl;
                std::cout << edge_dist << std::endl;
                std::cout << position << std::endl;
                std::cout << velocity << std::endl;
                std::cout << orientation << std::endl;
                std::cout << mean_density << std::endl;
                std::cout << radius << std::endl;
                std::cout << energy << std::endl;
                std::cout << ang_mom << std::endl;
                std::cout << closest_approach << std::endl;
                std::cout << excess_vel << std::endl;
                std::cout << impact_parameter << std::endl;
                return 0;
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

Vector3 Asteroid::get_rot_ang_mom() {
    Vector3 Omega = ang_mom / position.mag2();
    Matrix3 inv_mat = orientation.inverse().matrix();
    Matrix3 moiGlobal = orientation.matrix() * moi * inv_mat;
    return global_to_inertial() * (moiGlobal * (spin + Omega));
}

Vector3 Asteroid::get_com() const {
    Vector3 total_arm = Vector3::zero();
    for (Triangle const& t : triangles) {
        total_arm += t.get_lever_arm();
    }
    return 1 / mass * total_arm;
}

Vector3 Asteroid::get_torque() {
    Vector3 torque = Vector3::zero();
    for (Triangle t : triangles) {
        t *= orientation.matrix();
        torque += t.get_torque();
    }
    return (mu / pow(position.mag(), 3)) * torque;
}

void Asteroid::update_position(double dt) {
    accel = -mu / pow(position.mag(), 3) * position;
    velocity += accel * dt;
    position += velocity * dt;
}

void Asteroid::update_orientation(double dt) {
    Omega = ang_mom / position.mag2();

    inv_mat = orientation.inverse().matrix();
    moiGlobal = orientation.matrix() * moi * inv_mat;
    moiGlobalInverse = orientation.matrix() * moiInverse * inv_mat;

    //Vector3 torque = get_torque();
    torque = 3 * mu / pow(position.mag(), 3) *
        Vector3({-moiGlobal(1, 2), moiGlobal(0, 2), 0});

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
