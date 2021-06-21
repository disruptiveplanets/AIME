#include "backend.hpp"

Asteroid::Asteroid(std::string filename) {
    /* Expected format:
    L n m
    C00
    C1-1 C10 C11
    <...>
    rho00 rho01 rho02 rho03 rho04 rho05
    <...>
    spinx, spiny, spinz
    impact parameter
    velocity mps
    central mass
    */

    /* Note: since this code is intended to be used in a limited set of
     * scnarios, and since I will always have access to this source code
     * to debug any errors. I will not be writing a "safe parser". That is,
     * when I parse the file `filename`, I will not carefully print out the
     * error messages that could occur. The program will crash somehow and it
     * is up to the user to figure out the problem with the formatting of
     * `filename`. I'm doing this to save time.
    */


    std::ifstream f;
    f.open(filename);
    std::string line, word, x, y, z;
    std::istringstream ss;

    // Parameterization numbers
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x;
    ss >> y;
    ss >> z;
    L = std::stoi(x);// Max degree of spherical harmonics
    n = std::stoi(y);// Segmentation of cube face at outer shell
    m = std::stoi(z);// Number of shells

    int num_chunks = n * n * (m + 1) * (2 * m + 1) / m;
    chunks.reserve(num_chunks);
    make_chunks();

    std::vector<double> clms;

    // Clms
    for(int i = 0; i <= L; i++){
        std::getline(f, line);
        ss = std::istringstream (line);
        for(int m = -i; m <= i; m++) {
            ss >> word;
            clms.push_back(std::stod(word));
        }
    }
    radius = clms[0];

    // Density
    triangles.reserve(12 * num_chunks);
        // Over estimate the number of triangles about to be pushed
    std::getline(f, line);
    ss = std::istringstream (line);
    mean_density = 0;
    for(int i = 0; i <= num_chunks; i++) {
        ss >> word;
        mean_density += std::stod(word);
        chunks[i].shape(std::stod(word), L, clms, triangles);
    }
    mean_density /= chunks.size();

    // Spin
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x;
    ss >> y;
    ss >> z;
    spin = Vector3({std::stod(x), std::stod(y), std::stod(z)});

    // Impact parameter
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double impact_parameter = std::stod(word);

    // Velocity (units: times escape velocity)
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    velocity = Vector3({0, 0, std::stod(word)});

    // Central mass
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    mu = G * std::stod(word);

    calculate_mass();
    recenter();
    set_pos(impact_parameter);
    calculate_moi();

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

void Asteroid::draw(std::string filename, Vector3 axis) const {
    // Write data about how to draw an asteroid to a text file (.ast) to be read
    // by a python script and displayed with matplotlib.

    // Project everything along axis. Axis must be normalized
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
    double Ixx = 0, Iyy = 0, Izz = 0, Ixz = 0, Ixy = 0, Iyz = 0;
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
    std::array<double, 3> evals = moi.get_evals();
    std::array<Vector3, 3> evecs = moi.get_symmetric_evecs(evals);
    moiInverse = Matrix3::symmetric_invert(evals, evecs);

    Matrix3 new_moi = Matrix3::symmetric_reconstruct(evals, evecs);
        // Do if there are issues with moi * moiInverse != identity.

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

    #ifdef _DEBUG
    std::cout << "New COM: " << get_com() << " should be zero." << std::endl;
    #endif
}

Matrix3 Asteroid::global_to_inertial() const {
    double theta = atan2(position[1], -position[2]);
    return Matrix3::rotation_x(theta);
}

Matrix3 Asteroid::inertial_to_global() const {
    double theta = atan2(position[1], -position[2]);
    return Matrix3::rotation_x(-theta);
}

int Asteroid::simulate(std::ofstream&& resolved, std::ofstream&& unresolved) {
    // Motivation: Torque is proportional to 1/position^3. Angular acceleration
    // is now roughly constant every frame.
    const double scale_torque = mu / pow(closest_approach, 3) * mean_density *
        pow(radius, 5);
    const double min_delta_t = ONE_SECOND_TORQUE / scale_torque;
        // Toruqe at min delta t

    time = 0;
    int frames = 0;
    for (;position.mag() < edge_dist; frames++){
        double dt = min_delta_t * pow(position.mag() / closest_approach, 3);
        if (dt > MAX_DT) {
            dt = MAX_DT;
        }
        update_orientation(dt);
        update_position(dt);
        time += dt;
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
    Vector3 accel = -mu / pow(position.mag(), 3) * position;
    velocity += accel * dt;
    position += velocity * dt;
}

void Asteroid::update_orientation(double dt) {
    Vector3 Omega = ang_mom / position.mag2();
    Vector3 torque = Vector3::zero(); // get_torque();

    Matrix3 inv_mat = orientation.inverse().matrix();
    Matrix3 moiGlobal = orientation.matrix() * moi * inv_mat;
    Matrix3 moiGlobalInverse = orientation.matrix() * moiInverse * inv_mat;

    Vector3 omegaDot = moiGlobalInverse * (torque - Vector3::cross(
        Omega + spin, moiGlobal * (Omega + spin)))
        + 2 * Omega * Vector3::dot(position, velocity) / position.mag2()
        - Vector3::cross(Omega, spin);

    spin += dt * omegaDot;

    Quaternion d_quat = 0.5 * Quaternion(0, spin[0], spin[1], spin[2])
        * orientation;
    orientation += d_quat * dt;

    #ifdef _DEBUG
    max_quat_mag = max(max_quat_mag, d_quat.mag());
    #endif

    orientation /= orientation.mag();
}
