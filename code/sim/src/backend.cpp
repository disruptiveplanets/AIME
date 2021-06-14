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

    // Density
    triangles.reserve(12 * num_chunks);
        // Over estimate the number of triangles about to be pushed
    std::getline(f, line);
    ss = std::istringstream (line);
    for(int i = 0; i <= num_chunks; i++) {
        ss >> word;
        chunks[i].shape(std::stod(word), L, clms, triangles);
    }

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
    velocity = Vector3({-std::stod(word), 0, 0});

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

void Asteroid::make_chunks() {
    for (int f = 0; f < 6; f++) {
        for (int s = 0; s < m; s++){
            const int length = n * (s + 1) / m;
            const double alpha_above = float(s + 1) / m;
            const double alpha_below = float(s) / m;
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
    double dist = b * pow(INTEGRAL_LIMIT_FRAC, -1/3.0);
    position = Vector3({sqrt(dist * dist - b * b), b, 0});

    double a = mu / velocity.mag2(); // Positive
    closest_approach = sqrt(a * a + b * b) - a;
}

int Asteroid::simulate(std::ofstream&& resolved, std::ofstream&& unresolved){
    // Motivation: Torque is proportional to 1/position^3. Angular acceleration is now roughly constant every frame.

    int frames = 0;
    for (;;frames++){
        double dt = DELTA_T_MIN * pow(position.mag() / closest_approach, 3);
        update_orientation(dt);
        update_position(dt);
    }
    return frames;
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

Vector3 Asteroid::get_com() {
    Vector3 total_arm = Vector3({0, 0, 0});
    for (Triangle& t : triangles) {
        total_arm += t.get_lever_arm();
    }
    return 1 / mass * total_arm;
}

Vector3 Asteroid::get_torque() {
    Vector3 torque = Vector3({0, 0, 0});
    for (Triangle& t : triangles) {
        torque += t.get_torque();
    }
    return (mu / pow(position.mag(), 3)) * torque;
}

void Asteroid::update_position(double dt) {
    Vector3 accel = -mu / pow(position.mag(), 3) * position;
    velocity += accel * dt;
    position += accel * dt;
}

void Asteroid::update_orientation(double dt) {

}