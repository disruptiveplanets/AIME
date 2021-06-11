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
    central mass
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
    L = std::stoi(x);
    n = std::stoi(y);
    m = std::stoi(z);

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
    for(int s = 1; s <= m; s++){
        std::getline(f, line);
        ss = std::istringstream (line);
        for(int i = 0; i <= n/m * s; i++) {
            ss >> word;
            densities.push_back(std::stod(word));
        }
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
    position = Vector3({0, std::stod(word), 0});

    // Central mass
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    mu = G * std::stod(word);
}

void Asteroid::simulate(std::ofstream&& resolved, std::ofstream&& unresolved){

}

void Asteroid::generate_shape() {

}

void Asteroid::calculate_moi() {

}

void Asteroid::calculate_com() {

}
void Asteroid::calculate_mass() {

}
void Asteroid::update() {
    double dt = 1;// Change this
    update_orientation(dt);
    update_position(dt);
}
void Asteroid::update_position(double dt) {
    
}
void Asteroid::update_orientation(double dt) {

}