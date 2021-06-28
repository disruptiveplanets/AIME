#include <fstream>
#include <string>
#include <iostream>
#include <chrono>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"
#include "../sim/triangle.hpp"

#define CADENCE 3600.0 // Seconds between record
#define ROOT "/home/jtdinsmo/Dropbox (MIT)/12.420 project/"

/* Expected format for params.dat:
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

int main(int argc, char* argv[]) {
    // Open files
    std::string filename = "params.dat";
    if (argc >= 2) {
        filename = std::string(argv[1]);
    }

    if (filename.size() < 4 ||
        filename.substr(filename.size() - 4, 4) != ".dat") {
        std::cout << "The input file must be a .dat file." << std::endl;
        return 0;
    }
    std::string bare = filename.substr(0, filename.size() - 4);

    std::ifstream test;
    test.open(filename);
    if (!test.is_open()) {
        std::cout << "The file " << filename << " does not exist." << std::endl;
    }

    Vector3 axis = Vector3::z();
    std::string axis_name = "z";
    if (argc >= 3) {
        if (strcmp(argv[3], "x") == 0) {
            axis = Vector3::x();
            axis_name = "x";
        }
        if (strcmp(argv[3], "y") == 0) {
            axis = Vector3::y();
            axis_name = "y";
        }
    }

    // Load information
    /* Note: since this code is intended to be used in a limited set of
     * scnarios, and since I will always have access to this source code
     * to debug any errors. I will not be writing a "safe parser". That is,
     * when I parse the file `filename`, I will not carefully print out the
     * error messages that could occur. The program will crash somehow and it
     * is up to the user to figure out the problem with the formatting of
     * `filename`. I'm doing this to save time. */

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
    int L = std::stoi(x);// Max degree of spherical harmonics
    int n = std::stoi(y);// Segmentation of cube face at outer shell
    int m = std::stoi(z);// Number of shells

    int num_chunks = n * n * (m + 1) * (2 * m + 1) / m;

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

    std::vector<double> densities;
    std::getline(f, line);
    ss = std::istringstream (line);
    for(int i = 0; i <= num_chunks; i++) {
        ss >> word;
        densities.push_back(std::stod(word));
    }

    // Spin
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x;
    ss >> y;
    ss >> z;
    double spinx = std::stod(x);
    double spiny = std::stod(y);
    double spinz = std::stod(z);

    // Impact parameter
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double impact_parameter = std::stod(word);

    // Velocity (units: times escape velocity)
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double speed = std::stod(word);

    // Central mass
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double central_mass = std::stod(word);

    // Load asteroid
    Asteroid asteroid(L, n, m, clms, densities, spinx, spiny, spinz,
        impact_parameter, speed, central_mass);

    // Draw asteroid
    std::string ast_name = bare + "_" + axis_name + ".ast";
    asteroid.draw(ast_name, axis);

    std::string command = ("python3 \"" ROOT "code/draw/draw.py\" \"" ROOT
        "code/draw/") + ast_name + "\"";
    system(command.c_str());
    return 1;
}
