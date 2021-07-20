#include <fstream>
#include <string>
#include <iostream>
#include <chrono>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"
#include "../sim/triangle.hpp"

#define CADENCE 3600.0 // Seconds between record

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

    // Max ls
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x;
    ss >> y;
    int maxjl = std::stoi(x);
    int maxml = std::stoi(y);

    // Jlms
    std::vector<cdouble> jlms;
    double re = 0;
    double im = 0;
    for (int l = 0; l <= maxjl; l++){
        bool written = true;
        std::getline(f, line);
        ss = std::istringstream (line);
        for (int m = -l; m <= l; m++) {
            ss >> word;
            if (written) {
                re = std::stod(word);
            }
            else{
                im = std::stod(word);
                jlms.push_back({re, im});
            }
            written = !written;
        }
        jlms.push_back({re, 0});
    }

    // Mlms
    std::vector<cdouble> mlms;
    re = 0;
    im = 0;
    for (int l = 2; l <= maxml; l++){
        bool written = true;
        std::getline(f, line);
        ss = std::istringstream (line);
        for (int m = -l; m <= l; m++) {
            ss >> word;
            if (written) {
                re = std::stod(word);
            }
            else {
                im = std::stod(word);
                mlms.push_back({re, im});
            }
            written = !written;
        }
        mlms.push_back({re, 0});
    }

    // Spin
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x;
    double spin = std::stod(x);

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

    for(auto m : mlms) {
        std::cout << m << std::endl;
    }


    // Load asteroid
    Asteroid asteroid(jlms, mlms, spin,
        impact_parameter, speed);

    // Run asteroid
    std::vector<double> resolved_data;
    auto start = std::chrono::high_resolution_clock::now();

    int frames = asteroid.simulate(CADENCE, resolved_data);

    auto stop = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(
        stop - start).count() / 1000.0;

    std::cout << "Simulation took " << frames << " frames or "
        << time_taken << " s." << std::endl;

    std::ofstream resolved_file;
    resolved_file.open(bare + "-resolved.dat");
    for (uint i = 0; i < resolved_data.size(); i+=3) {
        resolved_file << resolved_data[i] << ' ';
        resolved_file << resolved_data[i + 1] << ' ';
        resolved_file << resolved_data[i + 2] << std::endl;
    }

    return 1;
}
