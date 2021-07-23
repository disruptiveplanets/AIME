#include <fstream>
#include <string>
#include <iostream>
#include <chrono>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"

#define CADENCE 3600.0 // Seconds between record

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
    uint maxjl = std::stoi(x);
    uint maxkl = std::stoi(y);

    // Jlms
    std::vector<cdouble> halfjlms;// 22, 21, 20, 33, 32, 31, 30
    for (uint l = 0; l <= maxjl; l++){
        std::getline(f, line);
        ss = std::istringstream (line);
        for (int m = -l; m < (int)l; m+=2) {
            ss >> x >> y;
            halfjlms.push_back({std::stod(x), std::stod(y)});
        }
        ss >> x;
        halfjlms.push_back({std::stod(x), 0});
    }

    // Mlms
    std::vector<cdouble> halfklms;
    for (uint l = 0; l <= maxkl; l++){
        std::getline(f, line);
        ss = std::istringstream (line);
        for (int m = -l; m < (int)l; m+=2) {
            ss >> x >> y;
            halfklms.push_back({std::stod(x), std::stod(y)});
        }
        ss >> x;
        halfklms.push_back({std::stod(x), 0});
    }

    // Spin
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x >> y >> z;
    Vector3 spin = Vector3({std::stod(x), std::stod(y), std::stod(z)});

    // Initial roll
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double initial_roll = std::stod(word);

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

    std::vector<cdouble> jlms, klms;
    for (uint l = 0; l <= maxjl; l++) {
        for (int m = -l; m <= (int)l; m++) {
            if (m < 0) {
                jlms.push_back(std::conj(halfjlms[l * (l + 1) / 2 + l - abs(m)])
                    * (double)parity(m));
            }
            else {
                jlms.push_back(halfjlms[l * (l + 1) / 2 + l - abs(m)]);
            }
        }
    }
    for (uint l = 0; l <= maxkl; l++) {
        for (int m = -l; m <= (int)l; m++) {
            if (m < 0) {
                klms.push_back(std::conj(halfklms[l * (l + 1) / 2 + l - abs(m)])
                    * (double)parity(m));
            }
            else {
                klms.push_back(halfklms[l * (l + 1) / 2 + l - abs(m)]);
            }
        }
    }


    // Load asteroid
    Asteroid asteroid(jlms, klms, spin, initial_roll,
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
