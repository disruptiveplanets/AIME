#include <fstream>
#include <string>
#include <iostream>
#include <chrono>
#include <sstream>

#include "../sim/backend.h"
#include "../sim/algebra.h"
#include "../sim/triangle.h"
#include "../sim/linkedlist.h"

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
std::ofstream resolved, unresolved;

int main(int argc, char* argv[]) {
    // Check the input file
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


    /* Note: since this code is intended to be used in a limited set of
     * scnarios, and since I will always have access to this source code
     * to debug any errors. I will not be writing a "safe parser". That is,
     * when I parse the file `filename`, I will not carefully print out the
     * error messages that could occur. The program will crash somehow and it
     * is up to the user to figure out the problem with the formatting of
     * `filename`. I'm doing this to save time.
    */

    // Get all the asteroid characteristics
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

    Cons* clms_head = NULL;
    Cons* clms = NULL;
    Cons* densities_head = NULL;
    Cons* densities = NULL;

    // Clms
    for(int i = 0; i <= L; i++){
        std::getline(f, line);
        ss = std::istringstream (line);
        for(int m = -i; m <= i; m++) {
            ss >> word;
            double* d = new double(std::stod(word));
            cons_push_inc(&clms, d);
            if (clms_head == NULL) {
                clms_head = clms;
            }
        }
    }

    // Densities
    int num_chunks = n * n * (m + 1) * (2 * m + 1) / m;
    std::getline(f, line);
    ss = std::istringstream (line);
    for(int i = 0; i <= num_chunks; i++) {
        ss >> word;
        double* d = new double(std::stod(word));
        cons_push_inc(&densities, d);
        if (densities_head == NULL) {
            densities_head = densities;
        }
    }

    // Spin
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> x;
    ss >> y;
    ss >> z;
    Vector3 spin = vector3_new(std::stod(x), std::stod(y), std::stod(z));

    // Impact parameter
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double impact_parameter = std::stod(word);

    // Velocity (units: times escape velocity)
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double velocity = std::stod(word);

    // Central mass
    std::getline(f, line);
    ss = std::istringstream (line);
    ss >> word;
    double central_mass = std::stod(word);

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

    Asteroid* asteroid = asteroid_new(L, n, m, clms_head, densities_head,
        spin, impact_parameter, velocity, central_mass);

    resolved.open(bare + "-resolved.dat");
    //unresolved.open(bare + "-unresolved.dat");

    auto start = std::chrono::high_resolution_clock::now();

    int frames = asteroid_simulate(asteroid, CADENCE);

    auto stop = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(
        stop - start).count() / 1000.0;

    std::cout << "Simulation took " << frames << " frames or "
        << time_taken << " s." << std::endl;

    cons_free_all(clms);
    cons_free_all(densities);
    asteroid_free(asteroid);
    free(asteroid);

    return 1;
}

void asteroid_record_resolved(Vector3 spin) {
    resolved << spin.x << ' ' << spin.y << ' ' << spin.z << std::endl;
}
void asteroid_record_unresolved() {}
