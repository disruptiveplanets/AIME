#include <fstream>
#include <string>
#include <iostream>
#include <chrono>

#include "backend.hpp"
#include "algebra.hpp"
#include "triangle.hpp"

#define CADENCE 3600.0 // Seconds between record

int main(int argc, char* argv[]) {
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

    Asteroid asteroid(filename);
    if (argc >= 4) {
        if (strcmp(argv[2], "draw") == 0) {
            Vector3 axis = Vector3::z();
            if (strcmp(argv[3], "x") == 0) {
                axis = Vector3::x();
            }
            if (strcmp(argv[3], "y") == 0) {
                axis = Vector3::y();
            }
            asteroid.draw(bare + ".ast", axis);
            return 1;
        }
        filename = std::string(argv[1]);
    }
    std::ofstream resolved, unresolved;
    resolved.open(bare + "-resolved.dat");
    unresolved.open(bare + "-unresolved.dat");

    auto start = std::chrono::high_resolution_clock::now();

    int frames = asteroid.simulate(CADENCE, std::move(resolved), std::move(unresolved));

    auto stop = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(
        stop - start).count() / 1000.0;

    std::cout << "Simulation took " << frames << " frames or "
        << time_taken << " s." << std::endl;


    return 1;
}
