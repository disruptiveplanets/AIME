#include <fstream>
#include <string>
#include <iostream>

#include "backend.cpp"
#include "algebra.cpp"

int main(int argc, char* argv[]) {
    std::string filename = "params.dat";
    if (argc == 2) {
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
    std::ofstream resolved, unresolved;
    resolved.open(bare + "-resolved.dat");
    unresolved.open(bare + "-unresolved.dat");

    asteroid.simulate(std::move(resolved), std::move(unresolved));
    
    return 1;
}