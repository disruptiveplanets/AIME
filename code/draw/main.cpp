#include <fstream>
#include <string>
#include <iostream>
#include <chrono>

#include "../sim/backend.hpp"
#include "../sim/algebra.hpp"
#include "../sim/triangle.hpp"

#define CADENCE 3600.0 // Seconds between record


void asteroid_draw(Asteroid* a, char* filename, Vector3 axis) {
    // Write data about how to draw an asteroid to a text file (.ast) to be read
    // by a python script and displayed with matplotlib.

    // Project everything along axis. Axis must be normalized
    Vector3 up = vector3_z();
    Vector3 cross = vector3_cross(&up, &axis);
    if (vector3_mag(&cross) < EPSILON) {
        up = vector3_x();
    }
    Vector3 ax_mul = vector3_mul(&axis, vector3_dot(&up, &axis));
    vector3_place_sub(&up, &ax_mul);
    vector3_place_div(&up, vector3_mag(&up));

    std::vector<double> depths;
    std::vector<std::string> lines;
    depths.reserve(cons_size(a->triangles_head));
    lines.reserve(cons_size(a->triangles_head));

    for (Cons* t = triangles_head; t != NULL; t = t->next) {
        if (!t->value->is_edge()){
            continue;
        }
        Vector3* corners = t->get_corners();
        Vector3 center = vector3_mul(vector3_add(
            vector3_add(corners[0], corners[1]), corners[2]), 3);
        Vector3 norm = Vector3::cross(vector3_sub(corners[1], corners[0]),
            vector3_sub(corners[2], corners[0]));
        vector3_place_div(norm, vector3_mag(norm);
        if (vector3_dot(norm, axis) < 0) {
            continue;
        }
        double depth = vector3_dot(center, axis);

        std::string line = "";
        for(int i = 0; i < 3; i++) {
            Vector3 offset = corners[i] -
                vector3_mul(axis, vector3_dot(corners[i], axis));
            if (vector3_mag(offset) == 0){
                line += "0 0 ";
                continue;
            }
            double theta = acos(vector3_dot(offset, up) / vector3_mag(offset));
            if (vector3_dot(vector3_cross(offset, up), axis) < 0) {
                theta = 2 * PI - theta;
            }
            double r = offset.mag();
            line += std::to_string(r * cos(theta)) + " "
                    + std::to_string(r * sin(theta)) + " ";
        }
        line += std::to_string(triangle_get_density(t->value)) + " "
                + std::to_string(vector3_dot(norm, axis));

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
        free(corners);
    }

    std::ofstream output;
    output.open(filename);
    for (std::string& line : lines) {
        output << line << '\n';
    }
    output.close();
}

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

    Vector3 axis = Vector3::z();
    if (argc >= 3) {
        if (strcmp(argv[3], "x") == 0) {
            axis = Vector3::x();
        }
        if (strcmp(argv[3], "y") == 0) {
            axis = Vector3::y();
        }
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
            clms.cons_push_inc(&clms, d);
            if (clms_head == NULL) {
                clms_head = clms;
            }
        }
    }

    // Densities
    std::getline(f, line);
    ss = std::istringstream (line);
    mean_density = 0;
    for(int i = 0; i <= num_chunks; i++) {
        ss >> word;
        mean_density += std::stod(word);
        double* d = new double(std::stod(word));
        densities.cons_push_inc(&densities, d);
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

    asteroid_draw(asteroid);

    cons_free_all(clms);
    cons_free_all(densities);
    asteroid_free(asteroid);

    return 1;
}
