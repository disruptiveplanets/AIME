#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <sstream>

#include "algebra.hpp"

#define G 6.67408e-11

using uint = unsigned int;

class Triangle {
public:
    Triangle(Vector3 corner, Vector3 l1, Vector3 l2) : corner(corner), l1(l1), 
    l2(l2) {
        norm = Vector3::cross(l1, l2);
    }

private:
    Vector3 corner;
    Vector3 l1;
    Vector3 l2;
    Vector3 norm;
};

class Asteroid {
public:
    Asteroid(std::string filename);

    void simulate(std::ofstream&& resolved, std::ofstream&& unresolved);

private:
    void generate_shape();
    void calculate_moi();
    void calculate_com();
    void calculate_mass();
    void update();
    void update_position(double dt);
    void update_orientation(double dt);

private:
    uint L; // Max degree of harmonics
    uint n; // Number of cube subdivisions on outermost shell
    uint m; // Number of shells
    std::vector<double> clms;
    std::vector<double> densities;
    double mass;
    Matrix3 moi;
    std::vector<Vector3> points;
    Vector3 position;
    Vector3 spin;
    double mu;
};