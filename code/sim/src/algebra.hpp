#pragma once

#include <array>
#include <math.h>
#include <iostream>
#include <complex>

#define PI 3.141592653589793438462643383
#define EPSILON 1e-10

class Vector3;

class Vector3 {
public:
    Vector3() {}
    Vector3(std::array<double, 3> e) : e(e), invalid_mag(true) {}

    Vector3 operator*(double d) const;
    Vector3 operator/(double d) const;
    Vector3 operator+(Vector3 v) const;
    Vector3 operator-(Vector3 v) const;
    Vector3 operator-() const;
    void operator*=(double v);
    void operator+=(Vector3 v);
    void operator-=(Vector3 v);
    static double dot(Vector3 v1, Vector3 v2);
    static Vector3 cross(Vector3 v1, Vector3 v2);
    double operator[](int i) const;
    double mag();
    double mag2() const;

    static Vector3 x() { return Vector3({1, 0, 0}); }
    static Vector3 y() { return Vector3({0, 1, 0}); }
    static Vector3 z() { return Vector3({0, 0, 1}); }
    static Vector3 zero() { return Vector3({0, 0, 0}); }

private:
    std::array<double, 3> e;
    double mag_;// Memoize mag
    bool invalid_mag;
};

class Matrix3 {
public:
    Matrix3() {}
    Matrix3(std::array<double, 9> e) : e(e) {}

    Matrix3 operator*(double d) const;
    Matrix3 operator/(double d) const;
    Vector3 operator*(Vector3 v) const;
    Matrix3 operator*(Matrix3 m) const;
    Matrix3 operator+(Matrix3 m) const;
    Matrix3 operator-(Matrix3 m) const;
    void operator+=(Matrix3 m);
    void operator-=(Matrix3 m);
    double operator()(int i, int j) const;
    Matrix3 transpose() const;
    void in_transpose();
    double det() const;
    double trace() const;
    std::array<double, 3> get_evals() const;
    std::array<Vector3, 3> get_symmetric_evecs(std::array<double, 3>const& evals) const;
    Matrix3 symmetric_invert();
    static Matrix3 symmetric_invert(std::array<double, 3>const& evals,
        std::array<Vector3, 3>const& evecs);
    static Matrix3 symmetric_reconstruct(std::array<double, 3> const& evals,
        std::array<Vector3, 3>const& evecs);

    static Matrix3 rotation_x(double theta);
    static Matrix3 rotation_y(double theta);
    static Matrix3 rotation_z(double theta);
    static Matrix3 identity() { return Matrix3({1, 0, 0, 0, 1, 0, 0, 0, 1}); }
    static Matrix3 diag(double a, double b, double c) { return Matrix3({
        a, 0, 0,
        0, b, 0,
        0, 0, c,
    }); }

private:
    std::array<double, 9> e;
};

std::ostream& operator<<(std::ostream& os, const Matrix3& dt);
std::ostream& operator<<(std::ostream& os, const Vector3& dt);
Vector3 operator*(double d, Vector3 const& v);
Matrix3 operator*(double d, Matrix3 const& m);
void force_diagonal(std::array<Vector3, 3>& vs);
