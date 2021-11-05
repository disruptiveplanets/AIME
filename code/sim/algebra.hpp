#pragma once

#include <array>
#include <math.h>
#include <iostream>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "wignerd.hpp"

#define PI 3.141592653589793438462643383
#define EPSILON 1e-10

using cdouble = std::complex<double>;

class Vector3;

class Vector3 {
public:
    Vector3() {}
    Vector3(std::array<double, 3> e) : e(e), invalid_mag(true) {}

    Vector3 operator*(double d) const;
    Vector3 operator/(double d) const;
    Vector3 operator+(Vector3 const& v) const;
    Vector3 operator-(Vector3 const& v) const;
    Vector3 operator-() const;
    void operator*=(double d);
    void operator+=(Vector3 const& v);
    void operator-=(Vector3 const& v);
    void operator/=(double d);
    static double dot(Vector3 const& v1, Vector3 const& v2);
    static Vector3 cross(Vector3 const& v1, Vector3 const& v2);
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
    Vector3 operator*(Vector3 const& v) const;
    Matrix3 operator*(Matrix3 const& m) const;
    Matrix3 operator+(Matrix3 const& m) const;
    Matrix3 operator-(Matrix3 const& m) const;
    void operator+=(Matrix3 const& m);
    void operator-=(Matrix3 const& m);
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


class Quaternion {
public:
    Quaternion() {};
    Quaternion(double r, double i, double j, double k) : r_(r), i_(i), j_(j),
        k_(k), invalid_mat(true) {}

    static Quaternion zero() {return Quaternion(0, 0, 0, 0); }
    static Quaternion identity() {return Quaternion(1, 0, 0, 0); }

    double r() const { return r_; }
    double i() const { return i_; }
    double j() const { return j_; }
    double k() const { return k_; }
    Quaternion operator+(Quaternion const& q) const;
    Quaternion operator-(Quaternion const& q) const;
    Quaternion operator*(double d) const;
    Quaternion operator/(double d) const;
    void operator+=(Quaternion const& q);
    void operator-=(Quaternion const& q);
    void operator*=(double d);
    void operator/=(double d);
    Quaternion operator*(Quaternion const& q) const;
    std::array<double, 3> euler_angles() const;
    bool is_nan() const {
        return isnan(r_) || isnan(i_) || isnan(j_) || isnan(k_);
    }

    Quaternion inverse() const;
    Matrix3 matrix();
    Vector3 vec() const;
    double mag() const;
    double mag2() const;
    Quaternion conjugate() const { return Quaternion(r_, -i_, -j_, -k_); }

private:
    double r_;
    double i_;
    double j_;
    double k_;
    Matrix3 mat;
    bool invalid_mat;
};


std::ostream& operator<<(std::ostream& os, const Matrix3& dt);
std::ostream& operator<<(std::ostream& os, const Vector3& dt);
std::ostream& operator<<(std::ostream& os, const Quaternion& dt);
Vector3 operator*(double d, Vector3 const& v);
Matrix3 operator*(double d, Matrix3 const& m);
Quaternion operator*(double d, Quaternion const& q);
void force_diagonal(std::array<Vector3, 3>& vs);

class cdouble {
    
}
