#pragma once

#include <array>
#include <math.h>

#define pi 3.141592653589793438462643383

class Vector3;

class Vector3 {
public:
    Vector3() {}
    Vector3(std::array<double, 3> entries) : entries(entries) {}

    Vector3 operator*(double d) const;
    Vector3 operator+(Vector3 v) const;
    Vector3 operator-(Vector3 v) const;
    void operator+=(Vector3 v);
    void operator-=(Vector3 v);
    static double dot(Vector3 v1, Vector3 v2);
    static Vector3 cross(Vector3 v1, Vector3 v2);
    double operator[](int i) const;
    double mag() const;
    double mag2() const;

private:
    std::array<double, 3> entries;
};

class Matrix3 {
public:
    Matrix3() {}
    Matrix3(std::array<double, 9> entries) : entries(entries) {}

    Matrix3 operator*(double d) const;
    Vector3 operator*(Vector3 v) const;
    Matrix3 operator*(Matrix3 m) const;
    Matrix3 operator+(Matrix3 m) const;
    Matrix3 operator-(Matrix3 m) const;
    double operator()(int i, int j) const;
    Matrix3 transpose() const;
    Matrix3 in_transpose();
    double det() const;

    static Matrix3 rotation_x(double theta);
    static Matrix3 rotation_y(double theta);
    static Matrix3 rotation_z(double theta);

private:
    std::array<double, 9> entries;
};

std::ostream& operator<<(std::ostream& os, const Matrix3& dt);
std::ostream& operator<<(std::ostream& os, const Vector3& dt);