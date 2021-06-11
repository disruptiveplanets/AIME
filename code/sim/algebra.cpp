#include "algebra.hpp"

Vector3 Vector3::operator*(double d) const {
    return Vector3({entries[0] * d, entries[1] * d, entries[2] * d});
}
Vector3 Vector3::operator+(Vector3 v) const {
    return Vector3({entries[0] + v[0], entries[1] + v[1], entries[2] + v[2]});
}
Vector3 Vector3::operator-(Vector3 v) const {
    return Vector3({entries[0] - v[0], entries[1] - v[1], entries[2] - v[2]});
}
void Vector3::operator+=(Vector3 v) {
    entries[0] += v[0];
    entries[1] += v[1];
    entries[2] += v[2];
}
void Vector3::operator-=(Vector3 v) {
    entries[0] -= v[0];
    entries[1] -= v[1];
    entries[2] -= v[2];
}
double Vector3::dot(Vector3 v1, Vector3 v2)  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
Vector3 Vector3::cross(Vector3 v1, Vector3 v2) {
    return Vector3({v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0]});
}
double Vector3::operator[](int i) const {
    return entries[i];
}
double Vector3::mag() const {
    return sqrt(entries[0] * entries[0] + entries[1] * entries[1] +
        entries[2] * entries[2]);
}
double Vector3::mag2() const {
    return entries[0] * entries[0] + entries[1] * entries[1] +
        entries[2] * entries[2];
}
std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    os << '[' << v[0] << " " << v[1] << " " << v[2] << ']';
    return os;
}


Matrix3 Matrix3::operator*(double d) const {
    return Matrix3({entries[0] * d, entries[1] * d, entries[2] * d,
                    entries[3] * d, entries[4] * d, entries[5] * d,
                    entries[6] * d, entries[7] * d, entries[8] * d});
}
Vector3 Matrix3::operator*(Vector3 v) const {
    return Vector3({entries[0] * v[0] + entries[1] * v[1] + entries[2] * v[2],
                    entries[3] * v[0] + entries[4] * v[1] + entries[5] * v[2],
                    entries[6] * v[0] + entries[7] * v[1] + entries[8] * v[2]});
}
Matrix3 Matrix3::operator*(Matrix3 m) const {
    return Matrix3(
        {entries[0] * m(0, 0) + entries[1] * m(1, 0) + entries[2] * m(2, 0),
         entries[0] * m(0, 1) + entries[1] * m(1, 1) + entries[2] * m(2, 1),
         entries[0] * m(0, 2) + entries[1] * m(1, 2) + entries[2] * m(2, 2),
         
         entries[3] * m(0, 0) + entries[4] * m(1, 0) + entries[5] * m(2, 0),
         entries[3] * m(0, 1) + entries[4] * m(1, 1) + entries[5] * m(2, 1),
         entries[3] * m(0, 2) + entries[4] * m(1, 2) + entries[5] * m(2, 2),
         
         entries[6] * m(0, 0) + entries[7] * m(1, 0) + entries[8] * m(2, 0),
         entries[6] * m(0, 1) + entries[7] * m(1, 1) + entries[8] * m(2, 1),
         entries[6] * m(0, 2) + entries[7] * m(1, 2) + entries[8] * m(2, 2),});
}
Matrix3 Matrix3::operator+(Matrix3 m) const {
    return Matrix3(
        {entries[0] + m(0, 0), entries[1] + m(0, 1), entries[2] + m(0, 2),
         entries[3] + m(1, 0), entries[4] + m(1, 1), entries[5] + m(1, 2),
         entries[6] + m(2, 0), entries[7] + m(2, 1), entries[8] + m(2, 2),});
}
Matrix3 Matrix3::operator-(Matrix3 m) const {
    return Matrix3(
        {entries[0] - m(0, 0), entries[1] - m(0, 1), entries[2] - m(0, 2),
         entries[3] - m(1, 0), entries[4] - m(1, 1), entries[5] - m(1, 2),
         entries[6] - m(2, 0), entries[7]- m(2, 1), entries[8] - m(2, 2),});

}
double Matrix3::operator()(int i, int j) const {
    return entries[3*i + j];
}
Matrix3 Matrix3::transpose() const {
    return Matrix3(
        {entries[0], entries[3], entries[6],
         entries[1], entries[4], entries[7],
         entries[2], entries[5], entries[8],});
}
Matrix3 Matrix3::in_transpose() {
    double temp[3] = {entries[3], entries[6], entries[7]};
    entries[3] = entries[1];
    entries[6] = entries[2];
    entries[7] = entries[5];
    entries[1] = temp[0];
    entries[2] = temp[1];
    entries[5] = temp[2];
}
double Matrix3::det() const {
    return entries[0] * entries[4] * entries[8] + 
        entries[1] * entries[5] + entries[6] + 
        entries[2] * entries[3] * entries[7] - 
        entries[6] * entries[4] * entries[2] - 
        entries[7] * entries[5] * entries[0] - 
        entries[8] * entries[1] * entries[3];
}
std::ostream& operator<<(std::ostream& os, const Matrix3& m) {
    os << '[' << m(0, 0) << " " << m(0, 1) << " " << m(0, 2) << "\n " << 
          m(1, 0) << " " << m(1, 1) << " " << m(1, 2) << "\n " << 
          m(2, 0) << " " << m(2, 1) << " " << m(2, 2) << ']';
    return os;
}
Matrix3 Matrix3::rotation_x(double theta) {
    return Matrix3({
        1, 0, 0,
        0, cos(theta), -sin(theta),
        0, sin(theta), cos(theta)
    });
}
Matrix3 Matrix3::rotation_y(double theta) {
    return Matrix3({
        cos(theta), 0, sin(theta),
        0, 1, 0,
        -sin(theta), 0, cos(theta)
    });
}
Matrix3 Matrix3::rotation_z(double theta) {
    return Matrix3({
        cos(theta), -sin(theta), 0,
        sin(theta), cos(theta), 0,
        0, 0, 1,
    });
}