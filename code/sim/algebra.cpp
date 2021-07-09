#include "algebra.hpp"

Vector3 Vector3::operator*(double d) const {
    return Vector3({e[0] * d, e[1] * d, e[2] * d});
}
Vector3 Vector3::operator/(double d) const {
    return Vector3({e[0] / d, e[1] / d, e[2] / d});
}
Vector3 Vector3::operator+(Vector3 const& v) const {
    return Vector3({e[0] + v[0], e[1] + v[1], e[2] + v[2]});
}
Vector3 Vector3::operator-(Vector3 const& v) const {
    return Vector3({e[0] - v[0], e[1] - v[1], e[2] - v[2]});
}
Vector3 Vector3::operator-() const {
    return Vector3({-e[0], -e[1], -e[2]});
}
void Vector3::operator*=(double d) {
    e[0] *= d;
    e[1] *= d;
    e[2] *= d;
    invalid_mag = true;
}
void Vector3::operator+=(Vector3 const& v) {
    e[0] += v[0];
    e[1] += v[1];
    e[2] += v[2];
    invalid_mag = true;
}
void Vector3::operator-=(Vector3 const& v) {
    e[0] -= v[0];
    e[1] -= v[1];
    e[2] -= v[2];
    invalid_mag = true;
}
void Vector3::operator/=(double d) {
    e[0] /= d;
    e[1] /= d;
    e[2] /= d;
    invalid_mag = true;
}
double Vector3::dot(Vector3 const& v1, Vector3 const& v2)  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
Vector3 Vector3::cross(Vector3 const& v1, Vector3 const& v2) {
    return Vector3({v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0]});
}
double Vector3::operator[](int i) const {
    return e[i];
}
double Vector3::mag() {
    if (invalid_mag) {
        mag_ = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        invalid_mag = false;
    }
    return mag_;
}
double Vector3::mag2() const {
    return e[0] * e[0] + e[1] * e[1] +
        e[2] * e[2];
}
std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    os << '[' << v[0] << " " << v[1] << " " << v[2] << ']';
    return os;
}


Matrix3 Matrix3::operator*(double d) const {
    return Matrix3({e[0] * d, e[1] * d, e[2] * d,
                    e[3] * d, e[4] * d, e[5] * d,
                    e[6] * d, e[7] * d, e[8] * d});
}
Matrix3 Matrix3::operator/(double d) const {
    return Matrix3({e[0] / d, e[1] / d, e[2] / d,
                    e[3] / d, e[4] / d, e[5] / d,
                    e[6] / d, e[7] / d, e[8] / d});
}
Vector3 Matrix3::operator*(Vector3 const& v) const {
    return Vector3({e[0] * v[0] + e[1] * v[1] + e[2] * v[2],
                    e[3] * v[0] + e[4] * v[1] + e[5] * v[2],
                    e[6] * v[0] + e[7] * v[1] + e[8] * v[2]});
}
Matrix3 Matrix3::operator*(Matrix3 const& m) const {
    return Matrix3(
        {e[0] * m(0, 0) + e[1] * m(1, 0) + e[2] * m(2, 0),
         e[0] * m(0, 1) + e[1] * m(1, 1) + e[2] * m(2, 1),
         e[0] * m(0, 2) + e[1] * m(1, 2) + e[2] * m(2, 2),

         e[3] * m(0, 0) + e[4] * m(1, 0) + e[5] * m(2, 0),
         e[3] * m(0, 1) + e[4] * m(1, 1) + e[5] * m(2, 1),
         e[3] * m(0, 2) + e[4] * m(1, 2) + e[5] * m(2, 2),

         e[6] * m(0, 0) + e[7] * m(1, 0) + e[8] * m(2, 0),
         e[6] * m(0, 1) + e[7] * m(1, 1) + e[8] * m(2, 1),
         e[6] * m(0, 2) + e[7] * m(1, 2) + e[8] * m(2, 2),});
}
Matrix3 Matrix3::operator+(Matrix3 const& m) const {
    return Matrix3(
        {e[0] + m(0, 0), e[1] + m(0, 1), e[2] + m(0, 2),
         e[3] + m(1, 0), e[4] + m(1, 1), e[5] + m(1, 2),
         e[6] + m(2, 0), e[7] + m(2, 1), e[8] + m(2, 2),});
}
Matrix3 Matrix3::operator-(Matrix3 const& m) const {
    return Matrix3(
        {e[0] - m(0, 0), e[1] - m(0, 1), e[2] - m(0, 2),
         e[3] - m(1, 0), e[4] - m(1, 1), e[5] - m(1, 2),
         e[6] - m(2, 0), e[7] - m(2, 1), e[8] - m(2, 2),});

}
void Matrix3::operator+=(Matrix3 const& m) {
    e[0] += m(0, 0); e[1] += m(0, 1); e[2] += m(0, 2);
    e[3] += m(1, 0); e[4] += m(1, 1); e[5] += m(1, 2);
    e[6] += m(2, 0); e[7] += m(2, 1); e[8] += m(2, 2);
}
void Matrix3::operator-=(Matrix3 const& m) {
    e[0] -= m(0, 0); e[1] -= m(0, 1); e[2] -= m(0, 2);
    e[3] -= m(1, 0); e[4] -= m(1, 1); e[5] -= m(1, 2);
    e[6] -= m(2, 0); e[7] -= m(2, 1); e[8] -= m(2, 2);
}
double Matrix3::operator()(int i, int j) const {
    return e[3*i + j];
}
Matrix3 Matrix3::transpose() const {
    return Matrix3(
        { e[0], e[3], e[6],
          e[1], e[4], e[7],
          e[2], e[5], e[8],});
}
void Matrix3::in_transpose() {
    double temp[3] = {e[3], e[6], e[7]};
    e[3] = e[1];
    e[6] = e[2];
    e[7] = e[5];
    e[1] = temp[0];
    e[2] = temp[1];
    e[5] = temp[2];
}
double Matrix3::det() const {
    return e[0] * e[4] * e[8] + e[1] * e[5] * e[6] + e[2] * e[3] * e[7]
         - e[6] * e[4] * e[2] - e[7] * e[5] * e[0] - e[8] * e[1] * e[3];
}
double Matrix3::trace() const {
    return e[0] + e[4] + e[8];
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
Matrix3 Matrix3::symmetric_invert() {
    std::array<double, 3> evals = get_evals();
    std::array<Vector3, 3> evecs = get_symmetric_evecs(evals);
    Matrix3 orthog = Matrix3({
        evecs[0][0], evecs[1][0], evecs[2][0],
        evecs[0][1], evecs[1][1], evecs[2][1],
        evecs[0][2], evecs[1][2], evecs[2][2],
    });
    Matrix3 diag = Matrix3::diag(1.0/evals[0], 1.0/evals[1], 1.0/evals[2]);
    return orthog *(diag * orthog.transpose());
}
Matrix3 Matrix3::symmetric_invert(std::array<double, 3>const& evals,
    std::array<Vector3, 3>const& evecs) {
    Matrix3 orthog = Matrix3({
        evecs[0][0], evecs[1][0], evecs[2][0],
        evecs[0][1], evecs[1][1], evecs[2][1],
        evecs[0][2], evecs[1][2], evecs[2][2],
    });
    Matrix3 diag = Matrix3::diag(1.0/evals[0], 1.0/evals[1], 1.0/evals[2]);
    return orthog *(diag * orthog.transpose());
}
Matrix3 Matrix3::symmetric_reconstruct(std::array<double, 3>const& evals,
    std::array<Vector3, 3>const& evecs) {
    Matrix3 orthog = Matrix3({
        evecs[0][0], evecs[1][0], evecs[2][0],
        evecs[0][1], evecs[1][1], evecs[2][1],
        evecs[0][2], evecs[1][2], evecs[2][2],
    });
    Matrix3 diag = Matrix3::diag(evals[0], evals[1], evals[2]);
    return orthog *(diag * orthog.transpose());
}
std::array<double, 3> Matrix3::get_evals() const {
    // Get the coefficients of the characteristic polynomial
    // a = -1
    double b = trace();
    double c = e[1] * e[3] + e[2] * e[6] + e[5] * e[7]
             - e[0] * e[4] - e[4] * e[8] - e[0] * e[8];
    double d = det();

    double p = b / 3.0;
    double q = p*p*p + (b*c + 3*d) / 6.0;
    std::complex<double> R = -c / 3.0;

    std::complex<double> r_c = pow(q + std::sqrt(q*q + pow(R - p*p, 3)), 1/3.0)
             + pow(q - std::sqrt(q*q + pow(R - p*p, 3)), 1/3.0) + p;

    double r = r_c.real();
    double disc = (b - r) * (b - r) - 4 * d / r;
    if (disc < EPSILON) { disc = 0; }
    double r1 = (b - r + sqrt(disc)) / 2;
    double r2 = (b - r - sqrt(disc)) / 2;
    return {r, r1, r2};
}
std::array<Vector3, 3> Matrix3::get_symmetric_evecs(std::array<double, 3>const& evals)
const {
    // Solve eigenvalue problem with Cramer's rule
    Vector3 evec1;
    int fix_column = -1;
    int a, b;// Not fixed columns
    double D = 0, Da = 0, Db = 0;
    while (abs(D) < EPSILON) {
        // We fix one entry of the eigenvalue, but only one which can be
        // chosen to be nonzero. If it must be zero, D will be zero.
        fix_column++;
        switch(fix_column) {
        case 0:
            a = 1;
            b = 2;
            break;
        case 1:
            a = 0;
            b = 2;
            break;
        case 2:
            a = 0;
            b = 1;
            break;
        default:
            // Matrix is diagonal
            return {Vector3::x(), Vector3::y(), Vector3::z()};
                // Return error result
        }
        D = ((*this)(a, a) - evals[0]) * ((*this)(b, b) - evals[0]) -
            (*this)(b, a) * (*this)(a, b);
        Da = -(*this)(a, fix_column) * ((*this)(b, b) - evals[0]) +
            (*this)(b, fix_column) * (*this)(a, b);
        Db = ((*this)(a, a) - evals[0]) * -(*this)(b, fix_column) +
            (*this)(b, a) * (*this)(a, fix_column);
    }

    // Write the eigenvalue to the array
    switch(fix_column) {
    case 0:
        evec1 = Vector3({1, Da/D, Db/D});
        break;
    case 1:
        evec1 = Vector3({Da/D, 1, Db/D});
        break;
    case 2:
        evec1 = Vector3({Da/D, Db/D, 1});
        break;
    }

    // Establish evec2 to be orthogonal.
    Vector3 evec2;
    fix_column = -1;
    D = 0;
    while (abs(D) < EPSILON) {
        // We fix one entry of the eigenvalue, but only one which can be
        // chosen to be nonzero. If it must be zero, D will be zero.
        fix_column++;
        switch(fix_column) {
        case 0:
            a = 1;
            b = 2;
            break;
        case 1:
            a = 0;
            b = 2;
            break;
        case 2:
            a = 0;
            b = 1;
            break;
        default:
            return {Vector3::x(), Vector3::y(), Vector3::z()};
        }
        D = ((*this)(a, a) - evals[1]) * ((*this)(b, b) - evals[1]) -
            (*this)(b, a) * (*this)(a, b);
        Da = -(*this)(a, fix_column) * ((*this)(b, b) - evals[1]) +
            (*this)(b, fix_column) * (*this)(a, b);
    }
    // Fix evec2 to be perpendicular
    switch(fix_column) {
    case 0:
        evec2 = Vector3({1, Da/D, 0});
        evec2 -= Vector3::dot(evec2, evec1) / evec1[2] * Vector3::z();
        break;
    case 1:
        evec2 = Vector3({Da/D, 1, 0});
        evec2 -= Vector3::dot(evec2, evec1) / evec1[2] * Vector3::z();
        break;
    case 2:
        evec2 = Vector3({Da/D, 0, 1});
        evec2 -= Vector3::dot(evec2, evec1) / evec1[1] * Vector3::y();
        break;
    }

    // Generate the third eigenvector
    Vector3 evec3 = Vector3::cross(evec1, evec2);

    // Normalize for cleanliness
    return {evec1 / evec1.mag(), evec2 / evec2.mag(), evec3 / evec3.mag()};
}


Vector3 operator*(double d, Vector3 const& v) {
    return Vector3({d * v[0], d * v[1], d * v[2]});
}
Matrix3 operator*(double d, Matrix3 const& m) {
    return Matrix3({d * m(0, 0), d * m(0, 1), d * m(0, 2),
                    d * m(1, 0), d * m(1, 1), d * m(1, 2),
                    d * m(2, 0), d * m(2, 1), d * m(2, 2)});
}



Matrix3 Quaternion::matrix() {
    if (invalid_mat) {
        double s = 1/mag2();
        mat = Matrix3({
            1 - 2*s*(j_*j_ + k_*k_), 2*s*(i_*j_ - k_*r_), 2*s*(i_*k_ + j_*r_),
            2*s*(i_*j_ + k_*r_), 1 - 2*s*(i_*i_ + k_*k_), 2*s*(j_*k_ - i_*r_),
            2*s*(i_*k_ - j_*r_), 2*s*(j_*k_ + i_*r_), 1 - 2*s*(i_*i_ + j_*j_)
        });
        invalid_mat = false;
    }
    return mat;
}
Quaternion Quaternion::inverse() const {
    return conjugate() / mag2();
}
Vector3 Quaternion::vec() const {
    return Vector3({i_, j_, k_});
}
Quaternion Quaternion::operator+(Quaternion const& q) const {
    return Quaternion(r_+q.r(), i_+q.i(), j_+q.j(), k_+q.k());
}
Quaternion Quaternion::operator-(Quaternion const& q) const {
    return Quaternion(r_-q.r(), i_-q.i(), j_-q.j(), k_-q.k());
}
Quaternion Quaternion::operator*(double d) const {
    return Quaternion(r_*d, i_*d, j_*d, k_*d);
}
Quaternion Quaternion::operator/(double d) const {
    return Quaternion(r_/d, i_/d, j_/d, k_/d);
}
void Quaternion::operator+=(Quaternion const& q) {
    r_ += q.r();
    i_ += q.i();
    j_ += q.j();
    k_ += q.k();
    invalid_mat = true;
}
void Quaternion::operator-=(Quaternion const& q) {
    r_ -= q.r();
    i_ -= q.i();
    j_ -= q.j();
    k_ -= q.k();
    invalid_mat = true;
}
void Quaternion::operator*=(double d) {
    r_ *= d;
    i_ *= d;
    j_ *= d;
    k_ *= d;
    invalid_mat = true;
}
void Quaternion::operator/=(double d) {
    r_ /= d;
    i_ /= d;
    j_ /= d;
    k_ /= d;
    invalid_mat = true;
}
Quaternion Quaternion::operator*(Quaternion const& q) const {
    Vector3 mv = vec();
    Vector3 qv = q.vec();
    double d = r_ * q.r() - Vector3::dot(mv, qv);
    Vector3 v = r_ * qv + q.r() * mv + Vector3::cross(mv, qv);
    return Quaternion(d, v[0], v[1], v[2]);
}
double Quaternion::mag() const {
    return sqrt(r_*r_+i_*i_+j_*j_+k_*k_);
}
double Quaternion::mag2() const {
    return r_*r_+i_*i_+j_*j_+k_*k_;
}

std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
    os << q.r() << " + " << q.i() << " i + " << q.j() << " j + " << q.k();
    return os;
}
Quaternion operator*(double d, Quaternion const& q) {
    return Quaternion(d*q.r(), d*q.i(), d*q.j(), d*q.k());
}

std::array<double, 3> Quaternion::euler_angles() const {
    // z-y-z
    return {atan2(j_*k_ - i_*r_, i_*k_ + j_*r_),
        acos(1 - 2*(i_*i_ + j_*j_)/mag2()),
        atan2(j_*k_ + i_*r_, -(i_*k_ - j_*r_))};
}
