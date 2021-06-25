#include "algebra.h"


Complex complex_sqrt(double s){
    if (s > 0) {
        return (Complex) {sqrt(s), 0};
    }
    return (Complex) {0, sqrt(-s)};
}
Complex complex_pow(Complex c, double p){
    double r = pow(c.r * c.r + c.i * c.i, p/2);
    double phi = atan2(c.i, c.r) * p;
    return (Complex) {r * cos(phi), r * sin(phi)};
}


void vector3_init(Vector3* v, double x, double y, double z) {
    v->x = x;
    v->y = y;
    v->z = z;
    v->invalid_mag = 1;// True
}
Vector3 vector3_new(double x, double y, double z) {
    return (Vector3) {x, y, z, 1, 1};
}
Vector3 vector3_zero() {
    return (Vector3) {0, 0, 0, 1, 1};
}
Vector3 vector3_x() {
    return (Vector3) {1, 0, 0, 1, 1};
}
Vector3 vector3_y() {
    return (Vector3) {0, 1, 0, 1, 1};
}
Vector3 vector3_z() {
    return (Vector3) {0, 0, 1, 1, 1};
}
Vector3 vector3_negative(Vector3 v) {
    return (Vector3) {-v.x, -v.y, -v.z, v.mag, v.invalid_mag};
}

Vector3 vector3_mul(const Vector3 v, double d) {
    return vector3_new(v.x * d, v.y * d, v.z * d);
}
void vector3_place_mul(Vector3* v, double d) {
    v->x *= d;
    v->y *= d;
    v->z *= d;
}
Vector3 vector3_div(const Vector3 v, double d) {
    return vector3_new(v.x / d, v.y / d, v.z / d);
}
void vector3_place_div(Vector3* v, double d) {
    v->x /= d;
    v->y /= d;
    v->z /= d;
}
Vector3 vector3_add(const Vector3 u, const Vector3 v) {
    return vector3_new(u.x + v.x, u.y + v.y, u.z + v.z);
}
void vector3_place_add(Vector3* u, const Vector3 v) {
    u->x += v.x;
    u->y += v.y;
    u->z += v.z;
}
Vector3 vector3_sub(const Vector3 u, const Vector3 v) {
    return vector3_new(u.x - v.x, u.y - v.y, u.z - v.z);
}
void vector3_place_sub(Vector3* u, const Vector3 v) {
    u->x -= v.x;
    u->y -= v.y;
    u->z -= v.z;
}

double vector3_dot(const Vector3 u, const Vector3 v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}
Vector3 vector3_cross(const Vector3 u, const Vector3 v) {
    return vector3_new(u.y * v.z - u.z * v.y,
                       u.z * v.x - u.x * v.z,
                       u.x * v.y - u.y * v.x);
}
Vector3 vector3_mat_mul(const Matrix3 m, const Vector3 v) {
    return vector3_new(
        m.a * v.x + m.b * v.y + m.c * v.z,
        m.d * v.x + m.e * v.y + m.f * v.z,
        m.g * v.x + m.h * v.y + m.i * v.z);
}

double vector3_mag(Vector3 v) {
    if (v.invalid_mag) {
        v.mag = (v.x * v.x + v.y * v.y + v.z * v.z);
        v.invalid_mag = 0;
    }
    return v.mag;
}
double vector3_mag2(const Vector3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

void matrix3_init(Matrix3* m, double a, double b, double c,
    double d, double e, double f,
    double g, double h, double i) {

    m->a = a;
    m->b = b;
    m->c = c;
    m->d = d;
    m->e = e;
    m->f = f;
    m->g = g;
    m->h = h;
    m->i = i;
}
Matrix3 matrix3_zero() {
    return (Matrix3) {0, 0, 0, 0, 0, 0, 0, 0, 0};
}
Matrix3 matrix3_rotation_x(double a) {
    return (Matrix3) {
        1, 0, 0,
        0, cos(a), -sin(a),
        0, sin(a), cos(a)
    };
}
Matrix3 matrix3_rotation_y(double a) {
    return (Matrix3) {
        cos(a), 0, sin(a),
        0, 1, 0,
        -sin(a), 0, cos(a)
    };
}
Matrix3 matrix3_rotation_z(double a) {
    return (Matrix3) {
        cos(a), -sin(a), 0,
        sin(a), cos(a), 0,
        0, 0, 1
    };
}
Matrix3 matrix3_identity() {
    return (Matrix3) {1, 0, 0, 0, 1, 0, 0, 0, 1};
}
Matrix3 matrix3_diag(double a, double b, double c) {
    return (Matrix3) {a, 0, 0, 0, b, 0, 0, 0, c};
}

Matrix3 matrix3_mul(const Matrix3 m, double d) {
    return (Matrix3) {
        m.a * d, m.b * d, m.c * d,
        m.d * d, m.e * d, m.f * d,
        m.g * d, m.h * d, m.i * d};
}
void matrix3_place_mul(Matrix3* m, double d) {
    m->a *= d;
    m->b *= d;
    m->c *= d;
    m->d *= d;
    m->e *= d;
    m->f *= d;
    m->g *= d;
    m->h *= d;
    m->i *= d;
}
Matrix3 matrix3_div(const Matrix3 m, double d) {
    return (Matrix3) {
        m.a / d, m.b / d, m.c / d,
        m.d / d, m.e / d, m.f / d,
        m.g / d, m.h / d, m.i / d};
}
void matrix3_place_div(Matrix3* m, double d) {
    m->a /= d;
    m->b /= d;
    m->c /= d;
    m->d /= d;
    m->e /= d;
    m->f /= d;
    m->g /= d;
    m->h /= d;
    m->i /= d;
}
Matrix3 matrix3_add(const Matrix3 m1, const Matrix3 m2) {
    return (Matrix3) {
        m1.a + m2.a, m1.b + m2.b, m1.c + m2.c,
        m1.d + m2.d, m1.e + m2.e, m1.f + m2.f,
        m1.g + m2.g, m1.h + m2.h, m1.i + m2.i};
}
void matrix3_place_add(Matrix3* m1, const Matrix3 m2) {
    m1->a += m2.a;
    m1->b += m2.b;
    m1->c += m2.c;
    m1->d += m2.d;
    m1->e += m2.e;
    m1->f += m2.f;
    m1->g += m2.g;
    m1->h += m2.h;
    m1->i += m2.i;
}
Matrix3 matrix3_sub(const Matrix3 m1, const Matrix3 m2) {
    return (Matrix3) {
        m1.a - m2.a, m1.b - m2.b, m1.c - m2.c,
        m1.d - m2.d, m1.e - m2.e, m1.f - m2.f,
        m1.g - m2.g, m1.h - m2.h, m1.i - m2.i};
}
void matrix3_place_sub(Matrix3* m1, const Matrix3 m2) {
    m1->a -= m2.a;
    m1->b -= m2.b;
    m1->c -= m2.c;
    m1->d -= m2.d;
    m1->e -= m2.e;
    m1->f -= m2.f;
    m1->g -= m2.g;
    m1->h -= m2.h;
    m1->i -= m2.i;
}
Matrix3 matrix3_transpose(const Matrix3 m) {
    return (Matrix3) {
        m.a, m.d, m.g,
        m.b, m.e, m.h,
        m.c, m.f, m.i};
}
void matrix3_place_transpose(Matrix3* m) {
    double temp1 = m->d;
    double temp2 = m->g;
    double temp3 = m->h;
    m->d = m->b;
    m->h = m->f;
    m->g = m->c;
    m->b = temp1;
    m->c = temp2;
    m->f = temp3;
}

Matrix3 matrix3_mat_mul(const Matrix3 m1, const Matrix3 m2) {
    return (Matrix3) {
        m1.a * m2.a + m1.b * m2.d + m1.c * m2.g,
        m1.a * m2.b + m1.b * m2.e + m1.c * m2.h,
        m1.a * m2.c + m1.b * m2.f + m1.c * m2.i,

        m1.d * m2.a + m1.e * m2.d + m1.f * m2.g,
        m1.d * m2.b + m1.e * m2.e + m1.f * m2.h,
        m1.d * m2.c + m1.e * m2.f + m1.f * m2.i,

        m1.g * m2.a + m1.h * m2.d + m1.i * m2.g,
        m1.g * m2.b + m1.h * m2.e + m1.i * m2.h,
        m1.g * m2.c + m1.h * m2.f + m1.i * m2.i
    };
}

double matrix3_get(const Matrix3 m, int a, int b) {
    switch(a) {
        case 0:
        switch(b) {
            case 0: return m.a;
            case 1: return m.b;
            case 2: return m.c;
        }
        case 1:
        switch(b) {
            case 0: return m.d;
            case 1: return m.e;
            case 2: return m.f;
        }
        case 2:
        switch(b) {
            case 0: return m.g;
            case 1: return m.h;
            case 2: return m.i;
        }
    }
    // Does not reach
    return 0;
}
double matrix3_det(const Matrix3 m) {
        return m.a * m.e * m.i + m.b * m.f * m.g + m.c * m.d * m.h
             - m.g * m.e * m.c - m.h * m.f * m.a - m.i * m.b * m.d;
}
double matrix3_trace(const Matrix3 m) {
    return m.a + m.e + m.i;
}

double* matrix3_get_evals(const Matrix3 m) {
    double* evals = (double*) malloc(sizeof(double) * 3);

    // Get the coefficients of the characteristic polynomial
    // a = -1
    double b = matrix3_trace(m);
    double c = m.b * m.d + m.c * m.g + m.f * m.h
             - m.a * m.d - m.d * m.i - m.a * m.i;
    double d = matrix3_det(m);

    double p = b / 3.0;
    double q = p*p*p + (b*c + 3*d) / 6.0;
    double R = -c / 3.0;

    Complex middle = complex_sqrt(q*q + pow(R - p*p, 3));
    Complex plus = {q + middle.r, middle.i};
    plus = complex_pow(plus, 1/3.0);
    Complex minus = {q - middle.r, middle.i};
    minus = complex_pow(minus, 1/3.0);

    double r = plus.r + minus.r + p;
    double disc = (b - r) * (b - r) - 4 * d / r;
    if (abs(disc) < EPSILON) { disc = 0; }

    evals[0] = r;
    evals[1] = (b - r + sqrt(disc)) / 2;
    evals[2] = (b - r - sqrt(disc)) / 2;

    return evals;
}
Vector3* matrix3_get_symmetric_evecs(const Matrix3 m, const double* evals) {
    // Solve eigenvalue problem with Cramer's rule
    Vector3* evecs = (Vector3*) malloc(3 * sizeof(Vector3));
    int fix_column = -1;
    int a, b;// Not fixed columns
    double D = 0;
    double Da = 0;
    double Db = 0;
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
            evecs[0] = vector3_x();
            evecs[1] = vector3_y();
            evecs[2] = vector3_z();
            return evecs;
        }
        D = (matrix3_get(m, a, a) - evals[0]) * (matrix3_get(m, b, b) - evals[0]) -
            matrix3_get(m, b, a) * matrix3_get(m, a, b);
        Da = -matrix3_get(m, a, fix_column) * (matrix3_get(m, b, b) - evals[0]) +
            matrix3_get(m, b, fix_column) * matrix3_get(m, a, b);
        Db = (matrix3_get(m, a, a) - evals[0]) * -matrix3_get(m, b, fix_column) +
            matrix3_get(m, b, a) * matrix3_get(m, a, fix_column);
    }

    // Write the eigenvalue to the array
    switch(fix_column) {
    case 0:
        evecs[0] = vector3_new(1, Da/D, Db/D);
        break;
    case 1:
        evecs[0] = vector3_new(Da/D, 1, Db/D);
        break;
    case 2:
        evecs[0] = vector3_new(Da/D, Db/D, 1);
        break;
    }

    // Establish evec2 to be orthogonal.
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
        }
        D = (matrix3_get(m, a, a) - evals[1]) * (matrix3_get(m, b, b) - evals[1]) -
            matrix3_get(m, b, a) * matrix3_get(m, a, b);
        Da = -matrix3_get(m, a, fix_column) * (matrix3_get(m, b, b) - evals[1]) +
            matrix3_get(m, b, fix_column) * matrix3_get(m, a, b);
    }
    // Fix evec2 to be perpendicular
    Vector3 dot_part;
    switch(fix_column) {
    case 0:
        evecs[1] = vector3_new(1, Da/D, 0);
        dot_part = vector3_mul(vector3_z(),
            vector3_dot(evecs[1], evecs[0]) / evecs[0].z);
        break;
    case 1:
        evecs[1] = vector3_new(Da/D, 1, 0);
        dot_part = vector3_mul(vector3_z(),
            vector3_dot(evecs[1], evecs[0]) / evecs[0].z);
        break;
    case 2:
        evecs[1] = vector3_new(Da/D, 0, 1);
        dot_part = vector3_mul(vector3_y(),
            vector3_dot(evecs[1], evecs[0]) / evecs[0].y);
        break;
    }
    vector3_place_sub(&evecs[1], dot_part);

    // Generate the third eigenvector
    evecs[2] = vector3_cross(evecs[0], evecs[1]);

    // Normalize for cleanliness
    vector3_place_div(&evecs[0], vector3_mag(evecs[0]));
    vector3_place_div(&evecs[1], vector3_mag(evecs[1]));
    vector3_place_div(&evecs[2], vector3_mag(evecs[2]));
    return evecs;
}
Matrix3 matrix3_symmetric_invert_mat(const Matrix3 m) {
    double* evals = matrix3_get_evals(m);
    Vector3* evecs = matrix3_get_symmetric_evecs(m, evals);
    Matrix3 orthog = {
        evecs[0].x, evecs[1].x, evecs[2].x,
        evecs[0].y, evecs[1].y, evecs[2].y,
        evecs[0].z, evecs[1].z, evecs[2].z,
    };
    Matrix3 diag = matrix3_diag(1.0/evals[0], 1.0/evals[1], 1.0/evals[2]);
    free(evals);
    free(evecs);
    return matrix3_mat_mul(orthog,
        matrix3_mat_mul(diag, matrix3_transpose(orthog)));
}
Matrix3 matrix3_symmetric_invert(const double* evals,  const Vector3* evecs) {
    Matrix3 orthog = {
        evecs[0].x, evecs[1].x, evecs[2].x,
        evecs[0].y, evecs[1].y, evecs[2].y,
        evecs[0].z, evecs[1].z, evecs[2].z,
    };
    Matrix3 diag = matrix3_diag(1.0/evals[0], 1.0/evals[1], 1.0/evals[2]);
    return matrix3_mat_mul(orthog,
        matrix3_mat_mul(diag, matrix3_transpose(orthog)));
}
Matrix3 matrix3_symmetric_reconstruct(const double* evals,  const Vector3* evecs) {
    Matrix3 orthog = {
        evecs[0].x, evecs[1].x, evecs[2].x,
        evecs[0].y, evecs[1].y, evecs[2].y,
        evecs[0].z, evecs[1].z, evecs[2].z,
    };
    Matrix3 diag = matrix3_diag(evals[0], evals[1], evals[2]);
    return matrix3_mat_mul(orthog,
        matrix3_mat_mul(diag, matrix3_transpose(orthog)));
}


Quaternion quaternion_new(double r, double i, double j, double k) {
    Quaternion q;
    q.r = r;
    q.i = i;
    q.j = j;
    q.k = k;
    q.invalid_mat = 1;
    return q;
}
Quaternion quaternion_zero() {
    Quaternion q;
    q.r = 0;
    q.i = 0;
    q.j = 0;
    q.k = 0;
    q.invalid_mat = 1;
    return q;
}
Quaternion quaternion_identity() {
    Quaternion q;
    q.r = 1;
    q.i = 0;
    q.j = 0;
    q.k = 0;
    q.invalid_mat = 1;
    return q;
}

Quaternion quaternion_mul(const Quaternion q, double d) {
    return quaternion_new(q.r * d, q.i * d, q.j * d, q.k * d);
}
void quaternion_place_mul(Quaternion* q, double d) {
    q->r *= d;
    q->i *= d;
    q->j *= d;
    q->k *= d;
}
Quaternion quaternion_div(const Quaternion q, double d) {
    return quaternion_new(q.r / d, q.i / d, q.j / d, q.k / d);
}
void quaternion_place_div(Quaternion* q, double d) {
    q->r /= d;
    q->i /= d;
    q->j /= d;
    q->k /= d;
}
Quaternion quaternion_add(const Quaternion q1, const Quaternion q2) {
    return quaternion_new(q1.r + q2.r, q1.i + q2.i, q1.j + q2.j, q1.k + q2.k);
}
void quaternion_place_add(Quaternion* q1, const Quaternion q2) {
    q1->r += q2.r;
    q1->i += q2.i;
    q1->j += q2.j;
    q1->k += q2.r;
}
Quaternion quaternion_sub(const Quaternion q1, const Quaternion q2) {
    return quaternion_new(q1.r - q2.r, q1.i - q2.i, q1.j - q2.j, q1.k - q2.k);
}
void quaternion_place_sub(Quaternion* q1, const Quaternion q2) {
    q1->r -= q2.r;
    q1->i -= q2.i;
    q1->j -= q2.j;
    q1->k -= q2.r;
}

Quaternion quaternion_mul_quat(const Quaternion q1, const Quaternion q2) {
    Vector3 v1 = quaternion_vec(q1);
    Vector3 v2 = quaternion_vec(q2);
    double d = q1.r * q2.r - vector3_dot(v1, v2);
    Vector3 p1 = vector3_mul(v2, q1.r);
    Vector3 p2 = vector3_mul(v1, q2.r);
    Vector3 v = vector3_add(vector3_add(vector3_cross(v1, v2), p2), p1);
    vector3_place_add(&v, p2);
    vector3_place_add(&v, p1);
    return quaternion_new(d, v.x, v.y, v.z);
}

Quaternion quaternion_inverse(const Quaternion q) {
    Quaternion conj = quaternion_conjugate(q);
    return quaternion_div(conj, quaternion_mag2(q));
}
Quaternion quaternion_conjugate(const Quaternion q) {
    return quaternion_new(q.r, -q.i, -q.j, -q.k);
}
void quaternion_calc_matrix(Quaternion* q) {
    if (q->invalid_mat) {
        double s = 1/quaternion_mag2(*q);
        q->mat = (Matrix3) {
            1 - 2*s*(q->j*q->j + q->k*q->k), 2*s*(q->i*q->j - q->k*q->r), 2*s*(q->i*q->k + q->j*q->r),
            2*s*(q->i*q->j + q->k*q->r), 1 - 2*s*(q->i*q->i + q->k*q->k), 2*s*(q->j*q->k - q->i*q->r),
            2*s*(q->i*q->k - q->j*q->r), 2*s*(q->j*q->k + q->i*q->r), 1 - 2*s*(q->i*q->i + q->j*q->j)
        };
        q->invalid_mat = 0;
    }
}
Vector3 quaternion_vec(const Quaternion q) {
    return vector3_new(q.i, q.j, q.k);
}
double quaternion_mag(const Quaternion q) {
    return sqrt(q.r * q.r + q.i * q.i + q.j * q.j + q.k * q.k);
}
double quaternion_mag2(const Quaternion q) {
    return q.r * q.r + q.i * q.i + q.j * q.j + q.k * q.k;
}
