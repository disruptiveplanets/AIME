#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793438462643383
#define EPSILON 1e-10

typedef struct complex {
    double r;
    double i;
} Complex;

Complex complex_sqrt(double);
Complex complex_pow(Complex, double);

typedef struct vector3 {
    double x;
    double y;
    double z;
    double mag;
    int invalid_mag;
} Vector3;

typedef struct matrix3 {
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;
    double h;
    double i;
} Matrix3;

typedef struct quaternion {
    double r;
    double i;
    double j;
    double k;
    Matrix3 mat;
    int invalid_mat;
} Quaternion;

// Vector functions

Vector3 vector3_new(double, double, double);
void vector3_init(Vector3*, double, double, double);
Vector3 vector3_zero();
Vector3 vector3_x();
Vector3 vector3_y();
Vector3 vector3_z();
Vector3 vector3_negative(Vector3);

Vector3 vector3_mul(const Vector3, double);
void vector3_place_mul(Vector3*, double);
Vector3 vector3_div(const Vector3, double);
void vector3_place_div(Vector3*, double);
Vector3 vector3_add(const Vector3, const Vector3);
void vector3_place_add(Vector3*, const Vector3);
Vector3 vector3_sub(const Vector3, const Vector3);
void vector3_place_sub(Vector3*, const Vector3);

double vector3_dot(const Vector3, const Vector3);
Vector3 vector3_cross(const Vector3, const Vector3);
Vector3 vector3_mat_mul(const Matrix3, const Vector3);

double vector3_mag(Vector3);
double vector3_mag2(const Vector3);


// Matrix functions
void matrix3_init(Matrix3*, double, double, double,
    double, double, double,
    double, double, double);
Matrix3 matrix3_zero();
Matrix3 matrix3_rotation_x(double);
Matrix3 matrix3_rotation_y(double);
Matrix3 matrix3_rotation_z(double);
Matrix3 matrix3_identity();
Matrix3 matrix3_diag(double, double, double);

Matrix3 matrix3_mul(const Matrix3, double);
void matrix3_place_mul(Matrix3*, double);
Matrix3 matrix3_div(const Matrix3, double);
void matrix3_place_div(Matrix3*, double);
Matrix3 matrix3_add(const Matrix3, const Matrix3);
void matrix3_place_add(Matrix3*, const Matrix3);
Matrix3 matrix3_sub(const Matrix3, const Matrix3);
void matrix3_place_sub(Matrix3*, const Matrix3);
Matrix3 matrix3_transpose(const Matrix3);
void matrix3_place_transpose(Matrix3*);

Matrix3 matrix3_mat_mul(const Matrix3, const Matrix3);

double matrix3_get(const Matrix3, int a, int b);
double matrix3_det(const Matrix3);
double matrix3_trace(const Matrix3);

double* matrix3_get_evals(const Matrix3);
Vector3* matrix3_get_symmetric_evecs(const Matrix3, const double* evals);
Matrix3 matrix3_symmetric_invert_mat(const Matrix3);
Matrix3 matrix3_symmetric_invert(const double* evals,  const Vector3* evecs);
Matrix3 matrix3_symmetric_reconstruct(const double* evals,  const Vector3* evecs);


// Quaternion functions

Quaternion quaternion_new(double, double, double, double);
Quaternion quaternion_zero();
Quaternion quaternion_identity();

Quaternion quaternion_mul(const Quaternion, double);
void quaternion_place_mul(Quaternion*, double);
Quaternion quaternion_div(const Quaternion, double);
void quaternion_place_div(Quaternion*, double);
Quaternion quaternion_add(const Quaternion, const Quaternion);
void quaternion_place_add(Quaternion*, const Quaternion);
Quaternion quaternion_sub(const Quaternion, const Quaternion);
void quaternion_place_sub(Quaternion*, const Quaternion);

Quaternion quaternion_mul_quat(const Quaternion, const Quaternion);

Quaternion quaternion_inverse(const Quaternion);
Quaternion quaternion_conjugate(const Quaternion);
void quaternion_calc_matrix(Quaternion*);
Vector3 quaternion_vec(const Quaternion);
double quaternion_mag(const Quaternion);
double quaternion_mag2(const Quaternion);
