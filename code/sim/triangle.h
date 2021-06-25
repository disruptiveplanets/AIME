#pragma once

#include "algebra.h"
#include "linkedlist.h"
#include "sphericalharmonics.h"

typedef struct theta_phi {
    double theta;
    double phi;
} ThetaPhi;

typedef struct chunk {
    double density;
    double ab;// alpha_above
    ThetaPhi ul, ur, ll, lr;
    double be;// alpha_below
} Chunk;

typedef struct triangle {
    Vector3 v;
    Vector3 l1;
    Vector3 l2;
    Vector3 norm;
    Vector3 premul;
    const Chunk* parent;
} Triangle;

// ThetaPhi functions
ThetaPhi theta_phi_from_cart(Vector3 c);
Vector3 theta_phi_to_vector(ThetaPhi* tp, double r);

// Chunk functions
Chunk* chunk_new(double alpha_above, double alpha_below,
    int n_here, int face, int a, int b);
void chunk_shape(Chunk* c, double density_, uint L,
    const Cons* clms, Cons** tris);


// Triangle functions
Triangle* triangle_new(const Chunk* parent, Vector3 v1, Vector3 v2, Vector3 v3);
void triangle_recenter(Triangle* t, const Vector3 l);
void triangle_rotate(Triangle* t, const Matrix3* m);

Vector3 triangle_get_lever_arm(const Triangle* t);
double triangle_get_mass(const Triangle* t);
Vector3 triangle_get_torque(const Triangle* t);
double triangle_get_Isame(const Triangle* t, int a);
double triangle_get_Idiff(const Triangle* t, int a, int b);
Vector3* triangle_get_corners(const Triangle* t);
int triangle_is_edge(const Triangle* t);
