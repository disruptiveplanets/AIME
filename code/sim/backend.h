#pragma once
#include "algebra.h"
#include "triangle.h"

#define DT_SLOPE 5e-6// 5e-6
#define ONE_SECOND_TORQUE 1e-8// 1e-8
#define INTEGRAL_LIMIT_FRAC 1.0e-5// 5e-6
#define NUM_THREADS 4
    // Torque at closest approach divided by torque at start of sim.
#define max(a, b) ((a) > (b) ? (a) : (b))
#define G 6.67408e-11
//#define _DEBUG

typedef unsigned int uint;

typedef struct asteroid {
    uint L; // Max degree of harmonics
    uint n; // Number of cube subdivisions on outermost shell
    uint m; // Number of shells
    double mass;
    Matrix3 moi;
    Matrix3 moiInverse;
    Cons* chunks_head;
    Cons* triangles_head;
    double edge_dist; // Limit of the integration region
    Vector3 position;
    Vector3 velocity;
    Vector3 spin;
    double mu;
    Quaternion orientation;
    double time;

    // Shape features
    double mean_density;// Used to generate rough parameters of the asteroid
    double radius;// Used to generate rough parameters of the asteroid

    // Orbital factors
    double energy;
    Vector3 ang_mom;
    double closest_approach;
    double excess_vel;
    double impact_parameter;

    double max_quat_mag;
} Asteroid;

// You must define these yourself
void asteroid_record_resolved(Vector3 spin);
void asteroid_record_unresolved();


// Asteroid functions
Asteroid* asteroid_new(int L, int n, int m, Cons* clms, Cons* densities,
    Vector3 spin, double impact_parameter, double velocity,
    double central_mass);
int asteroid_simulate(Asteroid* a, double cadence);
void asteroid_free(Asteroid* a);
