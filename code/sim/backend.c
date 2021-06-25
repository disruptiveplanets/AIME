#include "backend.h"

void make_chunks(Asteroid* a) {
    Cons* chunks = NULL;
    a->chunks_head = NULL;
    for (int f = 0; f < 6; f++) {
        for (int s = 0; s < a->m; s++){
            const int length = a->n * (s + 1) / a->m;
            const double alpha_above = (double)(s + 1) / a->m;
            const double alpha_below = (double)(s) / a->m;
            for (int ai = 0; ai < length; ai++){
                for (int bi = 0; bi < length; bi++) {
                    cons_push_inc(&chunks, chunk_new(alpha_above, alpha_below,
                        length, f, ai, bi));
                    if (a->chunks_head == NULL) {
                        a->chunks_head = chunks;
                    }
                }
            }
        }
    }
}

void set_pos(Asteroid* a, double b) {
    // Asteroid enters at +x and is traveling towards -x, with offset in +y
    // direction.
    a->edge_dist = b * pow(INTEGRAL_LIMIT_FRAC, -1/3.0);
    a->position = vector3_mul(vector3_new(0, b,
        -sqrt(a->edge_dist * a->edge_dist - b * b)), 1 - EPSILON);
        // So that it reads as just inside the allowed region

    a->energy = 0.5 * vector3_mag2(a->velocity) -
        a->mu / vector3_mag(a->position);
    a->ang_mom = vector3_cross(a->position, a->velocity);
    double velocity_periapsis = a->mu / vector3_mag(a->ang_mom) +
        sqrt(a->mu * a->mu / vector3_mag2(a->ang_mom) + 2 * a->energy);
    a->closest_approach = vector3_mag(a->ang_mom) / velocity_periapsis;
    a->excess_vel = sqrt(2 * a->energy);
    a->impact_parameter = vector3_mag(a->ang_mom) / a->excess_vel;

    Vector3 Omega = vector3_div(a->ang_mom, vector3_mag2(a->position));
    vector3_place_sub(&a->spin, Omega);
}

void calculate_moi(Asteroid* a) {
    double Ixx = 0, Iyy = 0, Izz = 0, Ixz = 0, Ixy = 0, Iyz = 0;
    for(Cons* t = a->triangles_head; t != NULL; t = t->next) {
        Ixx += triangle_get_Isame((Triangle*)t->value, 0);
        Iyy += triangle_get_Isame((Triangle*)t->value, 1);
        Izz += triangle_get_Isame((Triangle*)t->value, 2);
        Ixy += triangle_get_Idiff((Triangle*)t->value, 0, 1);
        Ixz += triangle_get_Idiff((Triangle*)t->value, 0, 2);
        Iyz += triangle_get_Idiff((Triangle*)t->value, 1, 2);
    }
    a->moi = (Matrix3) {Ixx, Ixy, Ixz,
                        Ixy, Iyy, Iyz,
                        Ixz, Iyz, Izz,};
    double* evals = matrix3_get_evals(a->moi);
    Vector3* evecs = matrix3_get_symmetric_evecs(a->moi, evals);
    a->moiInverse = matrix3_symmetric_invert(evals, evecs);
    a->moi = matrix3_symmetric_reconstruct(evals, evecs);
    a->orientation = quaternion_identity();

    free(evals);
    free(evecs);
}

void calculate_mass(Asteroid* a) {
    a->mass = 0;
    for (Cons* t = a->triangles_head; t != NULL; t = t->next) {
        a->mass += triangle_get_mass((Triangle*)t->value);
    }
}

Vector3 get_com(const Asteroid* a) {
    Vector3 total_arm = vector3_zero();
    for (Cons* t = a->triangles_head; t != NULL; t = t->next) {
        vector3_place_add(&total_arm,
            triangle_get_lever_arm((Triangle*)t->value));
    }
    return vector3_div(total_arm, 1 / a->mass);
}

void recenter(Asteroid* a) {
    Vector3 com = get_com(a);
    for (Cons* t = a->triangles_head; t != NULL; t = t->next) {
        triangle_recenter((Triangle*)t->value, vector3_negative(com));
    }

    #ifdef _DEBUG
    std::cout << "New COM: " << get_com() << " should be zero." << std::endl;
    #endif
}

Matrix3 global_to_inertial(const Asteroid* a) {
    double theta = atan2(a->position.y, -a->position.z);
    return matrix3_rotation_x(theta);
}

Matrix3 inertial_to_global(const Asteroid* a) {
    double theta = atan2(a->position.y, -a->position.z);
    return matrix3_rotation_x(-theta);
}


Vector3 get_rot_ang_mom(Asteroid* a) {
    Vector3 Omega = vector3_div(a->ang_mom, vector3_mag2(a->position));
    Quaternion inverse_quat = quaternion_inverse(a->orientation);
    quaternion_calc_matrix(&a->orientation);
    quaternion_calc_matrix(&inverse_quat);
    Matrix3 moiGlobal = matrix3_mat_mul(a->orientation.mat,
        matrix3_mat_mul(a->moi, inverse_quat.mat));
    return  vector3_mat_mul(global_to_inertial(a),
        vector3_mat_mul(moiGlobal, vector3_add(a->spin, Omega)));
}

Vector3 get_torque(Asteroid* a) {
    Vector3 torque = vector3_zero();
    quaternion_calc_matrix(&a->orientation);
    for (Cons* t = a->triangles_head; t != NULL; t = t->next) {
        triangle_rotate((Triangle*)t->value, &a->orientation.mat);
        vector3_place_add(&torque,  triangle_get_torque((Triangle*)t->value));
    }
    return vector3_mul(torque, (a->mu / pow(vector3_mag(a->position), 3)));
}

void update_position(Asteroid* a, double dt) {
    Vector3 accel = vector3_mul(a->position,
        -a->mu / pow(vector3_mag(a->position), 3));
    vector3_place_add(&a->velocity, vector3_mul(accel, dt));
    vector3_place_add(&a->position, vector3_mul(a->velocity, dt));
}

void update_orientation(Asteroid* a, double dt) {
    Vector3 Omega = vector3_div(a->ang_mom, vector3_mag2(a->position));
    Vector3 torque = get_torque(a);

    Quaternion inv_quat = quaternion_inverse(a->orientation);
    quaternion_calc_matrix(&a->orientation);
    quaternion_calc_matrix(&inv_quat);

    Matrix3 moiGlobal = matrix3_mat_mul(a->orientation.mat,
        matrix3_mat_mul(a->moi, inv_quat.mat));
    Matrix3 moiGlobalInverse = matrix3_mat_mul(a->orientation.mat,
        matrix3_mat_mul(a->moiInverse, inv_quat.mat));

    // Calculate change in omega
    Vector3 omega_dot = vector3_mat_mul(moiGlobalInverse, vector3_sub(torque,
        vector3_cross(vector3_add(Omega, a->spin),
            vector3_mat_mul(moiGlobal, vector3_add(Omega, a->spin)))));
    vector3_place_add(&omega_dot, vector3_mul(Omega, 2 * vector3_dot(a->position,
        a->velocity) / vector3_mag2(a->position)));
    vector3_place_sub(&omega_dot, vector3_cross(Omega, a->spin));

    // Change omega
    vector3_place_add(&a->spin, vector3_mul(omega_dot, dt));

    Quaternion d_quat = quaternion_new(
        0, 0.5 * a->spin.x, 0.5 * a->spin.y, 0.5 * a->spin.z);
    quaternion_place_add(&a->orientation,
        quaternion_mul(quaternion_mul_quat(d_quat, a->orientation), dt));

    #ifdef _DEBUG
    max_quat_mag = max(max_quat_mag, vector3_mag(d_quat));
    #endif

    quaternion_place_div(&a->orientation, quaternion_mag(a->orientation));
}

Asteroid* asteroid_new(int L, int n, int m, Cons* clms_head,
    Cons* densities_head, Vector3 spin, double impact_parameter,
    double velocity, double central_mass) {

    Asteroid* a = (Asteroid*) malloc(sizeof(Asteroid));
    a->L = L;
    a->n = n;
    a->m = m;
    a->radius = *(double*)clms_head->value;

    make_chunks(a);

    a->triangles_head = cons_new(NULL);// Start with a null value.
    Cons* triangles = a->triangles_head;
    Cons* density = densities_head;
    Cons* chunk = a->chunks_head;
    a->mean_density = 0;

    // Initialize all chunks and triangles
    int chunk_size = 0;
    while (chunk != NULL) {
        a->mean_density += *(double*)density->value;
        chunk_shape((Chunk*)chunk->value, *(double*)density->value, a->L,
            clms_head, &triangles);
        chunk = chunk->next;
        density = density->next;
        chunk_size++;
    }
    a->mean_density /= chunk_size;

    // Get rid of the head
    Cons* new_head = a->triangles_head->next;
    free(a->triangles_head);
    a->triangles_head = new_head;

    a->spin = spin;
    a->velocity = vector3_new(0, 0, velocity);
    a->mu = G * central_mass;

    calculate_mass(a);
    recenter(a);
    set_pos(a, impact_parameter);
    calculate_moi(a);

    return a;
}

int asteroid_simulate(Asteroid* a, double cadence) {
    // Motivation: Torque is proportional to 1/position^3. Angular acceleration
    // is now roughly constant every frame.
    const double scale_torque = a->mu / pow(a->closest_approach, 3)
        * a->mean_density * pow(a->radius, 5);
    const double min_delta_t = ONE_SECOND_TORQUE / scale_torque;
        // Toruqe at min delta t

    a->time = 0;
    int frames = 0;
    int cadence_index = -1;
    double dt;
    for (;vector3_mag(a->position) < a->edge_dist; frames++){
        dt = min_delta_t + DT_SLOPE * (pow(vector3_mag(a->position)
            / a->closest_approach, 3) - 1);
        update_orientation(a, dt);
        update_position(a, dt);
        a->time += dt;

        if ((int)(a->time / cadence) > cadence_index) {
            asteroid_record_resolved(a->spin);
            cadence_index = (int)(a->time / cadence);
        }
    }

    #ifdef _DEBUG
    printf("Simulation took %d seconds.", a->time);
    printf("Maximum dquaternion magnitude (want 0) %d", max_quat_mag);
    #endif

    return frames;
}

void asteroid_free(Asteroid* a) {
    cons_free_all(a->chunks_head);
    cons_free_all(a->triangles_head);
}
