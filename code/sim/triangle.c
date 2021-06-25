#include "triangle.h"

// Triangle auxilarry functions
double get_Isame_component(double l1b, double l2b, double vb) {
    return l1b * l1b * l1b + l1b * l1b * l2b + l1b * l2b * l2b + l2b * l2b * l2b
        +  5 * vb * (l1b * l1b + l1b * l2b + l2b * l2b)
        + 10 * vb * vb * (l1b + l2b) + 10 * vb * vb * vb;
}
double get_Idiff_component(double l1a, double l1b,
    double l2a, double l2b, double va, double vb) {

    return 6 * l1a * l1b + 2 * l2a * l2b + 2 * l1b * (l2a + 5 * va) +
        2 * l1a * (l2b + 5 * vb) + 5 * l2a * vb + 5 * l2b * va + 20 * va * vb;
}
double get_torque_component(const Vector3* l1, const Vector3* l2,
    const Vector3* v) {

    return l1->z * (
        (l2->x + 5 * v->x) * (l2->x + 2 * l1->x) +
        (l2->y + 5 * v->y) * (l2->y + 2 * l1->y) +
        3 * (l1->x * l1->x + l1->y * l1->y) +
        10 * (v->x * v->x + v->y * v->y));
}


// ThetaPhi functions
ThetaPhi theta_phi_from_cart(Vector3 c){
    return (ThetaPhi) { acos(c.z / vector3_mag(c)), atan2(c.y, c.x) };
}
Vector3 theta_phi_to_vector(ThetaPhi* tp, double r) {
    return vector3_new(r * sin(tp->theta) * cos(tp->phi),
                       r * sin(tp->theta) * sin(tp->phi),
                       r * cos(tp->theta));
}


// Chunk functions
Chunk* chunk_new(double alpha_above, double alpha_below,
    int n_here, int face, int a, int b) {

    Chunk* c = (Chunk*) malloc(sizeof(Chunk));

    Vector3 v, avec, bvec;
        // Vectors indicating the directions of a and b in cart. coords
        // avec cross bvec points out
    switch (face) {
    case 0:
        v = vector3_new(1, 1, 1);
        avec = vector3_new(0, 0, -2.0 / n_here);
        bvec = vector3_new(-2.0 / n_here, 0, 0);
        break;
    case 1:
        v = vector3_new(1, -1, 1);
        avec = vector3_new(0, 0, -2.0 / n_here);
        bvec = vector3_new(0, 2.0 / n_here, 0);
        break;
    case 2:
        v = vector3_new(-1, -1, 1);
        avec = vector3_new(0, 0, -2.0 / n_here);
        bvec = vector3_new(2.0 / n_here, 0, 0);
        break;
    case 3:
        v = vector3_new(-1, 1, 1);
        avec = vector3_new(0, 0, -2.0 / n_here);
        bvec = vector3_new(0, -2.0 / n_here, 0);
        break;
    case 4:
        v = vector3_new(-1, -1, 1);
        avec = vector3_new(2.0 / n_here, 0, 0);
        bvec = vector3_new(0, 2.0 / n_here, 0);
        break;
    case 5:
        v = vector3_new(1, 1, -1);
        avec = vector3_new(0, -2.0 / n_here, 0);
        bvec = vector3_new(-2.0 / n_here, 0, 0);
        break;
    }
    Vector3 amul = vector3_mul(avec, a);
    Vector3 aplusmul = vector3_mul(avec, a + 1);
    Vector3 bmul = vector3_mul(bvec, b);
    Vector3 bplusmul = vector3_mul(bvec, b + 1);

    c->density = 0;
    c->ab = alpha_above;
    c->be = alpha_below;
    c->ul = theta_phi_from_cart(vector3_add(vector3_add(v, amul), bmul));
    c->ur = theta_phi_from_cart(vector3_add(vector3_add(v, amul), bplusmul));
    c->ll = theta_phi_from_cart(vector3_add(vector3_add(v, aplusmul), bmul));
    c->lr = theta_phi_from_cart(vector3_add(vector3_add(v, aplusmul), bplusmul));
    return c;
}
void chunk_shape(Chunk* c, double density, uint L,
    const Cons* clms, Cons** tris) {

    c->density = density;
    double rul = 0, rur = 0, rll = 0, rlr = 0;

    const Cons* clm = clms;
    for (int l = 0; l <= L; l++) {
        for(int m = -l; m<= l; m++) {
            if (clm == NULL) {
                printf("ERROR: Ran out of clms to unpack!");
            }
            rul += *(double*)clm->value * real_spherical_harmonic(l, m,
                c->ul.theta, c->ul.phi);
            rur += *(double*)clm->value * real_spherical_harmonic(l, m,
                c->ur.theta, c->ur.phi);
            rll += *(double*)clm->value * real_spherical_harmonic(l, m,
                c->ll.theta, c->ll.phi);
            rlr += *(double*)clm->value * real_spherical_harmonic(l, m,
                c->lr.theta, c->lr.phi);
            clm = clm->next;
        }
    }

    Vector3 vul = theta_phi_to_vector(&c->ul, rul);
    Vector3 vur = theta_phi_to_vector(&c->ur, rur);
    Vector3 vll = theta_phi_to_vector(&c->ll, rll);
    Vector3 vlr = theta_phi_to_vector(&c->lr, rlr);

    if (c->be <= 0) {
        Vector3 o = vector3_zero(); // origin

        cons_push_inc(tris, triangle_new(c, vector3_mul(vur, c->ab),
            vector3_mul(vul, c->ab), vector3_mul(vlr, c->ab)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vll, c->ab),
            vector3_mul(vlr, c->ab), vector3_mul(vul, c->ab)));

        cons_push_inc(tris, triangle_new(c, vector3_mul(vur, c->ab),
            vector3_mul(vlr, c->ab), o));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vlr, c->ab),
            vector3_mul(vll, c->ab), o));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vll, c->ab),
            vector3_mul(vul, c->ab), o));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vul, c->ab),
            vector3_mul(vur, c->ab), o));
    }
    else {
        cons_push_inc(tris, triangle_new(c, vector3_mul(vur, c->ab),
            vector3_mul(vul, c->ab), vector3_mul(vlr, c->ab)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vll, c->ab),
            vector3_mul(vlr, c->ab), vector3_mul(vul, c->ab)));

        cons_push_inc(tris, triangle_new(c, vector3_mul(vur, c->ab),
            vector3_mul(vlr, c->ab), vector3_mul(vur, c->be)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vlr, c->ab),
            vector3_mul(vll, c->ab), vector3_mul(vlr, c->be)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vll, c->ab),
            vector3_mul(vul, c->ab), vector3_mul(vll, c->be)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vul, c->ab),
            vector3_mul(vur, c->ab), vector3_mul(vul, c->be)));

        cons_push_inc(tris, triangle_new(c, vector3_mul(vlr, c->be),
            vector3_mul(vur, c->be), vector3_mul(vlr, c->ab)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vll, c->be),
            vector3_mul(vlr, c->be), vector3_mul(vll, c->ab)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vul, c->be),
            vector3_mul(vll, c->be), vector3_mul(vul, c->ab)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vur, c->be),
            vector3_mul(vul, c->be), vector3_mul(vur, c->ab)));

        cons_push_inc(tris, triangle_new(c, vector3_mul(vul, c->be),
            vector3_mul(vur, c->be), vector3_mul(vll, c->be)));
        cons_push_inc(tris, triangle_new(c, vector3_mul(vlr, c->be),
            vector3_mul(vll, c->be), vector3_mul(vur, c->be)));
    }
}


// Triangle functions
Triangle* triangle_new(const Chunk* parent, Vector3 v1,
    Vector3 v2, Vector3 v3) {

    Triangle* t = (Triangle*) malloc(sizeof(Triangle));
    t->v = v1;
    t->l1 = vector3_sub(v2, v1);
    t->l2 = vector3_sub(v3, v1);
    t->norm = vector3_cross(t->l1, t->l2);
    t->parent = parent;
    return t;
}
void triangle_recenter(Triangle* t, const Vector3 l) {
    vector3_place_add(&t->v, l);
}
void triangle_rotate(Triangle* t, const Matrix3* m) {
    t->v = vector3_mat_mul(*m, t->v);
    t->l1 = vector3_mat_mul(*m, t->l1);
    t->l2 = vector3_mat_mul(*m, t->l2);
    t->norm = vector3_mat_mul(*m, t->norm);
}

Vector3 triangle_get_lever_arm(const Triangle* t) {
    return vector3_mul(t->norm, 1/24.0 * t->parent->density * (
        vector3_mag2(t->l1) + vector3_mag2(t->l2) + 6 * vector3_mag2(t->v)
        + vector3_dot(t->l1, t->l2)
        + 4 * vector3_dot(vector3_add(t->l1, t->l2), t->v)
    ));
}
double triangle_get_mass(const Triangle* t) {
    return 1/18.0 * t->parent->density * vector3_dot(t->norm,
        vector3_add(vector3_add(vector3_mul(t->v, 3), t->l1), t->l2));
}
Vector3 triangle_get_torque(const Triangle* t) {
    Vector3 cross = vector3_cross(t->norm, vector3_z());
    return vector3_mul(cross, 1/40.0 * t->parent->density * (
        get_torque_component(&t->l1, &t->l2, &t->v) +
        get_torque_component(&t->l2, &t->l1, &t->v) +
        5 * t->v.z * (
            t->l1.x * t->l1.x + t->l1.y * t->l1.y + t->l2.x * t->l2.x + t->l2.y * t->l2.y +
            t->l1.x * t->l2.x + t->l1.y * t->l2.y +
            2 * t->v.x * (2 * t->l1.x + 2 * t->l2.x + 3 * t->v.x) +
            2 * t->v.y * (2 * t->l1.y + 2 * t->l2.y + 3 * t->v.y)
        )));
}
double triangle_get_Isame(const Triangle* t, int a) {
    Vector3 b;
    Vector3 c;
    switch(a) {
    case 0:
        b = vector3_mul(vector3_y(), get_Isame_component(t->l1.y, t->l2.y, t->v.y));
        c = vector3_mul(vector3_z(), get_Isame_component(t->l1.z, t->l2.z, t->v.z));
        break;
    case 1:
        b = vector3_mul(vector3_x(), get_Isame_component(t->l1.x, t->l2.x, t->v.x));
        c = vector3_mul(vector3_z(), get_Isame_component(t->l1.z, t->l2.z, t->v.z));
        break;
    case 2:
        b = vector3_mul(vector3_y(), get_Isame_component(t->l1.y, t->l2.y, t->v.y));
        c = vector3_mul(vector3_x(), get_Isame_component(t->l1.x, t->l2.x, t->v.x));
        break;
    }
    return 1/60.0 * t->parent->density * vector3_dot(t->norm, vector3_add(b, c));
}
double triangle_get_Idiff(const Triangle* t, int a, int b) {
    Vector3 vc;
    int c;
    if ((a == 0 && b == 1) || (a == 1 && b == 0)) {
        vc = vector3_z();
        c = 2;
    }
    if ((a == 0 && b == 2) || (a == 2 && b == 0)) {
        vc = vector3_y();
        c = 1;
    }
    if ((a == 2 && b == 1) || (a == 1 && b == 2)) {
        vc = vector3_x();
        c = 0;
    }
    double l1a, l1b, l1c, l2a, l2b, l2c, va, vb, v_c;
    switch(a) {
    case 0:
        l1a = t->l1.x;
        l2a = t->l2.x;
        va = t->v.x;
        break;
    case 1:
        l1a = t->l1.y;
        l2a = t->l2.y;
        va = t->v.y;
        break;
    case 2:
        l1a = t->l1.z;
        l2a = t->l2.z;
        va = t->v.z;
        break;
    }
    switch(b) {
    case 0:
        l1b = t->l1.x;
        l2b = t->l2.x;
        vb = t->v.x;
        break;
    case 1:
        l1b = t->l1.y;
        l2b = t->l2.y;
        vb = t->v.y;
        break;
    case 2:
        l1b = t->l1.z;
        l2b = t->l2.z;
        vb = t->v.z;
        break;
    }
    switch(c) {
    case 0:
        l1c = t->l1.x;
        l2c = t->l2.x;
        v_c = t->v.x;
        break;
    case 1:
        l1c = t->l1.y;
        l2c = t->l2.y;
        v_c = t->v.y;
        break;
    case 2:
        l1c = t->l1.z;
        l2c = t->l2.z;
        v_c = t->v.z;
        break;
    }
    return -1/120.0 * t->parent->density * vector3_dot(vc, t->norm) * (
            l1c * get_Idiff_component(l1a, l1b, l2a, l2b, va, vb) +
            l2c * get_Idiff_component(l2a, l2b, l1a, l1b, va, vb) +
            5 * v_c * (l1a * (2 * l1b + l2b) + l2a * (2 * l2b + l1b)+
            4 * va * (l2b + l1b) + 4 * vb * (l2a + l1a) +
            12 * va * vb)
        );
}
Vector3* triangle_get_corners(const Triangle* t) {
    Vector3* corners = (Vector3*) malloc(sizeof(Vector3) * 3);
    corners[0] = t->v;
    corners[1] = vector3_add(t->v, t->l1);
    corners[2] = vector3_add(t->v, t->l2);
    return corners;
}
int triangle_is_edge(const Triangle* t)  {
    return abs(t->parent->ab - 1.0) < EPSILON;
}
