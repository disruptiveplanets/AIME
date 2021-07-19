#include "triangle.hpp"

double real_spherical_harmonic(uint l, int m, double theta, double phi) {
    if (m < 0) {
        return sqrt(2) * pow(-1, m) *
        boost::math::spherical_harmonic(l, -m, theta, phi).imag();
    }
    if (m == 0) {
        return boost::math::spherical_harmonic(l, 0, theta, phi).real();
    }
    else {
        return sqrt(2) * pow(-1, m) *
        boost::math::spherical_harmonic(l, m, theta, phi).real();
    }
}

Triangle::Triangle(Chunk& parent, Vector3 v1, Vector3 v2, Vector3 v3) :
    v(v1), l1(v2-v1), l2(v3-v1), norm(Vector3::cross(l1, l2)), parent(parent) {}
std::array<Vector3, 3> Triangle::get_corners() const {
    return {v, v + l1, v + l2};
}
void Triangle::operator*=(Matrix3 const& m) {
    v = m * v;
    l1 = m * l1;
    l2 = m * l2;
    norm = m * norm;// Assume this is a rotation matrix
}
double Triangle::get_density() const{
    return parent.density;
}
bool Triangle::is_edge() const {
    return abs(parent.ab - 1.0) < EPSILON;
}

double Triangle::get_mass() const {
    return 1/18.0 * parent.density * Vector3::dot(norm, l1 + l2 + 3 * v);
}
Vector3 Triangle::get_lever_arm() const {
    return 1/24.0 * parent.density * norm * (
        l1.mag2() + l2.mag2() + 6 * v.mag2() + Vector3::dot(l1, l2) +
        4 * Vector3::dot(l1 + l2, v)
    );
}
double Triangle::get_Isame(int a) const {
    Vector3 b, c;
    switch(a){
    case 0:
        b = Vector3::y() * get_Isame_component(l1[1], l2[1], v[1]);
        c = Vector3::z() * get_Isame_component(l1[2], l2[2], v[2]);
        break;
    case 1:
        b = Vector3::x() * get_Isame_component(l1[0], l2[0], v[0]);
        c = Vector3::z() * get_Isame_component(l1[2], l2[2], v[2]);
        break;
    case 2:
        b = Vector3::y() * get_Isame_component(l1[1], l2[1], v[1]);
        c = Vector3::x() * get_Isame_component(l1[0], l2[0], v[0]);
        break;
    }
    return 1/60.0 * parent.density *
        Vector3::dot(norm, b + c);
}
bool Triangle::is_nan() const {
    //Vector3 norm;
    //Vector3 premul;
    //Chunk* parent;
    return isnan(v[0]) || isnan(l1[0]) || isnan(l2[0]);
}
double Triangle::get_Isame_component(double l1b, double l2b, double vb) const {
    return l1b * l1b * l1b + l1b * l1b * l2b + l1b * l2b * l2b + l2b * l2b * l2b
        +  5 * vb * (l1b * l1b + l1b * l2b + l2b * l2b)
        + 10 * vb * vb * (l1b + l2b) + 10 * vb * vb * vb;
}
double Triangle::get_Idiff(int a, int b) const {
    Vector3 vc;
    int c;
    if ((a == 0 && b == 1) || (a == 1 && b == 0)) {
        vc = Vector3::z();
        c = 2;
    }
    if ((a == 0 && b == 2) || (a == 2 && b == 0)) {
        vc = Vector3::y();
        c = 1;
    }
    if ((a == 2 && b == 1) || (a == 1 && b == 2)) {
        vc = Vector3::x();
        c = 0;
    }
    return -1/120.0 * parent.density * Vector3::dot(vc, norm) * (
            l1[c] * get_Idiff_component(l1[a], l1[b], l2[a], l2[b], v[a], v[b]) +
            l2[c] * get_Idiff_component(l2[a], l2[b], l1[a], l1[b], v[a], v[b]) +
            5 * v[c] * (l1[a] * (2 * l1[b] + l2[b]) + l2[a] * (2 * l2[b] + l1[b])+
            4 * v[a] * (l2[b] + l1[b]) + 4 * v[b] * (l2[a] + l1[a]) +
            12 * v[a] * v[b])
        );
}
double Triangle::get_Idiff_component(double l1a, double l1b,
    double l2a, double l2b, double va, double vb) const {
    return 6 * l1a * l1b + 2 * l2a * l2b + 2 * l1b * (l2a + 5 * va) +
        2 * l1a * (l2b + 5 * vb) + 5 * l2a * vb + 5 * l2b * va + 20 * va * vb;
}
Vector3 Triangle::get_torque() const {
    return 1/40.0 * parent.density * (
        get_torque_component(l1, l2) + get_torque_component(l2, l1) +
        5 * v[2] * (
            l1[0] * l1[0] + l1[1] * l1[1] + l2[0] * l2[0] + l2[1] * l2[1] +
            l1[0] * l2[0] + l1[1] * l2[1] +
            2 * v[0] * (2 * l1[0] + 2 * l2[0] + 3 * v[0]) +
            2 * v[1] * (2 * l1[1] + 2 * l2[1] + 3 * v[1])
        )) * Vector3::cross(norm, Vector3::z());
}
double Triangle::get_torque_component(Vector3 const& l1_, Vector3 const& l2_)
const {
    return l1_[2] * (
        (l2_[0] + 5 * v[0]) * (l2_[0] + 2 * l1_[0]) +
        (l2_[1] + 5 * v[1]) * (l2_[1] + 2 * l1_[1]) +
        3 * (l1_[0] * l1_[0] + l1_[1] * l1_[1]) +
        10 * (v[0] * v[0] + v[1] * v[1]));
}

Chunk::Chunk(double alpha_above, double alpha_below,
    int n_here, int face, int a, int b) : ab(alpha_above), be(alpha_below) {
    Vector3 v, avec, bvec;
        // Vectors indicating the directions of a and b in cart. coords
        // avec cross bvec points out
    switch (face) {
    case 0:
        v = Vector3({1, 1, 1});
        avec = Vector3({0, 0, -2.0 / n_here});
        bvec = Vector3({-2.0 / n_here, 0, 0});
        break;
    case 1:
        v = Vector3({1, -1, 1});
        avec = Vector3({0, 0, -2.0 / n_here});
        bvec = Vector3({0, 2.0 / n_here, 0});
        break;
    case 2:
        v = Vector3({-1, -1, 1});
        avec = Vector3({0, 0, -2.0 / n_here});
        bvec = Vector3({2.0 / n_here, 0, 0});
        break;
    case 3:
        v = Vector3({-1, 1, 1});
        avec = Vector3({0, 0, -2.0 / n_here});
        bvec = Vector3({0, -2.0 / n_here, 0});
        break;
    case 4:
        v = Vector3({-1, -1, 1});
        avec = Vector3({2.0 / n_here, 0, 0});
        bvec = Vector3({0, 2.0 / n_here, 0});
        break;
    case 5:
        v = Vector3({1, 1, -1});
        avec = Vector3({0, -2.0 / n_here, 0});
        bvec = Vector3({-2.0 / n_here, 0, 0});
        break;
    }
    ul = ThetaPhi(v + avec * a + bvec * b);
    ll = ThetaPhi(v + avec * (a + 1) + bvec * b);
    ur = ThetaPhi(v + avec * a + bvec * (b + 1));
    lr = ThetaPhi(v + avec * (a + 1) + bvec * (b + 1));
}

void Chunk::shape(double density_, uint L, std::vector<double> const& clms,
    std::vector<Triangle>& tris) {

    density = density_;
    double rul = 0, rur = 0, rll = 0, rlr = 0;

    auto clm = clms.begin();
    for (uint l = 0; l <= L; l++) {
        for(int m = -l; m <= (int)l; m++) {
            if (clm == clms.end()) {
                std::cout << "Too few clms" << std::endl;
            }
            rul += *clm * real_spherical_harmonic(l, m, ul.theta, ul.phi);
            rur += *clm * real_spherical_harmonic(l, m, ur.theta, ur.phi);
            rll += *clm * real_spherical_harmonic(l, m, ll.theta, ll.phi);
            rlr += *clm * real_spherical_harmonic(l, m, lr.theta, lr.phi);
            clm++;
        }
    }

    Vector3 vul = ul.to_vector(rul);
    Vector3 vur = ur.to_vector(rur);
    Vector3 vll = ll.to_vector(rll);
    Vector3 vlr = lr.to_vector(rlr);

    if (be <= 0) {
        Vector3 o = Vector3::zero(); // origin

        tris.push_back(Triangle(*this, ab * vur, ab * vul, ab * vlr));
        tris.push_back(Triangle(*this, ab * vll, ab * vlr, ab * vul));

        tris.push_back(Triangle(*this, ab * vur, ab * vlr, o));
        tris.push_back(Triangle(*this, ab * vlr, ab * vll, o));
        tris.push_back(Triangle(*this, ab * vll, ab * vul, o));
        tris.push_back(Triangle(*this, ab * vul, ab * vur, o));
    }
    else{
        tris.push_back(Triangle(*this, ab * vur, ab * vul, ab * vlr));
        tris.push_back(Triangle(*this, ab * vll, ab * vlr, ab * vul));

        tris.push_back(Triangle(*this, ab * vur, ab * vlr, be * vur));
        tris.push_back(Triangle(*this, ab * vlr, ab * vll, be * vlr));
        tris.push_back(Triangle(*this, ab * vll, ab * vul, be * vll));
        tris.push_back(Triangle(*this, ab * vul, ab * vur, be * vul));

        tris.push_back(Triangle(*this, be * vlr, be * vur, ab * vlr));
        tris.push_back(Triangle(*this, be * vll, be * vlr, ab * vll));
        tris.push_back(Triangle(*this, be * vul, be * vll, ab * vul));
        tris.push_back(Triangle(*this, be * vur, be * vul, ab * vur));

        tris.push_back(Triangle(*this, be * vul, be * vur, be * vll));
        tris.push_back(Triangle(*this, be * vlr, be * vll, be * vur));
    }
}
