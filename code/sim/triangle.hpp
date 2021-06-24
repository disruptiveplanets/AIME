#pragma once

#include <vector>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "algebra.hpp"

class Chunk;

class ThetaPhi {
public:
    ThetaPhi() {}
    ThetaPhi(Vector3 cart) {
        theta = acos(cart[2] / cart.mag());
        phi = atan2(cart[1], cart[0]);
    }

    Vector3 to_vector(double r) const {
        return r * Vector3({ sin(theta) * cos(phi),
                             sin(theta) * sin(phi),
                             cos(theta) });
    }
    double theta;
    double phi;
};

class Triangle {
public:
    Triangle(Chunk* parent, Vector3 v1, Vector3 v2, Vector3 v3);

    void recenter(Vector3 const& l) {
        v += l;
    }

    void operator*=(Matrix3 const& m);

    Vector3 get_lever_arm() const; // Used for com
    double get_mass() const; // Used for com
    Vector3 get_torque() const; // Used tor torque
    double get_Isame(int a) const; // Used for MOI
    double get_Idiff(int a, int b) const; // Used for MOI
    std::array<Vector3, 3> get_corners() const;
    double get_density() const;
    bool is_edge() const;

private:
    double get_Isame_component(double l1b, double l2b, double vb) const;
    double get_Idiff_component(double l1a, double l1b,
        double l2a, double l2b, double va, double vb) const;
    double get_torque_component(Vector3 const& l1_, Vector3 const& l2_) const;

    Vector3 v;
    Vector3 l1;
    Vector3 l2;
    Vector3 norm;
    Vector3 premul;
    Chunk* parent;
};

class Chunk {
    friend class Triangle;
public:
    Chunk(double alpha_above, double alpha_below,
        int n_here, int face, int a, int b);

    void shape(double density_, uint L, std::vector<double> const& clms,
    std::vector<Triangle>& tris);

protected:
    double density;
    double ab;// alpha_above

private:
    ThetaPhi ul, ur, ll, lr;
    double be;// alpha_below
};
