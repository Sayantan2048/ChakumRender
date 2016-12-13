#if !defined(geometryPrimitives_h__)
#define geometryPrimitives_h__

#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>

enum PrimitiveType {sphere = 0, triangle};

//Axis aligned bounding box.
struct AABBox {
  Vec pMin, pMax;
  AABBox(const Vec &p): pMin(p), pMax(p){}
  AABBox(const Vec &p1, const Vec &p2);
  AABBox(const double &r, const Vec &p);
  static AABBox uNion(const AABBox &b, const Vec &p);
  static AABBox uNion(const AABBox &a, const AABBox &b);
  double intersect(const Ray &ray) const;
  int maximumExtent() const;
};

class BasePrimitive {
public:
  static const double timeStep = 1.0;
  // color
  Vec c;
  // reflectance
  double reflectance;
  // material type.
  MaterialType m;
  // Axis Aligned bounding box.
  AABBox box;
  // brdf
  double brdf(const Vec &n, const Vec &wo_ref, const Vec &wo, const Vec &wi, const BRDFApprox &brdfApprox = BRDFApprox()) const;
  BasePrimitive(Vec c_, double r_, MaterialType m_, AABBox b): c(c_), reflectance(r_), m(m_), box(b) {}
};

class Sphere: public BasePrimitive {
  static const double eps = 1e-4;
public:
  // sphere radius
  double r;
  // position
  Vec p;

  Sphere(double r_ = 1.0, Vec p_ = Vec(0., 0., 0.),
	 Vec c_ = Vec(0., 0., 0.), double ref_ = 0,
	 MaterialType m_ = MaterialType()): BasePrimitive(c_, ref_, m_, AABBox(r_, p_)) {
    r = r_;
    p = p_;
  }
  void updatePos(Vec vel = Vec(0., 0., 0.), Vec acc = Vec(0.,0.,0.)) {
    Vec d = vel * timeStep + acc * 0.5 * timeStep * timeStep;
    p = p + d;
    box = AABBox(p, r);
  }
  double intersect(const Ray &ray) const;
};

#define dV Vec(1., 1., 1.)
class Triangle: public BasePrimitive {
  static const double eps = 1e-10;
public:
  // Vertices
  Vec A, B, C;
  // normal
  Vec n;
  double area;
  Triangle(Vec A_ = dV, Vec B_ = dV, Vec C_ = dV, Vec c_ = dV, double r_ = 1, MaterialType m_ = MaterialType()):
    BasePrimitive(c_, r_, m_, AABBox::uNion(AABBox(A_, B_), C_)) {
	  Vec c;
	  A = A_;
	  B = B_;
	  C = C_;
	  n = (B - A);
	  c = (C - A);
	  n = c%n;
	  area = 0.5 * n.length();
	  n.norm();
  }
  void updatePos(Vec vel = Vec(0., 0., 0.), Vec acc = Vec(0.,0.,0.)) {
    Vec d = vel * timeStep + acc * 0.5 * timeStep * timeStep;

    A = A + d;
    B = B + d;
    C = C + d;

    Vec c = (B - A);
    n = (C - A);
    n = c%n;
    area = 0.5 * n.length();
    n.norm();

    box = AABBox::uNion(AABBox(A, B), C);
  }
  double intersect(const Ray &ray, Vec &N) const;
};
#undef dV
#endif