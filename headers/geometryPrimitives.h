#if !defined(geometryPrimitives_h__)
#define geometryPrimitives_h__
#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>

enum PrimitiveType {sphere = 0, triangle};

class BasePrimitive {
public:
  // color
  Vec c;
  // reflectance
  double reflectance;
  // material type.
  MaterialType m;
  // brdf
  double brdf(Vec n, Vec wo, Vec wi, Vec x) const;
  BasePrimitive(Vec c_, double r_, MaterialType m_): c(c_), reflectance(r_), m(m_) {}
};

class Sphere: public BasePrimitive {
  static const double eps = 1e-4;
public:
  // sphere radius
  double r;
  // position
  Vec p;

  Sphere(double r_, Vec p_, Vec c_, double ref_, MaterialType m_): BasePrimitive(c_, ref_, m_) {
    r = r_;
    p = p_;
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

  Triangle(Vec A_ = dV, Vec B_ = dV, Vec C_ = dV, Vec c_ = dV, double r_ = 1, MaterialType m_ = lambertian): BasePrimitive(c_, r_, m_) {
	  Vec c;
	  A = A_;
	  B = B_;
	  C = C_;
	  n = (B - A);
	  c = (C - A);
	  n = c%n;
	  n.norm();
  }
  double intersect(const Ray &ray, Vec &N) const;
};
#undef dV
#endif