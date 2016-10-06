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
  double brdf(Vec wi, Vec wr, Vec x) const;
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

class Triangle: public BasePrimitive {
  static const double eps = 1e-10;
public:
  // Vertices
  Vec A, B, C;
  // normal
  Vec n;

  Triangle(Vec A_, Vec B_, Vec C_, Vec c_, double r_, MaterialType m_): BasePrimitive(c_, r_, m_) {
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

#endif