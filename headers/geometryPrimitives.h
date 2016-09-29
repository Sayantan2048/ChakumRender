#if !defined(geometryPrimitives_h__)
#define geometryPrimitives_h__
#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>

enum PrimitiveType {sphere = 0, triangle};

class Sphere {
  MaterialType m;
  static const double eps = 1e-4;
public:
  // sphere radius
  double r;
  // position
  Vec p;
  // color
  Vec c;
  // reflectance of the sphere
  double reflectance;

  Sphere(double r_, Vec p_, Vec c_, double reflectance_, MaterialType m_) {
    r = r_;
    p = p_;
    c = c_;
    reflectance = reflectance_;
    m = m_;
  }

  double intersect(const Ray &ray) const;
  double brdf(Vec wi, Vec wr, Vec x) const;
};

class Triangle {
  MaterialType m;
  static const double eps = 1e-10;
public:
  // Vertices
  Vec A, B, C;
  // normal
  Vec n;
  // color
  Vec c;
  // reflectance of the triangle
  double reflectance;

  Triangle(Vec A_, Vec B_, Vec C_, Vec c_, double ref_, MaterialType m_) {
	  A = A_;
	  B = B_;
	  C = C_;
	  c = c_;
	  reflectance = ref_;
	  m = m_;
	  n = (B - A);
	  n = (C - A)%n;
	  n.norm();
	}
  double intersect(const Ray &ray, Vec &N) const;
  double brdf(Vec wi, Vec wr, Vec x) const;
};

#endif