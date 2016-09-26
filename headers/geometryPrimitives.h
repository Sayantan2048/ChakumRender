#if !defined(geometryPrimitives_h__)
#define geometryPrimitives_h__
#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>

enum PrimitiveType {sphere = 0, triangle};

class Sphere {
  MaterialType m;
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

#endif