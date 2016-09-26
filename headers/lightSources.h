#if !defined(lightSources_h__)
#define lightSources_h__
#include "mathPrimitives.h"
#include "geometryPrimitives.h"

struct PointSource {
  Vec p;
  Vec power;
  PointSource(Vec p_, Vec pw_): p(p_), power(pw_){}
  Vec radiance(Vec x) {
    double sqlen = (p - x).dot(p-x);
    return power / (4.0 * PI * sqlen);
  }
};

extern PointSource pSources[];
extern Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList);

#endif