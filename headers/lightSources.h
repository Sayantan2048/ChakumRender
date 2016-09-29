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

struct VolumeSource {
  Sphere p; //This sphere is not visible, if you want a visible sphere, put the sphere used in the constructor argument in the object list as well.
  Vec radiance;
  VolumeSource(Sphere p_, Vec r_): p(p_), radiance(r_) {
    double r = p_.r * (1.0 + eps);
    p = p_;
    p.r = r;
  }
private:
  static const double eps = 0.001;
};

extern Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList, Triangle *triangleList);
extern Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList, Triangle *triangleList);

#endif