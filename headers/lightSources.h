#if !defined(lightSources_h__)
#define lightSources_h__
#include "mathPrimitives.h"
#include "geometryPrimitives.h"
#include "objects.h"
#include <vector>

struct PointSource {
  Vec p;
  Vec power;
  PointSource(Vec p_, Vec pw_): p(p_), power(pw_){}
  // At any point x
  Vec radiance(Vec x) {
    double sqlen = (p - x).dot(p-x);
    return power / (4.0 * PI * sqlen);
  }
};

struct SphereSource {
  Vec radiance;
  Sphere p; //This sphere is not visible, if you want a visible sphere, put the sphere used in the constructor argument in the object list as well but with slightly smaller radii.
  Vec power;
  SphereSource(Sphere p_, Vec pw_): p(p_), power(pw_) {
    radiance = power / (4 * PI * PI * p.r * p.r);
  }

  SphereSource(Vec r_, Sphere p_): radiance(r_), p(p_) {
    power = radiance * ( 4 * PI * PI * p.r * p.r);
  }
  SphereSource() {}
};

class LightSource {
  std::vector<SphereSource> sList;
  std::vector<PointSource> pList;
  static const double eps = 0.01;
public:
  LightSource() {
    sList.reserve(5);
    pList.reserve(5);
  }
  void addSSource(const SphereSource &v) {
    sList.push_back(v);
    Sphere p = v.p;
    p.r /= (1 + eps);
    p.m.radiance = v.radiance;
    vSphereList.push_back(p); // For a visible sphere light.
  }
  void addPSource(const PointSource &p) {
    pList.push_back(p);
  }
  SphereSource* getSphereSourcePtr(uint32_t &nSources) {
    nSources = sList.size();
    SphereSource *ptr = new SphereSource[sList.size()];
    for (uint32_t i = 0; i < sList.size(); i++)
      ptr[i] = sList[i];

    return ptr;
  }
  Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
  Vec getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
};

extern LightSource *lSource;
extern void configureLightSources();

#endif