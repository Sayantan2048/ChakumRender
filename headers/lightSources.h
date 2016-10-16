#if !defined(lightSources_h__)
#define lightSources_h__
#include "mathPrimitives.h"
#include "geometryPrimitives.h"
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

struct VolumeSource {
  Vec radiance;
  Sphere p; //This sphere is not visible, if you want a visible sphere, put the sphere used in the constructor argument in the object list as well.
  Vec power;
  VolumeSource(Sphere p_, Vec pw_): p(p_), power(pw_) {
    double r = p_.r * (1.0 + eps);
    p = p_;
    p.r = r;
    radiance = power / (4 * PI * PI * p.r * p.r);
  }

  VolumeSource(Vec r_, Sphere p_): radiance(r_), p(p_) {
    double r = p_.r * (1.0 + eps);
    p = p_;
    p.r = r;
    power = radiance * ( 4 * PI * PI * p.r * p.r);
  }


private:
  static const double eps = 0.001;
};

class LightSource {
  std::vector<VolumeSource> vList;
  std::vector<PointSource> pList;

public:
  LightSource() {
    vList.reserve(5);
    pList.reserve(5);
  }
  void addVSource(const VolumeSource &v) {
    vList.push_back(v);
  }
  void addPSource(const PointSource &p) {
    pList.push_back(p);
  }
  Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
  Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
};

extern LightSource *lSource;

extern Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
extern Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);

#endif