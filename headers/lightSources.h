#if !defined(lightSources_h__)
#define lightSources_h__
#include "mathPrimitives.h"
#include "geometryPrimitives.h"
#include "objects.h"
#include <vector>
#include <iostream>
#include <cstring>

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

struct EnvSource {
  Vec *envMap; // takes two coordinate
  uint32_t width;
  uint32_t height;
  Vec radiance;
  uint32_t axis; // Axis around which image is wrapped
  const static uint32_t nSamples = 1000;
  const static uint32_t multiplier = 300;
  Vec impSamples[nSamples * multiplier];
  double normalizationArea;

  EnvSource(Vec *image, uint32_t w, uint32_t h, Vec rad, uint32_t ax): envMap(image), width(w), height(h), radiance(rad), axis(ax) {
    envMap = new Vec[h * w];
    memcpy(envMap, image, sizeof(Vec) * h * w);

    initSamples();
  }

  void getThetaPhi(double x, double y, double z, double &theta, double &phi) {
    theta = acos(z);
    double sintheta = sqrt(1 - z * z);
    if (sintheta < 1e-12) {
      phi = 0.0;
      return;
    }
    x /= sintheta;
    y /= sintheta;

    if (y >= 0 && x >= 0)
      phi = asin(y);
    else if (y >= 0 && x <= 0)
      phi = PI - asin(y);
    else if (y <= 0 && x <= 0)
      phi = PI + asin(-y);
    else
      phi = 2 * PI - asin(-y);
  }

  // Normalized direction
  Vec getRadiance(Vec direction) {
    double theta, phi, x, y, z;
    if (axis == 0) {
      x = direction.y;
      y = direction.z;
      z = direction.x;
    }
    else if (axis == 1) {
      x = direction.x;
      y = direction.z;
      z = direction.y;
    }
    else {
      x = direction.x;
      y = direction.y;
      z = direction.z;
    }
    getThetaPhi(x, y, z, theta, phi); // Wrapped around y axis.
    uint32_t col = (phi / (2 * PI)) * width;
    uint32_t row = (theta / PI) * height;

    if (col >= width && row >= height) {
      std::cout<<"Error:"<<col<<" "<<row<<"";
    }

    return envMap[row * width + col].mult(radiance);
  }
  void initSamples();
  double calcArea(); // Compute normalization factor for pdf.
  void calcMarginalDensity(double *pMarginal, uint32_t mSamples, double area);
  void calcConditionalDensity(double *pConditional, uint32_t cSamples, double *pMarginal, uint32_t mSamples, double area);

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

struct TriLight {
  Triangle t;
  Vec radiance;
  TriLight(const Triangle &t_, const Vec rad): t(t_), radiance(rad){}
};

class LightSource {
  std::vector<SphereSource> sList;
  std::vector<PointSource> pList;
  std::vector<EnvSource> eS;
  std::vector<TriLight> tList;
  static const double eps = 0.01;
public:
  LightSource() {
    sList.reserve(5);
    pList.reserve(5);
    eS.reserve(1);
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
  void addESource(const EnvSource &e) {
    eS.push_back(e);
  }
  void addTSource(const TriLight &t) {
    tList.push_back(t);
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
  Vec getLightFromTriSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
  Vec getLightFromEnvSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
};

extern LightSource *lSource;
extern void configureLightSources();

#endif