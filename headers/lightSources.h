#if !defined(lightSources_h__)
#define lightSources_h__
#include "mathPrimitives.h"
#include "geometryPrimitives.h"
#include "objects.h"
#include <vector>
#include <iostream>
#include <cstring>

struct BaseSource {
  Vec radiance;
  BaseSource(const Vec &col, const double &intensity) {
    Vec color = col;
    color.maxnorm();
    radiance = color * intensity;
  }
  BaseSource(double temperature, const double &intensity) {
    radiance = tempToColor(temperature) * intensity;
  }
};

struct PointSource: public BaseSource {
  Vec p;
  PointSource(Vec p_, const Vec &col, const double &intensity): BaseSource(col, intensity), p(p_){}
  PointSource(Vec p_, double temperature, const double &intensity): BaseSource(temperature, intensity), p(p_){}
  // At any point x
  Vec getRadiance(Vec x) {
    double sqlen = (p - x).dot(p-x);
    return radiance / (4.0 * PI * sqlen); // Power to radiance conversion. radiance is the name of variable that stores power!!
  }
};

struct EnvSource: public BaseSource {
  Vec *envMap; // takes two coordinate
  uint32_t width;
  uint32_t height;
  uint32_t axis; // Axis around which image is wrapped
  const static uint32_t nSamples = 2000;
  const static uint32_t multiplier = 100;
  Vec impSamples[nSamples * multiplier]; //bad idea!!
  double normalizationArea;

  EnvSource(Vec *image, uint32_t w, uint32_t h, uint32_t ax, const Vec &col, const double &intensity): BaseSource(col, intensity), envMap(image), width(w), height(h), axis(ax) {
    envMap = new Vec[h * w];
    memcpy(envMap, image, sizeof(Vec) * h * w);

    initSamples();
  }
  EnvSource(Vec *image, uint32_t w, uint32_t h, uint32_t ax, double temperature, const double &intensity): BaseSource(temperature, intensity), envMap(image), width(w), height(h), axis(ax) {
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
    getThetaPhi(x, y, z, theta, phi);
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

struct SphereSource : public BaseSource {
  Sphere p; //This sphere is not visible, if you want a visible sphere, put the sphere used in the constructor argument in the object list as well but with slightly smaller radii.
  SphereSource(const double &power, double radius, const Vec &center, const Vec &col): BaseSource(col, power) {
    radiance = radiance / (4 * PI * PI * p.r * p.r); // power to radiance conversion
    p = Sphere(radius, center, Vec(0,0,0), 0, MaterialType(VOLUME, radiance));
  }
  SphereSource(const double &power, double radius, const Vec &center, double temperature): BaseSource(temperature, power) {
    radiance = radiance / (4 * PI * PI * p.r * p.r); // power to radiance conversion
    p = Sphere(radius, center, Vec(0,0,0), 0, MaterialType(VOLUME, radiance));
  }
  SphereSource(double radius, const Vec &center, const Vec &col, const double &intensity): BaseSource(col, intensity) {
     p = Sphere(radius, center, Vec(0,0,0), 0, MaterialType(VOLUME, radiance));
  }
  SphereSource(double radius, const Vec &center, double temperature, const double &intensity): BaseSource(temperature, intensity) {
     p = Sphere(radius, center, Vec(0,0,0), 0, MaterialType(VOLUME, radiance));
  }
  SphereSource():BaseSource(0, 0) {}
};

struct TriLight : public BaseSource {
  Triangle t;
  TriLight(const Vec &A, const Vec &B, const Vec &C, const Vec &col, const double &intensity): BaseSource(col, intensity) {
    t = Triangle(A, B, C, Vec(0, 0, 0), 0.0,  MaterialType(PLANAR, radiance));
  }
  TriLight(const Vec &A, const Vec &B, const Vec &C, double temperature, const double &intensity): BaseSource(temperature, intensity) {
    t = Triangle(A, B, C, Vec(0, 0, 0), 0.0,  MaterialType(PLANAR, radiance));
  }
  TriLight():BaseSource(0,0){}
};

struct MeshLight {
  std::vector<TriLight> mesh;
  void add(const TriLight &t) {
    mesh.push_back(t);
  }
  void initMeshLight();
  double getSamples(uint32_t nSamples, Vec *store, uint32_t *ids);
  bool intersect(const Ray &r, uint32_t &id, double &t, Vec &N);
private:
  double area; //Total surface area of mesh light.
  const static uint32_t invCDFBins = 300;
  uint32_t invCDF[invCDFBins]; // Roll a dice between 0 to invCDFBins to get a Triangle.
};

class LightSource {
  std::vector<SphereSource> sList;
  std::vector<PointSource> pList;
  std::vector<EnvSource> eS;
  std::vector<TriLight> tList;
  std::vector<MeshLight> mList;
  static const double eps = 0.01;

  double integrateEdge(const Vec &v1, const Vec &v2, const Vec &n);
  double integrateEdge(const Vec &v1, const Vec &v2); // Simplified for n = Vec(0, 0, 1)
  bool clipTriangle(const Vec &A, const Vec &B,  const Vec &C, const Vec &n, const Vec &P, const double nA, const double nB, const double nC, Vec *clip);
  bool clipTriangle(const Vec &A, const Vec &B,  const Vec &C, const double nA, const double nB, const double nC, Vec *clip); //Simplified for n = Vec(0, 0, 1) and P = Vec(0, 0, 0)
  Vec getLightAnalyticDiffuse(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
  Vec getLightAnalytic(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const BRDFApprox &brdfApprox);

  //Experimental
  Vec getLightAnalyticDiffuseAdv(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, uint32_t &confidence);
public:
  LightSource() {
    sList.reserve(5);
    pList.reserve(5);
    eS.reserve(1);
  }
  void addSSource(const SphereSource &v) {
    sList.push_back(v);
    Sphere p = v.p;
    vSphereList.push_back(p); // For a visible sphere light.
  }
  void addPSource(const PointSource &p) {
    pList.push_back(p);
  }
  void addESource(const EnvSource &e) {
    if (eS.size() < 1)
      eS.push_back(e);
    else
      std::cout<<"Already added an Environment Source\n";
  }
  void addTSource(const TriLight &t) {
    tList.push_back(t);
    Triangle tt = t.t;
    vTriangleList.push_back(tt);
  }
  void addMSource(const MeshLight &m) {
    mList.push_back(m);
    for (uint32_t i = 0; i < m.mesh.size(); i++) {
      Triangle tt = m.mesh[i].t;
      vTriangleList.push_back(tt);
    }
  }
  SphereSource* getSphereSourcePtr(uint32_t &nSources) {
    nSources = sList.size();
    SphereSource *ptr = new SphereSource[sList.size()];
    for (uint32_t i = 0; i < sList.size(); i++)
      ptr[i] = sList[i];

    return ptr;
  }
  MeshLight* getMeshLightPtr(uint32_t &nSources) {
    nSources = mList.size();
    MeshLight *ptr = new MeshLight[mList.size()];
    for (uint32_t i = 0; i < mList.size(); i++)
      ptr[i] = mList[i];

    return ptr;
  }
  TriLight* getTriLightPtr(uint32_t &nSources) {
    nSources = tList.size();
    TriLight *ptr = new TriLight[tList.size()];
    for (uint32_t i = 0; i < tList.size(); i++)
      ptr[i] = tList[i];

    return ptr;
  }
  EnvSource* getEnvSourcePtr(){
    if (eS.size() == 1)
      return &eS[0];
    else
      return 0;
  }
  Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
  Vec getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromTriSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromMeshSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromEnvSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromMeshSource_AnalyticCV(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromMeshSource_Analytic(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);

  Vec getLightFromAllSources_MIS(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive);
  Vec getLightFromAllSources_BRDF(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromAllSources_UNIFORM(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
  Vec getLightFromAllSources_COS(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);

  //Experimental
  Vec getLightFromMeshSource_CVAdv(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples);
};

extern LightSource *lSource;

#endif