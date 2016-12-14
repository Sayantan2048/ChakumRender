#if !defined (mis_h__)
#define mis_h__
#include "lightSources.h"
#include "geometryPrimitives.h"

class MIS {
  uint32_t nSamples;
  static const uint32_t nSoAnIS = 5; // No. of samples for solid angle importance sampling.
  static const uint32_t nSuArIS = 5; // No. of samples for surface area importance sampling.
  static const uint32_t nCoIS = 2; // Cosine weighted importance sampling.
  static const uint32_t nPhBrIS = 5; // No. of phong brdf importance samples.

  SphereSource *sList; // List of Sphere Sources
  uint32_t nSS; // No. of SphereSources

  PointSource *pList;
  uint32_t nPS;

public:
  MIS(uint32_t nPointSource, PointSource *_pList, uint32_t nSphereSource, SphereSource *_sList) {
    nSamples = nSoAnIS * nSphereSource + nSuArIS * nSphereSource + nCoIS + nPhBrIS + nPointSource;
    sList = _sList;
    nSS = nSphereSource;
    pList = _pList;
    nPS = nPointSource;
  }
  uint32_t getSamples(const Vec &x, const Vec &n, const Vec &wo_ref, const Vec &wo, BasePrimitive * const bPtr, Vec *samples, double *weight);
};
extern MIS *lightSampler;
extern void loadLightSampler();
#endif
