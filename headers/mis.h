#if !defined (mis_h__)
#define mis_h__
#include "lightSources.h"
#include "geometryPrimitives.h"

class MIS {
  uint32_t nSamples;
  static const uint32_t nSoAnIS = 10; // No. of samples for solid angle importance sampling.
  static const uint32_t nSuArIS = 10; // No. of samples for surface area importance sampling.
  static const uint32_t nCoIS = 0; // Cosine weighted importance sampling.
  static const uint32_t nBrIS = 60; // No. of phong brdf importance samples.
  static const uint32_t nEnMpIS = 10; // No. of samples for environment map importance sampling.

  SphereSource *sList; // List of Sphere Sources
  uint32_t nSS; // No. of SphereSources

  MeshLight *mList; // List of Mesh Sources
  uint32_t nMS; // No. of Mesh Sources

  TriLight *tList; // List of Triangle Sources
  uint32_t nTS; // No. of Triangle Sources

  EnvSource *eS;
public:
  MIS(uint32_t nSphereSource, SphereSource *_sList, uint32_t nMeshSource, MeshLight *_mList, uint32_t nTriSource, TriLight *_tList, EnvSource *_eS) {
    sList = _sList;
    nSS = nSphereSource;

    mList = _mList;
    nMS = nMeshSource;

    tList = _tList;
    nTS = nTriSource;

    eS = _eS;
  }
  uint32_t getSamples(const Vec &x, const Vec &n, const Vec &wo_ref, const Vec &wo, BasePrimitive * const bPtr, Vec *samples, double *weight);
};
extern MIS *lightSampler;
extern void loadLightSampler();
#endif
