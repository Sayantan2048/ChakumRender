#include "lightSources.h"
#include "domainSampler.h"

class MIS {
  uint32_t nSamples;
  static const uint32_t nSoAnIS = 5; // No. of samples for solid angle importance sampling.
  static const uint32_t nSuArIS = 10; // No. of samples for surface area importance sampling.
  static const uint32_t nCoIS = 5; // Cosine weighted importance sampling.
  static const uint32_t nPhBrIS = 3; // No. of phong brdf importance samples.

  void getSolidAngleSamples(uint32_t &offset, Vec *storeSample, Vec *storePdf) {
  }

public:
  MIS(uint32_t nPointSource, uint32_t nSphereSource) {
    nSamples = nSoAnIS * nSphereSource + nSuArIS * nSphereSource + nCoIS + nPhBrIS + nPointSource;
  }

  // Return weight of the samples
  double * getSamples(Vec n, Vec wo, uint32_t &nSamples, double *store) {
  }

}
