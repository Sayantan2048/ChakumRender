#if !defined(domainSampler_h__)
#define domainSampler_h__

#include "mathPrimitives.h"

#define nSAMPLES 1000

class SphericalSampler {
  static Vec *sampleVolume;
  static Vec *sampleSurface;
  static int __initSamples;

  static int initSamples();
  SphericalSampler();

  public:
  // Vec n is direction of normal and Vec x is position of center.
  static double getSphericalVolumeSamples(Vec x, int nSamples, Vec *store);
  static double getSphericalSurfaceSamples(Vec x, int nSamples, Vec *store);
  static double getHemiSurfaceSamples(Vec n, Vec x, int nSamples, Vec *store);
  static double getHemiSurfaceSamplesTrue(Vec n, Vec x, int nSamples, Vec *store);
  static double getHemiVolumeSamples(Vec n, Vec x, int nSamples, Vec *store);
  static double getArcSurfaceSamples(Vec w, Vec x, double theta_max, int nSamples, Vec *store);
  static void getDistribution(Vec n, Vec x, int nSamples, Vec *samples);
};

#endif