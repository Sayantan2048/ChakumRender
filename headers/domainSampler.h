#if !defined(domainSampler_h__)
#define domainSampler_h__

#include "mathPrimitives.h"

#define nSAMPLES 1000

class SphericalSampler {
  static Vec *sampleVolume;
  static Vec *sampleSurface;
  static double *randoms; // A bag of random numbers between 0 and 1
  static Vec *cosineSurfaceSamples; // Cosine weighted samples along z-axis as normal and center origin.
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
  static double getSolidSurfaceSamples(Vec w, Vec x, double theta_max, int nSamples, Vec *store);
  static double getLightSurfaceSample(Vec c, double r, Vec x, int nSamples, Vec *store);
  static double getCosineSurfaceSamples(Vec n, Vec x0, int nSamples, Vec *store);
  static void getDistribution(Vec n, Vec x, int nSamples, Vec *samples);
};

#endif