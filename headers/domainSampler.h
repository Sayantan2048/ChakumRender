#if !defined(domainSampler_h__)
#define domainSampler_h__

#include "mathPrimitives.h"
#include "materialTypes.h"
#include <stdint.h>
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
  static double getSphericalVolumeSamples(const Vec &x, const int nSamples, Vec *store);
  static double getSphericalSurfaceSamples(const Vec &x, const int nSamples, Vec *store);
  static double getHemiDirectionSamples(const Vec &n, const int nSamples, Vec *store); //Assumes normalized n
  static double getHemiDirectionSamplesTrue(const Vec &n, const int nSamples, Vec *store); //Assumes normalized n
  static double getHemiVolumeSamples(const Vec &n, const Vec &x, int nSamples, Vec *store); //Assumes normalized n
  static double getSolidDirectionSamples(const Vec &w, const double theta_max, const int nSamples, Vec *store); // Assumes normalized w
  static double getLightDirectionSamples(const Vec &c, const double r, const Vec &x, const int nSamples, Vec *store);
  static double getCosineDirectionSamples(const Vec &n, const int nSamples, Vec *store); // Assumes normalized n.
  static double getClassicPhongDirectionSamples(const Vec &n,const Vec &w, const double e, const int nSamples, Vec *store); // Assumes normalized n, w
  static double getBrdfDirectionSamples(const Vec &n, const Vec &wo, const double alpha, const int nSamples, Vec *store, const BRDFType b); // Assumes normalized n, wo
  static void getTriLightSurfaceSamples(const Vec &p1, const Vec &p2, const Vec &p3, const uint32_t nSamples, Vec *store);//Area of the domain is simply area of triangle.
  static void getDistribution(const Vec &n, const Vec &x, const int nSamples, Vec *samples);
};

#endif