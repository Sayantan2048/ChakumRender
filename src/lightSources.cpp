#include "lightSources.h"
#include "mathPrimitives.h"
#include "shader.h"
#include "geometryPrimitives.h"
#include <iostream>
#include <cmath>
#include "domainSampler.h"

#define SAMPLES 		(10)
#define SAMPLING_TYPE	 	2 // 1 for Uniform hemispherical sampling, 2 Solid Angle Importance Sampling, 4 Light Surface Sampling, 8 cosine weighted sampling.

int nPointSources = 0;
PointSource pSources[] = {
  //PointSource(Vec(50, 10.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 68.6 - .27, 81.6), Vec(80000.0, 80000.0, 80000.0)),
  PointSource(Vec(50, 40.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 28.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(73, 16.6 - .27, 78), Vec(40000.0, 40000.0, 40000.0))
};

int nVolumeSources = 1;

VolumeSource vSources[] = {
 //VolumeSource(Sphere(10, Vec(50, 10.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian), Vec(0.0, 0.0, 40000.0))
 VolumeSource(Vec(10.0, 10.0, 10.0), Sphere(10, Vec(50, 68.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, MaterialType(1.0, 0.5))),
 //VolumeSource(Sphere(10.5, Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, lambertian), Vec(0.0, 0.0, 40000.0))
};

Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Vec refd; // direction from point toward light source
  Ray rr(0,0); // a ray from point toward light source.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;

  Vec sumLight = Vec(0.0, 0.0, 0.0);
  for (int i = 0; i < nPointSources; i++) {
    refd = (pSources[i].p - x).norm(); // direction from point towards source

    rr.d = refd;
    rr.o = x + n * eps;

    if (shadow(rr, (pSources[i].p - x).length()) < 0.5) {
      continue;
    }

    cosine = rr.d.dot(n);

    if (cosine < 0) {
       // This code colors the dark side of an object black.
	continue;
    }

    brdf = primitive->brdf(n, r.d * -1.0, refd, x);

    sumLight = sumLight + (pSources[i].radiance(x) * (cosine * brdf));
  }

  return sumLight;
}

#if SAMPLING_TYPE & 1

Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SAMPLES] = {0};

  // returns 1/pdf.
  double pdf = SphericalSampler::getHemiSurfaceSamples(n, x, SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (int j = 0; j < nVolumeSources; j++) {
    double sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vSources[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	continue;
      }

      brdf = primitive->brdf(n, r.d * -1.0, rr.d, x);

      sum += (cosine * brdf);
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#elif SAMPLING_TYPE & 2
Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[nSAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (int j = 0; j < nVolumeSources; j++) {
    double sum = 0;
    Vec w = vSources[j].p.p - x;
    // returns 1/pdf.
    double pdf = SphericalSampler::getSolidSurfaceSamples(w, x, asin(vSources[j].p.r/w.length()), SAMPLES - 1, samples);
    samples[SAMPLES - 1] = x + w.norm();
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vSources[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	  continue; // Color dark side of an object black.
      }

      brdf = primitive->brdf(n, r.d * -1.0, rr.d, x);

      sum += (cosine * brdf);
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#elif SAMPLING_TYPE & 4
Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[nSAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (int j = 0; j < nVolumeSources; j++) {
    double sum = 0;
    // returns 1/pdf.
    double pdf = SphericalSampler::getLightSurfaceSample(vSources[j].p.p, vSources[j].p.r, x, SAMPLES, samples);

    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x).norm();

      d = (samples[i] - x).length();

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot((samples[i] - vSources[j].p.p).norm());
      cosine1 = cosine1 < 0 ? -cosine1 : 0;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, r.d * -1.0, rr.d, x);

      sum += (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}
#elif SAMPLING_TYPE & 8

Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  double cosine = 0;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SAMPLES] = {0};

  // returns 1/pdf.
  double pdf = SphericalSampler::getCosineSurfaceSamples(n, x, SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (int j = 0; j < nVolumeSources; j++) {
    double sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vSources[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, r.d * -1.0, rr.d, x);

      sum += (brdf);
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#else
Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  double cosine = 0;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  Vec samples[SAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wr =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  // returns 1/pdf.
  double pdf = SphericalSampler::getPhongBRDFSamples(n, wr, x, primitive->m.phongExp, SAMPLES, samples);

  for (int j = 0; j < nVolumeSources; j++) {
    double sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vSources[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      sum += cosine;
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}
#endif