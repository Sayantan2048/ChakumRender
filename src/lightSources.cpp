#include "lightSources.h"
#include "mathPrimitives.h"
#include "shader.h"
#include "geometryPrimitives.h"
#include <cstdio>
#include <cmath>
#include "domainSampler.h"

#define SAMPLES 		(100)
#define SAMPLING_TYPE	 	2 // 1 for Uniform hemispherical sampling, 2 Solid Angle Importance Sampling, 4 Light Surface Sampling.

int nPointSources = 0;
PointSource pSources[] = {
  //PointSource(Vec(50, 10.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 68.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 40.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 28.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(73, 16.6 - .27, 78), Vec(40000.0, 40000.0, 40000.0))
};

int nVolumeSources = 1;

VolumeSource vSources[] = {
 //VolumeSource(Sphere(10, Vec(50, 10.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian), Vec(0.0, 0.0, 40000.0))
 VolumeSource(Sphere(10, Vec(50, 68.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian), Vec(40000.0, 40000.0, 40000.0)),
 //VolumeSource(Sphere(10.5, Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, lambertian), Vec(0.0, 0.0, 40000.0))
};

Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList, Triangle *triangleList) {
  int convex = 1; // check whether the surface is convex or not from the direction of viewing
  Vec refd; // direction from point toward light source
  Ray rr(0,0); // a ray from point toward light source.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;

  // test for convex or concave as seen from a point.
  if (r.d.dot(n) >= 0)
    convex = 0;

  Vec sumLight = Vec(0.0, 0.0, 0.0);
  for (int i = 0; i < nPointSources; i++) {
    refd = (pSources[i].p - x).norm(); // direction from point towards source

    rr.d = refd;
    rr.o = x + n * (convex?eps:-eps);

    if (shadow(rr, (pSources[i].p - x).length()) < 0.5) {
      continue;
    }

    cosine = rr.d.dot(n);

    if (cosine < 0) {
    cosine = -cosine;
      if (convex) { // This code colors the dark side of an object black.
	continue;
      }
    }

    if (m == sphere)
      brdf = sphereList[id].brdf(rr.d, rr.d, rr.d);
    else if (m == triangle)
      brdf = triangleList[id].brdf(rr.d, rr.d, rr.d);
    else
	fprintf(stderr, "PrimitiveType not found in getLightFromPointSources.\n");
    /* TODO: ADD BRDF FOR OTHER PrimitiveTypes */
    sumLight = sumLight + (pSources[i].radiance(x) * (cosine * brdf));
  }

  return sumLight;
}

#if SAMPLING_TYPE & 1

Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList, Triangle *triangleList) {
  int convex = 1; // check whether the surface is convex or not from the direction of viewing
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  static Vec samples[SAMPLES] = {0};

  // test for convex or concave as seen from a point.
  if (r.d.dot(n) >= 0) {
      convex = 0;
  }

  // returns 1/pdf.
  double pdf = SphericalSampler::getHemiSurfaceSamples(convex?n:(n*-1.0), x, SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * (convex?eps:-eps);

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
	cosine = -cosine;
	if (convex) { // This code colors the dark side of an object black.
	  continue;
	}
      }

      if (m == sphere)
	brdf = sphereList[id].brdf(rr.d, rr.d, rr.d);
      else if (m == triangle)
	brdf = triangleList[id].brdf(rr.d, rr.d, rr.d);
      else
	fprintf(stderr, "PrimitiveType not found in getLightFromPointSources.\n");
      /* TODO: ADD BRDF FOR OTHER PrimitiveTypes */
      sum += (cosine * brdf);
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#elif SAMPLING_TYPE & 2
Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList, Triangle *triangleList) {
  int convex = 1; // check whether the surface is convex or not from the direction of viewing
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  static Vec samples[nSAMPLES] = {0};

  // test for convex or concave as seen from a point.
  if (r.d.dot(n) >= 0) {
      convex = 0;
  }

  Vec sumLight = Vec();
  rr.o = x + n * (convex?eps:-eps);

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

      //double cosine1 = rr.d.dot((rr.o + rr.d * d - vSources[j].p.p).norm());
      //cosine1 = cosine1 < 0 ? -cosine1 : 0;

      if (cosine < 0) {
	cosine = -cosine;
	if (convex) { // This code colors the dark side of an object black.
	  continue;
	}
      }

      if (m == sphere)
	brdf = sphereList[id].brdf(rr.d, rr.d, rr.d);
      else if (m == triangle)
	brdf = triangleList[id].brdf(rr.d, rr.d, rr.d);
      else
	fprintf(stderr, "PrimitiveType not found in getLightFromPointSources.\n");
      /* TODO: ADD BRDF FOR OTHER PrimitiveTypes */
      sum += (cosine * brdf);
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#else
Vec getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList, Triangle *triangleList) {
  int convex = 1; // check whether the surface is convex or not from the direction of viewing
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  static Vec samples[nSAMPLES] = {0};

  // test for convex or concave as seen from a point.
  if (r.d.dot(n) >= 0) {
      convex = 0;
  }

  Vec sumLight = Vec();
  rr.o = x + n * (convex?eps:-eps);

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
	cosine = -cosine;
	if (convex) { // This code colors the dark side of an object black.
	  continue;
	}
      }

      if (m == sphere)
	brdf = sphereList[id].brdf(rr.d, rr.d, rr.d);
      else if (m == triangle)
	brdf = triangleList[id].brdf(rr.d, rr.d, rr.d);
      else
	fprintf(stderr, "PrimitiveType not found in getLightFromPointSources.\n");
      /* TODO: ADD BRDF FOR OTHER PrimitiveTypes */
      sum += (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + vSources[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#endif