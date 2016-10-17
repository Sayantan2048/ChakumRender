#include "lightSources.h"
#include "mathPrimitives.h"
#include "shader.h"
#include "objects.h"
#include "geometryPrimitives.h"
#include <iostream>
#include <cmath>
#include <stdint.h>
#include "domainSampler.h"

#define SAMPLES 		(10)
#define SAMPLING_TYPE	 	4 // 1 for Uniform hemispherical sampling, 2 Solid Angle Importance Sampling, 4 Light Surface Sampling, 8 cosine weighted sampling.
LightSource *lSource;
/*
int nPointSources = 0;
PointSource pSources[] = {
  //PointSource(Vec(50, 10.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 68.6 - .27, 81.6), Vec(80000.0, 80000.0, 80000.0)),
  PointSource(Vec(50, 40.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 28.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(73, 16.6 - .27, 78), Vec(40000.0, 40000.0, 40000.0))
};*/
/*
int nVolumeSources = 0;

VolumeSource vSources[] = {
 //VolumeSource(Sphere(10, Vec(50, 10.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian), Vec(0.0, 0.0, 40000.0))
 VolumeSource(Vec(10.0, 10.0, 10.0), Sphere(10, Vec(50, 68.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, MaterialType(1.0, 0.5))),
 //VolumeSource(Sphere(10.5, Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, lambertian), Vec(0.0, 0.0, 40000.0))
};*/

void configureLightSources() {
  lSource = new LightSource();
  //lSource->addPSource(PointSource(Vec(50, 68.6 - .27, 81.6), Vec(80000.0, 80000.0, 80000.0)));
  //lSource->addPSource(PointSource(Vec(50, 40.6 - .27, 81.6), Vec(80000.0, 80000.0, 80000.0)));
  lSource->addVSource(VolumeSource(Vec(10.0, 0.0, 0.0), Sphere(10, Vec(20, 40.6 - .27, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(1.0, 0.5, VOLUME, Vec(0., 0., 0.)))));
  lSource->addVSource(VolumeSource(Vec(0.0, 10.0, 0.0), Sphere(10, Vec(50, 40.6 - .27, 81.6), Vec(.0, 0.0, .0), 1.0, MaterialType(1.0, 0.5, VOLUME, Vec(0., 0., 0.)))));
  lSource->addVSource(VolumeSource(Vec(0.0, 0.0, 10.0), Sphere(10, Vec(35, 55.6 - .27, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(1.0, 0.5, VOLUME, Vec(0., 0., 0.)))));
}

Vec LightSource::getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Vec refd; // direction from point toward light source
  Ray rr(0,0); // a ray from point toward light source.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;

  Vec sumLight = Vec(0.0, 0.0, 0.0);
  for (uint32_t i = 0; i < pList.size(); i++) {
    refd = (pList[i].p - x).norm(); // direction from point towards source

    rr.d = refd;
    rr.o = x + n * eps;

    if (shadow(rr, (pList[i].p - x).length()) < 0.5) {
      continue;
    }

    cosine = rr.d.dot(n);

    if (cosine < 0) {
       // This code colors the dark side of an object black.
	continue;
    }

    brdf = primitive->brdf(n, r.d * -1.0, refd, x);

    sumLight = sumLight + (pList[i].radiance(x) * (cosine * brdf));
  }

  return sumLight;
}

#if SAMPLING_TYPE & 1

Vec LightSource::getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SAMPLES] = {0};

  // returns 1/pdf.
  double pdf = SphericalSampler::getHemiSurfaceSamples(n, x, SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (uint32_t j = 0; j < vList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vList[j].p.intersect(rr);

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
    sumLight = sumLight + vList[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#elif SAMPLING_TYPE & 2
Vec LightSource::getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[nSAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (uint32_t j = 0; j < vList.size(); j++) {
    double sum = 0;
    Vec w = vList[j].p.p - x;
    double sinThetaMax = vList[j].p.r/w.length();
    // returns 1/pdf.
    double pdf = SphericalSampler::getSolidSurfaceSamples(w, x, asin(sinThetaMax), SAMPLES - 1, samples);
    samples[SAMPLES - 1] = x + w.norm();
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vList[j].p.intersect(rr);

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
    sumLight = sumLight + vList[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#elif SAMPLING_TYPE & 4
Vec LightSource::getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[nSAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (uint32_t j = 0; j < vList.size(); j++) {
    double sum = 0;
    // returns 1/pdf.
    double pdf = SphericalSampler::getLightSurfaceSample(vList[j].p.p, vList[j].p.r, x, SAMPLES, samples);

    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x).norm();

      d = (samples[i] - x).length();

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot((samples[i] - vList[j].p.p).norm());
      cosine1 = cosine1 < 0 ? -cosine1 : 0;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, r.d * -1.0, rr.d, x);

      sum += (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + vList[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}
#elif SAMPLING_TYPE & 8

Vec LightSource::getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  double cosine = 0;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SAMPLES] = {0};

  // returns 1/pdf.
  double pdf = SphericalSampler::getCosineSurfaceSamples(n, x, SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (uint32_t j = 0; j < vList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vList[j].p.intersect(rr);

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
    sumLight = sumLight + vList[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}

#else
Vec LightSource::getLightFromVolumeSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  double cosine = 0;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  Vec samples[SAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wr =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  // returns 1/pdf.
  double pdf = SphericalSampler::getPhongBRDFSamples(n, wr, x, primitive->m.phongExp, SAMPLES, samples);

  for (uint32_t j = 0; j < vList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = vList[j].p.intersect(rr);

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
    sumLight = sumLight + vList[j].radiance * ( pdf * sum/SAMPLES);
  }
  return sumLight;
}
#endif