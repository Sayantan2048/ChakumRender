#include "lightSources.h"
#include "mathPrimitives.h"
#include "shader.h"
#include "objects.h"
#include "geometryPrimitives.h"
#include "mis.h"
#include "random.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <stdint.h>
#include "domainSampler.h"
#include "ppm.h"

#define SP_SAMPLES 		(20) // Should be less than the No. of samples in domainSampler.h
#define SP_SAMPLING_TYPE	 2// 1 for Uniform hemispherical sampling, 2 Solid Angle Importance Sampling, 4 Light Surface Sampling, 8 cosine weighted sampling.
//16 for phong BRDF sampling(works only if all objects are pure phong!!).
#define ENV_SAMPLES		1000 // Should be less than the No. of samples in domainSampler.h and nSamples in struct EnvSource
#define ENV_SAMPLING_TYPE	1// 1 for Uniform hemispherical sampling, 2 Importance sampling.
static double to_greyScale(Vec L);
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
  //lSource->addSSource(SphereSource(Vec(10.0, 0.0, 0.0), Sphere(10, Vec(20, 40.6 - .27, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(1.0, 0.5, VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(0.0, 10.0, 0.0), Sphere(10, Vec(50, 40.6 - .27, 81.6), Vec(.0, 0.0, .0), 1.0, MaterialType(1.0, 0.5, VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(10, Vec(50, 68.6 - .27, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(1.0, 0.0, VOLUME, Vec(0., 0., 0.)))));

  uint32_t w, h;
  Vec *image;
  readImage(image, w, h);
  lSource->addESource(EnvSource(image, w, h, Vec(1,1,1), 1));
  //delete image.
}

Vec LightSource::getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Vec refd; // direction from point toward light source
  Ray rr(0,0); // a ray from point toward light source.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

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

    brdf = primitive->brdf(n, wo_ref, refd, x);

    sumLight = sumLight + (pList[i].radiance(x) * (cosine * brdf));
  }

  return sumLight;
}
#if ENV_SAMPLING_TYPE & 1
Vec LightSource::getLightFromEnvSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  if (eS.size() < 1) return Vec();
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[ENV_SAMPLES] = {0};

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  // returns 1/pdf.
  double pdf = SphericalSampler::getHemiSurfaceSamples(n, x, ENV_SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (int i = 0; i < ENV_SAMPLES; i++) {
    rr.d = (samples[i] - x);

    if (shadow(rr, INF) < 0.5)
      continue;

    cosine = rr.d.dot(n);

    if (cosine < 0) {
	continue;
    }

    brdf = primitive->brdf(n, wo_ref, rr.d, x);

    sumLight = sumLight + eS[0].getRadiance(rr.d) * cosine * brdf;
  }
  sumLight = sumLight * (pdf/ENV_SAMPLES);

  return sumLight;
}
#else
Vec LightSource::getLightFromEnvSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  if (eS.size() < 1) return Vec();
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  Vec sumLight = Vec();
  rr.o = x;
  uint32_t offset = randomMTD(0, 299) * 1000;
  for (int i = 0; i < ENV_SAMPLES; i++) {
    rr.d = eS[0].impSamples[i + offset];

    if (shadow(rr, INF) < 0.5)
      continue;

    cosine = rr.d.dot(n);

    if (cosine < 0) {
	continue;
    }

    brdf = primitive->brdf(n, wo_ref, rr.d, x);
    Vec L = eS[0].getRadiance(rr.d);
    double pdf = to_greyScale(L);
    //std::cout<<pdf<<"\n";
    if(pdf > eps)
      sumLight = sumLight + (L * cosine * brdf)/pdf;
  }
  //std::cout<<sumLight.x<<" "<<sumLight.y<<" "<<sumLight.z<<" \n";
  sumLight = sumLight * (eS[0].normalizationArea/(ENV_SAMPLES));

  return sumLight;
}
#endif

#if SP_SAMPLING_TYPE & 1

Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SP_SAMPLES] = {0};

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  // returns 1/pdf.
  double pdf = SphericalSampler::getHemiSurfaceSamples(n, x, SP_SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < SP_SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = sList[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	continue;
      }

      brdf = primitive->brdf(n, wo_ref, rr.d, x);

      sum += (cosine * brdf);
    }
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/SP_SAMPLES);
  }
  return sumLight;
}

#elif SP_SAMPLING_TYPE & 2
Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SP_SAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    Vec w = sList[j].p.p - x;
    double sinThetaMax = sList[j].p.r/w.length();
    // returns 1/pdf.
    double pdf = SphericalSampler::getSolidSurfaceSamples(w, x, asin(sinThetaMax), SP_SAMPLES - 1, samples);
    samples[SP_SAMPLES - 1] = x + w.norm();
    for (int i = 0; i < SP_SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = sList[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	  continue; // Color dark side of an object black.
      }

      brdf = primitive->brdf(n, wo_ref, rr.d, x);

      sum += (cosine * brdf);
    }
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/SP_SAMPLES);
  }
  return sumLight;
}

#elif SP_SAMPLING_TYPE & 4
Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[nSP_SAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    // returns 1/pdf.
    double pdf = SphericalSampler::getLightSurfaceSample(sList[j].p.p, sList[j].p.r, x, SP_SAMPLES, samples);

    for (int i = 0; i < SP_SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x).norm();

      d = (samples[i] - x).length();

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot((samples[i] - sList[j].p.p).norm());
      cosine1 = cosine1 < 0 ? -cosine1 : 0;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, wo_ref, rr.d, x);

      sum += (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/SP_SAMPLES);
  }
  return sumLight;
}
#elif SP_SAMPLING_TYPE & 8

Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  double cosine = 0;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[SP_SAMPLES] = {0};

  // returns 1/pdf.
  double pdf = SphericalSampler::getCosineSurfaceSamples(n, x, SP_SAMPLES, samples);
  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < SP_SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = sList[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, wo_ref, rr.d, x);

      sum += (brdf);
    }
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/SP_SAMPLES);
  }
  return sumLight;
}

#elif SP_SAMPLING_TYPE & 16
Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  double cosine = 0;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  Vec samples[SP_SAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  // returns 1/pdf.
  double pdf = SphericalSampler::getPhongBRDFSamples(n, wo_ref, x, primitive->m.phongExp, SP_SAMPLES, samples);

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < SP_SAMPLES; i++) {
      double d;

      rr.d = (samples[i] - x);

      d = sList[j].p.intersect(rr);

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
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/SP_SAMPLES);
  }
  return sumLight;
}
#else
Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[100] = {0};
  double weight[100] = {0};

  uint32_t nSamples = lightSampler->getSamples(x, n, r.d, primitive->m.phongExp, samples, weight);

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    for (int i = 0; i < nSamples; i++) {
      double d;

      rr.d = (samples[i] - x).norm();

      d = sList[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	continue;
      }

      brdf = primitive->brdf(n, wo_ref, rr.d, x);

      sum += (cosine * brdf * weight[i]);
    }
    sumLight = sumLight + sList[j].radiance * sum;
  }
  return sumLight;
}
#endif
/*
double EnvSource::calcArea() {
  double d_theta = PI / height;
  double d_phi = 2 * PI / width;

  double area = 0;
  for (uint32_t i = 0; i < height; i++) {
    double sum = 0;
    for (uint32_t j = 0; j < width; j++) {
      Vec L = envMap[i * width + j].mult(radiance);
      sum += (L.x + L.y + L.z);
    }
    area +=  sum * sin(i * d_theta);
  }
  return area * d_theta * d_phi;
}

void EnvSource::calcMarginalDensity(double *pMarginal, double area) {
  double d_theta = PI / height;
  double d_phi = 2 * PI / width;

  for (uint32_t i = 0; i < height; i++) {
    double sum = 0;
    for (uint32_t j = 0; j < width; j++) {
      Vec L = envMap[i * width + j].mult(radiance);
      sum += (L.x + L.y + L.z);
    }
    pMarginal[i] = (sum * sin(i * d_theta) * d_phi) / area;
  }
}

void EnvSource::calcConditionalDensity(double *pConditional, double *pMarginal, double area) {
  double d_theta = PI / height;
  for (uint32_t i = 0; i < height; i++) {
     double pM = pMarginal[i];
     for (uint32_t j = 0; j < width; j++) {
       Vec L = envMap[i * width + j].mult(radiance);
       pConditional[i * width + j] = (L.x + L.y + L.z) * sin(i * d_theta) / (area * pM);
     }
  }
}*/

// Stores only the index, not the angles. Multiply by delta to convert index to angles.
void invertedCDF(double *pdf, uint32_t length, double delta, uint32_t *invCDF, uint32_t nBins) {
  double sum = 0;
  for (uint32_t j = 0; j < nBins; j++)
    invCDF[j] = 0;

  for (uint32_t i = 0; i < length; i++) {
    sum += pdf[i] * delta;
    //std::cout<<i<<" "<<sum<< " "<<pdf[i]<<"\n";
    uint32_t k = sum * nBins;
    invCDF[((k>=nBins)?nBins-1:k)] = i;
    for(uint32_t j = k - 1; j >= 0 && j < nBins && invCDF[j] == 0; j--)
      invCDF[j] = i;
  }

  for (uint32_t j = sum * nBins; j < nBins; j++)
    invCDF[j] = (length - 1);
}

double EnvSource::calcArea() {
  uint32_t resolutionTheta = 2000;
  uint32_t resolutionPhi = 4000;
  double d_theta = PI / resolutionTheta;
  double d_phi = 2.0 * PI / resolutionPhi;

  double area = 0;
  for (uint32_t i = 0; i < resolutionTheta; i++) {
    double sintheta = sin(d_theta * i);
    double costheta = cos(d_theta * i);
    double sum = 0;
    for (uint32_t j = 0; j < resolutionPhi; j++) {
      double cosphi = cos(d_phi * j);
      double sinphi = sin(d_phi * j);
      Vec dir = Vec(sintheta * cosphi, sintheta * sinphi, costheta);
      Vec L = getRadiance(dir);
      sum += to_greyScale(L);
    }
    area +=  sum * sintheta;
  }
  return area * d_theta * d_phi;
}

void EnvSource::calcMarginalDensity(double *pMarginal, uint32_t marginalSamples, double area) {
  uint32_t resolutionPhi = 4000;
  double d_theta = PI / marginalSamples;
  double d_phi = 2 * PI / resolutionPhi;

  for (uint32_t i = 0; i < marginalSamples; i++) {
    double sintheta = sin(d_theta * i);
    double costheta = cos(d_theta * i);
    double sum = 0;
    for (uint32_t j = 0; j < resolutionPhi; j++) {
      double cosphi = cos(d_phi * j);
      double sinphi = sin(d_phi * j);
      Vec dir = Vec(sintheta * cosphi, sintheta * sinphi, costheta);
      Vec L = getRadiance(dir);
      sum += to_greyScale(L);
    }
    pMarginal[i] = (sum * sintheta * d_phi) / area;
  }
}

void EnvSource::calcConditionalDensity(double *pConditional, uint32_t cSamples, double *pMarginal, uint32_t mSamples, double area) {
  double d_theta = PI / mSamples;
  double d_phi = 2.0 * PI / cSamples;

  for (uint32_t i = 0; i < mSamples; i++) {
     double pM = pMarginal[i];
     double sintheta = sin(d_theta * i);
     double costheta = cos(d_theta * i);
     for (uint32_t j = 0; j < cSamples; j++) {
      double cosphi = cos(d_phi * j);
      double sinphi = sin(d_phi * j);
      Vec dir = Vec(sintheta * cosphi, sintheta * sinphi, costheta);
      Vec L = getRadiance(dir);
       pConditional[i * cSamples + j] = to_greyScale(L) * sintheta / (area * pM);

    }
  }
}

static double to_greyScale(Vec L) {
  return (L.x * 0.3 + L.y * 0.59 + L.z * 0.11);
}

void EnvSource::initSamples() {
  /*double normalizationArea = calcArea();
  std::cout<<normalizationArea<<"\n";

  uint32_t pdfSamples = 3000;

  double *pMarginal = new double[height];
  calcMarginalDensity(pMarginal, normalizationArea);

  double *pConditional = new double[height * width];
  calcConditionalDensity(pConditional, pMarginal, normalizationArea);

  uint32_t *invMarginalCDF = new uint32_t[pdfSamples];
  invertedCDF(pMarginal, height, PI/height, invMarginalCDF, pdfSamples);

  uint32_t *invConditionalCDF = new uint32_t[height * pdfSamples];
  for (uint32_t i = 0; i < height; i++)
    invertedCDF(&pConditional[i * width], width, 2 * PI/width, &invConditionalCDF[i * pdfSamples], pdfSamples);

  seedMT(clock() & 0xFFFFFFFF);

  for (uint32_t j = 0; j < nSamples * multiplier; j++) {
    double e1 = randomMTD(0, 1);
    double e2 = randomMTD(0, 1);

    uint32_t thetaIdx = invMarginalCDF[(uint32_t)(e1 * pdfSamples)];
    double theta =  thetaIdx * (PI / height);
    uint32_t phiIdx = invConditionalCDF[thetaIdx * pdfSamples + (uint32_t)(e2 * pdfSamples)];
    double phi = phiIdx *  (2.0 * PI) / width;

    double x = sin(theta);
    double y = x;
    x *= cos(phi);
    y *= sin(phi);
    double z = cos(theta);

    if (axis == 0) {
      impSamples[j] = Vec(z, x, y);
    }
    else if (axis == 1) {
      impSamples[j] = Vec(x, z, y);
      //std::cout<<impSamples[j]<<" \n";
    }
    else {
      impSamples[j] = Vec(x, y, z);
    }
  }*/

  normalizationArea = calcArea();

  uint32_t marginalSamples = 3000;

  double *pMarginal = new double[marginalSamples];
  calcMarginalDensity(pMarginal, marginalSamples, normalizationArea);

  uint32_t conditionalSamples = 4000;
  double *pConditional = new double[marginalSamples * conditionalSamples];
  calcConditionalDensity(pConditional, conditionalSamples, pMarginal, marginalSamples, normalizationArea);

  uint32_t pdfSamples = 2000;
  uint32_t *invMarginalCDF = new uint32_t[pdfSamples];
  invertedCDF(pMarginal, marginalSamples, PI/marginalSamples, invMarginalCDF, pdfSamples);

  uint32_t *invConditionalCDF = new uint32_t[marginalSamples * pdfSamples];
  for (uint32_t i = 0; i < marginalSamples; i++)
    invertedCDF(&pConditional[i * conditionalSamples], conditionalSamples, 2 * PI/conditionalSamples, &invConditionalCDF[i * pdfSamples], pdfSamples);

  seedMT(clock() & 0xFFFFFFFF);

  for (uint32_t j = 0; j < nSamples * multiplier; j++) {
    double e1 = randomMTD(0, 1);
    double e2 = randomMTD(0, 1);

    uint32_t thetaIdx = invMarginalCDF[(uint32_t)(e1 * pdfSamples)];
    double theta =  thetaIdx * (PI / height);
    uint32_t phiIdx = invConditionalCDF[thetaIdx * pdfSamples + (uint32_t)(e2 * pdfSamples)];
    double phi = phiIdx *  (2.0 * PI) / width;

    double x = sin(theta);
    double y = x;
    x *= cos(phi);
    y *= sin(phi);
    double z = cos(theta);

    /*if (axis == 0) {
      impSamples[j] = Vec(z, x, y);
    }
    else if (axis == 1) {
      impSamples[j] = Vec(x, z, y);
      //std::cout<<impSamples[j]<<" \n";
    }
    else {*/
      impSamples[j] = Vec(x, y, z);
    //}
  }

}
