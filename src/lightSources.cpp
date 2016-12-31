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

#define MAX_SP_SAMPLES 		(1000) // Should be less than the No. of samples in domainSampler.h
#define SP_SAMPLING_TYPE	1// 1 Solid Angle Importance Sampling, 2 Light Surface Sampling

#define MAX_ENV_SAMPLES		1000 // Should be less than the No. of samples in domainSampler.h and nSamples in struct EnvSource
#define ENV_SAMPLING_TYPE	1// 1 for BRDF sampling, 2 Importance sampling.

#define MAX_TRI_SAMPLES 	1000 // Should be less than the No. of samples in domainSampler.h

#define MAX_MESH_SAMPLES	1000 // Should be less than the No. of samples in domainSampler.h


#define MAX_BRDF_SAMPLES	1000 // Should be less than the No. of samples in domainSampler.h
#define MAX_UNIFORM_SAMPLES	1000 // Should be less than the No. of samples in domainSampler.h
#define MAX_COSINE_SAMPLES	1000 // Should be less than the No. of samples in domainSampler.h

LightSource *lSource;

void configureLightSources() {
  lSource = new LightSource();
  //lSource->addPSource(PointSource(Vec(350, 100, -400), Vec(1, 1, 1), 40000));
  //lSource->addSSource(SphereSource(Vec(10.0, 0.0, 0.0), Sphere(10, Vec(20, 40.6 - .27, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(5.0, 0.0, 0.0), Sphere(10, Vec(50, 40.6 - .27, 81.6), Vec(.0, 0.0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(10, Vec(50, 68, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(30, Vec(200, 100, -400), 2000, 10));
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(10, Vec(100, 68, 81.6), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //Veach Scene
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(20, Vec(150, 68.6 - .27, 0), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(10, Vec(50, 68.6 - .27, 0), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(5, Vec(-50, 68.6 - .27, 0), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  //lSource->addSSource(SphereSource(Vec(10.0, 10.0, 10.0), Sphere(1, Vec(-150, 68.6 - .27, 0), Vec(.0, .0, .0), 1.0, MaterialType(VOLUME, Vec(0., 0., 0.)))));
  /*uint32_t w, h;
  Vec *image;
  readImage(image, w, h);
  lSource->addESource(EnvSource(image, w, h, 1, Vec(1.0, 1.0, 1.0), 2));
  delete []image;*/
  /*lSource->addTSource(TriLight(
    Triangle(Vec(30, 68, 60), Vec(70, 68, 100), Vec(30, 68, 100), Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 10)));
  lSource->addTSource(TriLight(
    Triangle(Vec(30, 68, 60), Vec(70, 68, 100), Vec(70, 80, 60), Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 10)));
*/
  
  double size = 20;
  Vec V1 = Vec(size, size, 0);
  Vec V2 = Vec(-size, -size, 0);
  Vec V3 = Vec(-size, size, 0);
  Vec V4 = Vec(size, -size, 0);
  Vec T = Vec(-600, 180, -1600);
  MeshLight mesh1;
  mesh1.add(TriLight(V1 + T, V3 + T, V2 + T, 40000, 640));
  mesh1.add(TriLight(V1 + T, V2 + T, V4 + T, 40000, 640));
  mesh1.initMeshLight();
  lSource->addMSource(mesh1);
  
  /*size = 14;
  V1 = Vec(size, size, 0);
  V2 = Vec(-size, -size, 0);
  V3 = Vec(-size, size, 0);
  V4 = Vec(size, -size, 0);
  T = Vec(-290, 250, -800);
  MeshLight mesh2;
  mesh2.add(TriLight(V1 + T, V3 + T, V2 + T, 20000, 320));
  mesh2.add(TriLight(V1 + T, V2 + T, V4 + T, 20000, 320));
  mesh2.initMeshLight();
  lSource->addMSource(mesh2);*/
  
  size = 40;
  V1 = Vec(size, size, 0);
  V2 = Vec(-size, -size, 0);
  V3 = Vec(-size, size, 0);
  V4 = Vec(size, -size, 0);
  T = Vec(-306, 180, -1600);
  MeshLight mesh3;
  mesh3.add(TriLight(V1 + T, V3 + T, V2 + T, 10000, 160));
  mesh3.add(TriLight(V1 + T, V2 + T, V4 + T, 10000, 160));
  mesh3.initMeshLight();
  lSource->addMSource(mesh3);
  
  /*size = 29;
  V1 = Vec(size, size, 0);
  V2 = Vec(-size, -size, 0);
  V3 = Vec(-size, size, 0);
  V4 = Vec(size, -size, 0);
  T = Vec(-56, 250, -800);
  MeshLight mesh4;
  mesh4.add(TriLight(V1 + T, V3 + T, V2 + T, 5000, 80));
  mesh4.add(TriLight(V1 + T, V2 + T, V4 + T, 5000, 80));
  mesh4.initMeshLight();
  lSource->addMSource(mesh4);*/
  
  size = 80;
  V1 = Vec(size, size, 0);
  V2 = Vec(-size, -size, 0);
  V3 = Vec(-size, size, 0);
  V4 = Vec(size, -size, 0);
  T = Vec(20, 180, -1600);
  MeshLight mesh5;
  mesh5.add(TriLight(V1 + T, V3 + T, V2 + T, 2500, 40));
  mesh5.add(TriLight(V1 + T, V2 + T, V4 + T, 2500, 40));
  mesh5.initMeshLight();
  lSource->addMSource(mesh5);
  
  /*size = 57;
  V1 = Vec(size, size, 0);
  V2 = Vec(-size, -size, 0);
  V3 = Vec(-size, size, 0);
  V4 = Vec(size, -size, 0);
  T = Vec(230, 250, -800);
  MeshLight mesh6;
  mesh6.add(TriLight(V1 + T, V3 + T, V2 + T, 1250, 20));
  mesh6.add(TriLight(V1 + T, V2 + T, V4 + T, 1250, 20));
  mesh6.initMeshLight();
  lSource->addMSource(mesh6);*/
  
  size = 160;
  V1 = Vec(size, size, 0);
  V2 = Vec(-size, -size, 0);
  V3 = Vec(-size, size, 0);
  V4 = Vec(size, -size, 0);
  T = Vec(430, 180, -1600);
  MeshLight mesh7;
  mesh7.add(TriLight(V1 + T, V3 + T, V2 + T, 1000, 10));
  mesh7.add(TriLight(V1 + T, V2 + T, V4 + T, 1000, 10));
  mesh7.initMeshLight();
  lSource->addMSource(mesh7);
  
  
/*
 MeshLight mesh1;
  mesh1.add(TriLight(
    Triangle(Vec(30, 68, 60), Vec(30, 68, 100), Vec(70, 68, 100), Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 0)));
  mesh1.add(TriLight(
    Triangle(Vec(30, 68, 60), Vec(70, 68, 100), Vec(70, 68, 60), Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 0, 10)));

  mesh1.initMeshLight();
  lSource->addMSource(mesh1);*/
/*
  MeshLight mesh2;
  Vec A = Vec(30, 68, 60);
  Vec B = Vec(70, 68, 60);
  Vec C = Vec(50, 68, 100);
  Vec D = Vec(50, 40, 80);
  mesh2.add(TriLight(
    Triangle(A, B, C, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 0)));
  mesh2.add(TriLight(
    Triangle(B, D, C, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(0, 10, 0)));
  mesh2.add(TriLight(
    Triangle(A, D, B, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 10)));
  mesh2.add(TriLight(
    Triangle(A, D, C, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(0, 0, 10)));
  mesh2.initMeshLight();
  lSource->addMSource(mesh2);*/
/*
  lSource->addTSource(TriLight(
    Triangle(A, B, C, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 0)));
  lSource->addTSource(TriLight(
    Triangle(B, D, C, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(0, 10, 0)));
  lSource->addTSource(TriLight(
    Triangle(A, D, B, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(10, 10, 10)));
  lSource->addTSource(TriLight(
    Triangle(A, D, C, Vec(0, 0, 0), 1.0, MaterialType(PLANAR, Vec(0., 0., 0.))), Vec(0, 0, 10)));
*/
//#undef dX
//#define dX (-100)
  //lSource->addTSource(TriLight(Vec(30 + dX, 200, -800), Vec(-130 + dX, 200, -800), Vec(-130 + dX, 40, -800), 10000, 10));
 // lSource->addTSource(TriLight(Vec(30 + dX, 200, -800), Vec(-130 + dX, 40, -800), Vec(30 + dX, 40, -800), 10000, 10));
}

Vec LightSource::getLightFromAllSources_MIS(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  if (mList.size() < 1 && sList.size() < 1 && tList.size() < 1 && eS.size() != 1 && pList.size() < 1) return Vec();
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double throughput = 0;
  Vec samples[200] = {0};
  double weights[200] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t i = 0; i < pList.size(); i++) {
    rr.d = (pList[i].p - x).norm(); // direction from point towards source

    cosine = rr.d.dot(n);
    if (cosine < 0) {
       // This code colors the dark side of an object black.
	continue;
    }
    if (shadow(rr, (pList[i].p - x).length()) < 0.5)
      continue;

    sumLight = sumLight + (pList[i].getRadiance(x) * (cosine * primitive->brdf(n, wo_ref, r.d * -1.0, rr.d)));
  }

  uint32_t nSamples = lightSampler->getSamples(x, n, wo_ref, r.d * -1.0, primitive, samples, weights);

  for (uint32_t i = 0; i < nSamples; i++) {
    rr.d = samples[i];
    cosine = rr.d.dot(n);
    if (cosine < 0)
	continue;
    double d;
    bool hitASource = false;
    throughput = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d) * cosine * weights[i];

    for (uint32_t j = 0; j < mList.size() && !hitASource; j++) {
      uint32_t id;
      Vec dummy;
      if (!mList[j].intersect(rr, id, d, dummy))
	continue;
      if (shadow(rr, d) < 0.5)
	continue;
      hitASource = true;
      sumLight = sumLight +  mList[j].mesh[id].radiance * throughput;
    }
    for (uint32_t j = 0; j < sList.size() && !hitASource; j++) {
      d = sList[j].p.intersect(rr);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  sList[j].radiance * throughput;
    }
    for (uint32_t j = 0; j < tList.size() && !hitASource; j++) {
      Vec dummy;
      d = tList[j].t.intersect(rr, dummy);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  tList[j].radiance * throughput;
    }
    if (!hitASource && eS.size() == 1)
      if (shadow(rr, INF) > 0.5)
	sumLight = sumLight + eS[0].getRadiance(rr.d) * throughput;
  }

  return sumLight;
}

Vec LightSource::getLightFromAllSources_BRDF(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (mList.size() < 1 && sList.size() < 1 && tList.size() < 1 && eS.size() != 1 && pList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_BRDF_SAMPLES ? MAX_BRDF_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double throughput = 0;
  Vec samples[sampleCount];
  double weights[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t i = 0; i < pList.size(); i++) {
    rr.d = (pList[i].p - x).norm(); // direction from point towards source

    cosine = rr.d.dot(n);
    if (cosine < 0) {
       // This code colors the dark side of an object black.
	continue;
    }
    if (shadow(rr, (pList[i].p - x).length()) < 0.5)
      continue;

    sumLight = sumLight + (pList[i].getRadiance(x) * (cosine * primitive->brdf(n, wo_ref, r.d * -1.0, rr.d)));
  }

  // Sample BRDF direction
  primitive->m.getBrdfDirectionSamples(n, wo_ref, r.d * -1.0, samples, weights, sampleCount);

  for (uint32_t i = 0; i < nSamples; i++) {
    rr.d = samples[i];
    cosine = rr.d.dot(n);
    if (cosine < 0)
	continue;
    double d;
    bool hitASource = false;
    throughput = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d) * cosine * weights[i] / sampleCount;

    for (uint32_t j = 0; j < mList.size() && !hitASource; j++) {
      uint32_t id;
      Vec dummy;
      if (!mList[j].intersect(rr, id, d, dummy))
	continue;
      if (shadow(rr, d) < 0.5)
	continue;
      hitASource = true;
      sumLight = sumLight +  mList[j].mesh[id].radiance * throughput;
    }
    for (uint32_t j = 0; j < sList.size() && !hitASource; j++) {
      d = sList[j].p.intersect(rr);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  sList[j].radiance * throughput;
    }
    for (uint32_t j = 0; j < tList.size() && !hitASource; j++) {
      Vec dummy;
      d = tList[j].t.intersect(rr, dummy);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  tList[j].radiance * throughput;
    }
    if (!hitASource && eS.size() == 1)
      if (shadow(rr, INF) > 0.5)
	sumLight = sumLight + eS[0].getRadiance(rr.d) * throughput;
  }

  return sumLight;
}

Vec LightSource::getLightFromAllSources_UNIFORM(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (mList.size() < 1 && sList.size() < 1 && tList.size() < 1 && eS.size() != 1 && pList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_UNIFORM_SAMPLES ? MAX_UNIFORM_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  Vec samples[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t i = 0; i < pList.size(); i++) {
    rr.d = (pList[i].p - x).norm(); // direction from point towards source

    cosine = rr.d.dot(n);
    if (cosine < 0) {
       // This code colors the dark side of an object black.
	continue;
    }
    if (shadow(rr, (pList[i].p - x).length()) < 0.5)
      continue;

    sumLight = sumLight + (pList[i].getRadiance(x) * (cosine * primitive->brdf(n, wo_ref, r.d * -1.0, rr.d)));
  }

  double pdfInv = SphericalSampler::getHemiDirectionSamples(n, sampleCount, samples) / sampleCount;

  for (uint32_t i = 0; i < nSamples; i++) {
    rr.d = samples[i];
    cosine = rr.d.dot(n);
    if (cosine < 0)
	continue;
    double d;
    bool hitASource = false;
    double throughput = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d) * cosine * pdfInv;

    for (uint32_t j = 0; j < mList.size() && !hitASource; j++) {
      uint32_t id;
      Vec dummy;
      if (!mList[j].intersect(rr, id, d, dummy))
	continue;
      if (shadow(rr, d) < 0.5)
	continue;
      hitASource = true;
      sumLight = sumLight +  mList[j].mesh[id].radiance * throughput;
    }
    for (uint32_t j = 0; j < sList.size() && !hitASource; j++) {
      d = sList[j].p.intersect(rr);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  sList[j].radiance * throughput;
    }
    for (uint32_t j = 0; j < tList.size() && !hitASource; j++) {
      Vec dummy;
      d = tList[j].t.intersect(rr, dummy);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  tList[j].radiance * throughput;
    }
    if (!hitASource && eS.size() == 1)
      if (shadow(rr, INF) > 0.5)
	sumLight = sumLight + eS[0].getRadiance(rr.d) * throughput;
  }

  return sumLight;
}

Vec LightSource::getLightFromAllSources_COS(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (mList.size() < 1 && sList.size() < 1 && tList.size() < 1 && eS.size() != 1 && pList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_COSINE_SAMPLES ? MAX_COSINE_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  Vec samples[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t i = 0; i < pList.size(); i++) {
    rr.d = (pList[i].p - x).norm(); // direction from point towards source

    cosine = rr.d.dot(n);
    if (cosine < 0) {
       // This code colors the dark side of an object black.
	continue;
    }
    if (shadow(rr, (pList[i].p - x).length()) < 0.5)
      continue;

    sumLight = sumLight + (pList[i].getRadiance(x) * (cosine * primitive->brdf(n, wo_ref, r.d * -1.0, rr.d)));
  }

  double pdfInv = SphericalSampler::getCosineDirectionSamples(n, sampleCount, samples) / sampleCount;

  for (uint32_t i = 0; i < nSamples; i++) {
    rr.d = samples[i];
    cosine = rr.d.dot(n);
    if (cosine < 0)
	continue;
    double d;
    bool hitASource = false;
    double throughput = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d) * pdfInv;

    for (uint32_t j = 0; j < mList.size() && !hitASource; j++) {
      uint32_t id;
      Vec dummy;
      if (!mList[j].intersect(rr, id, d, dummy))
	continue;
      if (shadow(rr, d) < 0.5)
	continue;
      hitASource = true;
      sumLight = sumLight +  mList[j].mesh[id].radiance * throughput;
    }
    for (uint32_t j = 0; j < sList.size() && !hitASource; j++) {
      d = sList[j].p.intersect(rr);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  sList[j].radiance * throughput;
    }
    for (uint32_t j = 0; j < tList.size() && !hitASource; j++) {
      Vec dummy;
      d = tList[j].t.intersect(rr, dummy);
      if (d >= INF)
	continue;
      if (shadow(rr, d) < 0.5)
	continue;

      hitASource = true;
      sumLight = sumLight +  tList[j].radiance * throughput;
    }
    if (!hitASource && eS.size() == 1)
      if (shadow(rr, INF) > 0.5)
	sumLight = sumLight + eS[0].getRadiance(rr.d) * throughput;
  }

  return sumLight;
}

Vec LightSource::getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  if (pList.size() < 1) return Vec();
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

    brdf = primitive->brdf(n, wo_ref, r.d * -1.0, refd);

    sumLight = sumLight + (pList[i].getRadiance(x) * (cosine * brdf));
  }

  return sumLight;
}

Vec LightSource::getLightFromTriSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (tList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_TRI_SAMPLES ? MAX_TRI_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < tList.size(); j++) {
    double sum = 0;
    SphericalSampler::getTriLightSurfaceSamples(tList[j].t.A, tList[j].t.B, tList[j].t.C, sampleCount, samples);
    // returns 1/pdf.
    double pdf = tList[j].t.area;
    for (uint32_t i = 0; i < sampleCount; i++) {
      double d;

      rr.d = (samples[i] - x).norm();

      d = (samples[i] - x).length();

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot(tList[j].t.n);
      cosine1 = cosine1 < 0 ? -cosine1 : cosine1;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d);

      sum += (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + tList[j].radiance * ( pdf * sum/sampleCount);
  }
  return sumLight;
}

Vec LightSource::getLightFromMeshSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (mList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_MESH_SAMPLES ? MAX_MESH_SAMPLES : nSamples;

  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[sampleCount];
  uint32_t ids[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  BRDFApprox brdfApprox;
  primitive->m.getBrdfApproxParam(n, r.d * -1.0, brdfApprox);

  for (uint32_t j = 0; j < mList.size(); j++) {
    Vec sum;
    // returns 1/pdf.
    double pdf = mList[j].getSamples(sampleCount, samples, ids);
    for (uint32_t i = 0; i < sampleCount; i++) {
      double d;

      rr.d = (samples[i] - x).norm();

      d = (samples[i] - x).length();

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot(mList[j].mesh[ids[i]].t.n);
      cosine1 = cosine1 < 0 ? -cosine1 : cosine1;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d, brdfApprox);

      sum = sum +  mList[j].mesh[ids[i]].radiance * (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + sum * (pdf/sampleCount);
  }
  return sumLight;
}

#if ENV_SAMPLING_TYPE & 1
Vec LightSource::getLightFromEnvSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (eS.size() != 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_ENV_SAMPLES ? MAX_ENV_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double eps = 1e-08;
  Vec sumLight;
  Vec samples[sampleCount];
  double weights[sampleCount];

  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  // Sample BRDF direction
  primitive->m.getBrdfDirectionSamples(n, wo_ref, r.d * -1.0, samples, weights, sampleCount);

  for (uint32_t i = 0; i < nSamples; i++) {
    rr.d = samples[i];
    double cosine = rr.d.dot(n);
    if (cosine < 0)
	continue;

    if (shadow(rr, INF) > 0.5)
      sumLight = sumLight + eS[0].getRadiance(rr.d) * primitive->brdf(n, wo_ref, r.d * -1.0, rr.d) * cosine * weights[i];
  }

  return sumLight / sampleCount;
}
#else
Vec LightSource::getLightFromEnvSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (eS.size() != 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_ENV_SAMPLES ? MAX_ENV_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-15;
  double brdf = 0;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();
  Random random(clock() & 0xFFFFFFFF);

  Vec sumLight = Vec();
  rr.o = x;
  uint32_t offset = random.randomMTD(0, (eS[0].multiplier -1) * sampleCount);
  for (uint32_t i = 0; i < sampleCount; i++) {
    rr.d = eS[0].impSamples[i + offset];

    if (shadow(rr, INF) < 0.5)
      continue;

    cosine = rr.d.dot(n);

    if (cosine < 0) {
	continue;
    }

    brdf = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d);
    Vec L = eS[0].getRadiance(rr.d);
    double pdf = to_greyScale(L);

    if(pdf > eps)
      sumLight = sumLight + (L * cosine * brdf)/pdf;
  }
  //std::cout<<sumLight.x<<" "<<sumLight.y<<" "<<sumLight.z<<" \n";
  sumLight = sumLight * (eS[0].normalizationArea/(sampleCount));

  return sumLight;
}
#endif

#if SP_SAMPLING_TYPE & 1
Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (sList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_SP_SAMPLES ? MAX_SP_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    Vec w = sList[j].p.p - x;
    double sinThetaMax = sList[j].p.r/w.length();
    // returns 1/pdf.
    double pdf = SphericalSampler::getSolidDirectionSamples(w.norm(), asin(sinThetaMax), sampleCount, samples);
    //samples[SP_SAMPLES - 1] = x + w.norm();
    for (uint32_t i = 0; i < sampleCount; i++) {
      double d;

      rr.d = samples[i];

      d = sList[j].p.intersect(rr);

      if (d >= INF)
	continue;

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      if (cosine < 0) {
	  continue; // Color dark side of an object black.
      }

      brdf = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d);

      sum += (cosine * brdf);
    }
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/sampleCount);
  }
  return sumLight;
}

#elif SP_SAMPLING_TYPE & 2
Vec LightSource::getLightFromSphereSources(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (sList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_SP_SAMPLES ? MAX_SP_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  for (uint32_t j = 0; j < sList.size(); j++) {
    double sum = 0;
    // returns 1/pdf.
    double pdf = SphericalSampler::getLightDirectionSamples(sList[j].p.p, sList[j].p.r, x, sampleCount, samples);

    for (uint32_t i = 0; i < sampleCount; i++) {
      double d;

      rr.d = samples[i];

      d = sList[j].p.intersect(rr);

      if (shadow(rr, d) < 0.5)
	continue;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot((x + rr.d * d - sList[j].p.p).norm());
      cosine1 = cosine1 < 0 ? -cosine1 : 0;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d);

      sum += (cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + sList[j].radiance * ( pdf * sum/sampleCount);
  }
  return sumLight;
}
#endif

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

void EnvSource::initSamples() {
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

  Random random(clock() & 0xFFFFFFFF);

  for (uint32_t j = 0; j < nSamples * multiplier; j++) {
    double e1 = random.randomMTD(0, 1);
    double e2 = random.randomMTD(0, 1);

    uint32_t thetaIdx = invMarginalCDF[(uint32_t)(e1 * pdfSamples)];
    double theta =  thetaIdx * (PI / marginalSamples);
    uint32_t phiIdx = invConditionalCDF[thetaIdx * pdfSamples + (uint32_t)(e2 * pdfSamples)];
    double phi = phiIdx *  (2.0 * PI) / conditionalSamples;

    double x = sin(theta);
    double y = x;
    x *= cos(phi);
    y *= sin(phi);
    double z = cos(theta);

    impSamples[j] = Vec(x, y, z);
  }

}

void MeshLight::initMeshLight() {
    area = 0;
    for (uint32_t i = 0; i < mesh.size(); i++)
      area += mesh[i].t.area;

    double cdf = 0;
    uint32_t lastIndex = 0;
    for (uint32_t i = 0; i < mesh.size(); i++) {
      cdf += mesh[i].t.area / area;

      uint32_t index = cdf * invCDFBins;

      if (index >= invCDFBins)
	index = invCDFBins - 1;

      for (uint32_t j = lastIndex; j <= index; j++)
        invCDF[j] = i;

      lastIndex = index + 1;
    }

    for (uint32_t i = lastIndex; i < invCDFBins; i++)
      invCDF[i] = mesh.size() - 1;
}

double MeshLight::getSamples(uint32_t nSamples, Vec *store, uint32_t *ids) {
    Random random(clock() & 0xFFFFFFFF);
    for (uint32_t i = 0; i < nSamples; i++) {
      double e = random.randomMTD(0, 1);
      double e1 = sqrt(random.randomMTD(0, 1));
      double e2 = random.randomMTD(0, 1);
      uint32_t index = e * invCDFBins;
      index = index >= invCDFBins ? invCDFBins - 1: index;

      index = invCDF[index];
      //index = e > 0.5 ? 0:1;
      if (ids)
	ids[i] = index;
      store[i] =  mesh[index].t.A * (1.0 - e1) +
		 mesh[index].t.B * e1 * (1.0 - e2) +
		 mesh[index].t.C * e2 * e1;
    }

    return area;
}

bool MeshLight::intersect(const Ray &r, uint32_t &id, double &t, Vec &N) {
    double d;
    Vec N_;
    t = INF;
    id = 0xFFFFFFFF;

    // Find the closest sphere with which the ray intersects.
    for (int i = mesh.size(); i--; ) {
      if((d = mesh[i].t.intersect(r, N_)) && d < t) {
	t = d; // set the distance of intersection.
	N = N_;
	id = i;
      }
    }

    // return true if the intersection distance is finite.
    return t < INF;
}
