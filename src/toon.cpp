//All experimental stuff!!
#include "lightSources.h"
#include "shader.h"
#include <cmath>

#define MAX_MESH_SAMPLES 1000 // Should be less than the No. of samples in domainSampler.h

// We set confidence to zero when the point x is in penubra region of a shadow. We compute shading at this point due to all light sources,
// irrespective of visibility. The excess light is taken care of in the control variate MC estimator.
Vec LightSource::getLightAnalyticDiffuseAdv(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, uint32_t &confidence) {
   if (mList.size() < 1) return Vec();
    confidence = 1;

   double isTotalShadow = 0; // zero means total shadow.
   for (uint32_t j = 0; j < mList.size(); j++) {
    for (uint32_t i = 0; i < mList[j].mesh.size(); i++) {
      Vec e1 = mList[j].mesh[i].t.A - x;
      Vec e2 = mList[j].mesh[i].t.B - x;
      Vec e3 = mList[j].mesh[i].t.C - x;

      double d1 = e1.length();
      double d2 = e2.length();
      double d3 = e3.length();

      e1.norm();
      e2.norm();
      e3.norm();

      Ray r1(x, e1);
      Ray r2(x, e2);
      Ray r3(x, e3);

      double s = shadow(r1, d1) + shadow(r2, d2) + shadow(r3, d3);
      if (s > 0.01 && s < 2.99)
	confidence = 0;
      isTotalShadow += s;
    }
  }

  if (isTotalShadow < 0.01) return Vec();

   Vec sumLight = Vec();
   for (uint32_t j = 0; j < mList.size(); j++) {

    for (uint32_t i = 0; i < mList[j].mesh.size(); i++) {
      Vec e1 = mList[j].mesh[i].t.A - x;
      Vec e2 = mList[j].mesh[i].t.B - x;
      Vec e3 = mList[j].mesh[i].t.C - x;

      e1.norm();
      e2.norm();
      e3.norm();

      double nA = e1.dot(n);
      double nB = e2.dot(n);
      double nC = e3.dot(n);

      double sum = 0;
      if (nA > 0 && nB > 0 && nC > 0) {
	sum = integrateEdge(e1, e2, n);
	sum += integrateEdge(e2, e3, n);
	sum += integrateEdge(e3, e1, n);
      }
      else if (nA < 0 && nB < 0 && nC < 0)
	continue;
      else {
	Vec clip[4];
	if (!clipTriangle(mList[j].mesh[i].t.A, mList[j].mesh[i].t.B, mList[j].mesh[i].t.C, n, x, nA, nB, nC, clip))
	  continue;

	e1 = (clip[0] - x).norm();
	e2 = (clip[1] - x).norm();
	e3 = (clip[2] - x).norm();
	Vec e4 = (clip[3] - x).norm();
	sum = integrateEdge(e1, e2, n);
	sum += integrateEdge(e2, e3, n);
	sum += integrateEdge(e3, e4, n);
	sum += integrateEdge(e4, e1, n);

	// Do not set confidence to 1 in here as there can be two visibility events 1. Partially visibile light source. 2. An object casting a penubra at the same location due to same or another light source.
	// Had it only been a scene case 1, confidence could have been set to 1.
	//confidence = 1;
      }
      sum = sum > 0 ? sum : sum * -1;

      sumLight = sumLight + mList[j].mesh[i].radiance * sum;
    }
  }

  return sumLight / (2.0 * PI);
}

Vec LightSource::getLightFromMeshSource_CVAdv(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (mList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_MESH_SAMPLES ? MAX_MESH_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[sampleCount] = {0};
  uint32_t ids[sampleCount] = {0};
  uint32_t cvConfidence = 1;

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();
  Vec G = getLightAnalyticDiffuseAdv(r, n, x, primitive, cvConfidence);

  if (cvConfidence > 0) return G;

  for (uint32_t j = 0; j < mList.size(); j++) {
    Vec sum;

    // returns 1/pdf.
    double pdf = mList[j].getSamples(sampleCount, samples, ids);
    for (uint32_t i = 0; i < sampleCount; i++) {
      double visibility = 1;
      double d;

      rr.d = (samples[i] - x).norm();

      d = (samples[i] - x).length();

      if (shadow(rr, d) < 0.5)
	visibility = 0;

      cosine = rr.d.dot(n);

      double cosine1 = rr.d.dot(mList[j].mesh[ids[i]].t.n);
      cosine1 = cosine1 < 0 ? -cosine1 : cosine1;

      if (cosine < 0) {
	 // This code colors the dark side of an object black.
	  continue;
      }

      brdf = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d);

      sum = sum +  mList[j].mesh[ids[i]].radiance * ((visibility - 1.0) * cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + (sum * (pdf/sampleCount));

    sumLight = sumLight + G;
  }
  return sumLight;
}