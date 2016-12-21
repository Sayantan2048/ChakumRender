#include "lightSources.h"
#include "shader.h"
#include <cmath>

#define MAX_MESH_SAMPLES 1000 // Should be less than the No. of samples in domainSampler.h

double LightSource::integrateEdge(const Vec &v1, const Vec &v2, const Vec &n) {
    double cosTheta = v1.dot(v2);
    double theta = acos(cosTheta);
    double d = (v1%v2).dot(n);
    double res = d * ((theta > 0.001) ? theta/sin(theta) : 1.0);

    return res;
}

// Simplified for n = Vec(0, 0, 1)
double LightSource::integrateEdge(const Vec &v1, const Vec &v2) {
    double cosTheta = v1.dot(v2);
    double theta = acos(cosTheta);
    double d = (v1%v2).z;
    double res = d * ((theta > 0.001) ? theta/sin(theta) : 1.0);

    return res;
}

//https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
bool LightSource::clipTriangle(const Vec &A, const Vec &B,  const Vec &C, const Vec &n, const Vec &P, const double nA, const double nB, const double nC, Vec *clip) {
  double eps = 1e-15;

  bool intersectAB = false;
  bool intersectBC = false;
  bool intersectCA = false;

  Vec dirAB = A - B;
  Vec dirBC = B - C;
  Vec dirCA = C - A;

  double nAB = n.dot(dirAB);
  double nBC = n.dot(dirBC);
  double nCA = n.dot(dirCA);

  double pA = n.dot(P - A);
  double pB = n.dot(P - B);
  double pC = n.dot(P - C);

  double d;
  int idx = 0;
  Vec iAB, iBC, iCA;

  if (nAB * nAB > eps) {
      d = pA/nAB;
      iAB = dirAB * d + A;
    /*
     * https://mathspace.co/learn/world-of-maths/coordinate-geometry/division-of-an-interval-in-a-given-ratio-18577/division-of-an-interval-in-a-given-ratio-1065/
     */
      Vec k1 = iAB - A;
      Vec k2 = B - iAB;

      if (k1.dot(k2) > 0) {
	intersectAB = true;
	if (nA < 0) {
	    clip[idx++] = iAB;
	    clip[idx++] = B;
	}
	else {
	    clip[idx++] = A;
	    clip[idx++] = iAB;
	}
      }
  }
  else if (pA * pA < eps) // plane grazes AB
    return false;

  if (nBC * nBC > eps) {
    d = pB/nBC;
    iBC = dirBC * d + B;

    Vec k1 = iBC - B;
    Vec k2 = C - iBC;

    if (k1.dot(k2) > 0) {
      intersectBC = true;
      if (nB < 0) {
	    clip[idx++] = iBC;
	    clip[idx++] = C;
	}
	else {
	    clip[idx++] = B;
	    clip[idx++] = iBC;
	}
    }
  }
  else if (pB * pB < eps ) // plane grazes BC
    return false;

  if (nCA * nCA > eps) {
    d = pC / nCA;
    iCA = dirCA * d + C;

    Vec k1 = iCA.x - C.x;
    Vec k2 = A.x - iCA.x;

    if (k1.dot(k2) > 0) {
      intersectCA = true;
      if (nC < 0) {
	    clip[idx++] = iCA;
	    clip[idx++] = A;
	}
	else {
	    clip[idx++] = C;
	    clip[idx++] = iCA;
	}
    }
  }
  else if (pC * pC < eps) // plane grazes CA
    return false;
  //std::cout<<intersectAB<<" "<<intersectBC<<" "<<intersectCA<<" "<<nA<<" "<<nB<<" "<<nC<<" N:("<<n.x<<","<<n.y<<","<<n.z<<") P:("<<P.x<<","<<P.y<<","<<P.z<<") iAB:(";
  //std::cout<<iAB.x<<","<<","<<iAB.y<<","<<iAB.z<<") iBC:("<<iBC.x<<","<<iBC.y<<","<<iBC.z<<") clip:";
  //clip[0].show();
  //clip[1].show();
  //clip[2].show();
  //clip[3].show();
  //std::cout<<"\n";
  return (intersectAB + intersectBC + intersectCA) > 1;
}

//Simplfied for n = Vec(0,0,1) and P = Vec(0, 0, 0)
bool LightSource::clipTriangle(const Vec &A, const Vec &B,  const Vec &C, const double nA, const double nB, const double nC, Vec *clip) {
  double eps = 1e-15;

  bool intersectAB = false;
  bool intersectBC = false;
  bool intersectCA = false;

  Vec dirAB = A - B;
  Vec dirBC = B - C;
  Vec dirCA = C - A;

  double d;
  int idx = 0;
  Vec iAB, iBC, iCA;

  if (dirAB.z * dirAB.z > eps) {
      d = -A.z/dirAB.z;
      iAB = dirAB * d + A;
    /*
     * https://mathspace.co/learn/world-of-maths/coordinate-geometry/division-of-an-interval-in-a-given-ratio-18577/division-of-an-interval-in-a-given-ratio-1065/
     */
      Vec k1 = iAB - A;
      Vec k2 = B - iAB;

      if (k1.dot(k2) > 0) {
	intersectAB = true;
	if (nA < 0) {
	    clip[idx++] = iAB;
	    clip[idx++] = B;
	}
	else {
	    clip[idx++] = A;
	    clip[idx++] = iAB;
	}
      }
  }
  else if (A.z * A.z < eps) // plane grazes AB
    return false;

  if (dirBC.z * dirBC.z > eps) {
    d = -B.z/dirBC.z;
    iBC = dirBC * d + B;

    Vec k1 = iBC - B;
    Vec k2 = C - iBC;

    if (k1.dot(k2) > 0) {
      intersectBC = true;
      if (nB < 0) {
	    clip[idx++] = iBC;
	    clip[idx++] = C;
	}
	else {
	    clip[idx++] = B;
	    clip[idx++] = iBC;
	}
    }
  }
  else if (B.z * B.z < eps ) // plane grazes BC
    return false;

  if (dirCA.z * dirCA.z > eps) {
    d = -C.z / dirCA.z;
    iCA = dirCA * d + C;

    Vec k1 = iCA.x - C.x;
    Vec k2 = A.x - iCA.x;

    if (k1.dot(k2) > 0) {
      intersectCA = true;
      if (nC < 0) {
	    clip[idx++] = iCA;
	    clip[idx++] = A;
	}
	else {
	    clip[idx++] = C;
	    clip[idx++] = iCA;
	}
    }
  }
  else if (C.z * C.z < eps) // plane grazes CA
    return false;
  //std::cout<<intersectAB<<" "<<intersectBC<<" "<<intersectCA<<" "<<nA<<" "<<nB<<" "<<nC<<" N:("<<n.x<<","<<n.y<<","<<n.z<<") P:("<<P.x<<","<<P.y<<","<<P.z<<") iAB:(";
  //std::cout<<iAB.x<<","<<","<<iAB.y<<","<<iAB.z<<") iBC:("<<iBC.x<<","<<iBC.y<<","<<iBC.z<<") clip:";
  //clip[0].show();
  //clip[1].show();
  //clip[2].show();
  //clip[3].show();
  //std::cout<<"\n";
  return (intersectAB + intersectBC + intersectCA) > 1;
}

Vec LightSource::getLightAnalyticDiffuse(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
   if (mList.size() < 1) return Vec();

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
      }
      sum = sum > 0 ? sum : sum * -1;

      sumLight = sumLight + mList[j].mesh[i].radiance * sum;
    }
  }

  return sumLight / (2.0 * PI);
}

Vec LightSource::getLightAnalytic(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const BRDFApprox &brdfApprox) {
   if (mList.size() < 1) return Vec();
   if (!brdfApprox.isSet)
     std::cerr<<"BRDF Approximations are not set.\n";

   Vec sumLight = Vec();
   for (uint32_t j = 0; j < mList.size(); j++) {

    for (uint32_t i = 0; i < mList[j].mesh.size(); i++) {
      Vec e1 = mList[j].mesh[i].t.A - x;
      Vec e2 = mList[j].mesh[i].t.B - x;
      Vec e3 = mList[j].mesh[i].t.C - x;

      e1 = brdfApprox.Minv.mul(e1); // Light source in canonical original distribution
      e2 = brdfApprox.Minv.mul(e2);
      e3 = brdfApprox.Minv.mul(e3);

      double sum = 0;
      if (e1.z > 0 && e2.z > 0 && e3.z > 0) {
	e1.norm();
	e2.norm();
	e3.norm();
	sum = integrateEdge(e1, e2);
	sum += integrateEdge(e2, e3);
	sum += integrateEdge(e3, e1);
      }
      else if (e1.z < 0 && e2.z < 0 && e3.z < 0)
	continue;
      else {
	Vec clip[4];
	if (!clipTriangle(e1, e2, e3, e1.z/e1.length(), e2.z/e2.length(), e3.z/e3.length(), clip))
	  continue;

	e1 = clip[0].norm();
	e2 = clip[1].norm();
	e3 = clip[2].norm();
	Vec e4 = clip[3].norm();

	sum = integrateEdge(e1, e2);
	sum += integrateEdge(e2, e3);
	sum += integrateEdge(e3, e4);
	sum += integrateEdge(e4, e1);
      }
      sum = sum > 0 ? sum : -sum;

      sumLight = sumLight + mList[j].mesh[i].radiance * sum;
    }
  }

  return sumLight / (2.0 * PI);
}

Vec LightSource::getLightFromMeshSource_Analytic(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  if (mList.size() < 1) return Vec();

  BRDFApprox brdfApprox;
  primitive->m.getBrdfApproxParam(n, r.d * -1.0, brdfApprox, true);

  Vec G = getLightAnalytic(r, n, x, primitive, brdfApprox) * brdfApprox.amplitude * primitive->m.specularCoef;
  G = G + getLightAnalyticDiffuse(r, n, x, primitive) * (1 - primitive->m.specularCoef);

  return G;
}

Vec LightSource::getLightFromMeshSource_AnalyticCV(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive, const uint32_t nSamples) {
  if (mList.size() < 1) return Vec();
  const uint32_t sampleCount = nSamples > MAX_MESH_SAMPLES ? MAX_MESH_SAMPLES : nSamples;
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;

  Vec samples[sampleCount];
  uint32_t ids[sampleCount];

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  BRDFApprox brdfApprox;
  primitive->m.getBrdfApproxParam(n, r.d * -1.0, brdfApprox, true);

  Vec G = getLightAnalytic(r, n, x, primitive, brdfApprox) * brdfApprox.amplitude * primitive->m.specularCoef;
  G = G + getLightAnalyticDiffuse(r, n, x, primitive) * (1 - primitive->m.specularCoef);

  for (uint32_t j = 0; j < mList.size(); j++) {
    Vec sum;
    // returns 1/pdf.
    double pdf = mList[j].getSamples(sampleCount, samples, ids);
    for (uint32_t i = 0; i < sampleCount; i++) {
      double visibility = 1;
      double brdfActual = 0;
      double brdfApproximate = 0;
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

      brdfActual = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d);
      brdfApproximate = primitive->brdf(n, wo_ref, r.d * -1.0, rr.d, brdfApprox);
      sum = sum +  mList[j].mesh[ids[i]].radiance * ((visibility * brdfActual- brdfApproximate) * cosine * cosine1 / (d * d));
    }
    sumLight = sumLight + (sum * (pdf/sampleCount));

    sumLight = sumLight + G;
  }
  return sumLight;
}