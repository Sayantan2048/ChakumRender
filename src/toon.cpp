//All experimental stuff!!
#include "lightSources.h"
#include "shader.h"
#include <cmath>

#define MESH_SAMPLES 80
static double integrateEdge(Vec v1, Vec v2, Vec n) {
    double cosTheta = v1.dot(v2);
    double theta = acos(cosTheta);
    double d = (v1%v2).dot(n);
    double res = d * ((theta > 0.001) ? theta/sin(theta) : 1.0);

    return res;
}

//https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
static bool clipTriangle(Vec A, Vec B,  Vec C, Vec n, Vec P, double nA, double nB, double nC, Vec *clip) {
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
/*  std::cout<<intersectAB<<" "<<intersectBC<<" "<<intersectCA<<" "<<nA<<" "<<nB<<" "<<nC<<" N:("<<n.x<<","<<n.y<<","<<n.z<<") P:("<<P.x<<","<<P.y<<","<<P.z<<") iAB:(";
  std::cout<<iAB.x<<","<<","<<iAB.y<<","<<iAB.z<<") iBC:("<<iBC.x<<","<<iBC.y<<","<<iBC.z<<") clip:";
  clip[0].show();
  clip[1].show();
  clip[2].show();
  clip[3].show();
  std::cout<<"\n";*/
  return (intersectAB + intersectBC + intersectCA) > 1;
}

Vec LightSource::getLightFromToonSource(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
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

	//sum = 0;
      }
      sum = sum > 0 ? sum : sum * -1;

      sumLight = sumLight + mList[j].mesh[i].radiance * sum;
    }
  }

  return sumLight / (2.0 * PI);
}

Vec LightSource::getLightFromToonSourceMesh(const Ray &r, const Vec &n, const Vec &x, BasePrimitive *primitive) {
  if (mList.size() < 1) return Vec();
  Ray rr(0,0); // a ray from point toward random direction in hemispherical domain.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;
  Vec samples[MESH_SAMPLES] = {0};
  uint32_t ids[MESH_SAMPLES] = {0};

  Vec sumLight = Vec();
  rr.o = x + n * eps;

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();
  Vec G = getLightFromToonSource(r, n, x, primitive);

  for (uint32_t j = 0; j < mList.size(); j++) {
    Vec sum;
    Vec sumG;

    // returns 1/pdf.
    double pdf = mList[j].getSamples(MESH_SAMPLES, samples, ids);
    for (int i = 0; i < MESH_SAMPLES; i++) {
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

      brdf = primitive->brdf(n, wo_ref, rr.d, x);

      sum = sum +  mList[j].mesh[ids[i]].radiance * ((visibility - 1.0) * cosine * cosine1 * brdf / (d * d));
    }
    sumLight = sumLight + (sum * (pdf/MESH_SAMPLES));

    sumLight = sumLight + G;
  }
  return sumLight;
}