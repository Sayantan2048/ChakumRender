#include "mathPrimitives.h"
#include "geometryPrimitives.h"
#include "materialTypes.h"
#include "objects.h"
#include "shader.h"
#include "lightSources.h"
#include "domainSampler.h"
#include "dummyAccel.h"
#include "bvhAccel.h"
#include "random.h"
#include <cmath>
#include <cstdio>
#include <iostream>

// return true if there is an intersection of a ray with any of the spheres.
static inline bool intersectSphere(const Ray &r, double &t, Vec &N, int &id, int nSpheres, Sphere *list) {
  /*double d = INF;
  id = 0xFFFFFFFF;
  t = INF;

  // Find the closest sphere with which the ray intersects.
  for(int i = nSpheres; i--; ) {
    if((d = list[i].intersect(r)) && d < t) {
      t = d; // set the distance of intersection.
      id = i; // Set the serial no. of the intersecting sphere
    }
  }

  N = list[id].p - (r.o + r.d * t);
  N.norm();

  // Use a normal towards the viewer.
  if (r.d.dot(N) > 0)
    N = N * -1.0;

  // return true if the intersection distance is finite.
  return t < INF;*/

  bool hit = bvhAccelS->intersectS(r, t, id);
  // if (hit) std::cerr<<"Oye Oye";
  if (!hit) return hit;

  N = list[id].p - (r.o + r.d * t);
  N.norm();

  // Use a normal towards the viewer.
  if (r.d.dot(N) > 0)
    N = N * -1.0;

  return hit;
}

// return true if there is an intersection of a ray with any of the spheres.
static inline bool intersectTriangle(const Ray &r, double &t, Vec &N, int &id, int nTriangles, Triangle *list) {
  /*double d;
  Vec N_;
  id = 0xFFFFFFFF;
  t = INF;

  // Find the closest sphere with which the ray intersects.
  for(int i = nTriangles; i--; ) {
    if((d = list[i].intersect(r, N_)) && d < t) {
      t = d; // set the distance of intersection.
      id = i; // Set the serial no. of the intersecting sphere
      N = N_;
    }
  }

  // return true if the intersection distance is finite.
  return t < INF;*/
  return bvhAccelT->intersectT(r, t, N, id);
  //return DummyAccel::intersect(r, t, N, id, nTriangles, list);
}

// Return 0 if shadowed else 1.
double shadow(const Ray &shadowRay, double distanceLightSource) {
  // This code casts shadow, Comment the next few lines to remove shadows.
  // If the ray from point to a light source hit any object.
  // WARNING: SHADWO FROM ANY MaterialType AND ANY PRIMITIVE
  double t;
  int id;
  double eps = 1e-5;//Avoid shadow due to self-intersection with light source.
  Vec N;
  if (intersectSphere(shadowRay, t, N, id, nSpheres, sphereList))
    if (t + eps < distanceLightSource) // Check whether the object is closer than light source.
      return 0.0;

  if (intersectTriangle(shadowRay, t, N, id, nTriangles, triangleList))
    if (t + eps < distanceLightSource) // Check whether the object is closer than light source.
      return 0.0;
  /*.
   * TODO.
   * .
   * . intersect Triangles and other geometryPrimitives
   */
  return 1.0;
}

// directIllumination shading
inline Vec shadeDI(const Ray &r,const Vec &x, const Vec &N, BasePrimitive *list) {
 /* Vec light = lSource->getLightFromPointSources(r, N, x, list) +
	lSource->getLightFromSphereSources(r, N, x, list, 20) + lSource->getLightFromEnvSource(r, N, x, list, 1000)
	+ lSource->getLightFromTriSources(r, N, x, list, 100) + lSource->getLightFromMeshSources(r, N, x, list, 160);*/
  //Vec light = lSource->getLightFromMeshSource_Analytic(r, N, x, list);
  Vec light = lSource->getLightFromAllSources_MIS(r, N, x, list);
  return list->c.mult(light); // Compute color of intersection.
}

Vec shadeExplicit(const Ray &r) {
  double terminationProbability = 0.3;
  double tS = INF, tT = INF;
  int idS = 0, idT = 0;
  Vec nT, nS;

  Random random(clock() & 0xFFFFFFFF);

  Vec col = Vec();
  Ray ri = r;
  std::vector<Vec> DI;
  std::vector<Vec> TP;

  DI.reserve(3);
  TP.reserve(3);

  // Returns normalized N.
  bool iTriangle = intersectTriangle(ri, tT, nT, idT, nTriangles, triangleList);
  bool iSphere = intersectSphere(ri, tS, nS, idS, nSpheres, sphereList);

  if (!(iSphere || iTriangle))
    return Vec();

  Vec x = tS < tT ? ri.o + ri.d * tS : ri.o + ri.d * tT;
  Vec n = tS < tT ? nS: nT;
  BasePrimitive *ptr = tS < tT ? (BasePrimitive *)&sphereList[idS]: (BasePrimitive *)&triangleList[idT];

  if (ptr->m.l != NONE)
    return ptr->m.getRadiance(n, ri.d * -1.0);

  Vec x_new;
  Vec n_new;
  BasePrimitive *ptr_new;

  for (int i = 0; i < 100; i++) {
    double e = random.randomMTD(0, 1);

    DI.push_back(shadeDI(ri, x, n, ptr));

    if (e < terminationProbability) {
      TP.push_back(Vec());
      break;
    }

    Vec wo_ref =  n * 2.0 * (n.dot(ri.d * -1.0)) + ri.d;
    wo_ref.norm();
    Ray secondaryRay(0,0);
    secondaryRay.o = x + n * 1e-6;
    Vec sample[1];
    double pdf = 0;
    int rayStatus = 1;

    do {
      rayStatus = 1;
      //ptr->m.getBrdfDirectionSampleRandomized(n, wo_ref, ri.d * -1.0, sample[0], pdf, random);
      //secondaryRay.d = sample[0];
      pdf = SphericalSampler::getHemiDirectionSamples(n, 1, sample);
      secondaryRay.d = sample[0];

      iTriangle = intersectTriangle(secondaryRay, tT, nT, idT, nTriangles, triangleList);
      iSphere = intersectSphere(secondaryRay, tS, nS, idS, nSpheres, sphereList);

      if (!(iSphere || iTriangle)) {
	rayStatus = 2; // Ray escapes scene
	break;
      }

      x_new = tS < tT ? secondaryRay.o + secondaryRay.d * tS : secondaryRay.o + secondaryRay.d * tT; // Ch
      n_new = tS < tT ? nS: nT; // Ch
      ptr_new = tS < tT ? (BasePrimitive *)&sphereList[idS]: (BasePrimitive *)&triangleList[idT];

      if (ptr_new->m.l != NONE)
	rayStatus = 0;

     } while (rayStatus == 0);

     if (rayStatus == 2) {
       	TP.push_back(Vec());
	break;
     }
     double cosine = secondaryRay.d.dot(n);
     double brdf = ptr->brdf(n, wo_ref, ri.d * -1.0, secondaryRay.d);
     double product = (pdf * cosine * brdf);
     TP.push_back(ptr->c*product);

     ri = secondaryRay;
     x = x_new;
     n = n_new;
     ptr = ptr_new;
  }

  col = DI[DI.size() -1];
  for (int i = DI.size() - 2; i >= 0; i--) {
    col = DI[i] + col.mult(TP[i]) / (1 - terminationProbability);
  }
  return col;
}
#define EXPLICIT_SAMPLES 50
inline Vec shadeDI2(const Ray &r,const Vec &x, const Vec &N, BasePrimitive *list) {
  /*Vec light = lSource->getLightFromPointSources(r, N, x, list) +
	lSource->getLightFromSphereSources(r, N, x, list, 20) + lSource->getLightFromEnvSource(r, N, x, list, 1000)
	+ lSource->getLightFromTriSources(r, N, x, list, 100) + lSource->getLightFromMeshSources(r, N, x, list, 1);*/
  Vec light = lSource->getLightFromMeshSource_Analytic(r, N, x, list);
  return list->c.mult(light); // Compute color of intersection.
}
Vec shadeExplicitFirstBounce(const Ray &r) {
  double tS = INF, tT = INF;
  int idS = 0, idT = 0;
  Vec nT, nS;
  const uint32_t maxSamples = EXPLICIT_SAMPLES * 3;
  Vec samples[maxSamples];

  // Returns normalized N.
  bool iTriangle = intersectTriangle(r, tT, nT, idT, nTriangles, triangleList);
  bool iSphere = intersectSphere(r, tS, nS, idS, nSpheres, sphereList);

  if (!(iSphere || iTriangle)) return Vec();

  Vec x = tS < tT ? r.o + r.d * tS : r.o + r.d * tT;
  Vec n = tS < tT ? nS: nT;
  BasePrimitive *ptr = tS < tT ? (BasePrimitive *)&sphereList[idS]: (BasePrimitive *)&triangleList[idT];

  Vec hitLightSource = ptr->m.getRadiance(n, r.d * -1.0);

  Vec directIllumination = (iSphere || iTriangle) ? ((ptr->m.l != NONE) ? hitLightSource : shadeDI2(r, x, n, ptr)) : 0;

  if (ptr->m.l != NONE)
	return hitLightSource;

  double pdfInv = SphericalSampler::getHemiDirectionSamples(n, maxSamples, samples);

  Vec wo_ref =  n * 2.0 * (n.dot(r.d * -1.0)) + r.d;
  wo_ref.norm();

  Ray secondaryRay(0,0);
  secondaryRay.o = x + n * 1e-6;

  uint32_t nIndirectRays = 0;
  Vec indirectIllumination;
  for (uint32_t i = 0; i < maxSamples && nIndirectRays < EXPLICIT_SAMPLES; i++) {
    secondaryRay.d = samples[i];

    iTriangle = intersectTriangle(secondaryRay, tT, nT, idT, nTriangles, triangleList);
    iSphere = intersectSphere(secondaryRay, tS, nS, idS, nSpheres, sphereList);

    if (!(iSphere || iTriangle))
	continue; // Ray escapes scene

    BasePrimitive *ptr_new = tS < tT ? (BasePrimitive *)&sphereList[idS]: (BasePrimitive *)&triangleList[idT];

    if (ptr_new->m.l != NONE)
	continue; // Hits a light source

    double cosine = secondaryRay.d.dot(n);
    double brdf = ptr->brdf(n, wo_ref, r.d * -1.0, secondaryRay.d);
    double product = (pdfInv * cosine * brdf);

    indirectIllumination = indirectIllumination + shadeExplicit(secondaryRay).mult((ptr->c * product));

    nIndirectRays++;
  }

  return directIllumination + indirectIllumination / nIndirectRays;
}

Vec shadeImplicit(const Ray &r) {
  double tS = INF, tT = INF;
  int idS = 0, idT = 0;
  Vec nT, nS;

  Random random(clock() & 0xFFFFFFFF);

  Vec product = Vec(1.0, 1.0, 1.0);
  Vec col = Vec();
  Ray ri = r;

  for (int i = 0; i < 4; i++) {
    // Returns normalized N.
    bool iTriangle = intersectTriangle(ri, tT, nT, idT, nTriangles, triangleList);
    bool iSphere = intersectSphere(ri, tS, nS, idS, nSpheres, sphereList);

    if (!(iSphere || iTriangle)) {
      col = Vec();
      break;
    }

    Vec x = tS < tT ? ri.o + ri.d * tS : ri.o + ri.d * tT;
    Vec n = tS < tT ? nS: nT;
    BasePrimitive *ptr = tS < tT ? (BasePrimitive *)&sphereList[idS]: (BasePrimitive *)&triangleList[idT];

    if ((iSphere || iTriangle)) {
      if (ptr->m.l != NONE) {
	//std::cout<<product.x <<" "<<product.y<<" "<<product.z<<"\n";
	Vec hitLightSource = ptr->m.getRadiance(n, ri.d * -1.0);
	col = hitLightSource.mult(product);
	//std::cout<<col.x <<" "<<col.y<<" "<<col.z<<"\n";
	break;
      }
      else {

	Vec wo_ref =  n * 2.0 * (n.dot(ri.d * -1.0)) + ri.d;
	wo_ref.norm();

	Ray secondaryRay(0,0);
	secondaryRay.o = x + n * 1e-6;
	Vec sample[1];
	double pdfInv = 1;
	double cosine = 1;
	if (1) {
	  pdfInv = SphericalSampler::getCosineDirectionSamples(n, 1, sample);
	  secondaryRay.d = sample[0];
	}
	else {
	  ptr->m.getBrdfDirectionSampleRandomized(n, wo_ref, ri.d * -1.0, sample[0], pdfInv, random);
	  secondaryRay.d = sample[0];
	  cosine = secondaryRay.d.dot(n);
	}

	double brdf = ptr->brdf(n, wo_ref, ri.d * -1.0, secondaryRay.d);
	product = product * (pdfInv * cosine * brdf);
	product = ptr->c.mult(product);
	ri = secondaryRay;
      }
    }
  }
  return col;
}

Vec shadeDirectOnly(const Ray &r) {
  double tS = INF, tT = INF;
  int idS = 0, idT = 0;
  Vec nT, nS;

  // Returns normalized N.
  bool iTriangle = intersectTriangle(r, tT, nT, idT, nTriangles, triangleList);
  bool iSphere = intersectSphere(r, tS, nS, idS, nSpheres, sphereList);

  if (!(iSphere || iTriangle)) return Vec();

  Vec x = tS < tT ? r.o + r.d * tS : r.o + r.d * tT;
  Vec n = tS < tT ? nS: nT;
  BasePrimitive *ptr = tS < tT ? (BasePrimitive *)&sphereList[idS]: (BasePrimitive *)&triangleList[idT];

  Vec hitLightSource = ptr->m.getRadiance(n, r.d * -1.0);

  Vec directIllumination = (iSphere || iTriangle) ? ((ptr->m.l != NONE) ? hitLightSource : shadeDI(r, x, n, ptr)) : 0;

  return directIllumination;
}

// R.d must be normalized before passing to shade.
Vec shade(const Ray &r, int &depth) {
#if 1
  return shadeDirectOnly(r);
#endif

#if 0
  return shadeImplicit(r);
#endif

#if 1
  return shadeExplicit(r);
#endif

#if 0
  return shadeExplicitFirstBounce(r);
#endif
}


