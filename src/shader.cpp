#include "mathPrimitives.h"
#include "geometryPrimitives.h"
#include "materialTypes.h"
#include "objects.h"
#include "shader.h"
#include "lightSources.h"
#include "domainSampler.h"
#include "dummyAccel.h"
#include "bvhAccel.h"
#include <cmath>
#include <cstdio>

// return true if there is an intersection of a ray with any of the spheres.
static inline bool intersectSphere(const Ray &r, double &t, int &id, int nSpheres, Sphere *list) {
  double d;
  id = 0xFFFFFFFF;
  t = INF;

  // Find the closest sphere with which the ray intersects.
  for(int i = nSpheres; i--; ) {
    if((d = list[i].intersect(r)) && d < t) {
      t = d; // set the distance of intersection.
      id = i; // Set the serial no. of the intersecting sphere
    }
  }

  // return true if the intersection distance is finite.
  return t < INF;
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
  //return bvhAccel->intersect(r, t, N, id);
  return DummyAccel::intersect(r, t, N, id, nTriangles, list);
}

// Return 0 if shadowed else 1.
double shadow(const Ray &shadowRay, double distanceLightSource) {
  // This code casts shadow, Comment the next few lines to remove shadows.
  // If the ray from point to a light source hit any object.
  // WARNING: SHADWO FROM ANY MaterialType AND ANY PRIMITIVE
  double t;
  int id;
  Vec N;
  if (intersectSphere(shadowRay, t, id, nSpheres, sphereList))
    if (t < distanceLightSource) // Check whether the object is closer than light source.
      return 0.0;

  if (intersectTriangle(shadowRay, t, N, id, nTriangles, triangleList))
    if (t < distanceLightSource) // Check whether the object is closer than light source.
      return 0.0;
  /*.
   * TODO.
   * .
   * . intersect Triangles and other geometryPrimitives
   */
  return 1.0;
}

//Return pixel color
inline Vec shadeSphere(const Ray &r,const double t, int id, Sphere *list) {
  const Sphere &obj = list[id];
  Vec x = r.o + r.d * t; // Calulate the point of intersection.
  Vec n = (x - obj.p).norm(); // Since this is a sphere, normal is a vector form sphere center to point of intersection.

  Vec light = getLightFromPointSources(r, n, x, &list[id]) +
	getLightFromVolumeSources(r, n, x, &list[id]);

  return obj.c.mult(light); // Compute color of intersection.
}

inline Vec shadeTriangle(const Ray &r,const double t, const Vec &N, int id, Triangle *list) {
  const Triangle &obj = list[id];
  Vec x = r.o + r.d * t; // Calulate the point of intersection.

  Vec light = getLightFromPointSources(r, N, x, &list[id]) +
	getLightFromVolumeSources(r, N, x, &list[id]);

  return obj.c.mult(light); // Compute color of intersection.
}

Vec shade(const Ray &r, int &depth) {
  double tS, tT;
  int idS, idT;
  Vec n;

  if (depth == 0)
    return Vec();

  depth--;

  bool iTriangle = intersectTriangle(r, tT, n, idT, nTriangles, triangleList);
  bool iSphere = intersectSphere(r, tS, idS, nSpheres, sphereList);

  /* TODO:Find the color of other primitive type and Chose the closest one.
  if (intersectType1 == 1 and intersectType2 == 0)
    return color of tyoe1.
  elseif (intersectType1 == 0 and intersectType2 == 1)
    return color of type2
  elseif ( intersectType1 == 1 and intersectType2 == 1){
    return color of the closest one.
  else
    return black.
    */
  Ray secondaryRay(0,0);

  if (iSphere || iTriangle) {
    Vec sum = Vec(0, 0, 0);
    if (tS < tT) {
      /*secondaryRay.o = r.o + r.d * tS;
      n = (secondaryRay.o - sphereList[idS].p).norm();
      secondaryRay.o = secondaryRay.o + n * 1e-06;
      secondaryRay.d = n;
     /* Vec *samples = new Vec[100];
      double pdf = SphericalSampler::getHemiSurfaceSamples(n, secondaryRay.o, 100, samples);
      for (int i = 0; i < 100; i++) {
	int depth2 = depth;
	secondaryRay.d = samples[i] - secondaryRay.o;
	sum  = sum + ((sphereList[idS].m == phong)? shade(secondaryRay, depth2) : Vec());
      }
      delete []samples;*/
     /*int depth2 = depth;
     sum = shade(secondaryRay, depth2);
     secondaryRay.d = (n + r.d * -1.0).norm();
     depth2 = depth;
     sum = sum + shade(secondaryRay, depth2);
     secondaryRay.d = (n + r.d * -1.5).norm();
     depth2 = depth;
     sum = sum + shade(secondaryRay, depth2);
     secondaryRay.d = (n + r.d * -2.5).norm();
     depth2 = depth;
     sum = sum + shade(secondaryRay, depth2);*/
     // fprintf(stderr, "%f %f %f\n", sum.x, sum.y, sum.z);*/
      return shadeSphere(r, tS, idS, sphereList) + sum * 0.25 * 0.2;
    }
    else {
     /* secondaryRay.o = r.o + r.d * tT;
      secondaryRay.o = secondaryRay.o + n * 1e-06;
      secondaryRay.d = n;
     int depth2 = depth;
     sum = shade(secondaryRay, depth2);
     secondaryRay.d = (n + r.d * -1.0).norm();
     depth2 = depth;
     sum = sum + shade(secondaryRay, depth2);
     secondaryRay.d = (n + r.d * -1.5).norm();
     depth2 = depth;
     sum = sum + shade(secondaryRay, depth2);
     secondaryRay.d = (n + r.d * -2.5).norm();
     depth2 = depth;
     sum = sum + shade(secondaryRay, depth2);
    // if (r.d.dot(n)>0)
    // fprintf(stderr, "%f %f %f\n", sum.x, sum.y, r.d.dot(n)    );*/
      return shadeTriangle(r, tT, n, idT, triangleList) + sum * 0.25 * 0.7;


    }

  }
  return Vec();
}


