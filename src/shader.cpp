#include "mathPrimitives.h"
#include "geometryPrimitives.h"
#include "materialTypes.h"
#include "objects.h"
#include "shader.h"
#include "lightSources.h"
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
  double d;
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
  return t < INF;
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

  Vec light = getLightFromPointSources(r, n, x, sphere, id, list, 0) +
	getLightFromVolumeSources(r, n, x, sphere, id, list, 0);

  return obj.c.mult(light); // Compute color of intersection.
}

inline Vec shadeTriangle(const Ray &r,const double t, const Vec &N, int id, Triangle *list) {
  const Triangle &obj = list[id];
  Vec x = r.o + r.d * t; // Calulate the point of intersection.

  Vec light = getLightFromPointSources(r, N, x, triangle, id, 0, list) +
	getLightFromVolumeSources(r, N, x, triangle, id, 0, list);

  return obj.c.mult(light); // Compute color of intersection.
}

Vec shade(const Ray &r) {
  double tS, tT;
  int idS, idT;
  Vec N;

  bool iTriangle = intersectTriangle(r, tT, N, idT, nTriangles, triangleList);
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

  if (iSphere || iTriangle)
    return tS<tT?shadeSphere(r, tS, idS, sphereList):shadeTriangle(r, tT, N, idT, triangleList);

  return Vec();
}


