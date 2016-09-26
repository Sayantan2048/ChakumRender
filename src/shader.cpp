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
  double d, inf = t = 1e20;

  // Find the closest sphere with which the ray intersects.
  for(int i = nSpheres; i--; ) {
    if((d = list[i].intersect(r)) && d < t) {
      t = d; // set the distance of intersection.
      id = i; // Set the serial no. of the intersecting sphere
    }
  }

  // return true if the intersection distance is finite.
  return t < inf;
}

// Return 0 if shadowed else 1.
double shadow(const Ray &shadowRay, double distanceLightSource) {
  // This code casts shadow, Comment the next few lines to remove shadows.
  // If the ray from point to a light source hit any object.
  // WARNING: SHADWO FROM ANY MaterialType AND ANY PRIMITIVE
  double t;
  int id;
  if (intersectSphere(shadowRay, t, id, nSpheres, sphereList))
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

  Vec light = getLightFromPointSources(r, n, x, sphere, id, list);

  return obj.c.mult(light); // Compute color of intersection.
}

Vec shade(const Ray &r) {
  double t;
  int id = 0;

  // return zero vector if no intersection is found.
  if (!intersectSphere(r, t, id, nSpheres, sphereList))
    return Vec();

  return shadeSphere(r, t, id, sphereList);
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
}


