#include "geometryPrimitives.h"
#include "mathPrimitives.h"
#include <cmath>
#include "materialTypes.h"
#include <cstdio>

double BasePrimitive::brdf(Vec n, Vec wo, Vec wi, Vec x) const {
  if (m == lambertian)
    return Lambertian::brdf() * reflectance;
  else if (m == phong)
    return Phong::brdf(n, wo, wi) * reflectance;
  else
      fprintf(stderr, "Invalid material type.\n");

  return 0;
}

// Always return positve t or Does not intersect!!
double Sphere::intersect(const Ray &ray) const {
  Vec op = p - ray.o; // Vector from center of circle to ray origin
  double t; // Small error value.
  double b = op.dot(ray.d); // length of projection of op on ray.d. Assuming d is normalized!! Otherwise this quantity b should be diveded by ray.d.norm().
  double det = b * b - op.dot(op) + r * r; // If the ray just grazes the circle, then det is zero, if it intersect then det is greater than zero.

  if (det < 0) return INF;
  else det = sqrt(det); // In case the ray intersects, det is the distance between point of intersection and foot of the perpendicular from center.

  // Considering only the first intersection between ray and surface:
  // For covex intersection, b > 0 and b > det, so outer ?: handles this condition.
  // For concave intersection, det > abs(b), so inner ?: handles this condition.
  // For concave intersection there are two possiblites, 1. op and ray.d are in same direction, hence b > 0.
  // 2. op and ray.d in opposite directions, hence b < 0. But for both cases, t = b + det.
  return (t = b - det) > eps ? t : ((t = b + det) > eps ? t: INF); // t is always positve.
}

double Triangle::intersect(const Ray &ray, Vec &N) const {
  // Ray plane intersection formula: t = (A.dot(n) - ray.o.dot(n))/(ray.d.dot(n))
  double raydotn = ray.d.dot(n);
  double t;

  N = n;

  //For trangle we have two direction for normals, I choose the one which pointing towards ray.o and in opposite direction to ray.d.
  if (raydotn > 0) {
    N = n * -1.0;
    raydotn = -raydotn;
  }

  if (raydotn*raydotn <= eps) return INF;

  t = ((A - ray.o).dot(N))/raydotn;

  if (t < 0) return INF;

  Vec P = ray.o + ray.d * t;

  Vec v0 = P - A;
  Vec v1 = B - A;
  Vec v2 = C - A;

  double k0 = v0.dot(v1);
  double k1 = v0.dot(v2);

  double b0 = v1.dot(v1);
  double b1 = v1.dot(v2);

  double c0 = v2.dot(v1); // c0 = b1?? WTF?
  double c1 = v2.dot(v2);

  double gamma = (k0*b1 - k1*b0)/(c0*b1 - c1*b0);
  double beta  = (k0 - gamma*c0)/b0;

  if (gamma + beta > 1) return INF;
  else if (0 > gamma || gamma > 1) return INF;
  else if (0 > beta || beta > 1) return INF;

  return t;
}