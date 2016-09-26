#include "geometryPrimitives.h"
#include "mathPrimitives.h"
#include <cmath>
#include "materialTypes.h"
#include <cstdio>

double Sphere::intersect(const Ray &ray) const {
  Vec op = p - ray.o; // Vector from center of circle to ray origin
  double t, eps = 1e-4; // Small error value.
  double b = op.dot(ray.d); // length of projection of op on ray.d. Assuming d is normalized!! Otherwise this quantity b should be diveded by ray.d.norm().
  double det = b * b - op.dot(op) + r * r; // If the ray just grazes the circle, then det is zero, if it intersect then det is greater than zero.

  if (det < 0) return 0;
  else det = sqrt(det); // In case the ray intersects, det is the distance between point of intersection and foot of the perpendicular from center.

  // Considering only the first intersection between ray and surface:
  // For covex intersection, b > 0 and b > det, so outer ?: handles this condition.
  // For concave intersection, det > abs(b), so inner ?: handles this condition.
  // For concave intersection there are two possiblites, 1. op and ray.d are in same direction, hence b > 0.
  // 2. op and ray.d in opposite directions, hence b < 0. But for both cases, t = b + det.
  return (t = b - det) > eps ? t : ((t = b + det) > eps ? t: 0); // t is always positve.
}

double Sphere::brdf(Vec wi, Vec wr, Vec x) const {
  if (m == lambertian)
    return Lambertian::brdf() * reflectance;
  else if (m == diffuse)
    return Diffuse::brdf() * reflectance;
  else
      fprintf(stderr, "INVALID BRDF IN SPHERE.\n");

  return 0;
}