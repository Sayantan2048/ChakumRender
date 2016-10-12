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

AABBox::AABBox(const Vec &p1, const Vec &p2) {
  pMin = Vec(miN(p1.x, p2.x), miN(p1.y, p2.y), miN(p1.z, p2.z));
  pMax = Vec(maX(p1.x, p2.x), maX(p1.y, p2.y), maX(p1.z, p2.z));
}

AABBox::AABBox(const Vec &p, double r) {
  pMin = Vec(p.x - r, p.y - r, p.z - r);
  pMax = Vec(p.x + r, p.y + r, p.z + r);
}

AABBox AABBox::uNion(const AABBox &b, const Vec &p) {
  AABBox ret = b;
  ret.pMin.x = miN(b.pMin.x, p.x);
  ret.pMin.y = miN(b.pMin.y, p.y);
  ret.pMin.z = miN(b.pMin.z, p.z);
  ret.pMax.x = maX(b.pMax.x, p.x);
  ret.pMax.y = maX(b.pMax.y, p.y);
  ret.pMax.z = maX(b.pMax.z, p.z);
  return ret;
}

AABBox AABBox::uNion(const AABBox &a, const AABBox &b) {
  AABBox ret = b;
  ret.pMin.x = miN(a.pMin.x, b.pMin.x);
  ret.pMin.y = miN(a.pMin.y, b.pMin.y);
  ret.pMin.z = miN(a.pMin.z, b.pMin.z);
  ret.pMax.x = maX(a.pMax.x, b.pMax.x);
  ret.pMax.y = maX(a.pMax.y, b.pMax.x);
  ret.pMax.z = maX(a.pMax.z, b.pMax.x);
  return ret;
}

#define swap(a,b) {tmp = a; a = b; b = tmp;}
double AABBox::intersect(const Ray &ray) const {
  double t0 = -INF, t1 = INF;
  double tmp;

  double invRayDir  = 1. / ray.d.x;
  double tNear = (pMin.x - ray.o.x) * invRayDir;
  double tFar = (pMax.x - ray.o.x) * invRayDir;
  
  if (tNear > tFar) swap(tNear, tFar);
  t0 = tNear > t0 ? tNear : t0;
  t1 = tFar < t1 ? tFar : t1;
  if (t0 > t1) return INF;

  invRayDir  = 1.f / ray.d.y;
  tNear = (pMin.y - ray.o.y) * invRayDir;
  tFar = (pMax.y - ray.o.y) * invRayDir;

  if (tNear > tFar) swap(tNear, tFar);
  t0 = tNear > t0 ? tNear : t0;
  t1 = tFar < t1 ? tFar : t1;
  if (t0 > t1) return INF;

  invRayDir  = 1.f / ray.d.z;
  tNear = (pMin.z - ray.o.z) * invRayDir;
  tFar = (pMax.z - ray.o.z) * invRayDir;

  if (tNear > tFar) swap(tNear, tFar);
  t0 = tNear > t0 ? tNear : t0;
  t1 = tFar < t1 ? tFar : t1;
  if (t0 > t1) return INF;

  return t0;
}
#undef swap

int AABBox::maximumExtent() const {
  Vec diag = pMax - pMin;
  if (diag.x > diag.y && diag.x > diag.z)
    return 0;
  else if (diag.y > diag.z)
    return 1;
  else
    return 2;
}