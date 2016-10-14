#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>
#include <cstdio>

double Phong::e = 2.0;

double Lambertian::brdf() {
  return 1/PI;
}

double Phong::brdf(Vec n, Vec wo, Vec wi) {
  Vec wr =  n * 2.0 * (n.dot(wi)) - wi;
  wr.norm();
  double cosine = wr.dot(wo);

  return (0.5 * (e + 2) * pow (cosine > 0?cosine:0, e) / (2 * PI)) + 0.5 / PI ;
}
