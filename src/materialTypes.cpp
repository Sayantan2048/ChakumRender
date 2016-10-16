#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>
#include <cstdio>

double MaterialType::brdf(Vec n, Vec wo, Vec wi) {
  Vec wr =  n * 2.0 * (n.dot(wi)) - wi;
  wr.norm();
  double cosine = wr.dot(wo);
  //cosine > 0?cosine:0
  return (specularCoef * (phongExp + 2) * pow (cosine > 0?cosine:0, phongExp) / (2 * PI)) + (1 - specularCoef) / PI;
}
