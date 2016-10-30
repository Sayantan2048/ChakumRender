#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>
#include <cstdio>

double MaterialType::brdf(Vec n, Vec wo, Vec wi) {
  double cosine = wo.dot(wi);
  //cosine > 0?cosine:0
  return (specularCoef * (phongExp + 1) * pow (cosine > 0?cosine:0, phongExp) / (2 * PI)) + (1 - specularCoef) / PI;
}
