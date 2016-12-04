#include "mathPrimitives.h"
#include "materialTypes.h"
#include <cmath>
#include <iostream>

double MaterialType::brdf(Vec n, Vec wo_ref, Vec wo, Vec wi) {
  if (b == CLASSIC_PHONG) {
    double cosine = wo_ref.dot(wi);
    //cosine > 0?cosine:0
    return (specularCoef * (phongExp + 1) * pow (cosine > 0?cosine:0, phongExp) / (2 * PI)) + (1 - specularCoef) / PI;
  }
  else if (b == GGX) {
      Vec h = (wo + wi).norm();

      double F = fresnel(h.dot(wi));
      double D = ggxDistribution(h, n);
      double G = ggxG1(wo, h, n) * ggxG1(wi, h, n);

      return (specularCoef * F * D * G) / (4.0 * wi.dot(n) * wo.dot(n)) + (1 - specularCoef) / PI;
  }
  else
    std::cerr<<"Undefined BRDF!!\n";

  return 0.0;
}
//See Microfacet BSDF paper, Bruce Walter
double inline MaterialType::fresnel(double c) const {
  double g = intIOR/extIOR;
  g *= g;
  g += c*c - 1.0;
  g = sqrt(g);

  double t1 = g - c;
  double t2 = g + c;

  double temp = c * t2 - 1.0;
  double F = temp * temp;

  temp = c * t1 + 1.0;
  F /= (temp * temp);
  F += 1.0;

  temp = t1 / t2;
  temp *= temp;

  F *= 0.5 * temp;

  return F;
}

double MaterialType::ggxDistribution(const Vec &h,const Vec &n) const {
  double cosine_sq = h.dot(n);

  cosine_sq *= cosine_sq;
  double tan_sq = 1.0 / cosine_sq - 1.0;
  double alpha_sq = ggxAlpha * ggxAlpha;

  tan_sq += alpha_sq;
  tan_sq *= tan_sq;
  double D = alpha_sq / (PI * cosine_sq * cosine_sq * tan_sq);

  return D;
}

double MaterialType::ggxG1(const Vec &v, const Vec &h, const Vec &n) const {
  double cosv_sq = v.dot(n);

  cosv_sq *= cosv_sq;
  double tanv_sq = 1.0 / cosv_sq - 1.0;

  if (tanv_sq <= 1e-15)
    return 1.0;
  else if (h.dot(v) <= 1e-15)
    return 0;

  double alpha_sq = ggxAlpha * ggxAlpha;

  tanv_sq *= alpha_sq;
  tanv_sq += 1.0;
  tanv_sq = 1.0 + sqrt(tanv_sq);

  return 2.0/tanv_sq;
}

