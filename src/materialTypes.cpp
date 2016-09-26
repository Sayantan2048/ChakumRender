#include "mathPrimitives.h"
#include "materialTypes.h"

double Lambertian::brdf() {
  return 1/PI;
}

double Diffuse::brdf() {
  return 1/PI;
}
