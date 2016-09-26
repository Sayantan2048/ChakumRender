#include "lightSources.h"
#include "mathPrimitives.h"
#include "shader.h"
#include "geometryPrimitives.h"
#include <cstdio>

int nPointSources = 1;
PointSource pSources[] = {
  PointSource(Vec(50, 68.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 40.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 28.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0)),
  PointSource(Vec(50, 10.6 - .27, 81.6), Vec(40000.0, 40000.0, 40000.0))
};


Vec getLightFromPointSources(const Ray &r, const Vec &n, const Vec &x, PrimitiveType m, int id, Sphere *sphereList) {
  int convex = 1; // check whether the surface is convex or not from the direction of viewing
  Vec refd; // direction from point toward light source
  Ray rr(0,0); // a ray from point toward light source.
  double cosine = 0;
  double eps = 1e-08;
  double brdf = 0;

  Vec sumLight = Vec(0.0, 0.0, 0.0);
  for (int i = 0; i < nPointSources; i++) {
    refd = (pSources[i].p - x).norm(); // direction from point towards source

    // test for convex or concave as seen from a point.
    if (r.d.dot(n) >= 0)
      convex = 0;

    rr.d = refd;
    rr.o = x + n * (convex?eps:-eps);

    if (shadow(rr, (pSources[i].p - x).length()) < 0.5) {
      continue;
    }
    cosine = rr.d.dot(n);

    if (cosine < 0) {
    cosine = -cosine;
      if (convex) { // This code colors the dark side of an object black.
	continue;
      }
    }

    if (m == sphere)
      brdf = sphereList[id].brdf(rr.d, rr.d, rr.d);
    else
	fprintf(stderr, "PrimitiveType not found in getLightFromPointSources.\n");
    /* TODO: ADD BRDF FOR OTHER PrimitiveTypes */
    sumLight = sumLight + (pSources[i].radiance(x) * (cosine * brdf));
  }

  return sumLight;
}