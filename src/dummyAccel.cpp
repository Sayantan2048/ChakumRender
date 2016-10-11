#include "geometryPrimitives.h"
#include "dummyAccel.h"
#include "objects.h"
#include <cmath>

int DummyAccel::nBoxT;
BoxContent * DummyAccel::boxT;
AABBox DummyAccel::box(Vec(1.,1.,1.));

void DummyAccel::initAccel(int nTriangles, Triangle *triangleList) {
  nBoxT = (int)sqrt(nTriangles);

  boxT = new BoxContent[nBoxT];

  int i;
  for (i = 0; i < nBoxT; i++) {
    boxT[i].start = i * nBoxT;
    boxT[i].end = (i + 1) * nBoxT;
    boxT[i].box = triangleList[i * nBoxT].box;
  }
  boxT[i - 1].end = nTriangles;

  for (int j = 0; j < nBoxT; j++)
  for (int i = boxT[j].start; i < boxT[j].end; i++)
    boxT[j].box = AABBox::uNion(boxT[j].box, triangleList[i].box);

  box = boxT[0].box;
  for (int i = 1; i < nBoxT; i++)
    box = AABBox::uNion(box, boxT[i].box);
}

bool DummyAccel::intersect(const Ray &r, double &t, Vec &N, int &id, int nTriangles, Triangle *list) {
  if (box.intersect(r) >= INF) return 0;

  double d;
  Vec N_;
  id = 0xFFFFFFFF;
  t = INF;

  for (int i = 0; i < nBoxT; i++) {
    if (boxT[i].box.intersect(r) < INF) {
      for(int j = boxT[i].start; j < boxT[i].end; j++) {
	if((d = list[j].intersect(r, N_)) && d < t) {
	  t = d; // set the distance of intersection.
	  id = i; // Set the serial no. of the intersecting sphere
	  N = N_;
	}
      }
    }
  }

  return t < INF;

}

