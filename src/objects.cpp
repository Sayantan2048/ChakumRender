#include "geometryPrimitives.h"
#include "materialTypes.h"
#include "objLoader.h"
#include "objects.h"
#include "mathPrimitives.h"
#include "transformations.h"
#include "dummyAccel.h"
#include "bvhAccel.h"
#include "lightSources.h"
#include <vector>
#include <iostream>

std::vector<Sphere> vSphereList;
std::vector<Triangle> vTriangleList;

uint32_t nSpheres = 0;
Sphere *sphereList;
/*Sphere sphereList[] = {
  //Sphere(1e5,  Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), 1.0, MaterialType(1.0, 0.0)),
  //Sphere(1e5,  Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), 1.0, MaterialType(1.0, 0.0)),
  //Sphere(1e5,  Vec(50, 40.8, 1e5), Vec(.75, .75, .75), 1.0, MaterialType(1.0, 0.0)),
  ,
  //Sphere(1e5,  Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), 1.0, MaterialType(1.0, 0.0)),
  Sphere(16.5,  Vec(27, 16.5, 47), Vec(.000, .999, .999), 1.0, MaterialType(1.0, 0.5)),
  Sphere(16.5,  Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, MaterialType(1.0, 0.5)),
  //Sphere(10,  Vec(50, 68.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, MaterialType(1.0, 0.0)),
  //Sphere(10,  Vec(50, 10.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, MaterialType(1.0, 0.0))
};*/

uint32_t nTriangles = 0;

#define VA Vec(27, 16.5, 47)
#define VB Vec(59, 16.5, 99)
#define VC Vec(55, 50.5, 66)
#define VD Vec(73, 16.5, 78)
Triangle * triangleList;

void loadAccels() {
  nTriangles = vTriangleList.size();
  triangleList = new Triangle[nTriangles];

  for (uint32_t i = 0; i < nTriangles; i++)
    triangleList[i] = vTriangleList[i];

  bvhAccelT = new BvhAccel((uint8_t *)triangleList, nTriangles, (std::size_t)sizeof(Triangle));
  bvhAccelT -> initAccel();

  nSpheres = vSphereList.size();
  sphereList = new Sphere[nSpheres];

  for (uint32_t i = 0; i < nSpheres; i++)
    sphereList[i] = vSphereList[i];

  bvhAccelS = new BvhAccel((uint8_t *)sphereList, nSpheres, (std::size_t)sizeof(Sphere));
  bvhAccelS -> initAccel();
}
