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
#define X (45)
#define Y 0
#define Z (-150)
#define VA Vec(27 + X, 16.5 + Y, 47 + Z)
#define VB Vec(59 + X, 16.5 + Y, 99 + Z)
#define VC Vec(55 + X, 50.5 + Y, 66 + Z)
#define VD Vec(73 + X, 16.5 + Y, 78 + Z)
Triangle * triangleList;

void loadObjects() {
  /*vSphereList.push_back(Sphere(1e5,  Vec(50, 1e5, 81.6), Vec(.75, .75, .75), 1.0, MaterialType(1.0, 0.5, NONE, Vec(0., 0., 0.))));
  vSphereList.push_back(Sphere(16.5,  Vec(27, 16.5, 47), Vec(.000, .999, .999), 1.0, MaterialType(1.0, 0.5, NONE, Vec(0., 0., 0.))));
  vSphereList.push_back(Sphere(16.5,  Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, MaterialType(1.0, 0.5, NONE, Vec(0., 0., 0.))));
  */
  /*vSphereList.push_back(Sphere(1e5,  Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25 * 0, .25 * 0), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));
  vSphereList.push_back(Sphere(1e5,  Vec(-1e5 + 99, 40.8, 81.6), Vec(.25 * 0, .25 * 0, .75), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));
  vSphereList.push_back(Sphere(1e5,  Vec(50, 40.8, 1e5), Vec(0*.75, .75, 0*.75), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));
  vSphereList.push_back(Sphere(1e5,  Vec(50, 1e5, 81.6), Vec(.75, .75, .75), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));
  vSphereList.push_back(Sphere(1e5,  Vec(50, -1e5 + 81.6, 81.6), Vec(0*.75, .75, .75), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));
  vSphereList.push_back(Sphere(16.5,  Vec(27, 16.5, 47), Vec(.999, .999, .999), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));
  vSphereList.push_back(Sphere(16.5,  Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, MaterialType(0.01, 10000.5, 1, 0.9)));*/
  //Dummy
  vSphereList.push_back(Sphere(1,  Vec(27, 1e5, 47), Vec(.999, .999, .999), 1.0, MaterialType(200.2, 1)));
  //Sphere(1e5,  Vec(1e5 + 1, 40.8, 81.6), Vec(0, 0, 0), Vec(.75, .25, .25)),
  //Sphere(1e5,  Vec(-1e5 + 99, 40.8, 81.6), Vec(0, 0, 0), Vec(.25, .25, .75)),
  //Sphere(1e5,  Vec(50, 40.8, 1e5), Vec(0, 0, 0), Vec(.75, .75, .75)),
  //Sphere(1e5,  Vec(50, 1e5, 81.6), Vec(0, 0, 0), Vec(.75, .75, .75)),
  //Sphere(1e5,  Vec(50, -1e5 + 81.6, 81.6), Vec(0, 0, 0), Vec(.75, .75, .75)),
  //Sphere(16.5,  Vec(27, 16.5, 47), Vec(0, 0, 0), Vec(.999, .999, .999)),
  //Sphere(16.5,  Vec(73, 16.5, 78), Vec(0, 0, 0), Vec(.999, .999, .999) ),
  //Sphere(10.5, Vec(50, 68.6 - .27, 81.6), Vec(400, 400, 400), Vec(1, 1, 1))

  //Dummy Triangle
  vTriangleList.push_back(Triangle(VA + 1e6, VB + 1e6, VC + 1e6, Vec(0.0, 0, 0.9), 0.9, MaterialType(1.0, 0.5)));
  //Eric Veach Scene
  //vTriangleList.push_back(Triangle(Vec(-5000, -30, 300 - 100), Vec(5000, -30, 300 - 100), Vec(5000, -30, 200 - 100), Vec(0.9, 0.9, 0.9), 0.9, MaterialType(200000.0, 1)));
  //vTriangleList.push_back(Triangle(Vec(-5000, -30, 300 - 100), Vec(-5000, -30, 200 - 100), Vec(5000, -30, 200 - 100), Vec(0.9, 0.9, 0.9), 0.9, MaterialType(200000.0, 1)));
  //vTriangleList.push_back(Triangle(Vec(-5000, -30, 300 - 150), Vec(5000, -30, 300 - 150), Vec(5000, -30, 200 - 150), Vec(0.9, 0.9, 0.9), 0.9, MaterialType(20.0, 1)));
  //vTriangleList.push_back(Triangle(Vec(-5000, -30, 300 - 150), Vec(-5000, -30, 200 - 150), Vec(5000, -30, 200 - 150), Vec(0.9, 0.9, 0.9), 0.9, MaterialType(20.0, 1)));
  //vTriangleList.push_back(Triangle(Vec(-5000, -30, 5000 + Z), Vec(5000, -30, 5000 + Z), Vec(5000, -30, -5000 + Z), Vec(0.9, 0.9, 0.9), 0.9, MaterialType(20.0, 1)));
  //vTriangleList.push_back(Triangle(Vec(-5000, -30, 5000 + Z), Vec(-5000, -30, -5000 + Z), Vec(5000, -30, -5000 + Z), Vec(0.9, 0.9, 0.9), 0.9, MaterialType(20.0, 1)));
  //vTriangleList.push_back(Triangle(VA, VB, VC, Vec(0.0, 0, 0.9), 0.9, MaterialType(1.0, 0.5)));
  //vTriangleList.push_back(Triangle(VB, VC, VD, Vec(0.9, 0.9, 0), 0.9, MaterialType(1.0, 0.5)));
  //vTriangleList.push_back(Triangle(VC, VD, VA, Vec(0, 0.9, 0.9), 0.9, MaterialType(1.0, 0.5)));
  //vTriangleList.push_back(Triangle(VD, VA, VB, Vec(0.9, 0.0, 0.9), 0.9, MaterialType(1.0, 0.5)));

  //vTriangleList.push_back(Triangle(Vec(-5000, 0, 5000), Vec(5000, 0, 5000), Vec(5000, 0, -5000), Vec(0.9, 0.9, 0.9), 1.0, MaterialType(400.0, 1)));
  //vTriangleList.push_back(Triangle(Vec(-5000, 0, 5000), Vec(-5000, 0, -5000), Vec(5000, 0, -5000), Vec(0.9, 0.9, 0.9), 1.0, MaterialType(400.0, 1)));

  vTriangleList.push_back(Triangle(Vec(-5000, 0, 5000), Vec(5000, 0, 5000), Vec(5000, 0, -5000), Vec(0.9, 0.9, 0.9), 1.0, MaterialType(0.15, 10000, 1, 1, GGX)));
  vTriangleList.push_back(Triangle(Vec(-5000, 0, 5000), Vec(-5000, 0, -5000), Vec(5000, 0, -5000), Vec(0.9, 0.9, 0.9), 1.0, MaterialType(0.15, 10000, 1, 1, GGX)));
  /*  objLoader *objData = new objLoader();
  objData->load("Aventador1.obj");

  nTriangles = objData->faceCount;

  triangleList = new Triangle[nTriangles];

  mat3 M = scale(35, 35, 35).mul(rotateY(PI/4 + PI));

  for(int i=0; i < nTriangles; i++) {
	obj_face *o = objData->faceList[i];

	obj_vector A, B, C;
	A = *objData->vertexList[o->vertex_index[0]];
	B = *objData->vertexList[o->vertex_index[1]];
	C = *objData->vertexList[o->vertex_index[2]];

	Vec Av, Bv, Cv;
	Av = M.mul(Vec(A.e[0], A.e[1], A.e[2])) + Vec(X, Y, Z);
	Bv = M.mul(Vec(B.e[0], B.e[1], B.e[2])) + Vec(X, Y, Z);
	Cv = M.mul(Vec(C.e[0], C.e[1], C.e[2])) + Vec(X, Y, Z);

	// Assignment operator, data of Triangle() will be copied to triangleList[i].
	// Triangle() object is however temporary.
	triangleList[i] = Triangle(Av, Bv, Cv, Vec(1.0, 0.01, 0.01), 1.0, MaterialType(2.2, 1.0));

  }*/


  //DummyAccel::initAccel(nTriangles, triangleList);
}

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