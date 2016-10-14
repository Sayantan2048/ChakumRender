#include "geometryPrimitives.h"
#include "materialTypes.h"
#include "objLoader.h"
#include "objects.h"
#include "mathPrimitives.h"
#include "transformations.h"
#include "dummyAccel.h"
#include "bvhAccel.h"
#include <iostream>

int nSpheres = 1;
Sphere sphereList[] = {
  //Sphere(1e5,  Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), 1.0, lambertian),
  //Sphere(1e5,  Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), 1.0, lambertian),
  //Sphere(1e5,  Vec(50, 40.8, 1e5), Vec(.75, .75, .75), 1.0, lambertian),
  Sphere(1e5,  Vec(50, 1e5, 81.6), Vec(.75, .75, .75), 1.0, lambertian),
  //Sphere(1e5,  Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), 1.0, lambertian),
  //Sphere(16.5,  Vec(27, 16.5, 47), Vec(.000, .999, .999), 0.2, phong),
  //Sphere(16.5,  Vec(73, 16.5, 78), Vec(.999, .999, .999), 0.3, phong),
  //Sphere(10,  Vec(50, 68.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian),
  //Sphere(10,  Vec(50, 10.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian)
};
int nTriangles = 0;
#define X (45)
#define Y 0
#define Z (-20)
#define VA Vec(27 + X, 16.5 + Y, 47 + Z)
#define VB Vec(59 + X, 16.5 + Y, 99 + Z)
#define VC Vec(55 + X, 50.5 + Y, 66 + Z)
#define VD Vec(73 + X, 16.5 + Y, 78 + Z)
Triangle *triangleList; // = {
 // Triangle(VA, VB, VC, Vec(0.0, 0, 0.9), 0.9, phong),
//  Triangle(VB, VC, VD, Vec(0.9, 0.9, 0), 0.9, phong),
//  Triangle(VC, VD, VA, Vec(0, 0.9, 0.9), 0.9, phong),
//  Triangle(VD, VA, VB, Vec(0.9, 0.0, 0.9), 0.9, phong),

//};

int load() {
  objLoader *objData = new objLoader();
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
	triangleList[i] = Triangle(Av, Bv, Cv, Vec(1.0, 0.01, 0.01), 1.0, phong);

  }
  bvhAccel = new BvhAccel((uint8_t *)triangleList, nTriangles, (std::size_t)sizeof(Triangle));
  bvhAccel -> initAccel();

  //DummyAccel::initAccel(nTriangles, triangleList);
}