#include "geometryPrimitives.h"
#include "materialTypes.h"
#include "objects.h"

int nSpheres = 6;
Sphere sphereList[] = {
  Sphere(1e5,  Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), 1.0, lambertian),
  Sphere(1e5,  Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), 1.0, lambertian),
  Sphere(1e5,  Vec(50, 40.8, 1e5), Vec(.75, .75, .75), 1.0, lambertian),
  Sphere(1e5,  Vec(50, 1e5, 81.6), Vec(.75, .75, .75), 1.0, lambertian),
  Sphere(1e5,  Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), 1.0, lambertian),
  //Sphere(16.5,  Vec(27, 16.5, 47), Vec(.999, .999, .999), 1.0, lambertian),
  Sphere(16.5,  Vec(73, 16.5, 78), Vec(.999, .999, .999), 1.0, lambertian),
  //Sphere(10,  Vec(50, 68.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian),
  //Sphere(10,  Vec(50, 10.6 - .27, 81.6), Vec(.999, .999, .999), 1.0, lambertian)
};

int nTriangles = 4;
#define X 0
#define Y -5
#define Z 0
Triangle triangleList[] = {
  Triangle(Vec(27 + X, 16.5 + Y, 47 + Z),  Vec(73 + X, 16.5 + Y, 78 + Z), Vec(55 + X, 50.5 + Y, 66 + Z), Vec(0.9, 0.9, 0.9), 0.9, lambertian),
  Triangle(Vec(59 + X, 16.5 + Y, 99 + Z),  Vec(73 + X, 16.5 + Y, 78 + Z), Vec(55 + X, 50.5 + Y, 66 + Z), Vec(0.9, 0.9, 0.9), 0.9, lambertian),
  Triangle(Vec(59 + X, 16.5 + Y, 99 + Z),  Vec(73 + X, 16.5 + Y, 78 + Z), Vec(27 + X, 16.5 + Y, 47 + Z), Vec(0.9, 0.9, 0.9), 0.9, lambertian),
  Triangle(Vec(59 + X, 16.5 + Y, 99 + Z),  Vec(55 + X, 50.5 + Y, 66 + Z), Vec(27 + X, 16.5 + Y, 47 + Z), Vec(0.9, 0.0, 0.9), 0.9, lambertian),

};