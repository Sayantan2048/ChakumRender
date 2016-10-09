#include "mathPrimitives.h"
#include <cmath>

// assumes previous axis ex, ey, ez(canonical).
// transform z axis to newZ.
mat3 transformAxis(Vec newZ) {
   Vec z = Vec(0., 0., 1.);
   newZ.norm();
   Vec newX = (newZ%z).norm();
   Vec newY = (newZ%newX).norm();
   return mat3(newX, newY, newZ);
}

mat3 rotateZ(double theta) {
    return mat3(Vec(cos(theta), sin(theta), 0.),
              Vec(-sin(theta), cos(theta), 0),
	      Vec(0., 0., 1.));
}

mat3 rotateY(double theta) {
    return mat3(Vec(cos(theta), 0., sin(theta)),
		Vec(0., 1., 0.),
		Vec(-sin(theta), 0., cos(theta)));
}

mat3 scale(double sx, double sy, double sz) {
  return mat3(Vec(sx, 0., 0.),
              Vec(0, sy, 0),
	      Vec(0., 0., sz));
}

