#if !defined (objects_h__)
#define objects_h__
#include <stdint.h>
#include "geometryPrimitives.h"
#include <vector>

extern uint32_t nSpheres;
extern Sphere *sphereList;

extern uint32_t nTriangles;
extern Triangle *triangleList;

// For ease of use during initial scene setup.
extern std::vector<Sphere> vSphereList;
extern std::vector<Triangle> vTriangleList;

extern void loadAccels();

#endif
