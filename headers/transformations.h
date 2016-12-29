#if !defined (transformations_h__)
#define transformations_h__

#include "mathPrimitives.h"

extern mat3 transformAxis(Vec newZ);
extern mat3 rotateZ(double theta);
extern mat3 rotateY(double theta);
extern mat3 rotate(const Vec &axis, double theta);
extern mat3 scale(double sx, double sy, double sz);

#endif