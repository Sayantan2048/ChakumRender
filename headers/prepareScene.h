#if !defined(prepareScene_h__)
#define prepareScene_h__
#include "sceneParser.h"
#include "mathPrimitives.h"
#include "lightSources.h"

extern bool initSceneParser(const char *fileName);
extern bool prepareScreenCamera(uint32_t &width, uint32_t &height, Vec &camPos, Vec &camDir, Vec &cx, Vec &cy, bool verbose);
extern bool prepareLightSources(bool verbose);
extern bool prepareObjects(bool verbose);
#endif
