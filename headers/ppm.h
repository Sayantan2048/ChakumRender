#if !defined (ppm_h__)
#define ppm_h__
#include "mathPrimitives.h"
#include <stdint.h>

int readImage(Vec * &image, uint32_t &w, uint32_t &h);
int writeImage(Vec *image, uint32_t w, uint32_t h);
#endif
