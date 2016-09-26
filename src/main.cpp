#include <iostream>
#include <cstdio>
#include "domainSampler.h"
#include "mathPrimitives.h"
#include "materialTypes.h"
#include "geometryPrimitives.h"
#include "objects.h"
#include "shader.h"
#include <cmath>

// clamp x between 0 and 1.
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }

// Clamp pixel values between 0 to 255. 1/2.2 is gamma correction factor!
inline int toDisplayValue(double x){ return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5); }

int main(int argc, char *argv[]) {
  int w = 512, h = 384;
  // camera location and direction of looking. Imagine right direction is x, up is y, and z is out of screen. Camera is mostly looking towards -z direction!!
  Ray camera( Vec(50, 50, 275.0), Vec(0, -0.05, -1).norm());
  //Define pixel pitch along width of the screen. Field of view is 30 + 30 or 60 degrees.
  Vec cx = Vec( w * 0.57735 / h, 0., 0.1); // hint : tan( 30 / 180.0 * M_PI ) == 0.57735
  //Define pixel pitch along height of the screen.
  Vec cy = (cx % camera.d).norm() * 0.57735;
  // 2D Array of pixels
  Vec pixelValue, *pixelColors = new Vec[w * h];

  #pragma omp parallel for schedule(dynamic, 1) private(pixelValue)
  for(int y = 0; y < h; y++) {
    // Percentage completion!!
    fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
    for(int x = 0; x < w; x++ ) {
      // Start from top left to bottm right of the screen.
      int idx = (h - y - 1) * w + x;
      pixelValue = Vec();
      // Shoot ray from camera thur each pixel: Computed as: camera direction +/- deviation from camera direction in terms of pixel pitch.
      Vec cameraRayDir = cx * ( double(x)/w - .5) + cy * ( double(y)/h - .5) + camera.d;
      // Find color of intersection. In case no intersection is found color the pixel black.
      pixelValue = shade(Ray(camera.o, cameraRayDir.norm()));
      // Clamp the rgb values.
      pixelColors[idx] = Vec(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
    }
  }

  fprintf(stderr,"\n");

  // Save pixelColors in ppm format.
  // hint: Google the PPM image format
  FILE *f = fopen("image.ppm", "w");
  // Format header.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int p = 0; p < w * h; p++) {
    //Save each pixel as rgb triplet.
    fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), toDisplayValue(pixelColors[p].z));
  }
  fclose(f);

  delete pixelColors;

  return 0;
}
