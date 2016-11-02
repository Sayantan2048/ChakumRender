#include "ppm.h"
#include "mathPrimitives.h"
#include <stdint.h>
#include <cstdio>

// Return 0 in case of failure.
/*int readImage(Vec * &image, uint32_t &width, uint32_t &height) {
  FILE *f = fopen("env1.ppm", "rb");
  uint32_t w, h, range;
  char string[100];
  if (f) {
    fscanf(f, "%s\n%d %d\n%d\n", string, &w, &h, &range);
    width = w;
    height = h;
    image = new Vec[width * height];

    printf("%d %d %d\n", w, h, range);

    uint32_t r, g, b;
    for (uint32_t i = 0; i <  width * height; i++) {
      fscanf(f, "%c %c %c ", &r, &g, &b);
      printf("%d %d %d\n", r, g, b);
      image[i].x = double(r)/255.0;
      image[i].y = double(g)/255.0;
      image[i].z = double(b)/255.0;
    }
    fclose(f);
  }
  else {
    fprintf(stderr, "Can't open requested file.\n");
    return 0;
  }

  return 1;
}*/
int readImage(Vec * &image, uint32_t &width, uint32_t &height) {
  char buff[16];
  FILE *f;
  int c, rgb_comp_color;
  char *filename = "env1.ppm";

  //open PPM file for reading
  f = fopen(filename, "rb");
  if (!f) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    return 0;
  }

  //read image format
  if (!fgets(buff, sizeof(buff), f)) {
    perror(filename);
    return 0;
  }

  //check the image format
  if (buff[0] != 'P' || buff[1] != '6') {
    fprintf(stderr, "Invalid image format (must be 'P6')\n");
    return 0;
  }

  //check for comments
  c = getc(f);
  while (c == '#') {
    while (getc(f) != '\n') ;
         c = getc(f);
  }

  ungetc(c, f);

  //read image size information
  if (fscanf(f, "%d %d", &width, &height) != 2) {
    fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
    return 0;
  }

  //read maxCol component
  if (fscanf(f, "%d", &rgb_comp_color) != 1) {
    fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
    return 0;
  }

    //check rgb component depth
  if (rgb_comp_color!= 255) {
    fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
    return 0;
  }

  while (fgetc(f) != '\n');

  //memory allocation for pixel data
  image = new Vec[width * height];
  unsigned char *pixelRawValues = new unsigned char[width * height * 3];

  if (!image) {
    fprintf(stderr, "Unable to allocate memory\n");
    return 0;
  }

  //read pixel data from file
  if (fread(pixelRawValues, 3 * width, height, f) != height) {
    fprintf(stderr, "Error loading image '%s'\n", filename);
    return 0;
  }

  for (uint32_t i = 0; i <  width * height; i++) {
      //printf("%d %d %d\n", pixelRawValues[i * 3], pixelRawValues[i * 3 + 1], pixelRawValues[i * 3 + 2]);
      image[i].x = double(pixelRawValues[i * 3])/255.0;
      image[i].y = double(pixelRawValues[i * 3 + 1])/255.0;
      image[i].z = double(pixelRawValues[i * 3 + 2])/255.0;
    }

  delete []pixelRawValues;
  fclose(f);
  return 1;
}

int writeImage(Vec *image, uint32_t width, uint32_t height) {
  FILE *f;
  char *filename = "env1_cpy.ppm";
  //open file for output
  f = fopen(filename, "wb");
  if (!f) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    return 0;
  }

  //write the header file
  //image format
  fprintf(f, "P6\n");

  //comments
  fprintf(f, "# Created by %s\n","SAYANTAN DATTA");

  //image size
  fprintf(f, "%d %d\n", width, height);

  // rgb component depth
  fprintf(f, "%d\n", 255);

  unsigned char *pixelRawValues = new unsigned char[width * height * 3];

  for (uint32_t i = 0; i < width * height; i++) {
    pixelRawValues[i * 3] = (image[i].x > 1 ? 1.0 : image[i].x) * 255;
    pixelRawValues[i * 3 + 1] = (image[i].y > 1 ? 1.0 : image[i].y) * 255;
    pixelRawValues[i * 3 + 2] = (image[i].z > 1 ? 1.0 : image[i].z) * 255;
  }
  // pixel data
  fwrite(pixelRawValues, 3 * width, height, f);
  fclose(f);

  return 1;
}
