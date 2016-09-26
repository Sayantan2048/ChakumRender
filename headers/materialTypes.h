#if !defined(materialTypes_h__)
#define materialTypes_h__

class Lambertian {
  Lambertian();
public:
  static double brdf();
};

class Diffuse {
  Diffuse();
public:
  static double brdf();
};

enum MaterialType {lambertian = 0, diffuse};

#endif
