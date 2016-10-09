#if !defined(materialTypes_h__)
#define materialTypes_h__

class Lambertian {
  Lambertian();
public:
  static double brdf();
};

class Phong {
  Phong();
  static double e;
public:
  static double brdf(Vec n, Vec wo, Vec wi);
};

enum MaterialType {lambertian = 0, phong};

#endif
