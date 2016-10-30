#if !defined(materialTypes_h__)
#define materialTypes_h__

enum LightType {NONE = 0, POINT, VOLUME};
class MaterialType {
public:
  double phongExp;
  double specularCoef;
  LightType l;
  Vec radiance;
  // Assumes w and n is normalized!!
  Vec getRadiance(const Vec &n, const Vec &w) {
    // This is true only for sphere volume sources.
    return radiance;
  }
  MaterialType(double pE, double sC, LightType ll, const Vec &rad) {
    phongExp = (pE >= 0) ? pE : -pE;
    specularCoef = (sC * ((sC < 0)? -1.0: 1.0));
    specularCoef = specularCoef <= 1 ? specularCoef : 1;
    l = ll;
    radiance = rad;
  }
  double brdf(Vec n, Vec wo, Vec wi);
};

#endif
