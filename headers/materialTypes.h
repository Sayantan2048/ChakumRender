#if !defined(materialTypes_h__)
#define materialTypes_h__

class MaterialType {
public:
  double phongExp;
  double specularCoef;
  MaterialType(double pE, double sC) {
    phongExp = (pE >= 0) ? pE : -pE;
    specularCoef = (specularCoef * ((specularCoef < 0)? -1.0: 1.0));
    specularCoef = specularCoef <= 1 ? specularCoef : 1;
  }
  double brdf(Vec n, Vec wo, Vec wi);
};

#endif
