#if !defined (shader_h__)
#define shader_h__

extern Vec shade(const Ray &r);
// Return 0 if shadowed else 1.
extern double shadow(const Ray &shadowRay, double distanceLightSource);

#endif