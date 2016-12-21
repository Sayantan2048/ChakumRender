#include "mathPrimitives.h"
#include <cmath>
#include <iostream>

struct mat33
{
    operator mat3() const
    {
        return mat3(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
    }

    double m[9];
};


#include "ltc.inc"

void M_GGX(const float theta, const float alpha, mat3 &M, mat3 &Minv, double &amplitude) {
  int t = maX(0, miN(size-1, (int)floor(theta / (0.5 * PI) * size)));
  int a = maX(0, miN(size-1, (int)floor(sqrt(alpha) * size)));
  int tPlusOne = miN(size - 1, t + 1);
    
  M = tabM[a + t*size];
  Minv = tabMinv[a + t * size];
  amplitude = (double)tabAmplitude[a + t*size];
  
  if (t == tPlusOne) return;
  
  // Perform linear interpolation, http://mathworld.wolfram.com/Two-PointForm.html 
  float mult = theta / (0.5 * PI) * size - (float)t;
  
  mat3 Mnext = tabM[a + tPlusOne*size];
  mat3 MinvNext =  tabMinv[a + tPlusOne * size];
  double ampNext = (double)tabAmplitude[a + tPlusOne*size];

  M.a[0] = M.a[0] + (Mnext.a[0] - M.a[0]) * mult;
  M.a[1] = M.a[1] + (Mnext.a[1] - M.a[1]) * mult;
  M.a[2] = M.a[2] + (Mnext.a[2] - M.a[2]) * mult;
  
  M.a[3] = M.a[3] + (Mnext.a[3] - M.a[3]) * mult;
  M.a[4] = M.a[4] + (Mnext.a[4] - M.a[4]) * mult;
  M.a[5] = M.a[5] + (Mnext.a[5] - M.a[5]) * mult;

  M.a[6] = M.a[6] + (Mnext.a[6] - M.a[6]) * mult;
  M.a[7] = M.a[7] + (Mnext.a[7] - M.a[7]) * mult;	
  M.a[8] = M.a[8] + (Mnext.a[8] - M.a[8]) * mult;

  Minv.a[0] = Minv.a[0] + (MinvNext.a[0] - Minv.a[0]) * mult;
  Minv.a[1] = Minv.a[1] + (MinvNext.a[1] - Minv.a[1]) * mult;
  Minv.a[2] = Minv.a[2] + (MinvNext.a[2] - Minv.a[2]) * mult;
  
  Minv.a[3] = Minv.a[3] + (MinvNext.a[3] - Minv.a[3]) * mult;
  Minv.a[4] = Minv.a[4] + (MinvNext.a[4] - Minv.a[4]) * mult;
  Minv.a[5] = Minv.a[5] + (MinvNext.a[5] - Minv.a[5]) * mult;

  Minv.a[6] = Minv.a[6] + (MinvNext.a[6] - Minv.a[6]) * mult;
  Minv.a[7] = Minv.a[7] + (MinvNext.a[7] - Minv.a[7]) * mult;	
  Minv.a[8] = Minv.a[8] + (MinvNext.a[8] - Minv.a[8]) * mult;
  
  amplitude = amplitude + (ampNext - amplitude) * mult;
}
/*
#include "ltc_ggx.inc"
const int size = 64;
void M_GGX(const float theta, const float alpha, mat3 &M, mat3 &Minv, double &amplitude)
{
	int t = maX(0, miN(size-1, (int)floor(theta / (0.5 * PI) * size)));
	int a = maX(0, miN(size-1, (int)floor(sqrt(alpha) * size)));

	//std::cout<<a<<" "<<t<<"\n";
	const double *mat = &g_ltc_mat[(a + t*size) * 4];
	Minv = mat3(Vec(1, 0, mat[1]),
		    Vec(0, mat[2], 0),
		    Vec(mat[3], 0, mat[0]));
	M = Minv.inv();

	amplitude = g_ltc_mag[a + t*size];
}*/
