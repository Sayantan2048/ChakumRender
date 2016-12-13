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

/*
#include "ltc.inc"

void M_GGX(const float theta, const float alpha, mat3 &M, mat3 &Minv, double &amplitude)
{
	int t = maX(0, miN(size-1, (int)floor(theta / (0.5 * PI) * size)));
	int a = maX(0, miN(size-1, (int)floor(sqrt(alpha) * size)));

	//std::cout<<a<<" "<<t<<"\n";
	M = tabM[a + t*size];
	Minv = tabMinv[a + t * size];
	amplitude = (double)tabAmplitude[a + t*size];
}*/

#include "ltc_ggx.inc"
const int size = 64;
void M_GGX2(const float theta, const float alpha, mat3 &M, mat3 &Minv, double &amplitude)
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
}
