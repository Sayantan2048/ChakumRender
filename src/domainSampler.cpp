#include "mathPrimitives.h"
#include "random.h"
#include "domainSampler.h"
#include <cmath>
#include <cstdio>
#include <ctime>

#define MULTIPLIER 300

Vec* SphericalSampler::sampleVolume;
Vec* SphericalSampler::sampleSurface;
int SphericalSampler::__initSamples = SphericalSampler::initSamples();

int SphericalSampler::initSamples() {
    seedMT(0x191123FB);
    double x, y, z;
    /* Optimize...Use aligned memory??*/
    sampleVolume = new Vec[nSAMPLES * MULTIPLIER];
    sampleSurface = new Vec[nSAMPLES * MULTIPLIER];
    for (int i = 0; i < nSAMPLES * MULTIPLIER;) {
      x = randomMTD(-1.0, 1.0);
      y = randomMTD(-1.0, 1.0);
      z = randomMTD(-1.0, 1.0);

      if ((x*x + y*y + z*z <= 1)) {
	sampleVolume[i] = Vec(x, y, z);
	sampleSurface[i] = Vec(x, y, z).norm();
	i++;
      }
    }
    return 0;
}

double SphericalSampler::getSphericalVolumeSamples(Vec x, int nSamples, Vec *store) {
    seedMT(clock() & 0xFFFFFFFF);
    int offset = randomMTD(0, (MULTIPLIER - 1) * nSAMPLES);
    // 0xFFFFFFF0 gives better alignment and performance!!
    for (int i = offset & 0xFFFFFFF0, j = 0; j < nSamples; i++, j++)
      store[j] = sampleVolume[i] + x;

    return 4.0 * PI /3.0;
}

double SphericalSampler::getSphericalSurfaceSamples(Vec x, int nSamples, Vec *store) {
    seedMT(clock() & 0xFFFFFFFF);
    int offset = randomMTD(0, (MULTIPLIER - 1) * nSAMPLES);
    // 0xFFFFFFF0 gives better alignment and performance!!
    for (int i = offset & 0xFFFFFFF0, j = 0; j < nSamples; i++, j++)
      store[j] = sampleSurface[i] + x;

    return 4 * PI;
}

double SphericalSampler::getHemiVolumeSamples(Vec n, Vec x, int nSamples, Vec *store) {
    seedMT(clock() & 0xFFFFFFFF);
    int offset = randomMTD(0, (MULTIPLIER - 1) * nSAMPLES);
    // 0xFFFFFFF0 gives better alignment and performance!!
    for (int i = offset & 0xFFFFFFF0, j = 0; j < nSamples; i++, j++)
      store[j] = n.dot(sampleVolume[i]) >= 0? x+sampleVolume[i] : x-sampleSurface[i];

    return 2.0 * PI / 3.0;
}

double SphericalSampler::getHemiSurfaceSamples(Vec n, Vec x, int nSamples, Vec *store) {
    seedMT(clock() & 0xFFFFFFFF);
    int offset = randomMTD(0, (MULTIPLIER - 1) * nSAMPLES);
    // 0xFFFFFFF0 gives better alignment and performance!!
    for (int i = offset & 0xFFFFFFF0, j = 0; j < nSamples; i++, j++)
      store[j] = n.dot(sampleSurface[i]) >= 0? x+sampleSurface[i] : x-sampleSurface[i];

    return 2.0 * PI;
}

double SphericalSampler::getHemiSurfaceSamplesTrue(Vec n, Vec x, int nSamples, Vec *store) {
    seedMT(clock() & 0xFFFFFFFF);

    double sx, sy, sz;
    Vec s;
    for (int i = 0; i < nSAMPLES;) {
      sx = randomMTD(-1.0, 1.0);
      sy = randomMTD(-1.0, 1.0);
      sz = randomMTD(-1.0, 1.0);

      if ((sx*sx + sy*sy + sz*sz <= 1)) {
	s = Vec(sx, sy, sz).norm();
	store[i] = n.dot(s) >= 0? x + s : x - s;
	i++;
      }
    }

    return 2.0 * PI;
}

// Importance sampling as per visible surface of light from a point x.
double SphericalSampler::getLightSurfaceSample(Vec c, double r, Vec x, int nSamples, Vec *store) {
    seedMT(clock() & 0xFFFFFFFF);
    int offset = randomMTD(0, (MULTIPLIER - 1) * nSAMPLES);
    Vec dir = (x - c).norm();
    // 0xFFFFFFF0 gives better alignment and performance!!
    for (int i = offset & 0xFFFFFFF0, j = 0; j < nSamples; i++, j++)
      store[j] = dir.dot(sampleSurface[i]) >= 0? c + sampleSurface[i] * r : c - sampleSurface[i] * r;

    return 2 * PI * r * r;
}

#define DEBUG_ARCSS 0
//Solid angle imortance sampling!!
//Sample around w as axis, centered at x with sample points making maximum angle of theta_max with w.
double SphericalSampler::getSolidSurfaceSamples(Vec w, Vec x, double theta_max, int nSamples, Vec *store) {
  seedMT(clock() & 0xFFFFFFFF);
  double eps = 1e-14;

  w.norm();
  double eita = w.x * w.x + w.y * w.y;
  double alpha = w.z * w.z + eita; // Parabola opens upwards since alpha >= 0.

  for (int i = 0; i < nSamples;) {
    double theta = randomMTD(theta_max/1000.0, theta_max);
    double delta1 = cos(theta);

    double beta = -2.0 * delta1 * w.z;
    double gamma = delta1 * delta1 - eita;

    double disc = beta * beta - 4 * alpha * gamma;

    if (disc * disc < eps)  {// discriminant almost zero, degenerate case
      double o = -beta / (2 * alpha);

      double r = sqrt(1 - o*o);

      double m = randomMTD(0.0, r);
      double n1 = sqrt (r * r - m * m);
      double n2 = -n1;

#if DEBUG_ARCSS
      fprintf(stdout, "C1:%f %f %f %f %f %f\n", m , n1, o, (Vec(m, n1, o).norm()).dot(w), delta1, theta);
      fprintf(stdout, "C2:%f %f %f %f %f %f\n", m , n2, o, (Vec(m, n2, o).norm()).dot(w), delta1, theta);
#endif
      store[i] = x + Vec(m, n1, o).norm();
      store[i + 1] = x + Vec(m, n2, o).norm();

      i += 2;
    }
    else if (disc > 0) {
      double o1, o2;
      if (beta >= 0) { // This is more numerically stable than directly computing two roots.
        double temp = -(beta + sqrt(disc));
        o1 = temp / (2.0 * alpha);
        o2 = (2 * gamma) / temp;
      }
      else {
        double temp = (-beta) + sqrt(disc);
        o1 = (2 * gamma) / temp;
        o2 = temp / (2 * alpha);
      }

      double eps2 = (o2 - o1) / 1000.0;
      double o = randomMTD(o1 + eps2, o2 - eps2); //if o = o1 or o2 then disc1 is zero. In this case the line is tangential to the circle.

      double sdelta1 = 1 - o * o;
      double sdelta2 = delta1 - o * w.z;
      double alpha1 = w.x * w.x + w.y * w.y;
      double beta1 = -2.0 * w.y * sdelta2;
      double gamma1 = sdelta2 * sdelta2 - w.x * w.x *sdelta1;

      double disc1 = beta1 * beta1 - 4 * alpha1 * gamma1;
      double n1, n2;
      if (disc1 * disc1 < eps) {
        n1 = n2 = -beta1 / (2 * alpha1);
        //if (w.x * w.x < eps * 100000) { // Plane parallel to X axis, i.e if w.x = 0, then disc1 is 0.
	  double m1, m2;
	  m1 = sqrt(sdelta1 - n1 * n1);
	  m2 = -m1;
#if DEBUG_ARCSS
	  fprintf(stdout, "D1:%f %f %f %f %f %f\n", m1 , n1, o, (Vec(m1, n1, o).norm()).dot(w), delta1, theta);
          fprintf(stdout, "D2:%f %f %f %f %f %f\n", m2 , n2, o, (Vec(m2, n2, o).norm()).dot(w), delta1, theta);
#endif
	  store[i] = x + Vec(m1, n1, o).norm();
	  store[i + 1] = x + Vec(m2, n2, o).norm();
	  i += 2;
        //}
        // Other reason for disc1 is 0 because o = o1 or o = o2. But the way o is generated, such condition will never happen!!
      }
      else if (beta1 >= 0) {
        double temp = -(beta1 + sqrt(disc1));
        double m1, m2;
        n1 = temp / (2.0 * alpha1);
        n2 = (2 * gamma1) / temp;

        m1 = (sdelta2 - n1 * w.y) / w.x;
        m2 = (sdelta2 - n2 * w.y) / w.x;
#if DEBUG_ARCSS
        fprintf(stdout, "A1:%f %f %f %f %f %f\n", m1 , n1, o, (Vec(m1, n1, o).norm()).dot(w), delta1, theta);
        fprintf(stdout, "A2:%f %f %f %f %f %f\n", m2 , n2, o, (Vec(m2, n2, o).norm()).dot(w), delta1, theta);
#endif
	store[i] = x + Vec(m1, n1, o).norm();
	store[i + 1] = x + Vec(m2, n2, o).norm();

	i += 2;
      }
      else {
        double temp = (-beta1) + sqrt(disc1);
        double m1, m2;
        n1 = (2 * gamma1) / temp;
        n2 = temp / (2 * alpha1);

        m1 = (sdelta2 - n1 * w.y) / w.x;
        m2 = (sdelta2 - n2 * w.y) / w.x;
#if DEBUG_ARCSS
        fprintf(stdout, "B1:%f %f %f %f %f %f\n", m1, n1, o, (Vec(m1, n1, o).norm()).dot(w), delta1, theta);
        fprintf(stdout, "B2:%f %f %f %f %f %f\n", m2, n2, o, (Vec(m2, n2, o).norm()).dot(w), delta1, theta);
#endif
	store[i] = x + Vec(m1, n1, o).norm();
	store[i + 1] = x + Vec(m2, n2, o).norm();

	i += 2;
      }
    }
#if DEBUG_ARCSS
    fprintf(stdout, "%f %f %f %f\n", disc, w.x, w.y, w.z);
#endif
  }

  // Use Cap-hat theorem to return 1/pdf which is essentially the area of spherical cap subtended by the light.
  return 2.0 * PI * (1 - cos(theta_max));

}

void SphericalSampler::getDistribution(Vec n, Vec x, int nSamples, Vec *samples) {
    int *histogramElevation = new int[181]; // 180 degrees
    int *histogramAzmuth = new int[361]; // 360 degrees
    int *histogramRadii = new int[101]; // 100 uniform samples along radius length.

    double meanElevation = 0, meanAzmuth = 0, meanRadii = 0;

    Vec t = Vec(1,1,1);
    Vec nOrtho = (n%t).norm(); // Vector orthogonal to n. This is our reference for calculating azmuth.
    n.norm();

    int i;
    for (i = 0; i < 181; i++)
      histogramElevation[i] = 0;

    for (i = 0; i < 361; i++)
      histogramAzmuth[i] = 0;

    for (i = 0; i < 101; i++)
      histogramRadii[i] = 0;

    for (int i = 0; i < nSamples; i++) {
	Vec dir = (samples[i] - x);
	double length = dir.length();
	dir.norm();
	double magnitudeProjN = n.dot(dir);

	Vec projN = n * magnitudeProjN;
	Vec projAzmuth = (dir - projN).norm(); // Projection of dir on plane defined by normal n.
	double elevationAngle = acos(magnitudeProjN) * 180.0 / PI;
	double azmuthAngle = acos(projAzmuth.dot(nOrtho)) * 180.0 / PI;
	if ((projAzmuth%nOrtho).dot(n) < 0)
	    azmuthAngle = 360 - azmuthAngle;

	histogramElevation[(int)elevationAngle]++;
	histogramAzmuth[(int)azmuthAngle]++;
	if (length > 1.01)
	  fprintf(stderr, "Not a unit sphere.\n");
	else
	  histogramRadii[int(length * 100.0)]++;

	meanElevation += ((int)elevationAngle) + 0.5;
	meanAzmuth += ((int)azmuthAngle) + 0.5;
	meanRadii += ((int)(length * 100)) + 0.5;
    }

    meanAzmuth /= nSamples;
    meanElevation /= nSamples;
    meanRadii /= nSamples;

    double stdElevation = 0, stdAzmuth = 0, stdRadii = 0;
    //fprintf(stderr, "Elevation Distribution.\n");
    for (i = 0; i < 181; i++) {
	//fprintf(stderr, "%d %d %d\n", i, i+1, histogramElevation[i]);
	stdElevation += (i + 0.5 - meanElevation) * (i + 0.5 - meanElevation) * histogramElevation[i];
    }

    //fprintf(stderr, "Azmuth Distribution.\n");
    for (i = 0; i < 361; i++) {
	//fprintf(stderr, "%d %d %d\n", i, i+1, histogramAzmuth[i]);
	stdAzmuth += (i + 0.5 - meanAzmuth) * (i + 0.5 - meanAzmuth) * histogramAzmuth[i];
    }

    //fprintf(stderr, "Radii Distribution.\n");
    for (i = 0; i < 101; i++) {
	//fprintf(stderr, "%d %d %d\n", i, i+1, histogramRadii[i]);
	stdRadii += (i + 0.5 - meanRadii) * (i + 0.5 - meanRadii) * histogramRadii[i] / 10000.0;
    }

    stdAzmuth = sqrt(stdAzmuth/nSamples);
    stdElevation = sqrt(stdElevation/nSamples);
    stdRadii = sqrt(stdRadii/nSamples);

    fprintf(stderr, "Mean Elevation=%f, Mean Azmuth=%f, Mean Radii=%f\n", meanElevation, meanAzmuth, meanRadii / 100.0);
    fprintf(stderr, "StdDev Elevation=%f, StdDev Azmuth=%f, StdDev Radii=%f\n", stdElevation, stdAzmuth, stdRadii);
}
