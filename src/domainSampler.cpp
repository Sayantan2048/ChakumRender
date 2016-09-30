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
