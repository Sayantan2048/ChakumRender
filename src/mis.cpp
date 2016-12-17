#include "lightSources.h"
#include "domainSampler.h"
#include <vector>
#include <iostream>
#include "mis.h"

MIS *lightSampler;

enum PdfType {SoAnISSP = 0, SuArISSP, SuArISME, CoIS, PhBrIS};

#define USE_SPHERE_SOURCE_LIGHT_SAMPLE 0
struct Pdf {
  double capAreaInv; // Solid Angle Sampling.
  double cosThetaMax; // Solid Angle Sampling
  Vec lightDir; // Solid Angle Sampling
#if USE_SPHERE_SOURCE_LIGHT_SAMPLE
  SphereSource *sS;
#endif
  MeshLight *mS;
  MaterialType *mPtr;
  Vec x;
  Vec n;
  Vec wo_ref;
  PdfType t;
  double N; // No. of samples in this pdf type
  double pdfVal(Vec sample) {
    if (t == SoAnISSP) {
      // if the sample lies withing thetaMax of lightDir, then return 1/capArea;
      if (sample.dot(lightDir) >= cosThetaMax)
	return capAreaInv;
      else
	return 0.;
    }
#if USE_SPHERE_SOURCE_LIGHT_SAMPLE
    else if (t == SuArISSP) {
      Ray r(x, sample);
      double d;
      if ((d = sS->p.intersect(r)) < INF) {
	double cosine1 = sample.dot((x + sample * d - sS->p.p).norm());
	cosine1 = cosine1 < 0 ? -cosine1 : 1e-5;
	return (capAreaInv * d * d)/cosine1;
      }
      else
	return 0;
    }
#endif
    else if (t == SuArISME) {
      Ray r(x, sample);
      double d;
      uint32_t dummy;
      Vec n;
      mS->intersect(r, dummy, d, n);
      if (d < INF) {
	double cosine1 = sample.dot(n);
	cosine1 = cosine1 < 0 ? -cosine1 : cosine1;
	return (capAreaInv * d * d) / cosine1;
      }
      else
	return 0;
    }
    else if (t == CoIS) {
      return sample.dot(n) * capAreaInv;
    }
    else if (t == PhBrIS) {

      return mPtr->pdfEval(n, wo_ref, x, sample);
    }
    else {
      std::cerr<<"No Such PDF!!\n";
      return 0;
    }
  }
};

uint32_t MIS::getSamples(const Vec &x, const Vec &n, const Vec &wo_ref, const Vec &wo, BasePrimitive* const bPtr, Vec *samples, double *weight) {
    std::vector<Pdf> pdf;

    uint32_t offset = 0;

    for (uint32_t i = 0; i < nSS && nSoAnIS > 0; i++, offset += nSoAnIS) {
      Pdf p;
      Vec w = sList[i].p.p - x;
      double sinThetaMax = sList[i].p.r/w.length();
      double capArea = SphericalSampler::getSolidDirectionSamples(w.norm(), asin(sinThetaMax), nSoAnIS, samples + offset);
      p.lightDir = w;
      p.capAreaInv = 1.0/capArea;
      p.cosThetaMax = sqrt(1 - sinThetaMax * sinThetaMax);
      p.x = x;
      p.t = SoAnISSP;
      p.N = nSoAnIS;
      pdf.push_back(p);
    }
#if USE_SPHERE_SOURCE_LIGHT_SAMPLE
    for (uint32_t i = 0; i < nSS && nSuArIS > 0; i++, offset += nSuArIS) {
      Pdf p;
      double capArea = SphericalSampler::getLightDirectionSamples(sList[i].p.p, sList[i].p.r, x, nSuArIS, samples + offset);
      p.sS = &sList[i];
      p.capAreaInv = 1.0/capArea;
      p.x = x;
      p.t = SuArISSP;
      p.N = nSuArIS;
      pdf.push_back(p);
    }
#endif
    for (uint32_t i = 0; i < nMS && nSuArIS > 0; i++, offset += nSuArIS) {
      Pdf p;
      double capArea = mList[i].getSamples(nSuArIS, samples + offset, NULL);
      for (uint32_t j = offset; j < offset + nSuArIS; j++)
	samples[j] = (samples[j] - x).norm();
      p.mS = &mList[i];
      p.capAreaInv = 1.0/capArea;
      p.x = x;
      p.t = SuArISME;
      p.N = nSuArIS;
      pdf.push_back(p);
    }

    if (nCoIS > 0) {
      Pdf p;
      double capArea = SphericalSampler::getCosineDirectionSamples(n, nCoIS, samples + offset);
      p.capAreaInv = 1.0/capArea;
      p.x = x;
      p.n = n;
      p.t = CoIS;
      p.N = nCoIS;
      pdf.push_back(p);
      offset += nCoIS;
    }

    if (nPhBrIS > 0) {
      Pdf p;
      bPtr->m.getBrdfDirectionSamples(n, wo_ref, wo, samples + offset, 0, nPhBrIS);
      p.wo_ref = wo_ref;
      p.n = n;
      p.x = wo;
      p.mPtr = &(bPtr->m);
      p.t = PhBrIS;
      p.N = nPhBrIS;
      pdf.push_back(p);
      offset += nPhBrIS;
    }

/*
    offset = 0;
    // After all pdfs are initialized and all samples are initialized.
    for (uint32_t i = 0; i < pdf.size(); i++) {
      for (uint32_t j = 0; j < pdf[i].N; j++) {
	Vec sample = samples[j + offset];
	double numerator = pdf[i].N * pdf[i].pdfVal(sample);
	double denom = 0;
	for (uint32_t pdfIter = 0; pdfIter < pdf.size(); pdfIter++) {
	  denom += pdf[pdfIter].N * pdf[pdfIter].pdfVal(sample);
	}
	weight[j + offset] = numerator / (denom * pdf[i].N * pdf[i].pdfVal(sample));
      }
      offset += pdf[i].N;
    }*/

    uint32_t allSamples = 0;
    for (uint32_t i = 0; i < pdf.size(); i++)
      allSamples += pdf[i].N;

    for (uint32_t i = 0; i < allSamples; i++) {
      double denom = 0;
      Vec sample = samples[i];
      for (uint32_t pdfIter = 0; pdfIter < pdf.size(); pdfIter++) {
	  denom += pdf[pdfIter].N * pdf[pdfIter].pdfVal(sample);
      }
      weight[i] = 1.0/(denom + 1e-20);
    }

    return nSamples;
}

void loadLightSampler() {
  uint32_t nSS;
  uint32_t nMS;
  SphereSource *sPtr = lSource->getSphereSourcePtr(nSS);
  MeshLight *mPtr = lSource->getMeshLightPtr(nMS);
  lightSampler = new MIS(nSS, sPtr, nMS, mPtr);
}
