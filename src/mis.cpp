#include "lightSources.h"
#include "domainSampler.h"
#include <vector>
#include <iostream>
#include "mis.h"

#define USE_SPHERE_SOURCE_LIGHT_SAMPLE 0

MIS *lightSampler;

enum PdfType {SoAnISSP = 0,
#if USE_SPHERE_SOURCE_LIGHT_SAMPLE
  SuArISSP,
#endif
  SuArISME, SuArISTR, EnMpIS, CoIS, BrIS};


struct Pdf {
  double capAreaInv; // Solid Angle Sampling.
  double cosThetaMax; // Solid Angle Sampling
  Vec lightDir; // Solid Angle Sampling
#if USE_SPHERE_SOURCE_LIGHT_SAMPLE
  SphereSource *sS;
#endif
  MeshLight *mS;
  TriLight *tS;
  EnvSource *eS;
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
      if (mS->intersect(r, dummy, d, n)) {
	double cosine1 = sample.dot(n);
	cosine1 = cosine1 < 0 ? -cosine1 : cosine1;
	return (capAreaInv * d * d) / cosine1;
      }
      else
	return 0;
    }
    else if (t == SuArISTR) {
      Ray r(x, sample);
      double d;
      Vec n;
      if ((d = tS->t.intersect(r, n)) < INF) {
	double cosine1 = sample.dot(n);
	cosine1 = cosine1 < 0 ? -cosine1 : cosine1;
	return (capAreaInv * d * d) / cosine1;
      }
      else
	return 0;
    }
    else if (t == CoIS)
      return sample.dot(n) * capAreaInv;

    else if (t == BrIS)
      return mPtr->pdfEval(n, wo_ref, x, sample);
    else if (t == EnMpIS) {
      Vec L = eS->getRadiance(sample);
      return to_greyScale(L) * capAreaInv;
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

    for (uint32_t i = 0; i < nTS && nSuArIS > 0; i++, offset += nSuArIS) {
      Pdf p;
      SphericalSampler::getTriLightSurfaceSamples(tList[i].t.A, tList[i].t.B, tList[i].t.C, nSuArIS, samples + offset);
      double capArea = tList[i].t.area;
      for (uint32_t j = offset; j < offset + nSuArIS; j++)
	samples[j] = (samples[j] - x).norm();
      p.tS = &tList[i];
      p.capAreaInv = 1.0/capArea;
      p.x = x;
      p.t = SuArISTR;
      p.N = nSuArIS;
      pdf.push_back(p);
    }

    if(eS != 0 && nEnMpIS > 0) {
      Pdf p;
      Random random(clock() & 0xFFFFFFFF);
      uint32_t e_offset = random.randomMTD(0, (eS[0].multiplier -1) * nEnMpIS);
      for (uint32_t i = 0; i < nEnMpIS; i++)
	samples[i + offset] = eS[0].impSamples[i + e_offset];
      p.eS = eS;
      p.capAreaInv = 1.0/eS[0].normalizationArea;
      p.t = EnMpIS;
      p.N = nEnMpIS;
      pdf.push_back(p);
      offset += nEnMpIS;
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

    if (nBrIS > 0) {
      Pdf p;
      bPtr->m.getBrdfDirectionSamples(n, wo_ref, wo, samples + offset, 0, nBrIS);
      p.wo_ref = wo_ref;
      p.n = n;
      p.x = wo;
      p.mPtr = &(bPtr->m);
      p.t = BrIS;
      p.N = nBrIS;
      pdf.push_back(p);
      offset += nBrIS;
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

    return allSamples;
}

void loadLightSampler() {
  uint32_t nSS;
  uint32_t nMS;
  uint32_t nTS;
  SphereSource *sPtr = lSource->getSphereSourcePtr(nSS);
  MeshLight *mPtr = lSource->getMeshLightPtr(nMS);
  TriLight *tPtr = lSource->getTriLightPtr(nTS);
  EnvSource *eS = lSource->getEnvSourcePtr(); // No memory allocated, direct pointer returned.
  lightSampler = new MIS(nSS, sPtr, nMS, mPtr, nTS, tPtr, eS);
}
