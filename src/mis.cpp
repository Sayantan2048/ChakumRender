#include "lightSources.h"
#include "domainSampler.h"
#include <vector>
#include <iostream>
#include "mis.h"

MIS *lightSampler;

enum PdfType {SoAnISSP = 0, SuArISSP, CoIS, PhBrIS};

struct Pdf {
  double capAreaInv; // Solid Angle Sampling.
  double cosThetaMax; // Solid Angle Sampling
  Vec lightDir; // Solid Angle Sampling
  SphereSource *sS;
  double e;
  Vec x;
  Vec n;
  Vec wr;
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
    else if (t == CoIS) {
      return sample.dot(n) * capAreaInv;
    }
    else if (t == PhBrIS) {
      double cosine = sample.dot(wr);
      cosine = cosine > 0 ? cosine : 0;
      return capAreaInv * pow(cosine, e);
    }
    else {
      std::cerr<<"No Such PDF!!\n";
      return 0;
    }
  }
};

uint32_t MIS::getSamples(Vec x, Vec n, Vec wo_ref, Vec wo, double e, Vec *samples, double *weight) {
    std::vector<Pdf> pdf;

    uint32_t offset = 0;

    for (uint32_t i = 0; i < nSS && nSoAnIS > 0; i++, offset += nSoAnIS) {
      Pdf p;
      Vec w = sList[i].p.p - x;
      double sinThetaMax = sList[i].p.r/w.length();
      double capArea = SphericalSampler::getSolidDirectionSamples(w, asin(sinThetaMax), nSoAnIS, samples + offset);
      p.lightDir = w.norm();
      p.capAreaInv = 1.0/capArea;
      p.cosThetaMax = sqrt(1 - sinThetaMax * sinThetaMax);
      p.x = x;
      p.t = SoAnISSP;
      p.N = nSoAnIS;
      pdf.push_back(p);
    }

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
      SphericalSampler::getClassicPhongDirectionSamples(n, wo_ref, e, nPhBrIS, samples + offset);
      p.capAreaInv = (e + 1.0)/(2 * PI);
      p.wr = wo_ref;
      p.x = x;
      p.e = e;
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
      weight[i] = 1.0/denom;
    }

    return nSamples;
}

void loadLightSampler() {
  uint32_t nSS;
  SphereSource *sPtr = lSource->getSphereSourcePtr(nSS);
  lightSampler = new MIS(0, 0, nSS, sPtr);
}
