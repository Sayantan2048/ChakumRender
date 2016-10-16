class MIS {
  uint32_t nSamples;
  static const uint32_t nLightSamples = 10; // No. of samples per volume/area/mesh light sources
  static const uint32_t nCosineSamples = 5; // Cosine weighted samples.
  static const uint32_t nBrdfSamples = 10; // No. of samples per brdf.
  uint32_t nVAMLights; // No. of volume/mesh/area lights.
  uint32_t nPointLight; // No. of point lights.

  // Return weight of the samples
  double * getSamples(Vec n, Vec wo, uint32_t &nSamples, double *store) {
  }

}
