#if !defined(sceneParser_h__)
#define sceneParser_h__

#include "tinyxml2.h"
#include <vector>
#include <stdint.h>
#include <cstdio>
#include <iostream>
#include <string>

using namespace tinyxml2;

struct SP_Vec {
  double x, y, z;
  SP_Vec():x(0.0),y(0.0),z(0.0) {}
  void show() {
    std::cout<<"("<<x<<", "<<y<<", "<<z<<")";
  }
};

struct SP_Base {
  double alpha;
  double intIOR;
  double extIOR;
  double specCoef;
  double reflectance;
  std::string brdf;
  SP_Base(): alpha(1.0), intIOR(10000.0), extIOR(1), specCoef(0), reflectance(1.0), brdf("DIFFUSE"){}
};

struct SP_TR: public SP_Base {
  SP_Vec v1;
  SP_Vec v2;
  SP_Vec v3;
  SP_Vec col;
  double radiance;
    
  SP_TR():radiance(0.0) {}
};

struct SP_SPHERE: public SP_Base {
  double radius;
  SP_Vec pos;
  SP_Vec col;
  double radiance;
    
  SP_SPHERE():radiance(0.0) {}
};

struct SP_ROTSEQ {
  std::vector<SP_Vec> rotAxisSeq;
  std::vector<double> rotDegSeq;
};

struct SP_MESH {
  std::vector<SP_TR> tList;
  SP_Vec translate;
  SP_Vec scale;
  SP_ROTSEQ rotations;
};

struct SP_MESH_SPHERE {
  std::vector<SP_SPHERE> sList;
  SP_Vec translate;
};


class SceneParser {
  tinyxml2::XMLDocument doc;
  bool parseVec3(const XMLElement *parent, const char *attribName, SP_Vec &v);
  bool parseDouble(const XMLElement *parent, const char *attribName, double &x);
  bool parseUint32(const XMLElement *parent, const char *attribName, uint32_t &x);
  bool parseString(const XMLElement *parent, const char *attribName, std::string &s);
  bool parseColor(const XMLElement *parent, SP_Vec &col);
  uint32_t countSibling(XMLNode *parent, const char *attribName);
  bool parsePointSource(const XMLElement *source, SP_Vec &pos, SP_Vec &col, double &power);
  bool parseSphereSource(const XMLElement *source, double &radius, SP_Vec &pos, SP_Vec &col, double &radiance);
  bool parseTriangle(const XMLElement *ptr, SP_TR &t, bool parseMaterial);
  bool parseSphere(const XMLElement *ptr, SP_SPHERE &s, bool parseMaterial);
  bool parseRotSeq(const XMLElement *parent, SP_ROTSEQ &rotations);
  bool parseMaterialProperties(const XMLElement *parent, SP_Base &m);
  bool parseMesh(XMLElement *parent, SP_MESH &mesh, bool parseMaterial, bool verbose);
  bool parseMesh_Sphere(XMLElement *parent, SP_MESH_SPHERE &mesh, bool parseMaterial, bool verbose);
public:
  SceneParser(const char *fName, bool &status);
  bool parseScreenSize(uint32_t &width, uint32_t &height, bool verbose);
  bool parseCamera(SP_Vec &pos, SP_Vec &dir, SP_Vec &up, double &fovx, double &fovy, bool verbose);
  uint32_t findPointSources(bool verbose);
  uint32_t findSphereSources(bool verbose);
  uint32_t findTriSources(bool verbose);
  uint32_t findMeshSources(bool verbose);
  uint32_t findMeshObjects(bool verbose);
  uint32_t findSphereObjects(bool verbose); 
  bool parsePointSources(SP_Vec *pos, SP_Vec *col, double *power, bool verbose);
  bool parseSphereSources(double *radius, SP_Vec *pos, SP_Vec *col, double *radiance, bool verbose);
  bool parseTriSources(SP_MESH *sources, bool verbose);
  bool parseMeshSources(SP_MESH *sources, bool verbose);
  bool parseMeshObjects(SP_MESH *objects, bool verbose);
  bool parseSphereObjects(SP_MESH_SPHERE *objects, bool verbose);
};

extern SceneParser *sceneParser;
#endif