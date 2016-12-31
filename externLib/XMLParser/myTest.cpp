#include "tinyxml2.h"
#include <cstdio>
#include <iostream>
#include <vector>

using namespace tinyxml2;
struct SP_Vec {
  double x, y, z;
  SP_Vec():x(0.0),y(0.0),z(0.0) {}
  void show() {
    std::cout<<"("<<x<<", "<<y<<", "<<z<<")";
  }
};

struct SP_TR {
  SP_Vec v1;
  SP_Vec v2;
  SP_Vec v3;
  SP_Vec col;
  double radiance;
  SP_TR():radiance(0.0) {}
};

struct SP_ROTSEQ {
  std::vector<SP_Vec> rotAxisSeq;
  std::vector<double> rotDegSeq;
};

struct SP_MESH {
  std::vector<SP_TR> tList;
  SP_Vec translate;
  SP_ROTSEQ rotations;
};

class SceneParser {
  tinyxml2::XMLDocument doc;
  bool parseVec3(const XMLElement *parent, const char *attribName, SP_Vec &v);
  bool parseDouble(const XMLElement *parent, const char *attribName, double &x);
  bool parseUint32(const XMLElement *parent, const char *attribName, uint32_t &x);
  bool parseColor(const XMLElement *parent, SP_Vec &col);
  uint32_t countSibling(XMLNode *parent, const char *attribName);
  bool parsePointSource(const XMLElement *source, SP_Vec &pos, SP_Vec &col, double &power);
  bool parseSphereSource(const XMLElement *source, double &radius, SP_Vec &pos, SP_Vec &col, double &radiance);
  bool parseTriangle(const XMLElement *ptr, SP_TR &t);
  bool parseRotSeq(const XMLElement *parent, SP_ROTSEQ &rotations);
public:
  SceneParser(const char *fName, bool &status);
  bool parseScreenSize(uint32_t &width, uint32_t &height, bool verbose);
  bool parseCamera(SP_Vec &pos, SP_Vec &dir, double &fovx, double &fovy, bool verbose);
  uint32_t findPointSources(bool verbose);
  uint32_t findSphereSources(bool verbose);
  uint32_t findTriSources(bool verbose);
  uint32_t findMeshSources(bool verbose);
  bool parsePointSources(SP_Vec *pos, SP_Vec *col, double *power, bool verbose);
  bool parseSphereSources(double *radius, SP_Vec *pos, SP_Vec *col, double *radiance, bool verbose);
  bool parseMesh(XMLElement *parent, SP_MESH &mesh, bool verbose);
  bool parseTriSources(SP_MESH *sources, bool verbose);
  bool parseMeshSources(SP_MESH *sources, bool verbose);
};

bool SceneParser::parseVec3(const XMLElement *ptr, const char *attribName, SP_Vec &v) {
  if (!ptr->FirstChildElement(attribName) || !ptr->FirstChildElement(attribName)->FirstChild())
    return false;
  
  const XMLText* textNode = ptr->FirstChildElement(attribName)->FirstChild()->ToText();
  if (sscanf(textNode->Value(), "%lf %lf %lf", &v.x,  &v.y,  &v.z) != 3)
    return false;
  
  return true;
}

bool SceneParser::parseDouble(const XMLElement *ptr, const char *attribName, double &x) {
  if (!ptr->FirstChildElement(attribName) || !ptr->FirstChildElement(attribName)->FirstChild())
    return false;
  
  const XMLText* textNode = ptr->FirstChildElement(attribName)->FirstChild()->ToText();
  if (sscanf(textNode->Value(), "%lf", &x) != 1)
    return false;
  
  return true;
}

bool SceneParser::parseUint32(const XMLElement *parent, const char *attribName, uint32_t &x) {
  if (!parent->FirstChildElement(attribName) || !parent->FirstChildElement(attribName)->FirstChild())
    return false;
  
  const XMLText* textNode = parent->FirstChildElement(attribName)->FirstChild()->ToText();
  if (sscanf(textNode->Value(), "%u", &x) != 1)
    return false;
  
  return true;
}

uint32_t SceneParser::countSibling(XMLNode *parent, const char *attribName) {
  XMLElement *iter = parent->FirstChildElement(attribName);
  if (!iter) return 0;
  uint32_t count = 1;
  
  while ((iter = iter->NextSiblingElement(attribName)))
    count++;
     
  return count;
}

bool SceneParser::parseColor(const XMLElement *parent, SP_Vec &col) {
  if (!parseVec3(parent, "COLOR", col)) {
    double temperature = 0;
    if (!parseDouble(parent, "TEMPERATURE", temperature))
      return false;
    else {
      col.x = temperature;
      col.y = temperature;
      col.z = temperature;
    }
  }
  
  return true;
}

bool SceneParser::parsePointSource(const XMLElement *source, SP_Vec &pos, SP_Vec &col, double &power) {
  if (!parseVec3(source, "POSITION", pos)) {
    std::cerr<<"No Point Source Position Information."<<std::endl;
    return false;
  }
  else if (!parseColor(source, col)) {
    std::cerr<<"No Point Source Color or Temperature Information."<<std::endl;
    return false;
  }
  else if (!parseDouble(source, "POWER", power)) {
    std::cerr<<"No Point Source Power Information."<<std::endl;
    return false;
  }
    
  return true;
}

bool SceneParser::parseSphereSource(const XMLElement *source, double &radius, SP_Vec &pos, SP_Vec &col, double &radiance) {
  if (!parseDouble(source, "RADIUS", radius)) {
    std::cerr<<"No Sphere Source Radius Information."<<std::endl;
    return false;
  }
  else if (!parseVec3(source, "POSITION", pos)) {
    std::cerr<<"No Sphere Source Position Information."<<std::endl;
    return false;
  }
  else if (!parseColor(source, col)) {
    std::cerr<<"No Sphere Source Color or Temperature Information."<<std::endl;
    return false;
  }
  else if (!parseDouble(source, "RADIANCE", radiance)) {
    double power = 0;
    if (!parseDouble(source, "POWER", power)) {
      std::cerr<<"No Sphere Source Power/Radiance Information."<<std::endl;
      return false;
    }
    else {
      radiance = power / (4.0 * 3.1415927 * 3.1415927 * radius * radius);
    }
  }
    
  return true;
}

bool SceneParser::parseTriangle(const XMLElement *ptr, SP_TR &t) {
  if (!parseVec3(ptr, "V1", t.v1)) {
    std::cerr<<"Vertex V1 couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseVec3(ptr, "V2", t.v2)) {
    std::cerr<<"Vertex V2 couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseVec3(ptr, "V3", t.v3)) {
    std::cerr<<"Vertex V3 couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseColor(ptr, t.col)) {
    std::cerr<<"No Triangle Object Color or Temperature Information."<<std::endl;
    return false;
  }
  if (ptr->FirstChildElement("RADIANCE"))
    if (!parseDouble(ptr, "RADIANCE", t.radiance)) {
      std::cerr<<"No Triangle Object Radiance Information."<<std::endl;
      return false;
    }
  
  return true;  
}

bool SceneParser::parseRotSeq(const XMLElement *parent, SP_ROTSEQ &rot) {
  const XMLElement *ptr = parent->FirstChildElement("ROTSEQ");
  if (!ptr) return false;
  
  XMLElement *iter = ((XMLElement *)ptr)->FirstChildElement("R");
  if (!iter) {
    std::cerr<<"No rotaions specified."<<std::endl;
    return false;
  }
  
  uint32_t count = 0;
  do {
    SP_Vec axis;
    double deg;
    if (!parseVec3(iter, "AXIS", axis)) {
      std::cerr<<"No axis information in rotaion: "<<count<<std::endl;
      return false;
    }
    if (!parseDouble(iter, "ROT", deg)) {
      std::cerr<<"No rotaion magnitude information in rotaion: "<<count<<std::endl;
      return false;  
    }
    rot.rotAxisSeq.push_back(axis);
    rot.rotDegSeq.push_back(deg);
   
    count++;
  } while ((iter = iter->NextSiblingElement("R")));
  
  if (rot.rotAxisSeq.size() != rot.rotDegSeq.size()) {
    std::cerr<<"Error parsing rotaion seqence."<<std::endl;
    return false;
  }
  
  return true;  
}

SceneParser::SceneParser(const char *fName, bool &status) {
  doc.LoadFile(fName);
  status = true;

  if (doc.Error()) {
    std::cerr<<"XML Document Syntax Error."<<std::endl;
    status = false;
  }
}

bool SceneParser::parseScreenSize(uint32_t &width, uint32_t &height, bool verbose = false) {
  XMLElement *screen = doc.FirstChildElement("SCREEN");
  if (!screen) {
    std::cerr<<"No Screen Information."<<std::endl;
    return false;
  }
  else if (!parseUint32(screen, "WIDTH", width)) {
    std::cerr<<"No Screen Width Information."<<std::endl;
    return false;
  }
  else if (!parseUint32(screen, "HEIGHT", height)) {
    std::cerr<<"No Screen Height Information."<<std::endl;
    return false;
  }
  
  if (verbose) {
    std::cout<<"Screen Size:("<<width<<", "<<height<<")"<<std::endl;
  }
  
  return true;
}

bool SceneParser::parseCamera(SP_Vec &pos, SP_Vec &dir, double &fovx, double &fovy, bool verbose = false) {
  
  XMLElement *camera = doc.FirstChildElement("CAMERA");
  if (!camera) {
    std::cerr<<"No Camera Information."<<std::endl;
    return false;
  }
  else if (!parseVec3(camera, "POSITION", pos)) {
    std::cerr<<"No Camera Position Information."<<std::endl;
    return false;
  }
  else if (!parseVec3(camera, "DIRECTION", dir)) {
    std::cerr<<"No Camera Direction Information."<<std::endl;
    return false;
  }
  else if (!parseDouble(camera, "FOVX", fovx)) {
    std::cerr<<"No Camera Fov X Information."<<std::endl;
    return false;
  }
  else if (!parseDouble(camera, "FOVY", fovy)) {
    std::cerr<<"No Camera Fov Y Information."<<std::endl;
    return false;
  }
  
  if (verbose) {
    std::cout<<"Camera Position:"; pos.show(); std::cout<<std::endl;
    std::cout<<"Camera Direction:"; dir.show(); std::cout<<std::endl;
    std::cout<<"Camera FOV X and Y:("<<fovx<<", "<<fovy<<")"<<std::endl;
  }
  
  return true;
}

uint32_t SceneParser::findPointSources(bool verbose = false) {
  uint32_t count = 0;
  count = countSibling(&doc, "POINTSOURCE");
   
  if (verbose)
    std::cout<<"Point Sources found: "<<count<<std::endl;
  
  return count;
}

uint32_t SceneParser::findSphereSources(bool verbose = false) {
  uint32_t count = 0;
  count = countSibling(&doc, "SPHERESOURCE");
   
  if (verbose)
    std::cout<<"Sphere Sources found: "<<count<<std::endl;
  
  return count;
}

uint32_t SceneParser::findTriSources(bool verbose = false) {
  uint32_t count = 0;
  count = countSibling(&doc, "TRISOURCE");
   
  if (verbose)
    std::cout<<"Triangle Sources found: "<<count<<std::endl;
  
  return count;
}

uint32_t SceneParser::findMeshSources(bool verbose = false) {
  uint32_t count = 0;
  count = countSibling(&doc, "MESHSOURCE");
   
  if (verbose)
    std::cout<<"MESH Sources found: "<<count<<std::endl;
  
  return count;
}


bool SceneParser::parsePointSources(SP_Vec *pos, SP_Vec *col, double *power, bool verbose = false) {
  XMLElement *iter = doc.FirstChildElement("POINTSOURCE");
  if (!iter) return false;
   
  uint32_t count = 0;
  do {
    if (!parsePointSource(iter, pos[count], col[count], power[count])) {
      std::cerr<<"Point Source: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    if (verbose) {
      std::cout<<"Point Source: "<<count<<" Position:";pos[count].show(); std::cout<<std::endl;
      std::cout<<"Point Source: "<<count<<" Color:";col[count].show(); std::cout<<std::endl;
      std::cout<<"Point Source: "<<count<<" Power: "<<power[count]<<std::endl;
    }
    count++;
  } while ((iter = iter->NextSiblingElement("POINTSOURCE")));
  
  return true;
}

bool SceneParser::parseSphereSources(double *radius, SP_Vec *pos, SP_Vec *col, double *radiance, bool verbose = false) {
  XMLElement *iter = doc.FirstChildElement("SPHERESOURCE");
  if (!iter) return false;
   
  uint32_t count = 0;
  do {
    if (!parseSphereSource(iter, radius[count], pos[count], col[count], radiance[count])) {
      std::cerr<<"Sphere Source: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    if (verbose) {
      std::cout<<"Sphere Source: "<<count<<" Radius: "<<radius[count]<<std::endl;
      std::cout<<"Sphere Source: "<<count<<" Position:";pos[count].show(); std::cout<<std::endl;
      std::cout<<"Sphere Source: "<<count<<" Color:";col[count].show(); std::cout<<std::endl;
      std::cout<<"Sphere Source: "<<count<<" Radiance: "<<radiance[count]<<std::endl;
    }
    count++;
  } while ((iter = iter->NextSiblingElement("SPHERESOURCE")));
  
  return true;
}

bool SceneParser::parseMesh(XMLElement *parent, SP_MESH &mesh, bool verbose = false) {
   XMLElement *iter = parent->FirstChildElement("TR");
   if (!iter) {
     std::cerr<<"No Triangle fins found in mesh."<<std::endl;
     return false;
   }
   uint32_t count = 0;
   do {
    SP_TR t; 
    if (!parseTriangle(iter, t)) {
      std::cerr<<"Triangle: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    mesh.tList.push_back(t);
    if (verbose) {
      std::cout<<"Triangle: "<<count<<" Vertices:"; mesh.tList[count].v1.show(); mesh.tList[count].v2.show(); mesh.tList[count].v3.show(); std::cout<<std::endl;
      std::cout<<"Triangle: "<<count<<" Color:"; mesh.tList[count].col.show(); std::cout<<std::endl;
      std::cout<<"Triangle: "<<count<<" Radiance:"<<mesh.tList[count].radiance<<std::endl;
    }
    count++;
  } while ((iter = iter->NextSiblingElement("TR")));
  
  if (parent->FirstChildElement("TRANSLATE"))
    if (!parseVec3(parent, "TRANSLATE", mesh.translate)) {
      std::cerr<<"No Mesh Translation Information."<<std::endl;
      return false;
    }
    
  if (parent->FirstChildElement("ROTSEQ"))
    if (!parseRotSeq(parent, mesh.rotations)) {
      std::cerr<<"No Mesh Rotaion Information."<<std::endl;
      return false;
    }
    
  if (verbose) {
    std::cout<<"Mesh Translation: "; mesh.translate.show(); std::cout<<std::endl;
    for (uint32_t i = 0; i < mesh.rotations.rotAxisSeq.size(); i++) {
      std::cout<<"Mesh Rotation "<<i<<" Axis: ";mesh.rotations.rotAxisSeq[i].show();std::cout<<" Degrees:"<<mesh.rotations.rotDegSeq[i]<<std::endl;
    }
  }
  
  return true;
}

bool SceneParser::parseTriSources(SP_MESH *sources, bool verbose = false) {
  XMLElement *iter = doc.FirstChildElement("TRISOURCE");
  
  if (!iter) return false;
  
  uint32_t count = 0;
  do {
    if (!parseMesh(iter, sources[count], verbose)) {
      std::cerr<<"Triangle Source: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    if (sources[count].tList.size() > 1) {
      std::cerr<<"Triangle Source: "<<count<<" has more than one triangle fins. Use Mesh Source instead."<<std::endl;
    }
    count++;
  } while ((iter = iter->NextSiblingElement("TRISOURCE")));
  
  return true;
}

bool SceneParser::parseMeshSources(SP_MESH *sources, bool verbose = false) {
  XMLElement *iter = doc.FirstChildElement("MESHSOURCE");
  
  if (!iter) return false;
  
  uint32_t count = 0;
  do {
    if (!parseMesh(iter, sources[count], verbose)) {
      std::cerr<<"MESH Source: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    count++;
  } while ((iter = iter->NextSiblingElement("MESHSOURCE")));
  
  return true;
}

int main() {
  bool status;
  SceneParser s("myxml.xml", status);
  if (!status)
    return 0;
  uint32_t width, height;
  s.parseScreenSize(width, height, true);
  SP_Vec cpos, cdir;
  double fovx, fovy;
  s.parseCamera(cpos, cdir, fovx, fovy, true);
  uint32_t count = s.findSphereSources(true);
  
  double *radius = new double[count];
  SP_Vec *pos = new SP_Vec[count];
  SP_Vec *col = new SP_Vec[count];
  double *power = new double[count];
  
  s.parseSphereSources(radius, pos, col, power, true);
  
  delete[] radius;
  delete[] pos;
  delete[] col;
  delete[] power;
  
  count = s.findPointSources(true);
  pos = new SP_Vec[count];
  col = new SP_Vec[count];
  power = new double[count];
  
  s.parsePointSources(pos, col, power, true);
 
  count = s.findMeshSources(true);
  SP_MESH *source = new SP_MESH[count];
  
  s.parseMeshSources(source, true);
  
  count = s.findTriSources(true);
  source = new SP_MESH[count];
  
  s.parseTriSources(source, true);
    
  return 0;
}