#include "sceneParser.h"
#include "mathPrimitives.h"

SceneParser *sceneParser;

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

bool SceneParser::parseString(const XMLElement *parent, const char *attribName, std::string &s) {
  if (!parent->FirstChildElement(attribName) || !parent->FirstChildElement(attribName)->FirstChild())
    return false;
  
  const XMLText* textNode = parent->FirstChildElement(attribName)->FirstChild()->ToText();
  
  s = textNode->Value();
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
      Vec _col = tempToColor(temperature);
      col.x = _col.x;
      col.y = _col.y;
      col.z = _col.z;
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

bool SceneParser::parseMaterialProperties(const XMLElement *parent, SP_Base &m) {
   if (!parseDouble(parent, "ALPHA", m.alpha)) {
    std::cerr<<"Material Property Alpha couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseDouble(parent, "SPECCOEF", m.specCoef)) {
    std::cerr<<"Material Property Specular Coefficient couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseDouble(parent, "REFLECTANCE", m.reflectance)) {
    std::cerr<<"Material Property Reflectance couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseString(parent, "BRDF", m.brdf)) {
    std::cerr<<"Material Property BRDF couldn't be parsed."<<std::endl;
    return false;
  }
  if (parent->FirstChildElement("INTIOR"))
    if (!parseDouble(parent, "INTIOR", m.intIOR)) {
      std::cerr<<"Material Property Internal IOR couldn't be parsed."<<std::endl;
      return false;
    }
  if (parent->FirstChildElement("EXTIOR"))
    if (!parseDouble(parent, "EXTIOR", m.extIOR)) {
      std::cerr<<"Material Property External IOR couldn't be parsed."<<std::endl;
      return false;
    }
    
  return true;  
}

bool SceneParser::parseTriangle(const XMLElement *ptr, SP_TR &t, bool parseMaterial = false) {
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
  else if (parseMaterial && !parseMaterialProperties(ptr, t))
    return false;
  if (ptr->FirstChildElement("RADIANCE"))
    if (!parseDouble(ptr, "RADIANCE", t.radiance)) {
      std::cerr<<"No Triangle Object Radiance Information."<<std::endl;
      return false;
    }
  
  return true;  
}

bool SceneParser::parseSphere(const XMLElement *ptr, SP_SPHERE &s, bool parseMaterial = false) {
  if (!parseDouble(ptr, "RADIUS", s.radius)) {
    std::cerr<<"Sphere radius couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseVec3(ptr, "POSITION", s.pos)) {
    std::cerr<<"Sphere position couldn't be parsed."<<std::endl;
    return false;
  }
  else if (!parseColor(ptr, s.col)) {
    std::cerr<<"No Sphere Color or Temperature Information."<<std::endl;
    return false;
  }
  else if (parseMaterial && !parseMaterialProperties(ptr, s))
    return false;
  
  if (ptr->FirstChildElement("RADIANCE"))
    if (!parseDouble(ptr, "RADIANCE", s.radiance)) {
      std::cerr<<"No Sphere Radiance Information."<<std::endl;
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

bool SceneParser::parseCamera(SP_Vec &pos, SP_Vec &dir, SP_Vec &up, double &fovx, double &fovy, bool verbose = false) {
  
  XMLElement *camera = doc.FirstChildElement("CAMERA");
  if (!camera) {
    std::cerr<<"No Camera Information."<<std::endl;
    return false;
  }
  else if (!parseVec3(camera, "POSITION", pos)) {
    std::cerr<<"No Camera Position Information."<<std::endl;
    return false;
  }
  else if (!parseVec3(camera, "LOOKAT", dir)) {
    std::cerr<<"No Camera Direction Information."<<std::endl;
    return false;
  }
  else if (!parseVec3(camera, "UP", up)) {
    std::cerr<<"No Camera Up Direction Information."<<std::endl;
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
    std::cout<<"Camera Up Direction:"; up.show(); std::cout<<std::endl;
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

uint32_t SceneParser::findMeshObjects(bool verbose = false) {
  uint32_t count = 0;
  count = countSibling(&doc, "MESHOBJECT");
   
  if (verbose)
    std::cout<<"MESH Objects found: "<<count<<std::endl;
  
  return count;
}

uint32_t SceneParser::findSphereObjects(bool verbose = false) {
  uint32_t count = 0;
  count = countSibling(&doc, "SPHEREOBJECT");
   
  if (verbose)
    std::cout<<"Sphere Mesh Objects found: "<<count<<std::endl;
  
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

bool SceneParser::parseMesh(XMLElement *parent, SP_MESH &mesh, bool parseMaterial, bool verbose = false) {
   XMLElement *iter = parent->FirstChildElement("TR");
   if (!iter) {
     std::cerr<<"No Triangle fins found in mesh."<<std::endl;
     return false;
   }
   uint32_t count = 0;
   do {
    SP_TR t; 
    if (!parseTriangle(iter, t, parseMaterial)) {
      std::cerr<<"Triangle: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    mesh.tList.push_back(t);
    if (verbose) {
      std::cout<<"Triangle: "<<count<<" Vertices:"; mesh.tList[count].v1.show(); mesh.tList[count].v2.show(); mesh.tList[count].v3.show(); std::cout<<std::endl;
      std::cout<<"Triangle: "<<count<<" Color:"; mesh.tList[count].col.show(); std::cout<<std::endl;
      std::cout<<"Triangle: "<<count<<" Radiance:"<<mesh.tList[count].radiance<<std::endl;
      if (parseMaterial) {
	std::cout<<"Triangle: "<<count<<" Alpha:"<<mesh.tList[count].alpha<<std::endl;
	std::cout<<"Triangle: "<<count<<" Specular Coefficient:"<<mesh.tList[count].specCoef<<std::endl;
	std::cout<<"Triangle: "<<count<<" Reflectance:"<<mesh.tList[count].reflectance<<std::endl;
	std::cout<<"Triangle: "<<count<<" BRDF:"<<mesh.tList[count].brdf<<std::endl;
	std::cout<<"Triangle: "<<count<<" Internal IOR:"<<mesh.tList[count].intIOR<<std::endl;
	std::cout<<"Triangle: "<<count<<" External IOR:"<<mesh.tList[count].extIOR<<std::endl;
      }
    }
    count++;
  } while ((iter = iter->NextSiblingElement("TR")));
  
  if (parent->FirstChildElement("TRANSLATE"))
    if (!parseVec3(parent, "TRANSLATE", mesh.translate)) {
      std::cerr<<"No Mesh Translation Information."<<std::endl;
      return false;
    }
    
  mesh.scale.x = mesh.scale.y = mesh.scale.z = 1;  
  if (parent->FirstChildElement("SCALE"))
    if (!parseVec3(parent, "SCALE", mesh.scale)) {
      std::cerr<<"No Mesh Scaling Information."<<std::endl;
      return false;
    }  
    
  if (parent->FirstChildElement("ROTSEQ"))
    if (!parseRotSeq(parent, mesh.rotations)) {
      std::cerr<<"No Mesh Rotaion Information."<<std::endl;
      return false;
    }
    
  if (verbose) {
    std::cout<<"Mesh Translation: "; mesh.translate.show(); std::cout<<std::endl;
    std::cout<<"Mesh Scaling: "; mesh.scale.show(); std::cout<<std::endl;
    for (uint32_t i = 0; i < mesh.rotations.rotAxisSeq.size(); i++) {
      std::cout<<"Mesh Rotation "<<i<<" Axis: ";mesh.rotations.rotAxisSeq[i].show();std::cout<<" Degrees:"<<mesh.rotations.rotDegSeq[i]<<std::endl;
    }
  }
  
  return true;
}

bool SceneParser::parseMesh_Sphere(XMLElement *parent, SP_MESH_SPHERE &mesh, bool parseMaterial, bool verbose) {
   XMLElement *iter = parent->FirstChildElement("SP");
   if (!iter) {
     std::cerr<<"No Sphere objects found in mesh."<<std::endl;
     return false;
   }
   uint32_t count = 0;
   do {
    SP_SPHERE s; 
    if (!parseSphere(iter, s, parseMaterial)) {
      std::cerr<<"Sphere: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    mesh.sList.push_back(s);
    if (verbose) {
      std::cout<<"Sphere: "<<count<<" Radius:"<<mesh.sList[count].radius<<std::endl;
      std::cout<<"Sphere: "<<count<<" Position:"; mesh.sList[count].pos.show(); std::cout<<std::endl;
      std::cout<<"Sphere: "<<count<<" Color:"; mesh.sList[count].col.show(); std::cout<<std::endl;
      std::cout<<"Sphere: "<<count<<" Radiance:"<<mesh.sList[count].radiance<<std::endl;
      if (parseMaterial) {
	std::cout<<"Sphere: "<<count<<" Alpha:"<<mesh.sList[count].alpha<<std::endl;
	std::cout<<"Sphere: "<<count<<" Specular Coefficient:"<<mesh.sList[count].specCoef<<std::endl;
	std::cout<<"Sphere: "<<count<<" Reflectance:"<<mesh.sList[count].reflectance<<std::endl;
	std::cout<<"Sphere: "<<count<<" BRDF:"<<mesh.sList[count].brdf<<std::endl;
	std::cout<<"Sphere: "<<count<<" Internal IOR:"<<mesh.sList[count].intIOR<<std::endl;
	std::cout<<"Sphere: "<<count<<" External IOR:"<<mesh.sList[count].extIOR<<std::endl;
      }
    }
    count++;
  } while ((iter = iter->NextSiblingElement("SP")));
  
  if (parent->FirstChildElement("TRANSLATE"))
     if (!parseVec3(parent, "TRANSLATE", mesh.translate)) {
      std::cerr<<"No Sphere mesh Translation Information."<<std::endl;
      return false;
     }
    
  if (parent->FirstChildElement("SCALE"))
    std::cout<<"Scaling is undefined for sphere objects."<<std::endl;
    
  if (parent->FirstChildElement("ROTSEQ"))
    std::cout<<"Rotaion is undefined for sphere objects."<<std::endl;
  
  if (verbose)
    std::cout<<"Sphere Mesh Translation: "; mesh.translate.show(); std::cout<<std::endl;
  
  return true;
}

bool SceneParser::parseTriSources(SP_MESH *sources, bool verbose = false) {
  XMLElement *iter = doc.FirstChildElement("TRISOURCE");
  
  if (!iter) return false;
  
  uint32_t count = 0;
  do {
    if (!parseMesh(iter, sources[count], false, verbose)) {
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
    if (!parseMesh(iter, sources[count], false, verbose)) {
      std::cerr<<"MESH Source: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    count++;
  } while ((iter = iter->NextSiblingElement("MESHSOURCE")));
  
  return true;
}

bool SceneParser::parseMeshObjects(SP_MESH *objects, bool verbose) {
  XMLElement *iter = doc.FirstChildElement("MESHOBJECT");
  
  if (!iter) return false;
  
  uint32_t count = 0;
  do {
    if (!parseMesh(iter, objects[count], true, verbose)) {
      std::cerr<<"MESH Object: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    count++;
  } while ((iter = iter->NextSiblingElement("MESHOBJECT")));
  
  return true;
}

bool SceneParser::parseSphereObjects(SP_MESH_SPHERE *objects, bool verbose) {
  XMLElement *iter = doc.FirstChildElement("SPHEREOBJECT");
  
  if (!iter) return false;
  
  uint32_t count = 0;
  do {
    if (!parseMesh_Sphere(iter, objects[count], true, verbose)) {
      std::cerr<<"Sphere Mesh Object: "<<count<<" couldn't be parsed."<<std::endl;
      return false;  
    }
    count++;
  } while ((iter = iter->NextSiblingElement("SPHEREOBJECT")));
  
  return true;
}
