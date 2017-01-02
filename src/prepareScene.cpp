#include "prepareScene.h"
#include "transformations.h"
#include "materialTypes.h"
#include <cmath>

bool initSceneParser(const char* fileName ) {
  bool status;
  sceneParser = new SceneParser(fileName, status);
  return status;
}

bool prepareScreenCamera(uint32_t &width, uint32_t &height, Vec &camPos, Vec &camDir, Vec &cx, Vec &cy, bool verbose = false) {
  if (!sceneParser->parseScreenSize(width, height, verbose))
      return false;
  
  double fovx, fovy;
  SP_Vec cpos, cdir, up; 
  if (!sceneParser->parseCamera(cpos, cdir, up, fovx, fovy, verbose))
    return false;
  
  camPos = Vec(cpos.x, cpos.y, cpos.z);
  camDir = Vec(cdir.x, cdir.y, cdir.z).norm();
  
  cx = camDir%Vec(up.x, up.y, up.z);
  if (cx.length() < 1e-20) {
    std::cerr<<"UP and LOOKAT cannot be parallel!"<<std::endl;
    return false;
  }
    
  cx.norm();
  cx = cx * ((double)width / height) * tan(fovx * PI / 360.0);
  
  cy = (cx % camDir).norm() * tan(fovy * PI / 360.0);
  
  return true;
}

static Vec SPVecToVec(const SP_Vec &x) {
  return Vec(x.x, x.y, x.z);
}

static mat3 SPVecToMat(const SP_Vec &s) {
  return mat3(s.x, 0, 0, 0, s.y, 0, 0, 0, s.z);
}

static mat3 RotSeqToMat(const SP_ROTSEQ &rs) {
  mat3 rot;
  for (uint32_t i = 0; i < rs.rotDegSeq.size(); i++) {
    rot = rot.mul(rotate(SPVecToVec(rs.rotAxisSeq[i]), rs.rotDegSeq[i] * PI / 180.0));
  }
  return rot;
}

bool prepareLightSources(bool verbose = false) {
  lSource = new LightSource();
  
  /* Mesh Sources. */
  uint32_t count = sceneParser->findMeshSources(verbose);
  if (count > 0) {
    SP_MESH *sources = new SP_MESH[count];
    if (!sceneParser->parseMeshSources(sources, verbose))
      return false;
    
    for (uint32_t i = 0; i < count; i++) {
	Vec trans = SPVecToVec(sources[i].translate);
	mat3 m = SPVecToMat(sources[i].scale);
	m = RotSeqToMat(sources[i].rotations).mul(m);
	MeshLight mesh;
	for (uint32_t j = 0; j < sources[i].tList.size(); j++) {
	  SP_TR t = sources[i].tList[j];
	  Vec v1 = SPVecToVec(t.v1);
	  Vec v2 = SPVecToVec(t.v2);
	  Vec v3 = SPVecToVec(t.v3);
	  Vec col = SPVecToVec(t.col);
	  v1 = m.mul(v1);
	  v2 = m.mul(v2);
	  v3 = m.mul(v3);
	  
	  mesh.add(TriLight(v1 + trans, v2 + trans, v3 + trans, col, t.radiance));
	}
	mesh.initMeshLight();
	lSource->addMSource(mesh);
    }
    delete []sources;
  }
  
  /* Triangle Sources*/
  count = sceneParser->findTriSources(verbose);
  if (count > 0) {
    SP_MESH *sources = new SP_MESH[count];
    if (!sceneParser->parseTriSources(sources, verbose))
      return false;
    
    for (uint32_t i = 0; i < count; i++) {
	Vec trans = SPVecToVec(sources[i].translate);
	mat3 m = SPVecToMat(sources[i].scale);
	m = RotSeqToMat(sources[i].rotations).mul(m);
	SP_TR t = sources[i].tList[0];
	
	Vec v1 = SPVecToVec(t.v1);
	Vec v2 = SPVecToVec(t.v2);
	Vec v3 = SPVecToVec(t.v3);
	Vec col = SPVecToVec(t.col);
	v1 = m.mul(v1);
	v2 = m.mul(v2);
	v3 = m.mul(v3);
	  
	lSource->addTSource(TriLight(v1 + trans, v2 + trans, v3 + trans, col, t.radiance));
    }
    delete []sources;
  }
  
  /*Point Sources*/
  count = sceneParser->findPointSources(verbose);
  if (count > 0) {
    SP_Vec *pos = new SP_Vec[count];
    SP_Vec *col = new SP_Vec[count];
    double *power = new double[count];
  
    if (!sceneParser->parsePointSources(pos, col, power, verbose))
      return false;
    
    for (uint32_t i = 0; i < count; i++)
      lSource->addPSource(PointSource(SPVecToVec(pos[i]), SPVecToVec(col[i]), power[i]));
  
    delete []pos;
    delete []col;
    delete []power;
  }
  
  count = sceneParser->findSphereSources(verbose);
  if (count > 0) {
    double *radius = new double[count];
    SP_Vec *pos = new SP_Vec[count];
    SP_Vec *col = new SP_Vec[count];
    double *radiance = new double[count];
  
    if (!sceneParser->parseSphereSources(radius, pos, col, radiance, verbose))
      return false;
    
    for (uint32_t i = 0; i < count; i++)
      lSource->addSSource(SphereSource(radius[i], SPVecToVec(pos[i]), SPVecToVec(col[i]), radiance[i]));
    
    delete []radius;
    delete []pos;
    delete []col;
    delete []radiance;
  }
  
  /*uint32_t w, h;
  Vec *image;
  readImage(image, w, h);
  lSource->addESource(EnvSource(image, w, h, 1, Vec(1.0, 1.0, 1.0), 2));
  delete []image;*/
   
  return true;
}

static BRDFType stringToBRDF(const std::string s) {
  BRDFType t;

  if (s.compare("GGX") == 0)
    t = GGX;
  else if (s.compare("BECKMANN") == 0)
    t = BECKMANN;
  else if (s.compare("PHONG") == 0)
    t = PHONG;
  else if (s.compare("GGXAPPROX") == 0)
    t = GGXAPPROX;
  else if (s.compare("CLASSIC_PHONG") == 0)
    t = CLASSIC_PHONG;
  else {
    std::cerr<<"Unrecognized BRDF!!"<<std::endl;
    t = BRDF_NONE;
  }
  
  return t;
}

bool prepareObjects(bool verbose = false) {
  /* Mesh Objects*/
  uint32_t count = sceneParser->findMeshObjects(verbose);
  if (count > 0) {
    SP_MESH *objects = new SP_MESH[count];
    if (!sceneParser->parseMeshObjects(objects, verbose))
      return false;
    
    for (uint32_t i = 0; i < count; i++) {
	Vec trans = SPVecToVec(objects[i].translate);
	mat3 m = SPVecToMat(objects[i].scale);
	m = RotSeqToMat(objects[i].rotations).mul(m);
	
	for (uint32_t j = 0; j < objects[i].tList.size(); j++) {
	  SP_TR t = objects[i].tList[j];
	  Vec v1 = SPVecToVec(t.v1);
	  Vec v2 = SPVecToVec(t.v2);
	  Vec v3 = SPVecToVec(t.v3);
	  Vec col = SPVecToVec(t.col);
	  v1 = m.mul(v1) + trans;
	  v2 = m.mul(v2) + trans;
	  v3 = m.mul(v3) + trans;
          
	  BRDFType b = stringToBRDF(t.brdf);
	  
	  vTriangleList.push_back(Triangle(v1, v2, v3, col, t.reflectance, 
		b == CLASSIC_PHONG ? MaterialType(t.alpha, t.specCoef) : MaterialType(t.alpha, t.intIOR, t.extIOR, t.specCoef, b)));
	}
    }
    
    delete []objects;
  }
  
  /* Sphere Mesh Objects*/
  count = sceneParser->findSphereObjects(verbose);
  if (count > 0) {
    SP_MESH_SPHERE *objects = new SP_MESH_SPHERE[count];
    if (!sceneParser->parseSphereObjects(objects, verbose))
      return false;
    
    for (uint32_t i = 0; i < count; i++) {
	Vec trans = SPVecToVec(objects[i].translate);
	for (uint32_t j = 0; j < objects[i].sList.size(); j++) {
	  SP_SPHERE s = objects[i].sList[j];
	  Vec pos = SPVecToVec(s.pos) + trans;
	  Vec col = SPVecToVec(s.col);
	           
	  BRDFType b = stringToBRDF(s.brdf);
	  
	  vSphereList.push_back(Sphere(s.radius, pos, col, s.reflectance,  
		b == CLASSIC_PHONG ? MaterialType(s.alpha, s.specCoef) : MaterialType(s.alpha, s.intIOR, s.extIOR, s.specCoef, b)));
	}
    }
    
    delete []objects;
  }

  //Dummy
  vSphereList.push_back(Sphere(1,  Vec(27, 1e5, 47), Vec(.999, .999, .999), 1.0, MaterialType(200.2, 1)));
  
  //Dummy Triangle
  vTriangleList.push_back(Triangle(1e6, 1e6, 1e6, Vec(0.0, 0, 0.9), 0.9, MaterialType(1.0, 0.5)));
  
  //vTriangleList.push_back(Triangle(VA, VB, VC, Vec(0.0, 0, 0.9), 0.9, MaterialType(1.0, 0.5)));
  //vTriangleList.push_back(Triangle(VB, VC, VD, Vec(0.9, 0.9, 0), 0.9, MaterialType(1.0, 0.5)));
  //vTriangleList.push_back(Triangle(VC, VD, VA, Vec(0, 0.9, 0.9), 0.9, MaterialType(1.0, 0.5)));
  //vTriangleList.push_back(Triangle(VD, VA, VB, Vec(0.9, 0.0, 0.9), 0.9, MaterialType(1.0, 0.5)));

  /*  objLoader *objData = new objLoader();
  objData->load("Aventador1.obj");

  nTriangles = objData->faceCount;

  triangleList = new Triangle[nTriangles];

  mat3 M = scale(35, 35, 35).mul(rotateY(PI/4 + PI));

  for(int i=0; i < nTriangles; i++) {
	obj_face *o = objData->faceList[i];

	obj_vector A, B, C;
	A = *objData->vertexList[o->vertex_index[0]];
	B = *objData->vertexList[o->vertex_index[1]];
	C = *objData->vertexList[o->vertex_index[2]];

	Vec Av, Bv, Cv;
	Av = M.mul(Vec(A.e[0], A.e[1], A.e[2])) + Vec(X, Y, Z);
	Bv = M.mul(Vec(B.e[0], B.e[1], B.e[2])) + Vec(X, Y, Z);
	Cv = M.mul(Vec(C.e[0], C.e[1], C.e[2])) + Vec(X, Y, Z);

	// Assignment operator, data of Triangle() will be copied to triangleList[i].
	// Triangle() object is however temporary.
	triangleList[i] = Triangle(Av, Bv, Cv, Vec(1.0, 0.01, 0.01), 1.0, MaterialType(2.2, 1.0));

  }*/


  //DummyAccel::initAccel(nTriangles, triangleList);
  
  return true;
}

