#if !defined(dummyAccel_h__)
#define dummyAccel_h__

// This is a dummy accelaration structure for speeding up single objects.
struct BoxContent {
  int start, end;
  AABBox box;
  BoxContent(int a = 0, int b = 0, AABBox bb = AABBox(Vec(1,1,1))): start(a), end(b), box(bb) {}
};

class DummyAccel {
  static int nBoxT;
  static BoxContent *boxT;
  static AABBox box;
public:
  // Not thread safe
  static void initAccel(int nTriangles, Triangle *triangleList);
  // Thread safe
  static bool intersect(const Ray &r, double &t, Vec &N, int &id, int nTriangles, Triangle *list);
};

#endif

