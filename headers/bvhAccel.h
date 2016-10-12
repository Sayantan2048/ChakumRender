#if !defined(bvhAccel_h__)
#define bvhAccel_h__
#include <vector>
#include <stdint.h>

struct PrimitiveInfo {
  AABBox box;
  uint32_t primitiveIdx;
  Vec centroid;
  PrimitiveInfo(unsigned int pn, const AABBox &b) : box(b), primitiveIdx(pn) {
    centroid = b.pMin * 0.5 + b.pMax * 0.5;
  }
};

struct BvhNode {
  AABBox bounds;
  BvhNode *children[2];
  uint32_t splitAxis; // for interior nodes, this information is used to optimize traversal algorithm.
  uint32_t firstPrimOffset, nPrimitives; // Used only for leaf nodes. Offset and number of primitives in the leaf node are linked with array orderedPrimitives.

  BvhNode(AABBox b = AABBox(Vec(1,1,1))) : bounds(b) { children[0] = children[1] = 0; }

  void initLeaf(uint32_t first, uint32_t n, const AABBox &b) {
    firstPrimOffset = first;
    nPrimitives = n;
    bounds = b;
  }
  void initInterior(uint32_t axis, BvhNode *c0, BvhNode *c1) {
    children[0] = c0;
    children[1] = c1;
    bounds = AABBox::uNion(c0->bounds, c1->bounds);
    splitAxis = axis;
    nPrimitives = 0;
  }
};

struct LinearBvhNode {
  AABBox bounds;
  union {
    uint32_t primitivesOffset; // leaf, Offset to re-ordered primitiveList array.
    uint32_t secondChildOffset; // interior, offset into flattended nodes array for second child.
  };
  uint32_t nPrimitives; // 0 -> interior node
  uint8_t axis; // interior node: xyz
  //uint8_t pad[2]; // ensure 32 byte total size
  LinearBvhNode(const AABBox &b = AABBox(Vec(1., 1., 1.))) : bounds(b) {}
};

class BvhAccel {
   uint8_t *primitiveList; //original list of primitives. This list is re-ordered as per orderedList.
   uint32_t nPrimitives; // number of primitives in original list.
   std::size_t primActualSize; // i.e sizeof(Triangle) or sizeof(Sphere)
   std::vector<uint32_t> orderedPrimitives; // Contains indexes to primitives in depth first order.
   std::vector<PrimitiveInfo> buildData;
   BvhNode *root;
   uint32_t totalNodes;
   void build();
   BvhNode* recursiveBuild(uint32_t start, uint32_t end, uint32_t *);
   struct CompareToMid {
      CompareToMid(int d, Vec m) { dim = d; mid = m; }
      int dim;
      Vec mid;
      bool operator()(const PrimitiveInfo &a) const {
	return (dim == 0 && a.centroid.x < mid.x) || (dim == 1 && a.centroid.y < mid.y) || (dim == 2 && a.centroid.z < mid.z);
      }
   };
   LinearBvhNode *nodes; //Array of flattended nodes in depth first order.
   void flatten();
   uint32_t recursiveFlatten(BvhNode *node, uint32_t *offset);
public:
   BvhAccel(uint8_t *list, uint32_t nP, std::size_t sz) : primitiveList(list), nPrimitives(nP), primActualSize(sz) {root = 0;}
   void initAccel(); // build the BVH tree and flatten it.
   bool intersect(const Ray &r, double &t, Vec &N, int &id);

};

extern BvhAccel *bvhAccel;

#endif