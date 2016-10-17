#include <vector>
#include <stdint.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include "geometryPrimitives.h"
#include "mathPrimitives.h"
#include "bvhAccel.h"
#define nPRIM_LEAF	1// See LinearBvhNode struct and check data-type of nPrimitives before changing this.

BvhAccel *bvhAccelT, *bvhAccelS;

void BvhAccel::initAccel() {
  totalNodes = 0;
  buildData.reserve(nPrimitives);
  for (unsigned int i = 0; i < nPrimitives; ++i) {
    AABBox bbox = ((BasePrimitive *)(primitiveList + i * primActualSize))->box;
    buildData.push_back(PrimitiveInfo(i, bbox));
  }
  build();
  flatten();
}

BvhNode* BvhAccel::recursiveBuild(uint32_t start, uint32_t end, uint32_t *totalNodes) {
  BvhNode *node = 0;
  try {
    node = new BvhNode();
  } catch(std::bad_alloc) {
    std::cerr<<"Failed to allocate memory in BvhAccel::recursiveBuild\n";
    exit(0);
  }
  (*totalNodes)++;
  if ((*totalNodes) == 0xFFFFFFFF) {
    std::cerr<<"Too many nodes.\n";
    exit(0);
  }
  AABBox bbox = buildData[start].box;
  for (uint32_t i = start + 1; i < end; ++i)
    bbox = AABBox::uNion(bbox, buildData[i].box);

  uint32_t nPrimitives = end - start;

  if (nPrimitives == nPRIM_LEAF) {
    uint32_t firstPrimOffset = orderedPrimitives.size();
    for (uint32_t i = start; i < end; ++i) {
      uint32_t primIdx = buildData[i].primitiveIdx;
      orderedPrimitives.push_back(primIdx);
    }
    node->initLeaf(firstPrimOffset, nPrimitives, bbox);
  }
  else {
    AABBox centroidBounds = AABBox(buildData[start].centroid);
    for (uint32_t i = start + 1; i < end; ++i)
      centroidBounds = AABBox::uNion(centroidBounds, buildData[i].centroid);

    int dim = centroidBounds.maximumExtent();

    // Next if is for the unlikey case when all centroids in the chosen dimension are coinciding.
    // Since we are taking the maximumExtent but still getting zero bound, that means bounding volume is actually zero.
    // So, it is likely that checking any one dimension should suffice. But still for the sake of completeness we check all dimension
    uint32_t mid = (start + end) / 2;
    if ((centroidBounds.pMax.x == centroidBounds.pMin.x && dim == 0) || // Next conditions are not needed.
	(centroidBounds.pMax.y == centroidBounds.pMin.y && dim == 1) ||
	(centroidBounds.pMax.z == centroidBounds.pMin.z && dim == 2)) {

      uint32_t firstPrimOffset = orderedPrimitives.size();
      for (uint32_t i = start; i < end; ++i) {
	uint32_t primIdx = buildData[i].primitiveIdx;
	orderedPrimitives.push_back(primIdx);
      }
      node->initLeaf(firstPrimOffset, nPrimitives, bbox);
      return node;
    }

    Vec pmid = (centroidBounds.pMin + centroidBounds.pMax) * 0.5;
    PrimitiveInfo *midPtr = std::partition(&buildData[start], &buildData[end - 1] + 1, CompareToMid(dim, pmid));
    mid = midPtr - &buildData[0];

    if (mid == end || mid == start) {
      uint32_t firstPrimOffset = orderedPrimitives.size();
      for (uint32_t i = start; i < end; ++i) {
	uint32_t primIdx = buildData[i].primitiveIdx;
	orderedPrimitives.push_back(primIdx);
      }
      node->initLeaf(firstPrimOffset, nPrimitives, bbox);
      return node;
    }
    node->initInterior(dim, recursiveBuild(start, mid, totalNodes),
	    recursiveBuild(mid, end, totalNodes));

  }
  return node;
}

void BvhAccel::build() {
  orderedPrimitives.reserve(nPrimitives);
  root = recursiveBuild(0, nPrimitives, &totalNodes);
  // reorder original primitiveList as per orderedPrimitives
  // code used from https://github.com/mission-peace/interview/blob/master/src/com/interview/array/ReorderArrayByIndex.java
  /*uint8_t *sVal = (uint8_t *) malloc(primActualSize);
  for (uint32_t i = 0 ; i < orderedPrimitives.size(); i++) {
    while (orderedPrimitives[i] != i) {
      uint32_t sIndex = orderedPrimitives[orderedPrimitives[i]];
      memcpy(sVal, &primitiveList[orderedPrimitives[i] * primActualSize], primActualSize);

       orderedPrimitives[orderedPrimitives[i]] = orderedPrimitives[i];
       memcpy(&primitiveList[orderedPrimitives[i] * primActualSize], &primitiveList[i * primActualSize], primActualSize);

       orderedPrimitives[i] = sIndex;
       memcpy(&primitiveList[i * primActualSize], sVal, primActualSize);
     }
  }
  free (sVal);*/
  uint8_t *tempArr = (uint8_t *) malloc(primActualSize * nPrimitives);

  for (uint32_t i = 0; i < nPrimitives; i++) {
    memcpy(&tempArr[i*primActualSize], &primitiveList[orderedPrimitives[i] * primActualSize], primActualSize);
  }

  for (uint32_t i = 0; i < nPrimitives; i++)
    memcpy(&primitiveList[i * primActualSize], &tempArr[i*primActualSize], primActualSize);

  free(tempArr);
}

uint32_t BvhAccel::recursiveFlatten(BvhNode *node, uint32_t *offset) {
  LinearBvhNode *linearNode = &nodes[*offset];
  linearNode->bounds = node->bounds;
  uint32_t myOffset = (*offset)++;

  if (node->nPrimitives > 0) {
    linearNode->primitivesOffset = node->firstPrimOffset;
    linearNode->nPrimitives = node->nPrimitives;
  }
  else {
    linearNode->axis = node->splitAxis;
    linearNode->nPrimitives = 0;
    recursiveFlatten(node->children[0], offset);
    linearNode->secondChildOffset = recursiveFlatten(node->children[1], offset);
 }

 return myOffset;
}

void BvhAccel::flatten() {
  nodes = new LinearBvhNode[totalNodes];
  uint32_t offset = 0;
  recursiveFlatten(root, &offset);
}

bool BvhAccel::intersectS(const Ray &r, double &t, int &id) {
  if (!nodes) return false;
  bool hit = false;

  Vec invDir(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
  uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };

  uint32_t todoOffset = 0, nodeNum = 0;
  uint32_t todo[64];

  double d;
  id = 0xFFFFFFFF;
  t = INF;

  while (true) {
    const LinearBvhNode *node = &nodes[nodeNum];
    if (node->bounds.intersect(r) < INF) {
	if (node->nPrimitives > 0) {
	  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
	    Sphere * ptr = (Sphere *)(primitiveList + (node->primitivesOffset + i) * primActualSize);
	    if ((d = ptr -> intersect(r)) && d < t) {
	      t = d; // set the distance of intersection.
	      id = node->primitivesOffset + i; // Set the serial no. of the intersecting sphere
	      hit = true;
	    }
	  }
	  if (todoOffset == 0) break;
	  nodeNum = todo[--todoOffset];
	}
	else {
	  if (dirIsNeg[node->axis]) {
	    todo[todoOffset++] = nodeNum + 1;
	    nodeNum = node->secondChildOffset;
	  }
	  else {
	    todo[todoOffset++] = node->secondChildOffset;
	    nodeNum = nodeNum + 1;
	  }
      }
    }
    else {
      if (todoOffset == 0) break;
      nodeNum = todo[--todoOffset];
    }

  }
  return hit;
}

bool BvhAccel::intersectT(const Ray &r, double &t, Vec &N, int &id) {
  if (!nodes) return false;
  bool hit = false;

  Vec invDir(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
  uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };

  uint32_t todoOffset = 0, nodeNum = 0;
  uint32_t todo[64];

  double d;
  Vec N_;
  id = 0xFFFFFFFF;
  t = INF;

  while (true) {
    const LinearBvhNode *node = &nodes[nodeNum];
    if (node->bounds.intersect(r) < INF) {
	if (node->nPrimitives > 0) {
	  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
	    Triangle * ptr = (Triangle *)(primitiveList + (node->primitivesOffset + i) * primActualSize);
	    if ((d = ptr -> intersect(r, N_)) && d < t) {
	      t = d; // set the distance of intersection.
	      id = node->primitivesOffset + i; // Set the serial no. of the intersecting sphere
	      N = N_;
	      hit = true;
	    }
	  }
	  if (todoOffset == 0) break;
	  nodeNum = todo[--todoOffset];
	}
	else {
	  if (dirIsNeg[node->axis]) {
	    todo[todoOffset++] = nodeNum + 1;
	    nodeNum = node->secondChildOffset;
	  }
	  else {
	    todo[todoOffset++] = node->secondChildOffset;
	    nodeNum = nodeNum + 1;
	  }
      }
    }
    else {
      if (todoOffset == 0) break;
      nodeNum = todo[--todoOffset];
    }

  }
  return hit;
}