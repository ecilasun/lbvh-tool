#pragma once

#include "emath.h"
#include <vector>

#define INVALID_NODE 0xFFFFFFFF

struct SRadixTreeNode
{
    SBoundingBox m_bounds;
    uint32_t m_spatialKey{0xFFFFFFFF}, m_primitiveIndex{0xFFFFFFFF};
    uint32_t m_left{0xFFFFFFFF}, m_right{0xFFFFFFFF};
};

struct SPackedRadixTreeNode
{
	float m_leftBoundsMinX_orV0X, m_leftBoundsMinY_orV0Y, m_leftBoundsMinZ_orV0Z;
	float m_leftBoundsMaxX_orV1X, m_leftBoundsMaxY_orV1Y, m_leftBoundsMaxZ_orV1Z;
	float m_rightBoundsMinX_orV2X, m_rightBoundsMinY_orV2Y, m_rightBoundsMinZ_orV2Z;
	float m_rightBoundsMaxX_orV3X, m_rightBoundsMaxY_orV3Y, m_rightBoundsMaxZ_orV3Z;
	uint32_t m_leftNode{0xFFFFFFFF}, m_rightNode{0xFFFFFFFF};
	uint32_t m_spatialKey{0xFFFFFFFF};
};

struct triangle final {
  SVec128 coords[3];
  SVec128 normals[3];
};

struct HitInfo
{
	SVec128 hitPos;
	uint32_t triIndex;
	uint32_t traversalCount;
	triangle *geometryIn;
	triangle *geometryOut;
};

typedef bool(*HitTestFunc)(const SRadixTreeNode &_self, const SVec128 &_rayStart, const SVec128 &_rayEnd, float &_t, const float _tmax, HitInfo &_hitinfo);

void GenerateLBVH(SRadixTreeNode *_nodes, std::vector<SRadixTreeNode> &_leafNodes, const int _numNodes);
void GeneratePackedLBVH(SPackedRadixTreeNode *_nodes, std::vector<SPackedRadixTreeNode> &_leafNodes, const int _numNodes);
void FindClosestHitLBVH(SRadixTreeNode *_nodes, const int _numNodes, const SVec128 &_rayStart, const SVec128 &_rayEnd, float &_t, uint32_t &_hitNode, HitInfo &_hitinfo, HitTestFunc _hitTestFunc);
void FindClosestHitLBVHPacked(SPackedRadixTreeNode *_nodes, const int _numNodes, const SVec128 &_rayStart, const SVec128 &_rayEnd, float &_t, uint32_t &_hitNode, HitInfo &_hitInfo);
