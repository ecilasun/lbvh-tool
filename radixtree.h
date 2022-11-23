#pragma once

#include "emath.h"
#include <vector>

struct SRadixTreeNode
{
    SBoundingBox m_bounds;
    uint32_t m_spatialKey{0xFFFFFFFF}, m_primitiveIndex{0xFFFFFFFF};
    uint32_t m_left{0xFFFFFFFF}, m_right{0xFFFFFFFF};
};

struct SRTGeometryNode
{
    SVec128 v0, v1, v2;
};

typedef bool(*HitTestFunc)(const SRadixTreeNode &_self, const SVec128 &_rayStart, const SVec128 &_rayDir, float &_t, const float _tmax);

void GenerateLBVH(SRadixTreeNode *_nodes, std::vector<SRadixTreeNode> &_leafNodes, const int _numNodes);
void FindClosestHitLBVH(SRadixTreeNode *_nodes, const int _numNodes, const SVec128 &_rayStart, const SVec128 &_rayEnd, float &_t, SVec128 &_hitPos, uint32_t &_hitNode, HitTestFunc _hitTestFunc);
void FindAnyHitLBVH(SRadixTreeNode *_nodes, const int _numNodes, const SVec128 &_rayStart, const SVec128 &_rayEnd, float &_t, SVec128 &_hitPos, uint32_t &_hitNode, HitTestFunc _hitTestFunc);
