#include "platform.h"
#include "math.h"
#include "locklesspipe.h"
#include "bvh8.h"
#include "objloader.h"

// Folder/file name of the scene to work with
//#define SCENENAME "sibenik"
#define SCENENAME "testscene"

// Define to use 1.f as cell size
//#define USE_UNIT_CELL

// Define this to get double-sided hits (i.e. no backface culling against incoming ray)
//#define DOUBLE_SIDED

// Define to enable reflections
#define ENABLE_REFLECTIONS

// Define to enable exponential fog (requires DOUBLE_SIDED)
//#define ENABLE_EXPFOG

// Define to enable shadows
#define ENABLE_SHADOWS

// Define to enable light occlusion or disable for NdotL
//#define ENABLE_OCCLUSION

// NOTE: Once there's a correct cell vs triangle test, this will be reduced
#define MAX_NODE_TRIS 64

// Depth of traversal stack
#define MAX_STACK_ENTRIES 16

// Define this to ignore triangles and work with BVH child nodes only, and see in x-ray vision
//#define IGNORE_CHILD_DATA

// Define to ignore hits and keep marching - use this together with SHOw_HEATMAP
//#define XRAY_MODE

// Define to see traversal count per tile
//#define SHOW_HEATMAP

// Define this to get worker tile debug view
//#define SHOW_WORKER_IDS

// Define this to show depth information only (no barycentrics/lighting etc)
//#define SHOW_DEPTH

// Number of worker threads
#define MAX_WORKERS 12

// Define to use Morton curve order instead of scanline-first
//#define USE_MORTON_ORDER

static bool g_done = false;
#if defined(USE_MORTON_ORDER)
static const uint32_t tilewidth = 4;
static const uint32_t tileheight = 4;
static const uint32_t width = 512;
static const uint32_t height = 512;
#else
// 320x240 but x2
static const uint32_t tilewidth = 4;
static const uint32_t tileheight = 8;
static const uint32_t width = 640;
static const uint32_t height = 480;
#endif
static const uint32_t tilecountx = width/tilewidth;
static const uint32_t tilecounty = height/tileheight;
static const float cameradistance = 20.f;
static const float raylength = cameradistance + 100.f;

struct SRenderContext
{
	float rotAng{0.f};
	float aspect{float(height) / float(width)};

	SVec128 rayOrigin;
	SVec128 lookAt;
	SVec128 upVec;
	SMatrix4x4 lookMat;
	SVec128 pzVec;
	SVec128 F;
	SVec128 nil{0.f, 0.f, 0.f, 0.f};
	SVec128 epsilon{0.02f, 0.02f, 0.02f, 0.f};
	SVec128 negepsilon{-0.02f, -0.02f, -0.02f, 0.f};
	uint8_t *pixels;
};

struct SWorkerContext
{
	uint32_t workerID{0};
	uint8_t rasterTile[4*tilewidth*tileheight]; // Internal rasterization tile (in hardware, to avoid arbitration need)
	SRenderContext *rc;
	CLocklessPipe<10> dispatchvector; // 1024 byte queue per worker
};

// NOTE: This structure is very expensive for E32E, a data reduction method has to be applied here.
// One approach could be to only stash connected triangle fans here.
struct BVH8LeafNode
{
	uint32_t m_numTriangles{0};
	uint32_t m_triangleIndices[MAX_NODE_TRIS]{};
};

struct triangle final {
  SVec128 coords[3];
  SVec128 normals[3];
};

struct hitinfo final {
	float hitT;
};

triangle *testtris;

SBVH8Database<BVH8LeafNode>* testBVH8;

static const uint32_t MAX_CLIPPED_VERTICES = 36;

void ClipEdgeToPlaneBegin(SVec128 _planeNormal, SVec128 _planeOrigin, SVec128& _prevPos, SVec128 _pos, float &_prevDist, int& _clipIndex, SVec128 _vClip[MAX_CLIPPED_VERTICES])
{
	SVec128 dotCurr = EVecDot3(_planeNormal, EVecSub(_pos, _planeOrigin));
	float dist = EVecGetFloatX(dotCurr);

	// No previous point, this is the first vertex [A>]------>B
	if (dist >= 0.f) // vertex is in front of or on the plane
		_vClip[_clipIndex++] = _pos;
	// else, vertex is behind the plane, ignore it

	// Remember where we came from
	_prevPos = _pos;
	_prevDist = dist;
}

void ClipEdgeToPlaneNext(SVec128 _planeNormal, SVec128 _planeOrigin, SVec128& _prevPos, SVec128 _pos, float &_prevDist, int& _clipIndex, SVec128 _vClip[MAX_CLIPPED_VERTICES])
{
	SVec128 dotCurr = EVecDot3(_planeNormal, EVecSub(_pos, _planeOrigin));
	float dist = EVecGetFloatX(dotCurr);

	float clipRatio = _prevDist / (_prevDist - dist);
	SVec128 clipRatioVec = { clipRatio, clipRatio, clipRatio, 1.f };

	uint32_t slotSelDist = dist < 0.f ? 1 : 0;
	uint32_t slotSelPrevDist = _prevDist < 0.f ? 1 : 0;
	uint32_t slot1 = (!slotSelDist && slotSelPrevDist) ? _clipIndex + 1 : MAX_CLIPPED_VERTICES - 1;

	SVec128 clippedPos = EVecAdd(_prevPos, EVecMul(EVecSub(_pos, _prevPos), clipRatioVec));
	_vClip[_clipIndex] = (slotSelDist ^ slotSelPrevDist) ? clippedPos : _pos;
	_vClip[slot1] = _pos;

	uint32_t stepSel = (slotSelDist && slotSelPrevDist) ? 0 : ((!slotSelDist && slotSelPrevDist) ? 2 : 1);
	_clipIndex += stepSel;

	// Remember where we came from
	_prevPos = _pos;
	_prevDist = dist;
}


void bvh8Builder(triangle* _triangles, uint32_t _numTriangles, SBVH8Database<BVH8LeafNode>* _bvh8)
{
	SVec128 sceneMin{FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
	SVec128 sceneMax{-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

	// Generate scene AABB
	for (uint32_t i = 0; i < _numTriangles; ++i)
	{
		sceneMin = EVecMin(sceneMin, _triangles[i].coords[0]);
		sceneMin = EVecMin(sceneMin, _triangles[i].coords[1]);
		sceneMin = EVecMin(sceneMin, _triangles[i].coords[2]);

		sceneMax = EVecMax(sceneMax, _triangles[i].coords[0]);
		sceneMax = EVecMax(sceneMax, _triangles[i].coords[1]);
		sceneMax = EVecMax(sceneMax, _triangles[i].coords[2]);
	}

	sceneMin = EVecFloor(sceneMin);
	sceneMax = EVecCeiling(sceneMax);

#ifdef USE_UNIT_CELL
	float minCell = 1.f;
#else
	// Find min ideal cell size
	SVec128 maxcelledge{0.0f,0.0f,0.0f,0.f};
	for (uint32_t i = 0; i < _numTriangles; ++i)
	{
		SVec128 triMin{FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
		SVec128 triMax{-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

		SVec128 v0 = _triangles[i].coords[0];
		SVec128 v1 = _triangles[i].coords[1];
		SVec128 v2 = _triangles[i].coords[2];

		triMin = EVecMin(triMin, v0);
		triMin = EVecMin(triMin, v1);
		triMin = EVecMin(triMin, v2);

		triMax = EVecMax(triMax, v0);
		triMax = EVecMax(triMax, v1);
		triMax = EVecMax(triMax, v2);

		// Always grab longest edge
		maxcelledge = EVecMax(EVecSub(triMax, triMin), maxcelledge);
	}
	float cellx = 2.f*EVecGetFloatX(maxcelledge);
	float celly = 2.f*EVecGetFloatY(maxcelledge);
	float cellz = 2.f*EVecGetFloatZ(maxcelledge);
	float minCell = EMinimum(cellx, EMinimum(celly, cellz)) / float(MAX_NODE_TRIS);
#endif

	printf("Using automatic cell size: %f\n", minCell);

	if (_bvh8->LoadBVH8(SCENENAME ".cache.bv8") == 0) // Cannot load from cache, rebuild and cache result for next time
	{
		// Set up grid bounds and cell scale
		_bvh8->m_GridAABBMin = sceneMin;
		_bvh8->m_GridCellSize = SVec128{minCell, minCell, minCell, 0.f};

		// Get scene span in cell units
		uint32_t gridMin[3], gridMax[3];
		_bvh8->QuantizePosition(sceneMin, gridMin);
		_bvh8->QuantizePosition(sceneMax, gridMax);

		// Scan the grid for intersecting triangles per cell
		for (int z = gridMin[2]; z <= gridMax[2]; ++z)
		for (int y = gridMin[1]; y <= gridMax[1]; ++y)
		for (int x = gridMin[0]; x <= gridMax[0]; ++x)
		{
			uint32_t spatialKey = _bvh8->EncodeKey(x, y, z);

			// World bounds of this cell
			SVec128 cellMin, cellMax;
			SVec128 cellV000 = EVecConvertIntToFloat(EVecConsti(x,y,z,0));
			SVec128 cellV111 = EVecConvertIntToFloat(EVecConsti(x+1,y+1,z+1,0));
			_bvh8->ToWorldUnits(cellV000, cellMin);
			_bvh8->ToWorldUnits(cellV111, cellMax);

			uint32_t keyIndex;
			uint32_t dataIndex = _bvh8->FindCell(spatialKey, keyIndex);

			for (uint32_t i = 0; i < _numTriangles; ++i)
			{
				SVec128 triMin{FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
				SVec128 triMax{-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

				SVec128 v0 = _triangles[i].coords[0];
				SVec128 v1 = _triangles[i].coords[1];
				SVec128 v2 = _triangles[i].coords[2];

				triMin = EVecMin(triMin, v0);
				triMin = EVecMin(triMin, v1);
				triMin = EVecMin(triMin, v2);

				triMax = EVecMax(triMax, v0);
				triMax = EVecMax(triMax, v1);
				triMax = EVecMax(triMax, v2);

				uint32_t elementMin[3], elementMax[3];
				_bvh8->QuantizePosition(triMin, elementMin);
				_bvh8->QuantizePosition(triMax, elementMax);

				// This primitive's aabb is crossed by the cell aabb, it's a good candidate to test
				if (x>=elementMin[0] && x<=elementMax[0] && y>=elementMin[1] && y<=elementMax[1] && z>=elementMin[2] && z<=elementMax[2])
				{
					// Do a more precise check by clipping the triangle into the cell bounds
					// This is used to discard empty cells that don't actually belong to the geometry
					SVec128 planeOrigin = cellV000;
					SVec128 nmx = EVecConst(-1.f,0.f,0.f,0.f);
					SVec128 nx = EVecConst(1.f,0.f,0.f,0.f);
					SVec128 nmy = EVecConst(0.f,-1.f,0.f,0.f);
					SVec128 ny = EVecConst(0.f,1.f,0.f,0.f);
					SVec128 nmz = EVecConst(0.f,0.f,-1.f,0.f);
					SVec128 nz = EVecConst(0.f,0.f,1.f,0.f);
					SVec128 planeNormals[6] = { nmx, nx, nmy, ny, nmz, nz };
					SVec128 planeOrigins[6] = {
						EVecAdd(planeOrigin,nx), planeOrigin,
						EVecAdd(planeOrigin,ny), planeOrigin,
						EVecAdd(planeOrigin,nz), planeOrigin };

					int clipSuccess = 0;
					int inputCount = 4;
					SVec128 vClip[MAX_CLIPPED_VERTICES], pos[MAX_CLIPPED_VERTICES];
					SVec128 l0, l1, l2;
					_bvh8->ToGridUnits(v0, pos[0]);
					_bvh8->ToGridUnits(v1, pos[1]);
					_bvh8->ToGridUnits(v2, pos[2]);
					_bvh8->ToGridUnits(v0, pos[3]);
					int clipIndex;
					for (int plane = 0; plane < 6; ++plane)
					{
						clipIndex = 0;
						SVec128 prevPos = pos[0];
						float prevDist = 0.f;
						ClipEdgeToPlaneBegin(planeNormals[plane], planeOrigins[plane], prevPos, pos[0], prevDist, clipIndex, vClip);
						for (int v = 0; v < inputCount; ++v)
							ClipEdgeToPlaneNext(planeNormals[plane], planeOrigins[plane], prevPos, pos[v], prevDist, clipIndex, vClip);

						SVec128i lastFirstSameMask = EVecCmpEQ(vClip[0], vClip[clipIndex-1]);
						int lastFirstEqual = EVecMoveMask(lastFirstSameMask);
						if (lastFirstEqual != 0x0000FFFF) // Ended with dissimilar vertex, close the polygon
							vClip[clipIndex++] = vClip[0];

						if (clipIndex <= 3) // One of the planes fully failed or this is a degenerate
							break;

						// Copy to input for next phase
						for (int p = 0; p < clipIndex; ++p)
							pos[p] = vClip[p];
						inputCount = clipIndex;

						clipSuccess++;
					}

					// Add a new data entry if clip test passed
					if (clipSuccess > 5)
					{
						// There's nothing here yet, or we have an existing leaf to work with
						BVH8LeafNode *leaf;
						if (dataIndex == BVH8IllegalIndex)
							leaf = _bvh8->Append(spatialKey, dataIndex, keyIndex);
						else
							leaf = &_bvh8->m_data[dataIndex];

						// Append new primitive to the cell only if we're not going overboard
						if (leaf->m_numTriangles < MAX_NODE_TRIS)
						{
							// Expand current min/max
							SVec128 exMin = EVecMin(_bvh8->m_dataLookup[keyIndex].m_BoundsMin, triMin);
							SVec128 exMax = EVecMax(_bvh8->m_dataLookup[keyIndex].m_BoundsMax, triMax);

							// Clip to cell bounds
							_bvh8->m_dataLookup[keyIndex].m_BoundsMin = EVecMax(cellMin, exMin);
							_bvh8->m_dataLookup[keyIndex].m_BoundsMax = EVecMin(cellMax, exMax);

							// Append new triangle index
							leaf->m_triangleIndices[leaf->m_numTriangles] = i;
							leaf->m_numTriangles++;
						}
					}
				}
			}
		}

		_bvh8->SortAscending(0, _bvh8->m_dataLookup.size());
		_bvh8->GenerateBVH8();
		_bvh8->SaveBVH8(SCENENAME ".cache.bv8");
	}
}

void Barycentrics(SVec128& P, SVec128& v1, SVec128& v2, SVec128& v3, SVec128& uvw)
{
    SVec128 e1 = EVecSub(v3, v1);
    SVec128 e2 = EVecSub(v2, v1);
    SVec128 e = EVecSub(P, v1);
    float d00 = EVecGetFloatX(EVecDot3(e1, e1));
    float d01 = EVecGetFloatX(EVecDot3(e1, e2));
    float d11 = EVecGetFloatX(EVecDot3(e2, e2));
    float d20 = EVecGetFloatX(EVecDot3(e, e1));
    float d21 = EVecGetFloatX(EVecDot3(e, e2));
    float invdenom = 1.0 / (d00 * d11 - d01 * d01);

	float bu = (d11 * d20 - d01 * d21) * invdenom;
	float bv = (d00 * d21 - d01 * d20) * invdenom;
	float bw = 1.f - bu - bv;

	uvw = EVecConst(bu, bv, bw, 0.f);
}

float TriHit(SVec128& origin, SVec128& direction, SVec128& v1, SVec128& v2, SVec128& v3, float max_t)
{
    SVec128 e1 = EVecSub(v3, v1);
    SVec128 e2 = EVecSub(v2, v1);
    SVec128 s1 = EVecCross3(direction, e2);
	SVec128 K = EVecDot3(s1, e1);

#if defined(DOUBLE_SIDED)
	// No facing check in this case
#else
	if (EVecGetFloatX(K) >= 0.f)
		return max_t; // Ignore backfacing (TODO: enable/disable this)
#endif

    SVec128 invd = EVecRcp(K);
    SVec128 d = EVecSub(origin, v1);
    SVec128 b1 = EVecMul(EVecDot3(d, s1), invd);
    SVec128 s2 = EVecCross3(d, e1);
    SVec128 b2 = EVecMul(EVecDot3(direction, s2), invd);
    SVec128 temp = EVecMul(EVecDot3(e2, s2), invd);

	float fb1 = EVecGetFloatX(b1);
	float fb2 = EVecGetFloatX(b2);
	float ftemp = EVecGetFloatX(temp);

    if (fb1 < 0.f || fb1 > 1.f ||
        fb2 < 0.f || fb1 + fb2 > 1.f ||
        ftemp < 0.f || ftemp > max_t)
		return max_t; // Missed
    else
		return ftemp; // Hit
}

bool SlabTestFast(SVec128& p0, SVec128& p1, SVec128& rayOrigin, SVec128& rayDir, SVec128& invRayDir)
{
	SVec128 t0 = EVecMul(EVecSub(p0, rayOrigin), invRayDir);
	SVec128 t1 = EVecMul(EVecSub(p1, rayOrigin), invRayDir);
	SVec128 tmin = EVecMin(t0, t1);
	SVec128 tmax = EVecMax(t0, t1);
	float enter = EVecMaxComponent3(tmin);
	float exit = EVecMinComponent3(tmax);
	return enter <= exit;
}

bool SlabTest(SVec128& p0, SVec128& p1, SVec128& rayOrigin, SVec128& rayDir, SVec128& invRayDir, SVec128& isect, SVec128& exitpos)
{
	SVec128 t0 = EVecMul(EVecSub(p0, rayOrigin), invRayDir);
	SVec128 t1 = EVecMul(EVecSub(p1, rayOrigin), invRayDir);
	SVec128 tmin = EVecMin(t0, t1);
	SVec128 tmax = EVecMax(t0, t1);
	float enter = EVecMaxComponent3(tmin);
	float exit = EVecMinComponent3(tmax);
	isect = EVecAdd(rayOrigin, EVecMul(EVecConst(enter,enter,enter,0.f), rayDir));
	exitpos = EVecAdd(rayOrigin, EVecMul(EVecConst(exit,exit,exit,0.f), rayDir));
	return enter <= exit;
}

uint32_t GatherChildNodes(SBVH8Database<BVH8LeafNode>* bvh, uint32_t rootNode, SVec128 startPos, uint32_t rayOctantMask, SVec128 deltaVec, SVec128 invDeltaVec, uint32_t* hitcells)
{
	uint32_t octantallocationmask = 0;

	// Generate octant allocation mask and store items in hit order (search on caller side will be in reverse bit scan order, i.e. from MSB to LSB)
	uint32_t childcount = bvh->m_dataLookup[rootNode].m_ChildCount;
	uint32_t firstChildNode = bvh->m_dataLookup[rootNode].m_DataIndex;
	for (uint32_t i = 0; i < childcount; ++i)
	{
		uint32_t idx = firstChildNode + i;
		SVec128 subminbounds = bvh->m_dataLookup[idx].m_BoundsMin;
		SVec128 submaxbounds = bvh->m_dataLookup[idx].m_BoundsMax;
		if (SlabTestFast(subminbounds, submaxbounds, startPos, deltaVec, invDeltaVec))
		{
			uint32_t byoctant = (bvh->m_dataLookup[idx].m_SpatialKey & 0x00000007) ^ rayOctantMask;
			hitcells[byoctant] = idx;
			octantallocationmask |= (1 << byoctant);
		}
	}

	return octantallocationmask;
}

EInline int32_t BitCountLeadingZeroes32(uint32_t x)
{
	#if defined(PLATFORM_WINDOWS)
		return __lzcnt(x);
	#else
		return x ? __builtin_clz(x) : 32;
	#endif
}

int traceBVH8NoTris(SBVH8Database<BVH8LeafNode>* bvh, uint32_t& marchCount, float& t, SVec128& startPos, uint32_t& hitID, SVec128& deltaVec, SVec128& invDeltaVec, SVec128& hitPos)
{
	uint32_t traversalStack[MAX_STACK_ENTRIES]{};
	int stackpointer = 0;

	// Store root node to start travelsal with
	traversalStack[stackpointer++] = bvh->m_LodStart[bvh->m_RootBVH8Node];

	// Generate ray octant mask
	uint32_t ray_octant = (EVecGetFloatX(deltaVec)<0.f ? 1:0) | (EVecGetFloatY(deltaVec)<0.f ? 2:0) | (EVecGetFloatZ(deltaVec)<0.f ? 4:0);

	hitID = 0xFFFFFFFF; // Miss
	hitPos = EVecAdd(startPos, deltaVec); // End of ray
	t = 1.f;

	// Trace until stack underflows
	while (stackpointer > 0)
	{
		// Number of loops through until we find a hit
		++marchCount;

		--stackpointer;
		uint32_t currentNode = traversalStack[stackpointer];

		if (bvh->m_dataLookup[currentNode].m_ChildCount != 0)
		{
			uint32_t hitcells[8];
			uint32_t octantallocationmask = GatherChildNodes(bvh, currentNode, startPos, ray_octant, deltaVec, invDeltaVec, hitcells);

			uint32_t idx = BitCountLeadingZeroes32(octantallocationmask);
			while (idx != 32)
			{
				// Convert to cell index
				uint32_t cellindex = 31 - idx;

				// Stack overflow
				if (stackpointer > MAX_STACK_ENTRIES)
					return -1;

				// Push valid child nodes onto stack
				traversalStack[stackpointer++] = hitcells[cellindex];

				// Clear this bit for next iteration and get new index
				octantallocationmask ^= (1 << cellindex);
				idx = BitCountLeadingZeroes32(octantallocationmask);
			}
		}
		else // Leaf node reached
		{
			// Time to invoke a 'hit test' callback and stop if we have an actual hit
			// It is up to the hit test callback to determine which primitive in this cell is closest etc
			/*if (leafNodeHitTest(currentNode, m_dataLookup[currentNode].m_DataIndex, this, ray, hit))
				return 1;*/

			SVec128 subminbounds = bvh->m_dataLookup[currentNode].m_BoundsMin;
			SVec128 submaxbounds = bvh->m_dataLookup[currentNode].m_BoundsMax;
			SVec128 exitpos;
			hitID=currentNode;
			SlabTest(subminbounds, submaxbounds, startPos, deltaVec, invDeltaVec, hitPos, exitpos);
			t = EVecGetFloatX(EVecLen3(EVecSub(hitPos, startPos)));
		#ifndef XRAY_MODE
			return 1; // NOTE: do not return to generate an x-ray view
		#endif
		}
	}

	return 0;
}

int traceBVH8(SBVH8Database<BVH8LeafNode>* bvh, uint32_t& marchCount, float& t, SVec128& startPos, uint32_t& hitID, SVec128& deltaVec, SVec128& invDeltaVec, SVec128& hitPos)
{
	uint32_t traversalStack[MAX_STACK_ENTRIES]{};
	int stackpointer = 0;

	// Store root node to start travelsal with
	traversalStack[stackpointer++] = bvh->m_LodStart[bvh->m_RootBVH8Node];

	// Generate ray octant mask
	uint32_t ray_octant = (EVecGetFloatX(deltaVec)<0.f ? 1:0) | (EVecGetFloatY(deltaVec)<0.f ? 2:0) | (EVecGetFloatZ(deltaVec)<0.f ? 4:0);

	hitID = 0xFFFFFFFF; // Miss
	hitPos = EVecAdd(startPos, deltaVec); // End of ray
	t = 1.f;

	// Trace until stack underflows
	while (stackpointer > 0)
	{
		// Number of loops through until we find a hit
		++marchCount;

		--stackpointer;
		uint32_t currentNode = traversalStack[stackpointer];

		if (bvh->m_dataLookup[currentNode].m_ChildCount != 0)
		{
			uint32_t hitcells[8];
			uint32_t octantallocationmask = GatherChildNodes(bvh, currentNode, startPos, ray_octant, deltaVec, invDeltaVec, hitcells);

			uint32_t idx = BitCountLeadingZeroes32(octantallocationmask);
			while (idx != 32)
			{
				// Convert to cell index
				uint32_t cellindex = 31 - idx;

				// Stack overflow
				if (stackpointer > MAX_STACK_ENTRIES)
					return -1;

				// Push valid child nodes onto stack
				traversalStack[stackpointer++] = hitcells[cellindex];

				// Clear this bit for next iteration and get new index
				octantallocationmask ^= (1 << cellindex);
				idx = BitCountLeadingZeroes32(octantallocationmask);
			}
		}
		else // Leaf node reached
		{
			// Time to invoke a 'hit test' callback and stop if we have an actual hit
			// It is up to the hit test callback to determine which primitive in this cell is closest etc
			/*if (leafNodeHitTest(currentNode, m_dataLookup[currentNode].m_DataIndex, this, ray, hit))
				return 1;*/

			// Only hit the BVH8 leaf node with no data access
#ifdef IGNORE_CHILD_DATA
			SVec128 subminbounds = bvh->m_dataLookup[currentNode].m_BoundsMin;
			SVec128 submaxbounds = bvh->m_dataLookup[currentNode].m_BoundsMax;
			SVec128 exitpos;
			hitID=currentNode;
			SlabTest(subminbounds, submaxbounds, startPos, deltaVec, invDeltaVec, hitPos, exitpos);
			t = EVecGetFloatX(EVecLen3(EVecSub(hitPos, startPos)));
			#ifndef XRAY_MODE
			return 1; // NOTE: do not return to generate an x-ray view
			#endif
#else

			// Default inline hit test returning hit position
			uint32_t modelNode = bvh->m_dataLookup[currentNode].m_DataIndex;
			float last_t = 1.f;
			uint32_t hitTriangleIndex = 0xFFFFFFFF;
			for (uint32_t tri=0; tri<bvh->m_data[modelNode].m_numTriangles; ++tri)
			{
				uint32_t triangleIndex = bvh->m_data[modelNode].m_triangleIndices[tri];

				SVec128 v1 = testtris[triangleIndex].coords[0];
				SVec128 v2 = testtris[triangleIndex].coords[1];
				SVec128 v3 = testtris[triangleIndex].coords[2];

				// Hit if closer than before
				float curr_t = TriHit(startPos, deltaVec, v1, v2, v3, last_t);
				if (curr_t < last_t)
				{
					last_t = curr_t;
					hitTriangleIndex = triangleIndex;
				}
			}

			if (hitTriangleIndex == 0xFFFFFFFF) // Nothing hit
				continue;

			// TODO: hit position, hit normal etc

			t = last_t;
			hitID = hitTriangleIndex;
			hitPos = EVecAdd(startPos, EVecMul(EVecConst(t, t, t, 0.f), deltaVec));
			#ifndef XRAY_MODE
			// Positive hit
			return 1;
			#endif
#endif
		}
	}

	return 0;
}

static int DispatcherThread(void *data)
{
	SWorkerContext *vec = (SWorkerContext *)(data);
	while(!g_done)
	{
		uint32_t workItem = 0xFFFFFFFF;
		if (vec->dispatchvector.Read(&workItem, sizeof(uint32_t)))
		{
			// x/y tile indices
			uint32_t tx = workItem%tilecountx;
			uint32_t ty = workItem/tilecountx;

			// Upper left pixel position
			uint32_t ox = tx*tilewidth;
			uint32_t oy = ty*tileheight;

			uint8_t *p = vec->rasterTile;
			for (uint32_t iy = 0; iy<tileheight; ++iy)
			{
				float py = vec->rc->aspect * (float(height)/2.f-float(iy+oy))/float(height);
				for (uint32_t ix = 0; ix<tilewidth; ++ix)
				{
					float px = (float(ix+ox) - float(width)/2.f)/float(width);

					float t = raylength;

					SVec128 pyVec{py,py,py,0.f};
					SVec128 U = EVecMul(vec->rc->lookMat.r[1], pyVec);
					SVec128 raylenvec{t,t,t,0.f};
					SVec128 pxVec{px,px,px,0.f};
					SVec128 L = EVecMul(vec->rc->lookMat.r[0], pxVec);
					SVec128 traceRay = EVecMul(EVecAdd(EVecAdd(L, U), vec->rc->F), raylenvec);
					SVec128 invRay = EVecRcp(traceRay);

					uint32_t hitID = 0xFFFFFFFF;
					SVec128 hitpos = EVecAdd(vec->rc->rayOrigin, traceRay);
					uint32_t marchCount = 0;

					traceBVH8(testBVH8, marchCount, t, vec->rc->rayOrigin, hitID, traceRay, invRay, hitpos);

					float final = 0.f;

					if (hitID != 0xFFFFFFFF)
					{
						SVec128 sunPos{20.f,35.f,20.f,1.f};
						SVec128 sunRay = EVecSub(sunPos, hitpos);
						SVec128 invSunRay = EVecRcp(sunRay);
						SVec128 nrm;

#ifdef SHOW_DEPTH
						// Depth
						final = t;
#else
						SVec128 viewRay = EVecNorm3(traceRay);

#if defined(ENABLE_EXPFOG)
						final = 0.f;
						int through = 0;
						for (uint32_t i=0; i<4; ++i)
						{
							// Advance forward
							SVec128 marchpos = EVecAdd(hitpos, EVecMul(viewRay, vec->rc->epsilon));

							SVec128 exitPos;
							float t2 = raylength;
							hitID = 0xFFFFFFFF;
							float rx = float(rand()%100-50);
							float ry = float(rand()%100-50);
							float rz = float(rand()%100-50);
							SVec128 rray = EVecConst(rx, ry, rz, 0.f);
							SVec128 randomRay = EVecMul(EVecNorm3(rray), raylenvec);
							SVec128 invRandomRay = EVecRcp(randomRay);
							traceBVH8(testBVH8, marchCount, t2, marchpos, hitID, traceRay, invRay, exitPos);

							if (hitID == 0xFFFFFFFF)
								break;

							if ((through%2)==0)
							{
								float L = EVecGetFloatX(EVecLen3(EVecSub(exitPos, hitpos)));
								float F = EMaximum(0.f, 1.f - expf(-L*0.04f));
								final += F;
							}
							hitpos = exitPos;
							through++;
						}
						final *= 2.f;
#else

#if defined(ENABLE_OCCLUSION)
						{
							SVec128 uvw;
							Barycentrics(hitpos,
								testtris[hitID].coords[0],
								testtris[hitID].coords[1],
								testtris[hitID].coords[2], uvw);
							SVec128 uvwx = EVecSplatX(uvw);
							SVec128 uvwy = EVecSplatY(uvw);
							SVec128 uvwz = EVecSplatZ(uvw);
							SVec128 uvwzA = EVecMul(uvwz, testtris[hitID].normals[0]); // A*uvw.zzz
							SVec128 uvwyB = EVecMul(uvwy, testtris[hitID].normals[1]); // B*uvw.yyy
							SVec128 uvwxC = EVecMul(uvwx, testtris[hitID].normals[2]); // C*uvw.xxx
							nrm = EVecNorm3(EVecAdd(uvwzA, EVecAdd(uvwyB, uvwxC))); // A*uvw.zzz + B*uvw.yyy + C*uvw.xxx
						}

						hitpos = EVecAdd(EVecMul(nrm, EVecConst(1.2f,1.2f,1.2f,0.f)), hitpos);

						float occsum = 0.f;
						for(uint32_t i=0;i<8;++i)
						{
							SVec128 occHitPos;
							float t2 = raylength;
							hitID = 0xFFFFFFFF;
							float rx = float(rand()%100-50);
							float ry = float(rand()%100-50);
							float rz = float(rand()%100-50);
							SVec128 rray = EVecConst(rx, ry, rz, 0.f);
							SVec128 randomRay = EVecMul(EVecNorm3(rray), raylenvec);
							SVec128 invRandomRay = EVecRcp(randomRay);
							occsum += float(traceBVH8NoTris(testBVH8, marchCount, t2, hitpos, hitID, randomRay, invRandomRay, occHitPos));
						}
						occsum /= 8.f;
						occsum = 1.f-occsum;
						final  = occsum;
#else
						// Global + NdotL
						{
							SVec128 uvw;
							Barycentrics(hitpos,
								testtris[hitID].coords[0],
								testtris[hitID].coords[1],
								testtris[hitID].coords[2], uvw);
							SVec128 uvwx = EVecSplatX(uvw);
							SVec128 uvwy = EVecSplatY(uvw);
							SVec128 uvwz = EVecSplatZ(uvw);
							SVec128 uvwzA = EVecMul(uvwz, testtris[hitID].normals[0]); // A*uvw.zzz
							SVec128 uvwyB = EVecMul(uvwy, testtris[hitID].normals[1]); // B*uvw.yyy
							SVec128 uvwxC = EVecMul(uvwx, testtris[hitID].normals[2]); // C*uvw.xxx
							nrm = EVecNorm3(EVecAdd(uvwzA, EVecAdd(uvwyB, uvwxC))); // A*uvw.zzz + B*uvw.yyy + C*uvw.xxx
							float L = fabs(EVecGetFloatX(EVecDot3(nrm, EVecNorm3(sunRay))));
							final += L;
						}

						// Hit position bias/offset
						//hitpos = EVecAdd(hitpos, EVecMul(viewRay, vec->rc->negepsilon));
						hitpos = EVecAdd(EVecMul(nrm, vec->rc->epsilon), hitpos);
#endif // ENABLE_OCCLUSION
#endif // ENABLE_EXPFOG

						// Reflections
#if defined(ENABLE_REFLECTIONS)
						{
							float tr = raylength;
							hitID = 0xFFFFFFFF;
							SVec128 reflHitPos;
							//r = viewRay-2*dot(viewRay, n)*n;
							SVec128 two{2.f,2.f,2.f,0.f};
							//SVec128 negone{-1.f,-1.f,-1.f,0.f};
							SVec128 negViewRay = viewRay;//EVecMul(negone, viewRay);
							SVec128 reflRay = EVecMul(EVecSub(negViewRay, EVecMul(EVecDot3(negViewRay, nrm), EVecMul(two, nrm))), raylenvec);
							SVec128 invReflRay = EVecRcp(reflRay);
							traceBVH8(testBVH8, marchCount, tr, hitpos, hitID, reflRay, invReflRay, reflHitPos);
							if (hitID != 0xFFFFFFFF)
							{
								SVec128 uvw;
								Barycentrics(reflHitPos,
									testtris[hitID].coords[0],
									testtris[hitID].coords[1],
									testtris[hitID].coords[2], uvw);
								SVec128 uvwx = EVecSplatX(uvw);
								SVec128 uvwy = EVecSplatY(uvw);
								SVec128 uvwz = EVecSplatZ(uvw);
								SVec128 uvwzA = EVecMul(uvwz, testtris[hitID].normals[0]); // A*uvw.zzz
								SVec128 uvwyB = EVecMul(uvwy, testtris[hitID].normals[1]); // B*uvw.yyy
								SVec128 uvwxC = EVecMul(uvwx, testtris[hitID].normals[2]); // C*uvw.xxx
								SVec128 rnrm = EVecNorm3(EVecAdd(uvwzA, EVecAdd(uvwyB, uvwxC))); // A*uvw.zzz + B*uvw.yyy + C*uvw.xxx
								float rL = fabs(EVecGetFloatX(EVecDot3(rnrm, EVecNorm3(sunRay))));
								float diminish = fabs(1.f/(2.f*raylength*tr+0.01f));
								diminish = EMinimum(1.f, EMaximum(0.f, diminish));
								final += EMinimum(final, rL*diminish);
							}
						}
#endif // ENABLE_REFLECTIONS

						// Shadow
#if defined(ENABLE_SHADOWS)
						{
							float t2 = raylength;
							hitID = 0xFFFFFFFF;
							SVec128 shadowHitPos;
							traceBVH8(testBVH8, marchCount, t2, hitpos, hitID, sunRay, invSunRay, shadowHitPos);
							if (t2<1.f)
								final *= 0.65f;
						}
#endif // ENABLE_SHADOWS
#endif // SHOW_DEPTH
					}

#if defined(SHOW_HEATMAP)
					uint8_t C = marchCount*4;
#else
					uint8_t C = uint8_t(final*255.f);
#endif // SHOW_HEATMAP

#if defined(SHOW_WORKER_IDS)
					p[(ix+iy*tilewidth)*4+0] = (vec->workerID&4) ? 128:C; // B
					p[(ix+iy*tilewidth)*4+1] = (vec->workerID&2) ? 128:C; // G
					p[(ix+iy*tilewidth)*4+2] = (vec->workerID&1) ? 128:C; // R
					p[(ix+iy*tilewidth)*4+3] = 0xFF;
#else
					p[(ix+iy*tilewidth)*4+0] = C; // B
					p[(ix+iy*tilewidth)*4+1] = C; // G
					p[(ix+iy*tilewidth)*4+2] = C; // R
					p[(ix+iy*tilewidth)*4+3] = 0xFF;
#endif // SHOW_WORKER_IDS
				}
			}

			// Copy to final surface
			// In hardware, this is 32 bytes of data, so it's much less than a
			// cache line size (which is 512 bits, i.e. 64 bytes)
			for (uint32_t h=0;h<tileheight;++h)
				memcpy(&vec->rc->pixels[4*(ox + (oy+h)*width)], &vec->rasterTile[h*tilewidth*4], tilewidth*4);
		}
	}

	return 0;
}

#if defined(USE_MORTON_ORDER)
void EMorton2DDecode(const uint32_t morton, uint32_t &x, uint32_t &y)
{
  uint32_t res = morton&0x5555555555555555;
  res=(res|(res>>1)) & 0x3333333333333333;
  res=(res|(res>>2)) & 0x0f0f0f0f0f0f0f0f;
  res=(res|(res>>4)) & 0x00ff00ff00ff00ff;
  res=res|(res>>8);
  x = res;
  res = (morton>>1)&0x5555555555555555;
  res=(res|(res>>1)) & 0x3333333333333333;
  res=(res|(res>>2)) & 0x0f0f0f0f0f0f0f0f;
  res=(res|(res>>4)) & 0x00ff00ff00ff00ff;
  res=res|(res>>8);
  y = res;
}
#endif

#if defined(PLATFORM_LINUX)
int main(int _argc, char** _argv)
#else
int SDL_main(int _argc, char** _argv)
#endif
{
#if defined(PLATFORM_LINUX)
	// Console always there for Linux
#else
	HWND hWnd = GetConsoleWindow();
	ShowWindow( hWnd, SW_SHOW ); // In case it's not shown at startup
#endif

	objl::Loader objloader;
#if defined(PLATFORM_WINDOWS)
	if (!objloader.LoadFile(SCENENAME "\\" SCENENAME ".obj"))
#else
	if (!objloader.LoadFile(SCENENAME "/" SCENENAME ".obj"))
#endif
	{
		printf("Failed to load OBJ file\n");
		return 1;
	}

	// Set up triangle data
	int t=0;
	int totaltriangles = 0;

	for (auto &mesh : objloader.LoadedMeshes)
		totaltriangles += int(mesh.Indices.size()/3);

	if (totaltriangles==0)
	{
		printf("No triangles in model\n");
		return 1;
	}

	testtris = new triangle[totaltriangles];

	//auto &mesh = objloader.LoadedMeshes[1]; // Just the floor plane..
	for (auto &mesh : objloader.LoadedMeshes) // ..or, the entire scene
	{
		int triCount = mesh.Indices.size()/3;
		for (int i=0;i<triCount;++i)
		{
			int tri = i*3;
			unsigned int i0 = mesh.Indices[tri+0];
			unsigned int i1 = mesh.Indices[tri+1];
			unsigned int i2 = mesh.Indices[tri+2];

			testtris[t].coords[0] = EVecConst( mesh.Vertices[i0].Position.X, mesh.Vertices[i0].Position.Y, mesh.Vertices[i0].Position.Z, 0.f);
			testtris[t].coords[1] = EVecConst( mesh.Vertices[i1].Position.X, mesh.Vertices[i1].Position.Y, mesh.Vertices[i1].Position.Z, 0.f);
			testtris[t].coords[2] = EVecConst( mesh.Vertices[i2].Position.X, mesh.Vertices[i2].Position.Y, mesh.Vertices[i2].Position.Z, 0.f);

			testtris[t].normals[0] = EVecConst( mesh.Vertices[i0].Normal.X, mesh.Vertices[i0].Normal.Y, mesh.Vertices[i0].Normal.Z, 0.f);
			testtris[t].normals[1] = EVecConst( mesh.Vertices[i1].Normal.X, mesh.Vertices[i1].Normal.Y, mesh.Vertices[i1].Normal.Z, 0.f);
			testtris[t].normals[2] = EVecConst( mesh.Vertices[i2].Normal.X, mesh.Vertices[i2].Normal.Y, mesh.Vertices[i2].Normal.Z, 0.f);

			++t;
		}
	}

	// BVH8
	testBVH8 = new SBVH8Database<BVH8LeafNode>;
	bvh8Builder(testtris, totaltriangles, testBVH8);

	// Trace the BVH
	if(SDL_Init(SDL_INIT_VIDEO) != 0)
	{
		fprintf(stderr, "Could not init SDL2: %s\n", SDL_GetError());
		return 1;
	}

    SDL_Window *screen = SDL_CreateWindow("BVH8Tool",
            SDL_WINDOWPOS_UNDEFINED,
            SDL_WINDOWPOS_UNDEFINED,
            width, height, 0);

    if(!screen)
	{
        fprintf(stderr, "Could not create window\n");
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(screen, -1, SDL_RENDERER_SOFTWARE);
    if(!renderer)
	{
        fprintf(stderr, "Could not create renderer\n");
        return 1;
    }

	SDL_Surface *surface = SDL_GetWindowSurface(screen);

    SDL_SetRenderDrawColor(renderer, 32, 32, 128, 255);
    SDL_RenderClear(renderer);
    SDL_RenderPresent(renderer);

	bool done = false;

	SRenderContext rc;
	SWorkerContext wc[MAX_WORKERS];
	SDL_Thread *thrd[MAX_WORKERS];
	for (uint32_t i=0; i<MAX_WORKERS; ++i)
	{
		wc[i].workerID = i;
		wc[i].rc = &rc;
		thrd[i] = SDL_CreateThread(DispatcherThread, "DispatcherThread", (void*)&wc[i]);
	}

	rc.rotAng = 0.f;
	rc.aspect = float(height) / float(width);
	rc.nil = SVec128{0.f, 0.f, 0.f, 0.f};
	rc.negepsilon = SVec128{-0.02f, -0.02f, -0.02f, 0.f};
	rc.epsilon = SVec128{0.02f, 0.02f, 0.02f, 0.f};

	do{
		SDL_Event event;
		while(SDL_PollEvent(&event))
		{
			if(event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE))
				done = true;
		}

		if (SDL_MUSTLOCK(surface))
			SDL_LockSurface(surface);

		// Set up camera data
		rc.rayOrigin = SVec128{sinf(rc.rotAng)*cameradistance, (1.f+sinf(rc.rotAng*0.5f))*cameradistance*0.5f, cosf(rc.rotAng)*cameradistance, 1.f};
		rc.lookAt = SVec128{0.f,0.f,0.f,1.f};
		rc.upVec = SVec128{0.f,1.f,0.f,0.f};
		rc.lookMat = EMatLookAtRightHanded(rc.rayOrigin, rc.lookAt, rc.upVec);
		rc.pzVec = SVec128{-1.f,-1.f,-1.f,0.f};
		rc.F = EVecMul(rc.lookMat.r[2], rc.pzVec);
		rc.pixels = (uint8_t*)surface->pixels;

		int distributedAll = 0;
		uint32_t workunit = 0;
		do {
			// Distribute all tiles across all work queues
			for (uint32_t i=0; i<MAX_WORKERS; ++i)
			{
				if (wc[i].dispatchvector.FreeSpace()!=0) // We have space in this worker's queue
				{
#if defined(USE_MORTON_ORDER)
					uint32_t x, y;
					EMorton2DDecode(workunit, x, y);
					uint32_t mortonindex = x+y*tilecountx;
					if (mortonindex < tilecountx*tilecounty) // Ran out of tiles yet?
						wc[i].dispatchvector.Write(&mortonindex, sizeof(uint32_t));
					else
						distributedAll = 1; // Done with all tiles
#else
					if (workunit < tilecountx*tilecounty) // Ran out of tiles yet?
						wc[i].dispatchvector.Write(&workunit, sizeof(uint32_t));
					else
						distributedAll = 1; // Done with all tiles
#endif
					// Next tile
					workunit++;
				}
			}
		} while (!distributedAll); // We're done handing out jobs

		// Rotate
		rc.rotAng += 0.01f;

		// Wait for all threads to be done with locked image pointer before updating window image
		/*int tdone;
		do
		{
			tdone = 0;
			for (uint32_t i=0; i<MAX_WORKERS; ++i)
				tdone += wc[i].dispatchvector.BytesAvailable() ? 0 : 1;
		} while(tdone != MAX_WORKERS);*/

		// TODO: Copy tiles to their respective positions
		//(uint8_t*)surface->pixels;

		if (SDL_MUSTLOCK(surface))
			SDL_UnlockSurface(surface);
		SDL_UpdateWindowSurface(screen);

	} while (!done);

	delete [] testtris;
	delete testBVH8;

	g_done = true;

	for (uint32_t i=0; i<MAX_WORKERS; ++i)
	{
		int threadReturnValue;
		SDL_WaitThread(thrd[i], &threadReturnValue);
	}

	// Done
	SDL_FreeSurface(surface);
    SDL_DestroyWindow(screen);
	SDL_Quit();

	return 0;
}
