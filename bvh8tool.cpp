#if defined(PLATFORM_LINUX)
#include "platformlinux.h"
#include "SDL2/SDL.h"
#else // PLATFORM_WINDOWS
#include "platformwin.h"
#include "SDL.h"
#endif

#include <stdio.h>
#include "math.h"
#include "bvh8.h"
#include "objloader.h"

// Fog prototype
//#define FOG_WORK

// Size of each cell in world units
#define NODE_DIMENSION 0.6f
// NOTE: Once there's a correct cell vs triangle test, this will be reduced
#define MAX_NODE_TRIS 40
// Depth of traversal stack
#define MAX_STACK_ENTRIES 16

// Define this to ignore triangles and work with BVH child nodes only, and see in x-ray vision
//#define IGNORE_CHILD_DATA

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

uint32_t width = 640;
uint32_t height = 480;

uint8_t* pixels;
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

	// Set up grid bounds and cell scale
	_bvh8->m_GridAABBMin = sceneMin;
	_bvh8->m_GridCellSize = SVec128{NODE_DIMENSION, NODE_DIMENSION, NODE_DIMENSION, 0.f};

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

#ifdef FOG_WORK
	// Fog requires double-sided hits
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
			//return 1; // NOTE: no 'hit' returned since we want x-ray vision
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
			// Positive hit
			return 1;
#endif
		}
	}

	return 0;
}

void block(int x, int y, uint8_t B, uint8_t G, uint8_t R)
{
	for (int oy=y;oy<y+2;++oy)
	for (int ox=x;ox<x+2;++ox)
	{
		pixels[(ox+oy*width)*4+0] = B;
		pixels[(ox+oy*width)*4+1] = G;
		pixels[(ox+oy*width)*4+2] = R;
		pixels[(ox+oy*width)*4+3] = 0xFF;
	}
}

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
	ShowWindow( hWnd, SW_HIDE );
#endif

	objl::Loader objloader;
	if (!objloader.LoadFile("test.obj"))
	{
		printf("Failed to load test.obj\n");
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
	const float cameradistance = 30.f;
	uint32_t lowTraces = 0xFFFFFFFF;
	uint32_t highTraces = 0x00000000;

	do{
		SDL_Event event;
		while(SDL_PollEvent(&event))
		{
			if(event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE))
				done = true;
		}

		if (SDL_MUSTLOCK(surface))
			SDL_LockSurface(surface);
		
		static float rotAng = 0.f;
		float aspect = float(height) / float(width);

		SVec128 rayOrigin{sinf(rotAng)*cameradistance, (1.f+sinf(rotAng*0.5f))*cameradistance*0.5f, cosf(rotAng)*cameradistance, 1.f};
		SVec128 lookAt{0.f,0.f,0.f,1.f};
		SVec128 upVec{0.f,1.f,0.f,0.f};
		SMatrix4x4 lookMat = EMatLookAtRightHanded(rayOrigin, lookAt, upVec);
		SVec128 pzVec{-1.f,-1.f,-1.f,0.f};
		SVec128 F = EVecMul(lookMat.r[2], pzVec);

		pixels = (uint8_t*)surface->pixels;
		uint32_t maxTraces = 0;
		SVec128 nil{0.f, 0.f, 0.f, 0.f};
		SVec128 epsilon{-0.02f, -0.02f, -0.02f, 0.f};
		static int EVENODD = 0;
		EVENODD = (EVENODD+1)%4;
		// Checkerboard
		for (int y=(EVENODD%2)*2; y<height; y+=4)
		{
			float py = aspect * (float(height)/2.f-float(y))/float(height);
			SVec128 pyVec{py,py,py,0.f};
			SVec128 U = EVecMul(lookMat.r[1], pyVec);
			for (int x=(EVENODD/2)*2; x<width; x+=4)
			{
				float t = cameradistance*2.f;
				SVec128 raylenvec{t,t,t,0.f};

				// Rotating camera
				float px = (float(x) - float(width)/2.f)/float(width);
				SVec128 pxVec{px,px,px,0.f};
				SVec128 L = EVecMul(lookMat.r[0], pxVec);
				//SVec128 traceRay = EVecAdd(EVecAdd(L, U), F);
				SVec128 traceRay = EVecMul(EVecAdd(EVecAdd(L, U), F), raylenvec);
				SVec128 invRay = EVecRcp(traceRay);

				// Scene
				uint32_t marchCount = 0;
				uint32_t hitID = 0xFFFFFFFF;
				SVec128 hitpos = rayOrigin;

				float final = 0.f;
				SVec128 nrm;
				int color = 0;

#ifdef FOG_WORK
				// Mesh fog attempt
				int intersectioncount = 0;
				SVec128 forwardepsilon{0.02f, 0.02f, 0.02f, 0.f};
				SVec128 viewRay = EVecNorm3(traceRay);
				float F = 0.001f;
				for (int i=0; i<12; ++i)
				{
					// Find ray exit
					t = 512.f;
					SVec128 exitpos;
					traceBVH8(testBVH8, marchCount, t, hitpos, hitID, traceRay, invRay, exitpos);

					if (hitID != 0xFFFFFFFF)
					{
						if (intersectioncount%2==1)
							F += 1.f-expf(-EVecGetFloatX(EVecLen3(EVecSub(hitpos, exitpos))));
						hitpos = EVecAdd(exitpos, EVecMul(viewRay, forwardepsilon));
						++intersectioncount;
					}
					else
					{
						F += 0.01f;
						break;
					}
				}
				final += 1.f-expf(-F*F*1.f);
				color = int(final*255.f);

#else

				traceBVH8(testBVH8, marchCount, t, rayOrigin, hitID, traceRay, invRay, hitpos);
				maxTraces = EMaximum(maxTraces, marchCount);

				/*if (hitID!=0xFFFFFFFF) // HAVE_HIT: Ray hit a primitive in a leaf nodes
				{
					// Hit position
					//hitpos = EVecAdd(rayOrigin, EVecMul(traceRay,  EVecConst(t-0.01f,t-0.01f,t-0.01f,1.f)));
					// Barycentric coortinates for attribute interpolation
					//SVec128 uvw;
					//Barycentrics(hitpos, testtris[hitID].coords[0], testtris[hitID].coords[1], testtris[hitID].coords[2], uvw);
					//block(x,y, ((hitID>>1)%2)*128, ((hitID>>1)%4)*64, ((hitID>>2)%8)*32);
					float D = EVecGetFloatX(EVecLen3(hitpos));
					int C = int(D*16.f);
					block(x,y, C, C, C);
				}
				else
				{
					block(x,y, 0,0,0);
				}*/

				// Shadow, only if we hit some geometry
#ifdef IGNORE_CHILD_DATA
				block(x,y, marchCount*6, marchCount*6, marchCount*6);
#else // !IGNORE_CHILD_DATA
				if (hitID != 0xFFFFFFFF)
				{
					SVec128 sunPos{20.f,35.f,20.f,1.f};
					SVec128 sunRay = EVecSub(sunPos, hitpos);
					SVec128 invSunRay = EVecRcp(sunRay);

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
						final += fabs(EVecGetFloatY(nrm)*0.5f+0.5f); // global
					}

					// Hit position bias/offset
					SVec128 viewRay = EVecNorm3(traceRay);
					hitpos = EVecAdd(hitpos, EVecMul(viewRay, epsilon));

					// Reflections
					float tr = cameradistance*2.f;
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
						float diminish = fabs(1.f/(64.f*tr+0.01f));
						diminish = EMinimum(1.f, EMaximum(0.f, diminish));
						final += rL*diminish;
					}

					// Shadow
					float t2 = cameradistance*2.f;
					hitID = 0xFFFFFFFF;
					SVec128 shadowHitPos;
					traceBVH8(testBVH8, marchCount, t2, hitpos, hitID, sunRay, invSunRay, shadowHitPos);
					if (t2<1.f)
						final *= 0.85f;

					color = int(final*64.f);
				}
#endif // IGNORE_CHILD_DATA
#endif // FOG_WORK
				block(x,y, color, color, color);
			}
		}

		lowTraces = EMinimum(maxTraces, lowTraces);
		highTraces = EMaximum(maxTraces, highTraces);
		printf("MaxTraces: %d Highest: %d Lowest: %d\n", maxTraces, highTraces, lowTraces);

		rotAng += 0.02f;

		if (SDL_MUSTLOCK(surface))
			SDL_UnlockSurface(surface);
		SDL_UpdateWindowSurface(screen);

	} while (!done);

	delete [] testtris;
	delete testBVH8;

	// Done
	SDL_FreeSurface(surface);
    SDL_DestroyWindow(screen);
	SDL_Quit();

	return 0;
}
