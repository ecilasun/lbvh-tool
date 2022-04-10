#include <stdio.h>

#include "bvh8.h"

typedef SBVH8Database<BVH8LeafNode> TBVH8Scene;

// Note: use -exec set output-radix 16 on debug console to set output to hex format

#define VOXEL_DIMENSION 3.f

static const uint32 MAX_CLIPPED_VERTICES = 36;

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

	uint32 slotSelDist = dist < 0.f ? 1 : 0;
	uint32 slotSelPrevDist = _prevDist < 0.f ? 1 : 0;
	uint32 slot1 = (!slotSelDist && slotSelPrevDist) ? _clipIndex + 1 : MAX_CLIPPED_VERTICES - 1;

	SVec128 clippedPos = EVecAdd(_prevPos, EVecMul(EVecSub(_pos, _prevPos), clipRatioVec));
	_vClip[_clipIndex] = (slotSelDist ^ slotSelPrevDist) ? clippedPos : _pos;
	_vClip[slot1] = _pos;

	uint32 stepSel = (slotSelDist && slotSelPrevDist) ? 0 : ((!slotSelDist && slotSelPrevDist) ? 2 : 1);
	_clipIndex += stepSel;

	// Remember where we came from
	_prevPos = _pos;
	_prevDist = dist;
}

void dumpTree(TBVH8Scene &tree, uint32 root, uint32& nodesvisited)
{
	uint32 childcount = tree.m_dataLookup[root].m_ChildCount;
	if (childcount==0)
		printf("%d ", root);
	else
		printf("%d[%d] -> ", root, childcount);
	nodesvisited++;
	for (uint32 c=0; c<childcount; ++c)
		dumpTree(tree, tree.m_dataLookup[root].m_DataIndex+c, nodesvisited);
}

int main(int _argc, char** _argv)
{
	//if (_argc<=1)
	{
		printf("BVH8 Tool v0.0\nUsage:\nbvh8tool inputfile outputfile\n");
		//return 0;
	}

	// Assume scene size = max vertex bounds
	float minx = -1024.f;
	float maxx = 1024.f;
	float miny = -1024.f;
	float maxy = 1024.f;
	float minz = -1024.f;
	float maxz = 1024.f;

	// Test scene with 4 triangles (base and front face of a 3x3x3 cube)
	float vertices[]={
		// XZ plane bottom face
		-3.f,-3.f,-3.f, 3.f,-3.f,-3.f, 3.f,-3.f,3.f,
		3.f,-3.f,3.f, -3.f,-3.f,3.f, -3.f,-3.f,-3.f,
		// XY plane front face
		-3.f,-3.f,3.f, 3.f,-3.f,3.f, 3.f,3.f,3.f,
		3.f,3.f,3.f, -3.f,3.f,3.f, -3.f,-3.f,3.f };

	SVec128 cellScaleVec = EVecConst(VOXEL_DIMENSION, VOXEL_DIMENSION, VOXEL_DIMENSION, 1.f);
	SVec128 aabbrange = EVecDiv(EVecConst(maxx-minx,maxy-miny,maxz-minz,0.f), cellScaleVec);
	SVec128i qmax = EVecConvertFloatToInt(EVecFloor(aabbrange));
	SVec128i dim = EVecAddi(qmax, EVecConsti(1, 1, 1, 0));

	TBVH8Scene spatialdatabase;
	spatialdatabase.SetGridAABB(EVecConst(minx, miny, minz, 1.f), EVecConst(maxx, maxy, maxz, 2.f));
	spatialdatabase.Clear();

	for (uint32 t=0; t<4; ++t)
	{
		float v0x = vertices[t*9+0];
		float v0y = vertices[t*9+1];
		float v0z = vertices[t*9+2];
		float v1x = vertices[t*9+3];
		float v1y = vertices[t*9+4];
		float v1z = vertices[t*9+5];
		float v2x = vertices[t*9+6];
		float v2y = vertices[t*9+7];
		float v2z = vertices[t*9+8];
		SVec128 v0 = EVecConst(v0x, v0y, v0z, 1.f);
		SVec128 v1 = EVecConst(v1x, v1y, v1z, 1.f);
		SVec128 v2 = EVecConst(v2x, v2y, v2z, 1.f);

		SVec128 l0, l1, l2;
		spatialdatabase.ToGridUnits(v0, l0);
		spatialdatabase.ToGridUnits(v1, l1);
		spatialdatabase.ToGridUnits(v2, l2);
		SVec128i q0 = EVecConvertFloatToInt(EVecFloor(l0));
		SVec128i q1 = EVecConvertFloatToInt(EVecFloor(l1));
		SVec128i q2 = EVecConvertFloatToInt(EVecFloor(l2));
		SVec128i maxBound = EVecMaxi(q0, EVecMaxi(q1, q2));
		SVec128i minBound = EVecMini(q0, EVecMini(q1, q2));

		// Scan grid and generate data for cells this polygon crosses
		for (int z = EVecGetIntZ(minBound); z <= EVecGetIntZ(maxBound); ++z)
		{
			for (int y = EVecGetIntY(minBound); y <= EVecGetIntY(maxBound); ++y)
			{
				for (int x = EVecGetIntX(minBound); x <= EVecGetIntX(maxBound); ++x)
				{
					if (x<0 || y<0 || z<0 || x>1023 || y>1023 || z>1023)
						continue;

					SVec128i cellOrdinal = EVecConsti(x, y, z, 0);

					SVec128 planeOrigin = EVecConvertIntToFloat(cellOrdinal);
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

					// Start with vertices of the triangle in grid space
					int inputCount = 4;
					SVec128 vClip[MAX_CLIPPED_VERTICES], pos[MAX_CLIPPED_VERTICES];
					pos[0] = l0;
					pos[1] = l1;
					pos[2] = l2;
					pos[3] = l0;
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

					if (clipSuccess > 5)
					{
						uint32 keyIndex = 0;
						// Instead of triangled centroid, use the cell center which is more stable
						uint32 cellKey = spatialdatabase.EncodeKey(cellOrdinal);
						uint32 dataIndex = spatialdatabase.FindCell(cellKey, keyIndex);
						BVH8LeafNode* tile = nullptr;
						if (dataIndex == BVH8IllegalIndex)
						{
							tile = spatialdatabase.Append(cellKey, dataIndex, keyIndex);
							tile->m_NumTriangles = 0;
						}
						else
						{
							tile = spatialdatabase.GetDataItem(dataIndex);
						}

						// Use clipped triangles to generate cell bounds
						for (int v = 1; v <= clipIndex-3; ++v)
						{
							// Clipped triangle vertices in world space
							SVec128 tv0, tv1, tv2;
							spatialdatabase.ToWorldUnits(pos[0], tv0);
							spatialdatabase.ToWorldUnits(pos[v], tv1);
							spatialdatabase.ToWorldUnits(pos[v+1], tv2);
							// Generate bounds, which will probably poke out of the cell bounds
							spatialdatabase.m_dataLookup[keyIndex].m_BoundsMin = EVecMin(spatialdatabase.m_dataLookup[keyIndex].m_BoundsMin, EVecMin(EVecMin(tv0, tv1), tv2));
							spatialdatabase.m_dataLookup[keyIndex].m_BoundsMax = EVecMax(spatialdatabase.m_dataLookup[keyIndex].m_BoundsMax, EVecMax(EVecMax(tv0, tv1), tv2));

							// Clip the bounds to the cell to avoid overlaps
							SVec128 gridwpos, gridwposplus;
							spatialdatabase.ToWorldUnits(EVecConst(float(x),float(y),float(z),0.f), gridwpos);
							gridwposplus = EVecAdd(gridwpos, spatialdatabase.m_GridCellSize);
							spatialdatabase.m_dataLookup[keyIndex].m_BoundsMin = EVecMax(spatialdatabase.m_dataLookup[keyIndex].m_BoundsMin, gridwpos);
							spatialdatabase.m_dataLookup[keyIndex].m_BoundsMax = EVecMin(spatialdatabase.m_dataLookup[keyIndex].m_BoundsMax, gridwposplus);
						}

						// NOTE: Since we're adding the non-clipped triangle (i.e. one triangle only)
						// we don't need to use the loop above to stash the clipped triangle fan.
						uint32 tricount = tile->m_NumTriangles;
						if (tricount < MAX_TRIS)
						{
							tile->m_Vertices[tricount*3+0] = v0;
							tile->m_Vertices[tricount*3+1] = v1;
							tile->m_Vertices[tricount*3+2] = v2;
							++tricount;
						}
						tile->m_NumTriangles = tricount;
					}
				}
			}
		}

	}

	// We're done, radix sort the end result (ascending)
	spatialdatabase.SortAscending(0, (uint32)spatialdatabase.m_dataLookup.size(), true);

	// Generate the final BVH8
	spatialdatabase.GenerateBVH8();

	uint32 rdi = spatialdatabase.m_LodStart[spatialdatabase.m_RootBVH8Node];
	printf("root node %d points at data index %d\n", spatialdatabase.m_RootBVH8Node, rdi);
	if (rdi == spatialdatabase.m_dataLookup.size()-1)
		printf("OK: correct layout since root node is last entry in table\n");
	else
		printf("ERROR: incorrect layout (last table entry at %d but root node is %d).\n", uint32(spatialdatabase.m_dataLookup.size()-1), rdi);

	if (rdi+1 == spatialdatabase.m_LodEnd[spatialdatabase.m_RootBVH8Node])
		printf("OK: root node spans only one entry\n");
	else
		printf("ERROR: root node spans more than one entry\n");

	printf("LOD table:\n");
	for (uint32 i=0; i<16; ++i)
		printf("LOD %d start: %d end: %d\n", i, spatialdatabase.m_LodStart[i], spatialdatabase.m_LodEnd[i]);

	// Scan the tree
	printf("Tree visit\n");
	uint32 currentnode = rdi;
	uint32 nodesvisited = 0;
	dumpTree(spatialdatabase, currentnode, nodesvisited);
	printf("\n");
	if (nodesvisited == spatialdatabase.m_dataLookup.size())
		printf("OK: can visit all nodes (visited %d out of %d)\n", nodesvisited, uint32(spatialdatabase.m_dataLookup.size()));
	else
		printf("ERROR: cannot visit all nodes, detached nodes found (visited %d out of %d)\n", nodesvisited, uint32(spatialdatabase.m_dataLookup.size()));

	return 0;
}
