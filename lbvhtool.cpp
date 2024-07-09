#include "platform.h"
#include "locklesspipe.h"
#include "radixtree.h"
#include "objloader.h"

// Folder/file name of the scene to work with
//#define SCENENAME "sibenik"
#define SCENENAME "testscene"

// NOTE: Once there's a correct cell vs triangle test, this will be reduced
#define MAX_NODE_TRIS 64

// Depth of traversal stack
#define MAX_STACK_ENTRIES 16

// Define to see traversal count per tile
//#define SHOW_HEATMAP

// Define this to get worker tile debug view
//#define SHOW_WORKER_IDS

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
static const uint32_t tilewidth = 8;
static const uint32_t tileheight = 8;
static const uint32_t width = 512;
static const uint32_t height = 512;
#endif
static const uint32_t tilecountx = width/tilewidth;
static const uint32_t tilecounty = height/tileheight;
static const float cameradistance = 20.f;
static const float raylength = cameradistance + 1000.f;

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
	uint32_t heat{0};
	uint8_t rasterTile[4*tilewidth*tileheight]; // Internal rasterization tile (in hardware, to avoid arbitration need)
	SRenderContext *rc;
	CLocklessPipe<12> dispatchvector;
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

triangle *sceneGeometry;
SRadixTreeNode* testLBVH;
std::vector<SRadixTreeNode> testLeafnodes;
uint32_t lbvhLeafCount = 0;

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


void bvhBuilder(triangle* _triangles, uint32_t _numTriangles, std::vector<SRadixTreeNode> &_leafnodes, SRadixTreeNode **_lbvh)
{
	// Generate bounds AABB
	SVec128 sceneMin{FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
	SVec128 sceneMax{-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

	for (uint32_t i = 0; i < _numTriangles; ++i)
	{
		sceneMin = EVecMin(sceneMin, _triangles[i].coords[0]);
		sceneMin = EVecMin(sceneMin, _triangles[i].coords[1]);
		sceneMin = EVecMin(sceneMin, _triangles[i].coords[2]);

		sceneMax = EVecMax(sceneMax, _triangles[i].coords[0]);
		sceneMax = EVecMax(sceneMax, _triangles[i].coords[1]);
		sceneMax = EVecMax(sceneMax, _triangles[i].coords[2]);
	}

	SVec128 tentwentythree{1023.f, 1023.f, 1023.f, 1.f};
	SVec128 sceneBounds = EVecSub(sceneMax, sceneMin);
	SVec128 gridCellSize = EVecDiv(sceneBounds, tentwentythree);

	// Generate geometry pool

	_leafnodes.clear();
	SVec128 onethird = EVecConst(0.3333333f, 0.3333333f, 0.3333333f, 0.f);

	for (uint32_t tri = 0; tri < _numTriangles; ++tri)
	{
		SBoundingBox primitiveAABB;
		EResetBounds(primitiveAABB);

		SVec128 v0 = EVecSetW(_triangles[tri].coords[0], 1.f);
		SVec128 v1 = EVecSetW(_triangles[tri].coords[1], 1.f);
		SVec128 v2 = EVecSetW(_triangles[tri].coords[2], 1.f);

		EExpandBounds(primitiveAABB, v0);
		EExpandBounds(primitiveAABB, v1);
		EExpandBounds(primitiveAABB, v2);

		SVec128 origin = EVecMul(EVecAdd(v0, EVecAdd(v1, v2)), onethird);

		uint32_t qXYZ[3];
		EQuantizePosition(origin, qXYZ, sceneMin, gridCellSize);
		uint32_t mortonCode = EMortonEncode(qXYZ[0], qXYZ[1], qXYZ[2]);

		SRadixTreeNode node;
		node.m_bounds.m_Min = primitiveAABB.m_Min;
		node.m_bounds.m_Max = primitiveAABB.m_Max;
		node.m_primitiveIndex = tri;
		node.m_spatialKey = mortonCode;
		_leafnodes.emplace_back(node);
	}

	lbvhLeafCount = (uint32_t)_leafnodes.size();
	*_lbvh = new SRadixTreeNode[2*lbvhLeafCount-1];

	GenerateLBVH(*_lbvh, _leafnodes, lbvhLeafCount);

	printf("LBVH data generated. Leaf node count:%d\n", lbvhLeafCount);
}

bool ClosestHitLBVH(const SRadixTreeNode &_self, const SVec128 &_rayStart, const SVec128 &_rayDir, float &_t, const float _tmax, uint32_t &_heat)
{
	uint32_t tri = _self.m_primitiveIndex;
	if (tri == 0xFFFFFFFF)
		return false;
	
	bool isHit = HitTriangle (
		sceneGeometry[tri].coords[0],
		sceneGeometry[tri].coords[1],
		sceneGeometry[tri].coords[2],
		_rayStart, _rayDir,
		_t, _tmax );

	_heat += isHit ? 1:0;

	return isHit;
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
					SVec128 rayDir = EVecMul(EVecAdd(EVecAdd(L, U), vec->rc->F), raylenvec);
					SVec128 rayEnd = EVecAdd(vec->rc->rayOrigin, rayDir);

					uint32_t hitNode = 0xFFFFFFFF;
					SVec128 hitpos;

					vec->heat = 0;
					FindClosestHitLBVH(testLBVH, lbvhLeafCount, vec->rc->rayOrigin, rayEnd, t, hitpos, hitNode, vec->heat, ClosestHitLBVH);

					float final = 0.f;

					if (hitNode != 0xFFFFFFFF)
					{
						auto &self = testLBVH[hitNode];
						uint32_t tri = self.m_primitiveIndex;

						SVec128 sunPos{20.f,35.f,20.f,1.f};
						SVec128 sunRay = EVecSub(sunPos, hitpos);
						SVec128 invSunRay = EVecRcp(sunRay);
						SVec128 nrm;

						SVec128 viewRay = EVecNorm3(rayDir);

						// Global + NdotL
						{
							SVec128 uvw;
							float fuvw[3];
							CalculateBarycentrics(hitpos,
								sceneGeometry[tri].coords[0],
								sceneGeometry[tri].coords[1],
								sceneGeometry[tri].coords[2], fuvw);
							SVec128 uvwx = EVecSplatX(EVecConst(fuvw[0]));
							SVec128 uvwy = EVecSplatX(EVecConst(fuvw[1]));
							SVec128 uvwz = EVecSplatX(EVecConst(fuvw[2]));
							SVec128 uvwzA = EVecMul(uvwz, sceneGeometry[tri].normals[0]); // A*uvw.zzz
							SVec128 uvwyB = EVecMul(uvwy, sceneGeometry[tri].normals[1]); // B*uvw.yyy
							SVec128 uvwxC = EVecMul(uvwx, sceneGeometry[tri].normals[2]); // C*uvw.xxx
							nrm = EVecNorm3(EVecAdd(uvwzA, EVecAdd(uvwyB, uvwxC))); // A*uvw.zzz + B*uvw.yyy + C*uvw.xxx
							float L = fabs(EVecGetFloatX(EVecDot3(nrm, EVecNorm3(sunRay))));
							final += L;
						}
					}

					uint8_t C = uint8_t(final*255.f);
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

#if defined(SHOW_HEATMAP)
					p[(ix+iy*tilewidth)*4+2] = (vec->heat*64)%255; // R
#endif
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

	sceneGeometry = new triangle[totaltriangles];

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

			sceneGeometry[t].coords[0] = EVecConst( mesh.Vertices[i0].Position.X, mesh.Vertices[i0].Position.Y, mesh.Vertices[i0].Position.Z, 0.f);
			sceneGeometry[t].coords[1] = EVecConst( mesh.Vertices[i1].Position.X, mesh.Vertices[i1].Position.Y, mesh.Vertices[i1].Position.Z, 0.f);
			sceneGeometry[t].coords[2] = EVecConst( mesh.Vertices[i2].Position.X, mesh.Vertices[i2].Position.Y, mesh.Vertices[i2].Position.Z, 0.f);

			sceneGeometry[t].normals[0] = EVecConst( mesh.Vertices[i0].Normal.X, mesh.Vertices[i0].Normal.Y, mesh.Vertices[i0].Normal.Z, 0.f);
			sceneGeometry[t].normals[1] = EVecConst( mesh.Vertices[i1].Normal.X, mesh.Vertices[i1].Normal.Y, mesh.Vertices[i1].Normal.Z, 0.f);
			sceneGeometry[t].normals[2] = EVecConst( mesh.Vertices[i2].Normal.X, mesh.Vertices[i2].Normal.Y, mesh.Vertices[i2].Normal.Z, 0.f);

			++t;
		}
	}

	// Build LBVH
	bvhBuilder(sceneGeometry, totaltriangles, testLeafnodes, &testLBVH);

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

	rc.rotAng = 3.141592f;
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
		rc.rayOrigin = SVec128{sinf(rc.rotAng)*cameradistance, sinf(rc.rotAng*0.1f)*12.f, cosf(rc.rotAng)*cameradistance, 1.f};
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
				if (wc[i].dispatchvector.FreeSpace() != 0) // We have space in this worker's queue
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
		int tdone;
		do
		{
			tdone = 0;
			for (uint32_t i=0; i<MAX_WORKERS; ++i)
				tdone += wc[i].dispatchvector.BytesAvailable() ? 0 : 1;
		} while(tdone != MAX_WORKERS);

		// TODO: Copy tiles to their respective positions
		//(uint8_t*)surface->pixels;

		if (SDL_MUSTLOCK(surface))
			SDL_UnlockSurface(surface);
		SDL_UpdateWindowSurface(screen);

	} while (!done);

	delete [] sceneGeometry;
	delete testLBVH;

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
