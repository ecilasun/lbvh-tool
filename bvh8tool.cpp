#include <stdio.h>
#include "math.h"
#include "lbvh.h"
#include "objloader.h"

#include "SDL2/SDL.h"

// Utilizes code from https://github.com/jimbok8/lbvh-1 by Taylor Holberton and contributors

struct triangle final {
  float coords[9];
};

struct hitinfo final {
	float hitT;
};

uint32 width = 1280;
uint32 height = 960;

uint8_t* pixels;
triangle *testtris;

/*triangle testtris[3] = {
	// Floor plane
	{-3,-3,-3, 3,-3,5, 3,-3,-3},
	{-3,-3,-3, -3,-3,5, 3,-3,5},
	// Single tri facing forward
	{0,4,2, 0,-9,2, 2,-9,2},
};*/

void dumpNodes(lbvh::bvh<float> &bvh, uint32 node)
{
	auto nb = bvh[node];

	printf("(%d)->", node);

	if (nb.left&0x80000000)
		printf("(L v0:%f,%f,%f v1:%f,%f,%f v2:%f,%f,%f)\n",
			testtris[node].coords[0],testtris[node].coords[1],testtris[node].coords[2],
			testtris[node].coords[3],testtris[node].coords[4],testtris[node].coords[5],
			testtris[node].coords[6],testtris[node].coords[7],testtris[node].coords[8]); // left leads to leaf node
	else
		dumpNodes(bvh, nb.left); // left leads to an internal node

	if (nb.right&0x80000000)
		printf("(R v0:%f,%f,%f v1:%f,%f,%f v2:%f,%f,%f)\n",
			testtris[node].coords[0],testtris[node].coords[1],testtris[node].coords[2],
			testtris[node].coords[3],testtris[node].coords[4],testtris[node].coords[5],
			testtris[node].coords[6],testtris[node].coords[7],testtris[node].coords[8]); // right leads to leaf node
	else
		dumpNodes(bvh, nb.right); // right leads to an internal node
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

	/*if (EVecGetFloatX(K) >= 0.f)
		return max_t; // Ignore backfacing (TODO: enable/disable this)*/

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

void traceBVH(lbvh::bvh<float> &bvh, uint32 rootnode, uint32& marchCount, float& t, SVec128& rayOrigin, uint32& hitID, SVec128& rayDir, SVec128& invRayDir)
{
	uint nodeid = rootnode&0x7FFFFFFF;
	marchCount++;

	if (rootnode&0x80000000) // This is a leaf node, do a triangle test
	{
		SVec128 v0{testtris[nodeid].coords[0], testtris[nodeid].coords[1], testtris[nodeid].coords[2],0.f};
		SVec128 v1{testtris[nodeid].coords[3], testtris[nodeid].coords[4], testtris[nodeid].coords[5],0.f};
		SVec128 v2{testtris[nodeid].coords[6], testtris[nodeid].coords[7], testtris[nodeid].coords[8],0.f};
		float hitT = TriHit(rayOrigin, rayDir, v0, v1, v2, 512.f);

		// Closer than closest hit?
		if (hitT < t && hitT > 0.f)
		{
			t = hitT;
			hitID = nodeid;
		}
	}
	else
	{
		SVec128 p0{bvh[nodeid].box.min.x, bvh[nodeid].box.min.y, bvh[nodeid].box.min.z, 0.f};
		SVec128 p1{bvh[nodeid].box.max.x, bvh[nodeid].box.max.y, bvh[nodeid].box.max.z, 0.f};

		SVec128 isect, exitpos;
		bool hit = SlabTest(p0, p1, rayOrigin, rayDir, invRayDir, isect, exitpos);
		if (hit) // We hit the bounds
		{
			// NOTE: This is easily solvable for BVH8 with a ray octant mask
			// In this case we'll need to check for the 'nearest' hit somehow
			// or know which node is the best hit candidate ahead of time

			traceBVH(bvh, bvh[nodeid].right, marchCount, t, rayOrigin, hitID, rayDir, invRayDir);
			traceBVH(bvh, bvh[nodeid].left, marchCount, t, rayOrigin, hitID, rayDir, invRayDir);

			// If either are missed, bring 't' closer
			//t = EMinimum(t, EVecGetFloatX(EVecLen3(EVecSub(rayOrigin, exitpos))));
		}
	}
}

void block(int x, int y, uint8_t B, uint8_t G, uint8_t R)
{
	for (int oy=y;oy<y+4;++oy)
	for (int ox=x;ox<x+4;++ox)
	{
		pixels[(ox+oy*width)*4+0] = B;
		pixels[(ox+oy*width)*4+1] = G;
		pixels[(ox+oy*width)*4+2] = R;
		pixels[(ox+oy*width)*4+3] = 0xFF;
	}
}

int main(int _argc, char** _argv)
{
	objl::Loader objloader;
	if (!objloader.LoadFile("test.obj"))
	{
		printf("Failed to load test.obj\n");
		return 1;
	}

	// Set up triangle data
	int t=0;
	int totaltriangles = 0; // 3

	for (auto &mesh : objloader.LoadedMeshes)
		totaltriangles += int(mesh.Indices.size()/3);

	if (totaltriangles==0)
	{
		printf("No triangles in model\n");
		return 1;
	}

	testtris = new triangle[totaltriangles];

	for (auto &mesh : objloader.LoadedMeshes)
	{
		for (int i=0;i<mesh.Indices.size()/3;++i)
		{
			unsigned int i0 = mesh.Indices[i*3+0];
			unsigned int i1 = mesh.Indices[i*3+1];
			unsigned int i2 = mesh.Indices[i*3+2];

			testtris[t].coords[0] = mesh.Vertices[i0].Position.X;
			testtris[t].coords[1] = mesh.Vertices[i0].Position.Y;
			testtris[t].coords[2] = mesh.Vertices[i0].Position.Z;

			testtris[t].coords[3] = mesh.Vertices[i1].Position.X;
			testtris[t].coords[4] = mesh.Vertices[i1].Position.Y;
			testtris[t].coords[5] = mesh.Vertices[i1].Position.Z;

			testtris[t].coords[6] = mesh.Vertices[i2].Position.X;
			testtris[t].coords[7] = mesh.Vertices[i2].Position.Y;
			testtris[t].coords[8] = mesh.Vertices[i2].Position.Z;

			++t;
		}
	}

	auto tri_to_box = [](const triangle& s) -> lbvh::aabb<float> {
		float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
		float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
		for (int i=0;i<3;++i)
		{
			float X = s.coords[i*3+0];
			float Y = s.coords[i*3+1];
			float Z = s.coords[i*3+2];
			minx = EMinimum(minx, X);
			miny = EMinimum(miny, Y);
			minz = EMinimum(minz, Z);
			maxx = EMaximum(maxx, X);
			maxy = EMaximum(maxy, Y);
			maxz = EMaximum(maxz, Z);
		}
		return lbvh::aabb<float> { {minx, miny, minz}, {maxx, maxy, maxz} };
	};

	lbvh::builder<float> builder;
	auto bvh = builder(testtris, totaltriangles, tri_to_box);

	/*dumpNodes(bvh, 0);
	printf("\n");*/

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
	const float cameradistance = 20.f;

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

		SVec128 rayOrigin{sinf(rotAng)*cameradistance, cameradistance*0.5f, cosf(rotAng)*cameradistance, 1.f};
		SVec128 lookAt{0.f,0.f,0.f,1.f};
		SVec128 upVec{0.f,1.f,0.f,0.f};
		SMatrix4x4 lookMat = EMatLookAtRightHanded(rayOrigin, lookAt, upVec);
		SVec128 pzVec{-1.f,-1.f,-1.f,0.f};
		SVec128 F = EVecMul(lookMat.r[2], pzVec);

		pixels = (uint8_t*)surface->pixels;
		for (int y=0; y<height; y+=4)
		{
			float py = aspect * (float(height)/2.f-float(y))/float(height);
			SVec128 pyVec{py,py,py,0.f};
			SVec128 U = EVecMul(lookMat.r[1], pyVec);
			for (int x=0; x<width; x+=4)
			{
				// Rotating camera
				float px = (float(x) - float(width)/2.f)/float(width);
				SVec128 pxVec{px,px,px,0.f};
				SVec128 L = EVecMul(lookMat.r[0], pxVec);
				SVec128 traceRay = EVecAdd(EVecAdd(L, U), F);
				SVec128 invRay = EVecRcp(traceRay);

				// Trace and return hit depth
				float t = 64.f;
				uint32 marchCount = 0;
				uint32 hitID = 0;
				traceBVH(bvh, 0, marchCount, t, rayOrigin, hitID, traceRay, invRay);

				SVec128 hitpos = EVecAdd(rayOrigin, EVecMul(traceRay,  EVecConst(t,t,t,1.f)));

				float D = t<64.f ? EVecGetFloatX(EVecLen3(EVecSub(hitpos, rayOrigin)))/16.f : 0.f;//(t<64.f?t:0.f)/16.f;
				int C = int(D*255.f);

				block(x,y, C, ((hitID>>1)%4)*32, ((hitID>>2)%8)*32);
				//block(x,y, hitID, 0, marchCount);
			}
		}

		rotAng += 0.02f;

		if (SDL_MUSTLOCK(surface))
			SDL_UnlockSurface(surface);
		SDL_UpdateWindowSurface(screen);

	} while (!done);

	//delete [] testtris;

	// Done
	SDL_FreeSurface(surface);
    SDL_DestroyWindow(screen);
	SDL_Quit();

	return 0;
}
