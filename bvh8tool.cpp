#include <stdio.h>
#include "math.h"
#include "lbvh.h"

#include "SDL2/SDL.h"

// Utilizes code from https://github.com/jimbok8/lbvh-1 by Taylor Holberton and contributors

struct triangle final {
  float coords[9];
};

struct hitinfo final {
	float hitT;
};

triangle testtris[3] = {
	// Floor plane
	{-3,-3,-3, 3,-3,13, 3,-3,-3},
	{-3,-3,-3, -3,-3,13, 3,-3,13},
	// Single tri facing forward
	{0,0,5, 0,-6,5, 2,-6,5},
};

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

	if (EVecGetFloatX(K) >= 0.f)
		return max_t; // Ignore backfacing (TODO: enable/disable this)

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


bool SlabTest(SVec128& p0, SVec128& p1, SVec128& rayOrigin, SVec128& rayDir, SVec128& invRayDir/*, SVec128& isect, SVec128& exitpos*/)
{
	SVec128 t0 = EVecMul(EVecSub(p0, rayOrigin), invRayDir);
	SVec128 t1 = EVecMul(EVecSub(p1, rayOrigin), invRayDir);
	SVec128 tmin = EVecMin(t0, t1);
	SVec128 tmax = EVecMax(t0, t1);
	float enter = EVecMaxComponent3(tmin);
	float exit = EVecMinComponent3(tmax);
	//SVec128 isect = EVecAdd(rayOrigin, EVecMul(EVecConst(enter,enter,enter,0.f), rayDir));
	//SVec128 exitpos = EVecAdd(rayOrigin, EVecMul(EVecConst(exit,exit,exit,0.f), rayDir));
	return enter <= exit;
}

void traceBVH(lbvh::bvh<float> &bvh, uint32 rootnode, uint32& marchCount, float& closestHit, SVec128& rayOrigin, SVec128& rayDir, SVec128& invRayDir)
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
		if (hitT < closestHit)
			closestHit = hitT;
		return;
	}

	SVec128 p0{bvh[nodeid].box.min.x, bvh[nodeid].box.min.y, bvh[nodeid].box.min.z, 0.f};
	SVec128 p1{bvh[nodeid].box.max.x, bvh[nodeid].box.max.y, bvh[nodeid].box.max.z, 0.f};

	/*SVec128 isect, exitpos;*/
	bool hit = SlabTest(p0, p1, rayOrigin, rayDir, invRayDir/*, isect, exitpos*/);
	if (hit) // We hit the bounds
	{
		// NOTE: This is easily solvable for BVH8 with a ray octant mask
		// In this case we'll need to check for the 'nearest' hit somehow
		// or know which node is the best hit candidate ahead of time

		traceBVH(bvh, bvh[nodeid].right, marchCount, closestHit, rayOrigin, rayDir, invRayDir);
		traceBVH(bvh, bvh[nodeid].left, marchCount, closestHit, rayOrigin, rayDir, invRayDir);
	}
}

SVec128 RayDirection(float fieldOfView, float fragCoordX, float fragCoordY)
{
	float width = 1920.f;
	float height = 1080.f;
    float x = fragCoordX - width/2.f;
    float y = fragCoordY - height/2.f;
    float z = height / tanf((fieldOfView*M_PI/180.f)/2.0);
    return EVecMul(EVecNorm3(EVecConst(x,y,-z,1.f)), EVecConst(64.f,64.f,64.f,1.f));
}

int main(int _argc, char** _argv)
{
	auto tri_to_box = [](const triangle& s) -> lbvh::aabb<float> {
		float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
		float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
		for (int i=0;i<3;++i)
		{
			minx = EMinimum(minx, s.coords[i*3+0]);
			miny = EMinimum(miny, s.coords[i*3+1]);
			minz = EMinimum(minz, s.coords[i*3+2]);
			maxx = EMaximum(maxx, s.coords[i*3+0]);
			maxy = EMaximum(maxy, s.coords[i*3+1]);
			maxz = EMaximum(maxz, s.coords[i*3+2]);
		}
		return lbvh::aabb<float> { {minx, miny, minz}, {maxx, maxy, maxz} };
	};

	lbvh::builder<float> builder;
	auto bvh = builder(testtris, 3, tri_to_box);

	dumpNodes(bvh, 0);
	printf("\n");

	// Trace the BVH
	if(SDL_Init(SDL_INIT_VIDEO) != 0)
	{
		fprintf(stderr, "Could not init SDL2: %s\n", SDL_GetError());
		return 1;
	}

    SDL_Window *screen = SDL_CreateWindow("BVH8Tool",
            SDL_WINDOWPOS_UNDEFINED,
            SDL_WINDOWPOS_UNDEFINED,
            1920, 1080, 0);

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

	float fov = 60.f;
	float imageAspectRatio = 1920.f / 1080.f; // width > height

	do{
		SDL_UpdateWindowSurface(screen);
		SDL_Event event;
		while(SDL_PollEvent(&event))
		{
			if(event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE))
				done = true;
		}

		if (SDL_MUSTLOCK(surface))
			SDL_LockSurface(surface);
		
		static float rotAng = 0.f;
		rotAng += 0.01f;

		uint8_t* pixels = (uint8_t*)surface->pixels;
		for (int y=0; y<1080; ++y)
		{
			for (int x=0; x<1920; ++x)
			{
				// Stationary camera
				SVec128 rayOrigin{0.f,0.f,30.f,0.f};

				// Rotate around world origin at 30 units distance
				//SVec128 rayOrigin{cosf(rotAng)*30.f,0.f,sinf(rotAng)*30.f,0.f};
				//SMatrix4x4 lookAt = EMatLookAtRightHanded(rayOrigin, EVecConst(0.f,0.f,0.f,0.f), EVecConst(0,1,0,0));

				// Look towards the scene
				SVec128 rayDir = RayDirection(60.f, float(x), float(1080-y));
				SVec128 invRayDir = EVecRcp(rayDir);

				// Rotate look direction to align with camera
				//SVec128 rotDir = EVecTransform3(lookAt, rayDir);
				//SVec128 invRayDir = EVecRcp(rotDir);

				// Trace and return hit depth
				float closestHit = 64.f;
				uint32 marchCount = 0;
				traceBVH(bvh, 0, marchCount, closestHit, rayOrigin, /*rotDir*/rayDir, invRayDir);

				if (closestHit<64.f)
				{
					pixels[(x+y*1920)*4+0] = int(closestHit*255.f);
					pixels[(x+y*1920)*4+1] = marchCount*8;
					pixels[(x+y*1920)*4+2] = 0;
				}
				else
				{
					pixels[(x+y*1920)*4+0] = 0;
					pixels[(x+y*1920)*4+1] = marchCount*8;
					pixels[(x+y*1920)*4+2] = 0;
				}
				pixels[(x+y*1920)*4+3] = 0xFF;
			}
		}

		if (SDL_MUSTLOCK(surface))
			SDL_UnlockSurface(surface);

	} while (!done);

	// Done
	SDL_FreeSurface(surface);
    SDL_DestroyWindow(screen);
	SDL_Quit();

	return 0;
}
