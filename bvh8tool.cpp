#include <stdio.h>
#include "math.h"
#include "lbvh.h"
#include "objloader.h"

#include "SDL2/SDL.h"

//#define USE_SMALL_DATA_SET

// Utilizes code from https://github.com/jimbok8/lbvh-1 by Taylor Holberton and contributors

struct triangle final {
  SVec128 coords[3];
  SVec128 normals[3];
};

struct hitinfo final {
	float hitT;
};

uint32 width = 1280;
uint32 height = 960;

uint8_t* pixels;
#ifdef USE_SMALL_DATA_SET
triangle testtris[3] = {
	// Floor plane
	SVec128{-3,-3,-3,0.f}, SVec128{3,-3,5,0.f}, SVec128{3,-3,-3,0.f},
	SVec128{-3,-3,-3,0.f}, SVec128{-3,-3,5,0.f}, SVec128{3,-3,5,0.f},
	// Single tri facing forward
	SVec128{0,4,2,0.f}, SVec128{0,-9,2,0.f}, SVec128{2,-9,2,0.f} };
#else
triangle *testtris;
#endif

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
	marchCount++;
	uint nodeid = rootnode&0x7FFFFFFF;

	if (rootnode&0x80000000) // This is a leaf node, do a triangle test
	{
		float hitT = TriHit(rayOrigin, rayDir, testtris[nodeid].coords[0], testtris[nodeid].coords[1], testtris[nodeid].coords[2], t+1.f);

		// Closer than closest hit and not missed
		if (hitT < t)
		{
			t = hitT;
			hitID = nodeid;
		}
	}
	else
	{
		/*if (hitID!=0xFFFFFFFF) // ANY_HIT: return as soon as we have a triangle
			return;*/

		SVec128 p0{bvh[nodeid].box.min.x, bvh[nodeid].box.min.y, bvh[nodeid].box.min.z, 0.f};
		SVec128 p1{bvh[nodeid].box.max.x, bvh[nodeid].box.max.y, bvh[nodeid].box.max.z, 0.f};

		SVec128 isect, exitpos;
		bool hit = SlabTest(p0, p1, rayOrigin, rayDir, invRayDir, isect, exitpos);
		if (hit) // We hit the bounds
		{
			// NOTE: This is easily solvable for BVH8 with a ray octant mask
			// In this case we'll need to check for the 'nearest' hit somehow
			// or know which node is the best hit candidate ahead of time
			SVec128 delta = EVecSub(isect, rayOrigin);
			float vlen = EVecGetFloatX(EVecLen3(delta));

			// Current ray is further than this cell, take the branch
			if(t > vlen) // CLOSEST_HIT: find the hit position closest to the ray origin (i.e. minimum t)
			{
				traceBVH(bvh, bvh[nodeid].left, marchCount, t, rayOrigin, hitID, rayDir, invRayDir);
				traceBVH(bvh, bvh[nodeid].right, marchCount, t, rayOrigin, hitID, rayDir, invRayDir);
			}
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
#ifndef USE_SMALL_DATA_SET
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

	for (auto &mesh : objloader.LoadedMeshes)
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
#endif

	auto tri_to_box = [](const triangle& s) -> lbvh::aabb<float> {
		SVec128 fmin = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
		SVec128 fmax = {-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

		fmin = EVecMin(fmin, s.coords[0]);
		fmin = EVecMin(fmin, s.coords[1]);
		fmin = EVecMin(fmin, s.coords[2]);
		fmax = EVecMax(fmax, s.coords[0]);
		fmax = EVecMax(fmax, s.coords[1]);
		fmax = EVecMax(fmax, s.coords[2]);

		float minx = EVecGetFloatX(fmin);
		float miny = EVecGetFloatY(fmin);
		float minz = EVecGetFloatZ(fmin);
		float maxx = EVecGetFloatX(fmax);
		float maxy = EVecGetFloatY(fmax);
		float maxz = EVecGetFloatZ(fmax);

		return lbvh::aabb<float> { {minx, miny, minz}, {maxx, maxy, maxz} };
	};

#ifdef USE_SMALL_DATA_SET
	int totaltriangles = 3;
#endif

	lbvh::builder<float> builder;
	auto bvh = builder(testtris, totaltriangles, tri_to_box);

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
	uint32 lowTraces = 0xFFFFFFFF;
	uint32 highTraces = 0x00000000;

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
		uint32 maxTraces = 0;
		SVec128 nil{0.f,0.f,0.f,0.f};
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

				// Scene
				float t = cameradistance*2.f;
				uint32 marchCount = 0;
				uint32 hitID = 0xFFFFFFFF;
				traceBVH(bvh, 0, marchCount, t, rayOrigin, hitID, traceRay, invRay);
				maxTraces = EMaximum(maxTraces, marchCount);

				if (hitID!=0xFFFFFFFF) // HAVE_HIT: Ray hit a primitive in a leaf nodes
				{
					// Hit position
					SVec128 hitpos = EVecAdd(rayOrigin, EVecMul(traceRay,  EVecConst(t-0.01f,t-0.01f,t-0.01f,1.f)));
					// Barycentric coortinates for attribute interpolation
					SVec128 uvw;
					Barycentrics(hitpos, testtris[hitID].coords[0], testtris[hitID].coords[1], testtris[hitID].coords[2], uvw);
					block(x,y, ((hitID>>1)%2)*128, ((hitID>>1)%4)*64, ((hitID>>2)%8)*32);
				}
				else
				{
					block(x,y, 0,0,0);
				}
				// X-Ray view
				//block(x,y, marchCount, marchCount, marchCount);

				// Shadow
				/*if (hitID!=0xFFFFFFFF)
				{
					float D = (t<128.f?t:0.f)/16.f;
					int C = int(D*255.f);*/
					//int hX = int(EVecGetFloatX(hitpos)*100.f);
					//int hY = int(EVecGetFloatY(hitpos)*100.f);
					//int hZ = int(EVecGetFloatZ(hitpos)*100.f);
					//block(x,y, C, C, C);

					/*SVec128 sunPos{20.f,35.f,sinf(rotAng*4.f)*20.f,1.f};
					SVec128 sunRay = EVecSub(sunPos, hitpos);
					SVec128 invSunRay = EVecRcp(sunRay);
					float t2 = 64.f;
					hitID = 0xFFFFFFFF;
					traceBVH(bvh, 0, marchCount, t2, hitpos, hitID, sunRay, invSunRay);
					float sunlen = EVecGetFloatX(EVecLen3(sunRay));
					if (t2<sunlen)
					{
						float D = fabs(t2-sunlen);
						int C = int(D);
						block(x,y, C, C, C);
					}
					else
					{
						float D = EVecGetFloatX(EVecLen3(hitpos));
						int C = int(D*16.f);
						block(x,y, C, C, C);
					}
				}*/
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

	//delete [] testtris;

	// Done
	SDL_FreeSurface(surface);
    SDL_DestroyWindow(screen);
	SDL_Quit();

	return 0;
}
