#include <stdio.h>
#include "math.h"
#include "lbvh.h"

// Utilizes code from https://github.com/jimbok8/lbvh-1 by Taylor Holberton and contributors

struct triangle final {
  float coords[9];
};

struct hitinfo final {
	float hitT;
};

triangle testtris[3] = {
	// Floor plane
	{3,-3,-3, 3,-3,3, -3,-3,-3},
	{-3,-3,-3, -3,-3,3, 3,-3,3},
	// Single tri facing forward
	{2,-3,0, -2,-3,0, 0,0,0},
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

	return 0;
}
