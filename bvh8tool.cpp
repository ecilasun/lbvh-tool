#include <stdio.h>

#include "bvh8.h"

int main(int _argc, char** _argv)
{
	if (_argc<=1)
	{
		printf("BVH8 Tool v0.0\nUsage:\nbvh8tool inputfile outputfile\n");
		return 0;
	}

	SBVH8Database<BVH8LeafNode> spatialdatabase;

	SVec128 half = EVecConst(0.5f,0.5f,0.5f,0.0f);

	spatialdatabase.SetGridAABB(EVecConst(-16, -16, -16, 0), EVecConst(16, 16, 16, 0));
	spatialdatabase.Clear();

	// Test scene with 3 triangles
	float vertices[]={
		-3.f,0.f,0.f, 3.f,0.f,0.f, 1.5f,3.f,0.f,
		-3.f,0.f,1.f, 3.f,0.f,1.f, 1.5f,3.f,1.f,
		-3.f,0.f,2.f, 3.f,0.f,2.f, 1.5f,3.f,2.f };

	// Append triangles
	for (uint32 t=0; t<3; ++t)
	{
		// Start with empty AABB
		SBoundingBox primitiveAABB;
		EResetBounds(primitiveAABB);

		// Build vertices
		SVec128 v0 = EVecConst(vertices[t*9+0], vertices[t*9+1], vertices[t*9+2], 0.f);
		SVec128 v1 = EVecConst(vertices[t*9+3], vertices[t*9+4], vertices[t*9+5], 0.f);
		SVec128 v2 = EVecConst(vertices[t*9+6], vertices[t*9+7], vertices[t*9+8], 0.f);

		// Generate the AABB for this primitive
		EExpandBounds(primitiveAABB, v0);
		EExpandBounds(primitiveAABB, v1);
		EExpandBounds(primitiveAABB, v2);

		// Grab center of this AABB
		SVec128 origin = EVecMul(EVecAdd(primitiveAABB.m_Min, primitiveAABB.m_Max), half);

		// Quantize positon and encode as spatial key
		uint32 qXYZ[3];
		spatialdatabase.QuantizePosition(origin, qXYZ);
		uint32 mortonCode = EMortonEncode(qXYZ[0],qXYZ[1],qXYZ[2]);

		// Append this triangle as either a new node or part of existing node
		// NOTE: The next loop will iterate over contents to fit the AABB
		uint32 dataIndex, keyIndex;
		BVH8LeafNode* existingNode = spatialdatabase.Append(mortonCode, dataIndex, keyIndex);
		existingNode->m_NumTriangles++;
		existingNode->m_Vertices[0] = v0;
		existingNode->m_Vertices[1] = v1;
		existingNode->m_Vertices[2] = v2;
	}

	// Generate bounding box information for each node
	for( auto& scenedata : spatialdatabase.m_dataLookup)
	{
		// Get to the data at this node if there's one
		if (scenedata.m_DataIndex != BVH8NodeInfo::SSDIllegalIndex)
		{
			SBoundingBox nodeAABB;
			EResetBounds(nodeAABB);

			// Expand the AABB for this node to match its contents
			BVH8LeafNode& leafData = spatialdatabase.m_data[scenedata.m_DataIndex];
			for (int i=0;i<leafData.m_NumTriangles*3;++i)
			{
				EExpandBounds(nodeAABB, leafData.m_Vertices[i*3+0]);
				EExpandBounds(nodeAABB, leafData.m_Vertices[i*3+1]);
				EExpandBounds(nodeAABB, leafData.m_Vertices[i*3+2]);
			}

			// Update
			scenedata.m_BoundsMin = nodeAABB.m_Min;
			scenedata.m_BoundsMax = nodeAABB.m_Max;
		}
	}

	return 0;
}
