#include <math.h>
#include <float.h>

#include <array>
#include <vector>
#include <algorithm>

#include "math.h"

#define MAX_TRIS 8

#define BVH8IllegalIndex 0xFFFFFFFF

// NOTE: This structure is very expensive for E32E, a data reduction method has to be applied here.
// One approach could be to only stash connected triangle fans here.
struct BVH8LeafNode
{
	uint32 m_NumTriangles{0};
	SVec128 m_Vertices[3*MAX_TRIS]; // X-Y-Z-A
};

// NOTE: It is best if this structure fits into E32E's 64 byte cache line
// Currently it's 48 bytes which means 
struct BVH8NodeInfo
{
	BVH8NodeInfo(){ }
	explicit BVH8NodeInfo(const uint32 _spatialKey, const uint32 _dataindex) : m_SpatialKey(_spatialKey), m_DataIndex(_dataindex) { }

	// LCM of 48 (current byte size) and 64 (E32E cache size) is 192
	// This means 4 of this structure will fit into 3 cache lines
	uint32 m_SpatialKey{0};
	uint32 m_DataIndex{BVH8IllegalIndex};
	uint32 m_ChildCount{0};
	uint32 m_unused_padding_0;
	SVec128 m_BoundsMin{FLT_MAX,FLT_MAX,FLT_MAX,0.f};
	SVec128 m_BoundsMax{-FLT_MAX,-FLT_MAX,-FLT_MAX,0.f}; // There are 3 unused DWORDS in this structure
};

// --------------------------------------------------------------------------------
// Radix sort for 4 byte 8bit unsigned integers
// Avoids compare, only reads value and directly clusters by radix
// --------------------------------------------------------------------------------

inline void ERadixSortClear(unsigned int count[256])
{
	// Clear the count table
	for (int t = 0; t < 256; ++t)
		count[t] = 0;
}

inline void ERadixSortScanAscending(const unsigned int _phase, BVH8NodeInfo *_sortArray, const size_t _numElements, unsigned int count[256])
{
	// Scan count each element's existance in the input array
	unsigned int shiftCount = _phase * 8;
	for (size_t elem = 0; elem < _numElements; ++elem)
	{
		unsigned int rdx = (_sortArray[elem].m_SpatialKey >> shiftCount) & 0xFF;
		count[rdx]++;
	}
}

inline void ERadixSortScanDescending(const unsigned int _phase, BVH8NodeInfo *_sortArray, const size_t _numElements, unsigned int count[256])
{
	// Scan count each element's existance in the input array
	unsigned int shiftCount = _phase * 8;
	for (size_t elem = 0; elem < _numElements; ++elem)
	{
		unsigned int rdx = (_sortArray[elem].m_SpatialKey >> shiftCount) & 0xFF;
		rdx = rdx ^ 0xFF;
		count[rdx]++;
	}
}

inline void ERadixSortGenerateOffsets(unsigned int count[256], unsigned int offsets[256])
{
	// Generate prefix sum
	offsets[0] = 0;
	for (int t = 1; t < 256; ++t)
		offsets[t] = offsets[t-1] + count[t-1];
}

inline void ERadixSortRemapAscending(const unsigned int _phase, BVH8NodeInfo *_sortArray, BVH8NodeInfo *_tmpArray, const size_t _numElements, unsigned int offsets[256])
{
	// Remap elements to new indices
	unsigned int shiftCount = _phase * 8;
	for (int elem = 0; elem <_numElements; ++elem)
	{
		unsigned int rdx = (_sortArray[elem].m_SpatialKey >> shiftCount) & 0xFF;
		int newEntry = offsets[rdx];
		_tmpArray[newEntry] = _sortArray[elem];
		offsets[rdx]++;
	}
}

inline void ERadixSortRemapDescending(const unsigned int _phase, BVH8NodeInfo *_sortArray, BVH8NodeInfo *_tmpArray, const size_t _numElements, unsigned int offsets[256])
{
	// Remap elements to new indices
	unsigned int shiftCount = _phase * 8;
	for (int elem = 0; elem <_numElements; ++elem)
	{
		unsigned int rdx = (_sortArray[elem].m_SpatialKey >> shiftCount) & 0xFF;
		rdx = rdx ^ 0xFF;
		int newEntry = offsets[rdx];
		_tmpArray[newEntry] = _sortArray[elem];
		offsets[rdx]++;
	}
}

inline void ERadixSortSpatialDatabaseAscending(BVH8NodeInfo *_sortArray, BVH8NodeInfo *_tmpArray, size_t _numElements)
{
	// Count table (per-byte sort, therefore only 2^8 elements long)
	unsigned int offsets[256];
	unsigned int count[256];

	// Start with src=source and tgt=target, swapping at each phase.
	BVH8NodeInfo *srce = _sortArray;
	BVH8NodeInfo *trgt = _tmpArray;

	// There are 4 phases in our sort (using 4 bytes of the unsiged integer being sorted)
	for (int phase = 0; phase < 4; ++phase)
	{
		// 0 - Reset count table (dispatch("ERadixSortClear.cp",1,1,1), only 256 entries)
		ERadixSortClear(count);

		// 1 - Scan left to right, grab current byte, increment its counter (dispatch("radixSortScan.cp",roundup(_numElements/64),1,1))
		ERadixSortScanAscending(phase, srce, _numElements, count);

		// 2 - Generate placement offset table (dispatch("radixSortGenerateOffsets.cp",1,1,1), only 256 entries)
		ERadixSortGenerateOffsets(count, offsets);

		// 3 - Reorder using count table (dispatch("radixSortRemap.cp",roundup(_numElements/64),1,1))
		ERadixSortRemapAscending(phase, srce, trgt, _numElements, offsets);

		// Swap structure pointers (on 4th byte, this will be back to original 'src')
		BVH8NodeInfo *temp = srce;
		srce = trgt;
		trgt = temp;
	}
}

inline void ERadixSortSpatialDatabaseDescending(BVH8NodeInfo *_sortArray, BVH8NodeInfo *_tmpArray, size_t _numElements)
{
	// Count table (per-byte sort, therefore only 2^8 elements long)
	unsigned int offsets[256];
	unsigned int count[256];

	// Start with src=source and tgt=target, swapping at each phase.
	BVH8NodeInfo *srce = _sortArray;
	BVH8NodeInfo *trgt = _tmpArray;

	// There are 4 phases in our sort (using 4 bytes of the unsiged integer being sorted)
	for (int phase = 0; phase < 4; ++phase)
	{
		// 0 - Reset count table (dispatch("ERadixSortClear.cp",1,1,1), only 256 entries)
		ERadixSortClear(count);

		// 1 - Scan left to right, grab current byte, increment its counter (dispatch("radixSortScan.cp",roundup(_numElements/64),1,1))
		ERadixSortScanDescending(phase, srce, _numElements, count);

		// 2 - Generate placement offset table (dispatch("radixSortGenerateOffsets.cp",1,1,1), only 256 entries)
		ERadixSortGenerateOffsets(count, offsets);

		// 3 - Reorder using count table (dispatch("radixSortRemap.cp",roundup(_numElements/64),1,1))
		ERadixSortRemapDescending(phase, srce, trgt, _numElements, offsets);

		// Swap structure pointers (on 4th byte, this will be back to original 'src')
		BVH8NodeInfo *temp = srce;
		srce = trgt;
		trgt = temp;
	}
}

// --------------------------------------------------------------------------------------------------------------------------
// BVH related structures
// --------------------------------------------------------------------------------------------------------------------------

struct SBVHLeafNode
{
	SVec128 m_AABBMin{0,0,0,0};
	SVec128 m_AABBMax{0,0,0,0};
	uint32 m_ParentIndex{0};
	uint32 m_LeftChildIndex{0};
	uint32 m_RightChildIndex{0};
	uint32 m_Level{0};
};

// --------------------------------------------------------------------------------------------------------------------------
// Spatial encoding helpers
// --------------------------------------------------------------------------------------------------------------------------

EInline uint32 EMortonEncode(uint32 _x, uint32 _y, uint32 _z)
{
	// Pack 3 10-bit indices into a 30-bit Morton code
	// Logic below is HLSL compatible
	_x &= 0x000003ff;	_y &= 0x000003ff;	_z &= 0x000003ff;
	_x |= (_x << 16);	_y |= (_y << 16);	_z |= (_z << 16);
	_x &= 0xff0000ff;	_y &= 0xff0000ff;	_z &= 0xff0000ff;
	_x |= (_x << 8);	_y |= (_y << 8);	_z |= (_z << 8);
	_x &= 0x0300f00f;	_y &= 0x0300f00f;	_z &= 0x0300f00f;
	_x |= (_x << 4);	_y |= (_y << 4);	_z |= (_z << 4);
	_x &= 0x030c30c3;	_y &= 0x030c30c3;	_z &= 0x030c30c3;
	_x |= (_x << 2);	_y |= (_y << 2);	_z |= (_z << 2);
	_x &= 0x09249249;	_y &= 0x09249249;	_z &= 0x09249249;
	return (_x | (_y << 1) | (_z << 2));
}

EInline void EMortonDecode(const uint32 _morton, uint32& _index1, uint32& _index2, uint32& _index3)
{
	// Unpack 3 10-bit indices from a 30-bit Morton code
	// Logic below is HLSL compatible
	uint32 value1 = _morton;
	uint32 value2 = (value1 >> 1);
	uint32 value3 = (value1 >> 2);
	value1 &= 0x09249249;		value2 &= 0x09249249;		value3 &= 0x09249249;
	value1 |= (value1 >> 2);	value2 |= (value2 >> 2);	value3 |= (value3 >> 2);
	value1 &= 0x030c30c3;		value2 &= 0x030c30c3;		value3 &= 0x030c30c3;
	value1 |= (value1 >> 4);	value2 |= (value2 >> 4);	value3 |= (value3 >> 4);
	value1 &= 0x0300f00f;		value2 &= 0x0300f00f;		value3 &= 0x0300f00f;
	value1 |= (value1 >> 8);	value2 |= (value2 >> 8);	value3 |= (value3 >> 8);
	value1 &= 0xff0000ff;		value2 &= 0xff0000ff;		value3 &= 0xff0000ff;
	value1 |= (value1 >> 16);	value2 |= (value2 >> 16);	value3 |= (value3 >> 16);
	value1 &= 0x000003ff;		value2 &= 0x000003ff;		value3 &= 0x000003ff;
	_index1 = value1;
	_index2 = value2;
	_index3 = value3;
}

EInline void EQuantizePosition(const SVec128 &_worldpos, uint32 _qXYZ[3], SVec128 _gridAABBMin, SVec128 _gridCellSize)
{
	SVec128 gridLocalPosition = EVecDiv(EVecSub(_worldpos, _gridAABBMin), _gridCellSize);

	// Clamp within 0-SSDGridEntryCountPerAxis-1 range (inclusive)
	static const SVec128 cellClampMax{ 1023.f, 1023.f, 1023.f, 0.f};
	static const SVec128 cellClampMin{ 0.f, 0.f, 0.f, 0.f };
	gridLocalPosition = EVecMin(EVecMax(gridLocalPosition, cellClampMin), cellClampMax);

	_qXYZ[0] = uint32(EFloor(EVecGetFloatX(gridLocalPosition)));
	_qXYZ[1] = uint32(EFloor(EVecGetFloatY(gridLocalPosition)));
	_qXYZ[2] = uint32(EFloor(EVecGetFloatZ(gridLocalPosition)));
}

// --------------------------------------------------------------------------------------------------------------------------
// Raytracing structures
// --------------------------------------------------------------------------------------------------------------------------
struct SBVH8Ray
{
	SVec128 m_RayStart;
	SVec128 m_RayEnd;
};

struct SBVH8Hit
{
	uint32 m_IsHit{0};
	uint32 m_HitNodeIndex{0xFFFFFFFF};			// Index of the current node hit, usually caller just needs to pass this out as it was passed it without modification
	uint32 m_HitDataIndex{0xFFFFFFFF};			// Index of the current data cell hit, usually caller just needs to pass this out as it was passed it without modification
	uint32 m_HitLeafItemIndex{ 0xFFFFFFFF };	// If the hit test implements precise data hits with cell contents, this is the index of the item in the cell that was hit
	SVec128 m_HitPosition{0,0,0,0};				// Hit position either on the slab (as returned by IntersectSlab()) or the point on planar/procedural geometry contained within
	SVec128 m_ExitPosition{0,0,0,0};			// Exit point of ray can be used to 'resume' tracing if this hit was not a viable candidate. Set it equal to m_HitPosition for infinitely thin geometry, or exit point returned by IntersectSlab()
	SVec128 m_HitNormal{ 0,0,0,0 };				// Normal at hit position, if applicable
	float m_UVWT[4];							// Barycentric coordinates at intersection. Contents can be interpreted as a triangle, quad, line etc based on data cell contents
};

// --------------------------------------------------------------------------------------------------------------------------
// Sparse spatial database
// --------------------------------------------------------------------------------------------------------------------------

template <typename T>
struct SBVH8Database
{
	static const uint32 MaxOctreeLevels = 16;
	// List of occupied tiles, encoded using the Morton_3D_Encode functions
	// Make sure the input is quantized properly to ensure 0..1023 range is not exceeded per axis in your local grid.
	std::vector<BVH8NodeInfo> m_dataLookup;
	// Data element per key. Optional, place your data items here if data access is needed.
	std::vector<T> m_data;
	float m_cellScale{ 1.f };
	SVec128 m_GridAABBMin, m_GridCellSize;
	uint32 m_LodStart[MaxOctreeLevels]{};
	uint32 m_LodEnd[MaxOctreeLevels]{};
	uint32 m_RootBVH8Node{};
	bool m_isSorted;

	typedef bool(*TLeafNodeHitTest) (const uint32 currentNode, const uint32 currentDataIndex, const SBVH8Database<T>* self, const SBVH8Ray& ray, SBVH8Hit& hit);

	explicit SBVH8Database() : m_isSorted(false), m_cellScale(1.f)
	{
		m_data.reserve(16384);
		m_dataLookup.reserve(16384);
	}

	// Appends a key at the end of the database, without any matching data.
	size_t AppendKey(const uint32 _cellKey)
	{
		size_t keyIndex = m_dataLookup.size();
		// Store illegal index to hint "no data" together with key
		m_dataLookup.emplace_back(BVH8NodeInfo(_cellKey, BVH8IllegalIndex));
		m_isSorted = false;
		return keyIndex; // Return index of appended entry
	}

	// Appends a data item and returns its index for the following 'bind' operation.
	size_t AppendData(T& _data)
	{
		uint32 dataindex = m_data.size();
		m_data.emplace_back(_data);
		m_isSorted = false;
		return dataindex;
	}

	// Attaches a manually appended key to a manually appended or existing data item.
	void BindKeyToData(const size_t _keyindex, const size_t _dataindex)
	{
		m_dataLookup[_keyindex].m_DataIndex = _dataindex;
		m_isSorted = false;
	}

	void SetGridAABB(const SVec128 &_aabbMin, const SVec128 &_aabbMax)
	{
		SVec128 tentwentythree = EVecConst(1023.f, 1023.f, 1023.f, 1.f);
		m_GridAABBMin = _aabbMin;
		m_GridCellSize = EVecDiv(EVecSub(_aabbMax,_aabbMin), tentwentythree);
	}

	EInline void ToGridUnits(const SVec128 &_worldPosition, SVec128 &_gridLocalPosition)
	{
		_gridLocalPosition = EVecDiv(EVecSub(_worldPosition, m_GridAABBMin), m_GridCellSize);
	}

	EInline void ToGridUnits(const SVec128 &_worldPosition, SVec128 &_gridLocalPosition, SVec128 _gridScale)
	{
		_gridLocalPosition = EVecDiv(EVecSub(_worldPosition, m_GridAABBMin), _gridScale);
	}

	EInline void ToWorldUnits(const SVec128 &_gridLocalPosition, SVec128 &_worldPosition)
	{
		_worldPosition = EVecAdd(EVecMul(_gridLocalPosition, m_GridCellSize), m_GridAABBMin);
	}

	EInline void QuantizePosition(const SVec128 &_worldpos, uint32 _qXYZ[3])
	{
		SVec128 gridLocalPosition;
		ToGridUnits(_worldpos, gridLocalPosition);

		// Clamp within 0-SSDGridEntryCountPerAxis-1 range (inclusive)
		static const SVec128 cellClampMax{ 1023.f, 1023.f, 1023.f, 0.f};
		static const SVec128 cellClampMin{ 0.f, 0.f, 0.f, 0.f };
		gridLocalPosition = EVecMin(EVecMax(gridLocalPosition, cellClampMin), cellClampMax);

		_qXYZ[0] = uint32(EFloor(EVecGetFloatX(gridLocalPosition)));
		_qXYZ[1] = uint32(EFloor(EVecGetFloatY(gridLocalPosition)));
		_qXYZ[2] = uint32(EFloor(EVecGetFloatZ(gridLocalPosition)));
	}

	EInline void QuantizePosition(const SVec128 &_worldpos, uint32 _qXYZ[3], SVec128 _gridScale)
	{
		SVec128 gridLocalPosition;
		ToGridUnits(_worldpos, gridLocalPosition, _gridScale);

		// Clamp within 0-SSDGridEntryCountPerAxis-1 range (inclusive)
		static const SVec128 cellClampMax{ 1023.f, 1023.f, 1023.f, 0.f};
		static const SVec128 cellClampMin{ 0.f, 0.f, 0.f, 0.f };
		gridLocalPosition = EVecMin(EVecMax(gridLocalPosition, cellClampMin), cellClampMax);

		_qXYZ[0] = uint32(EFloor(EVecGetFloatX(gridLocalPosition)));
		_qXYZ[1] = uint32(EFloor(EVecGetFloatY(gridLocalPosition)));
		_qXYZ[2] = uint32(EFloor(EVecGetFloatZ(gridLocalPosition)));
	}

	uint32 FindCell(const uint32 _cellKey, uint32 &_keyIndex)
	{
		auto foundItem = std::find_if(m_dataLookup.begin(), m_dataLookup.end(),
			[&_cellKey](BVH8NodeInfo const& item) { return item.m_SpatialKey == _cellKey; }
		);
		
		if (foundItem != m_dataLookup.end())
		{
			_keyIndex = foundItem - m_dataLookup.begin();
			return uint32(foundItem->m_DataIndex);
		}
		else
			_keyIndex = 0xFFFFFFFF;
		return BVH8IllegalIndex;
	}

	// Returns pointer to appended item.
	T* Append(const uint32 _cellKey, uint32& _dataIndex, uint32& _keyIndex)
	{
		// Store index at which 'data' is stored
		// (i.e. when storing a key with 2 bit metadata only, instead of actual data)
		_dataIndex = m_data.size();
		_keyIndex = m_dataLookup.size();
		m_dataLookup.emplace_back(BVH8NodeInfo(_cellKey, _dataIndex));
		m_data.emplace_back(T()); // Add empty item
		m_isSorted = false;
		return &m_data[_dataIndex];
	}

	// Must be called to be able to query the database.
	void SortAscending(const uint32 start_offset, const uint32 sort_length, bool update_lod = false)
	{
		if (sort_length<=1)
		{
			m_isSorted = true;
			return;
		}
		// Sort data lookup in the same order that keys would be sorted.
		// Ignores user bit and sorts only using the spatial location part.
		BVH8NodeInfo *m_tmpLookup = new BVH8NodeInfo[sort_length];
		ERadixSortSpatialDatabaseAscending(m_dataLookup.data()+start_offset, m_tmpLookup, sort_length);
		delete [] m_tmpLookup;

		if (update_lod)
		{
			m_LodStart[0] = 0;
			m_LodEnd[0] = m_dataLookup.size();
		}

		m_isSorted = true;
	}

	// Get data item at matching index
	T* GetDataItem(const uint32 _dataIndex) { return &m_data[_dataIndex]; }
	size_t GetDataItemCount() { return m_data.size(); }
	size_t GetDataLookupItemCount() { return m_dataLookup.size(); }

	// Encode key from position
	uint32 EncodeKey(SVec128i &_xyz) { return EMortonEncode(EVecGetIntX(_xyz), EVecGetIntY(_xyz), EVecGetIntZ(_xyz)); }
	uint32 EncodeKey(const uint32 _x, const uint32 _y, const uint32 _z) { return EMortonEncode(_x, _y, _z); }

	#define BIT(_x_) (32-(1 << (_x_)))
	uint32 EncodeKeyReverse(const uint32 _x, const uint32 _y, const uint32 _z)
	{
		uint32 mortonCode = 0;
		const uint32 numBits = 10;
		const uint32 numAxis = 3;
		uint32 coords[numAxis] = { _y, _x, _z };
		for (uint32 bitIndex = 0; bitIndex < numBits; bitIndex++)
		{
			for (uint32 axis = 0; axis < numAxis; axis++)
			{
				uint32 bit = BIT(bitIndex) & coords[axis];
				if (bit)
				{
					mortonCode |= BIT(bitIndex * numAxis + axis);
				}
			}
		}
		return mortonCode;
	}

	// Decode position from key
	void DecodeKey(const uint32 _cellKey, uint32 &_x, uint32 &_y, uint32 &_z)
	{
		EMortonDecode(_cellKey, _x, _y, _z);
	}

	void DecodeKey(const uint32 _cellKey, SVec128 &_xyz)
	{
		uint32 x, y, z;
		EMortonDecode(_cellKey, x, y, z);
		_xyz = EVecConst((float)x, (float)y, (float)z, 0);
	}

	void Clear()
	{
		m_data.clear();
		m_dataLookup.clear();
		m_isSorted = false;
	}

	void RemoveDuplicates()
	{
		m_dataLookup.erase( std::unique( m_dataLookup.begin(), m_dataLookup.end(), [](const BVH8NodeInfo& v1, const BVH8NodeInfo& v2) { return v1.m_SpatialKey == v2.m_SpatialKey; } ), m_dataLookup.end());
	}

	void GenerateBVH8()
	{
		uint32 startOffset = 0, fromLevel = 0;
		uint32 lowest_quality_reached = 0;
		while (!lowest_quality_reached)
		{
			startOffset = GenerateOctreeParentNodes(startOffset, fromLevel, lowest_quality_reached);
			++fromLevel;
		}
		m_RootBVH8Node = fromLevel;
	}

	uint32 GenerateOctreeParentNodes(const uint32 start_offset, const uint32 from_level, uint32 &lowest_quality_reached)
	{
		uint32 lod_start_offset = m_dataLookup.size();
		uint32 to_level = from_level + 1;
		uint32 i = start_offset; // Pass the start of LOD at 'from_level' here (higher nodes get pushed towards the end)

		// Mark range start
		m_LodStart[to_level] = lod_start_offset;

		while (i < lod_start_offset)
		{
			uint32 current_key = m_dataLookup[i].m_SpatialKey;

			// Strip low 3 bits to generate octree key for one level up
			uint32 octree_shared_key = current_key&0xFFFFFFF8;

			// Scan neighbourhood for keys that share these exact bits (they will belong to the same 'parent' node)
			// This is to figure out how many 'keys' to skip over to get to the next node that's not inside the parent cell
			uint32 last_matching_key = i;
			// Fit the bounding box to the first child node
			BVH8NodeInfo datalookup;
			datalookup.m_BoundsMin = EVecMin(datalookup.m_BoundsMin, m_dataLookup[i].m_BoundsMin);
			datalookup.m_BoundsMax = EVecMax(datalookup.m_BoundsMax, m_dataLookup[i].m_BoundsMax);
			uint32 rangeTop = i+1;
			uint32 rangeEnd = EMinimum(i+8,lod_start_offset);
			if (rangeTop == rangeEnd)
				break;
			for(uint32 j=rangeTop; j<rangeEnd; ++j)
			{
				uint32 sibling_key = m_dataLookup[j].m_SpatialKey;
				uint32 octree_sibling_key = sibling_key&0xFFFFFFF8;

				// This is a sibling as it shares the same exact parent node bits
				if (octree_shared_key == octree_sibling_key)
				{
					last_matching_key = j;
					datalookup.m_ChildCount++;
					// Fit the bounding box to the rest of the child nodes
					datalookup.m_BoundsMin = EVecMin(datalookup.m_BoundsMin, m_dataLookup[j].m_BoundsMin);
					datalookup.m_BoundsMax = EVecMax(datalookup.m_BoundsMax, m_dataLookup[j].m_BoundsMax);
				}
			}
			datalookup.m_ChildCount++; // One more for the initial child

			SVec128 centroid = EVecMul(EVecAdd(datalookup.m_BoundsMax, datalookup.m_BoundsMin), EVecConst(0.5f,0.5f,0.5f,0.f));
			uint32 qXYZ[3];
			float target_scale = powf(2.f, (float)to_level);
			QuantizePosition(centroid, qXYZ, EVecMul(m_GridCellSize, EVecConst(target_scale,target_scale,target_scale,0.f)));
			uint32 parent_key = EncodeKey(qXYZ[0], qXYZ[1], qXYZ[2]);

			// Generate a shared key for siblings
			// Here we shift out the last 3 bits to effectively reduce size by half compared to previous level
			//uint32 parent_key = (octree_shared_key>>3)&0x3FFFFFFF;
			datalookup.m_DataIndex = i; // Point at child node instead of data
			datalookup.m_SpatialKey = parent_key;
			m_dataLookup.emplace_back(datalookup);

			i = last_matching_key+1;
		}

		// Mark range end
		uint32 lod_end_offset = m_dataLookup.size();
		m_LodEnd[to_level] = lod_end_offset;

		// Sort to ensure higher keys move to the end
		SortAscending(lod_start_offset, lod_end_offset-lod_start_offset);

		lowest_quality_reached = (lod_end_offset - lod_start_offset) <= 1;

		// Return start offset of generated nodes
		return lod_start_offset;
	}
};
