#pragma once
#include<vector>
#include<string>
#include<sstream>
#define DIMENSION 2
#define METADATASPACE (2*DIMENSION+1)
typedef int CELL_TYPE;
extern int fnCallHierarchySummary[3][3];
inline CELL_TYPE* GetDPTableCell(int i, int j, CELL_TYPE* data)
{
	CELL_TYPE* cell = data+METADATASPACE;
	int side = data[DIMENSION];
	int iOffset=i - data[1];
	cell += (iOffset*1* side);
	int jOffset=j - data[0];
	cell += jOffset;
	return cell;
}

