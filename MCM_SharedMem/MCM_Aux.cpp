#include"HelperFunctions.h"
long int inputSizex[DIMENSION];
vector<int> matrixChain;
void _ReadInput(char* fileName, vector<int>& matrixChain, int maxChainLength)
{
	char matrixDimStr[8];
	int matrixDim;
	ifstream inFile(fileName, ifstream::in);
	if(!fileName || !inFile.is_open())
	{
		printf("ERROR: Unable to open input file\n");
		exit(0);
	}
	while(!inFile.eof())
	{
		inFile.getline(matrixDimStr,8);
		if(inFile.eof())
			break;
		matrixDim = atoi(matrixDimStr);
		memset(matrixDimStr,0,8);
		matrixChain.push_back(matrixDim);
		if(matrixChain.size() == maxChainLength)
			break;
	}
}

int ReadInput(int argc, char* argv[])
{
	if((strcmp(argv[1],"-h")==0)||(argc==1))
	{
		printf("Usage: ./<exe> <input> <chainLength>\n");
		exit(0);
	}
	_ReadInput(argv[1],matrixChain, atoi(argv[2]));
	inputSizex[0] = matrixChain.size()-1;inputSizex[1] = matrixChain.size()-1;
	return inputSizex[0];
}

void InitializeDPTable(Box* b, CELL_TYPE* data)
{
	int numCells = b->GetBoxSize();
	unsigned int relIndx = 0;
	for(unsigned int i=b->coords[1];i<=b->coords[3];i++)
		for(unsigned int j=b->coords[0];j<=b->coords[2];j++)
			data[relIndx++] = (i==j)?0:INT_MAX;
	return;
}

void PrintResults()
{
	vector<int> coords;
	coords.push_back(inputSizex[0]);
	coords.push_back(1);
	CELL_TYPE* data = GetTileOfDPTable(inputSizex[0],1);
	if(data != NULL)
	{
		CELL_TYPE* result = GetDPTableCell(1,inputSizex[0],data);
		if(result!=NULL) printf("Optimum cost:%d\n",*result);
	}
}

/*
 * Terminating case for the MCM problem. Copy paste the code between == into the terminating case of a recursive method in RecursiveFunctions.cpp
 *
 		
==================Terminating case for method A (no dependencies. cell is computed from initial conditions.)=============================
		int i=x->coords[1];
		int j=x->coords[0];
		int k=x->coords[0];
		if(i == j)
			return;

		CELL_TYPE w_ikj=matrixChain[i-1]*matrixChain[k]*matrixChain[j];
		CELL_TYPE* cost_ij = GetDPTableCell(i,j,xData);
		CELL_TYPE* cost_ik = GetDPTableCell(i,k,xData);
		CELL_TYPE* cost_kj = GetDPTableCell(k+1,j,xData);
		CELL_TYPE newCost = *cost_ik + *cost_kj + w_ikj;
		if(newCost < *cost_ij)
			*cost_ij =  newCost;
==================Terminating case for methods B and C (reads two cells that yData and zData point to.)=========================================
		int i=x->coords[1];
		int j=x->coords[0];
		int k=y->coords[0];
		if(i == j)
			return;
		int* kjData = NULL;
		if(((k+1)>z->coords[3]) || ((k+1)<z->coords[1]) || (j>z->coords[2]) || (j<z->coords[0])) 
		{
			if(!(((k+1)>a->coords[3]) || ((k+1)<a->coords[1]) || (j>a->coords[2]) || (j<a->coords[0]))) 
				kjData = aData;
			else if(!(((k+1)>x->coords[3]) || ((k+1)<x->coords[1]) || (j>x->coords[2]) || (j<x->coords[0]))) 
				kjData = xData;
		}
		else
			kjData = zData;
		if(!kjData)
			return;

		CELL_TYPE w_ikj=matrixChain[i-1]*matrixChain[k]*matrixChain[j];
		CELL_TYPE* cost_ij = GetDPTableCell(i,j,xData);
		CELL_TYPE* cost_ik = GetDPTableCell(i,k,yData);
		CELL_TYPE* cost_kj = GetDPTableCell(k+1,j,kjData);
		CELL_TYPE newCost = *cost_ik + *cost_kj + w_ikj;
		if(newCost < *cost_ij)
			*cost_ij =  newCost;
========================================================================================================================================


 */
