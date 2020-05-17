#include "HelperFunctions.h"
int fnCallHierarchySummary[3][3]={{2,1,0},{0,4,4},{0,0,8}};
map<int, vector<FunctionCall*> > fnCalls;
extern int recursionDepth;
int fnCounter=0;
map<int,int> tileUpdateLog;
int computeGrid = COMPUTE_UTM;
extern long int inputSizex[2];

#define _paramListFunctionA Box*, CELL_TYPE*, int
#define _paramListFunctionB Box*, CELL_TYPE*, Box*, CELL_TYPE*, Box*, CELL_TYPE*, Box*, CELL_TYPE*, int
#define _paramListFunctionC Box*, CELL_TYPE*, Box*, CELL_TYPE*, Box*, CELL_TYPE*, Box*, CELL_TYPE*, int

void A(Box* x, CELL_TYPE* xData, int callStackDepth);
void B(Box* x, CELL_TYPE* xData, Box* y, CELL_TYPE* yData, Box* z, CELL_TYPE* zData, Box* a, CELL_TYPE* aData, int callStackDepth);
void C(Box* x, CELL_TYPE* xData, Box* y, CELL_TYPE* yData, Box* a, CELL_TYPE* aData, Box* z, CELL_TYPE* zData, int callStackDepth);
void A_unroll(Box* x, CELL_TYPE* xData, int parentTileIDx, int callStackDepth);
void B_unroll(Box* x, CELL_TYPE* xData, int parentTileIDx, Box* y, CELL_TYPE* yData, int parentTileIDy, Box* z, CELL_TYPE* zData, int parentTileIDz, Box* a, CELL_TYPE* aData, int parentTileIDa, int callStackDepth);
void C_unroll(Box* x, CELL_TYPE* xData, int parentTileIDx, Box* y, CELL_TYPE* yData, int parentTileIDy, Box* a, CELL_TYPE* aData, int parentTileIDa, Box* z, CELL_TYPE* zData, int parentTileIDz, int callStackDepth);
void A(Box* x, CELL_TYPE* xData, int callStackDepth)
{
	if((x->coords[0]==x->coords[2]) && (x->coords[1]==x->coords[3]))
	{
		//Write the code for terminating case here.
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

		return;
	}

	Box x00(x->coords[0],x->coords[1],x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x01(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1],x->coords[2],x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x11(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2],x->coords[3]);

	if(callStackDepth > recursionDepth+4)
	{


	A(&x00, xData, callStackDepth+1);	A(&x11, xData, callStackDepth+1);

	B(&x01, xData, &x00, xData, &x11, xData, &x11, xData, callStackDepth+1);

	return;

	}


#ifdef PARALLEL
	cilk_spawn
#endif
	A(&x00, xData, callStackDepth+1);	A(&x11, xData, callStackDepth+1);
#ifdef PARALLEL
	cilk_sync;
#endif

	B(&x01, xData, &x00, xData, &x11, xData, &x11, xData, callStackDepth+1);

	return;
}

void B(Box* x, CELL_TYPE* xData, Box* y, CELL_TYPE* yData, Box* z, CELL_TYPE* zData, Box* a, CELL_TYPE* aData, int callStackDepth)
{
	if((x->coords[0]==x->coords[2]) && (x->coords[1]==x->coords[3]))
	{
		//Write the code for terminating case here.
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

		return;
	}

	Box a00(a->coords[0],a->coords[1],a->coords[2]-(a->coords[2]-a->coords[0]+1)/2,a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box a01(a->coords[0]+(a->coords[2]-a->coords[0]+1)/2,a->coords[1],a->coords[2],a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box x00(x->coords[0],x->coords[1],x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x01(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1],x->coords[2],x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x10(x->coords[0],x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]);
	Box x11(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2],x->coords[3]);
	Box y00(y->coords[0],y->coords[1],y->coords[2]-(y->coords[2]-y->coords[0]+1)/2,y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y01(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1],y->coords[2],y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y11(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1]+(y->coords[3]-y->coords[1]+1)/2,y->coords[2],y->coords[3]);
	Box z00(z->coords[0],z->coords[1],z->coords[2]-(z->coords[2]-z->coords[0]+1)/2,z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z01(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1],z->coords[2],z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z11(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1]+(z->coords[3]-z->coords[1]+1)/2,z->coords[2],z->coords[3]);

	if(callStackDepth > recursionDepth+4)
	{

	B(&x10, xData, &y11, yData, &z00, zData, &a00, aData, callStackDepth+1);

	C(&x00, xData, &y01, yData, &a00, aData, &x10, xData, callStackDepth+1);	C(&x11, xData, &x10, xData, &z11, zData, &z01, zData, callStackDepth+1);


	B(&x00, xData, &y00, yData, &z00, zData, &x10, xData, callStackDepth+1);	B(&x11, xData, &y11, yData, &z11, zData, &a01, aData, callStackDepth+1);

	C(&x01, xData, &y01, yData, &a01, aData, &x11, xData, callStackDepth+1);
	C(&x01, xData, &x00, xData, &z11, zData, &z01, zData, callStackDepth+1);
	B(&x01, xData, &y00, yData, &z11, zData, &x11, xData, callStackDepth+1);

	return;

	}

	B(&x10, xData, &y11, yData, &z00, zData, &a00, aData, callStackDepth+1);

#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x00, xData, &y01, yData, &a00, aData, &x10, xData, callStackDepth+1);	C(&x11, xData, &x10, xData, &z11, zData, &z01, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_sync;
#endif


#ifdef PARALLEL
	cilk_spawn
#endif
	B(&x00, xData, &y00, yData, &z00, zData, &x10, xData, callStackDepth+1);	B(&x11, xData, &y11, yData, &z11, zData, &a01, aData, callStackDepth+1);
#ifdef PARALLEL
	cilk_sync;
#endif

	C(&x01, xData, &y01, yData, &a01, aData, &x11, xData, callStackDepth+1);
	C(&x01, xData, &x00, xData, &z11, zData, &z01, zData, callStackDepth+1);
	B(&x01, xData, &y00, yData, &z11, zData, &x11, xData, callStackDepth+1);

	return;
}

void C(Box* x, CELL_TYPE* xData, Box* y, CELL_TYPE* yData, Box* a, CELL_TYPE* aData, Box* z, CELL_TYPE* zData, int callStackDepth)
{
	if((x->coords[0]==x->coords[2]) && (x->coords[1]==x->coords[3]))
	{
		//Write the code for terminating case here.
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

		return;
	}

	Box a00(a->coords[0],a->coords[1],a->coords[2]-(a->coords[2]-a->coords[0]+1)/2,a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box a01(a->coords[0]+(a->coords[2]-a->coords[0]+1)/2,a->coords[1],a->coords[2],a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box x00(x->coords[0],x->coords[1],x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x01(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1],x->coords[2],x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x10(x->coords[0],x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]);
	Box x11(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2],x->coords[3]);
	Box y00(y->coords[0],y->coords[1],y->coords[2]-(y->coords[2]-y->coords[0]+1)/2,y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y01(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1],y->coords[2],y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y10(y->coords[0],y->coords[1]+(y->coords[3]-y->coords[1]+1)/2,y->coords[2]-(y->coords[2]-y->coords[0]+1)/2,y->coords[3]);
	Box y11(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1]+(y->coords[3]-y->coords[1]+1)/2,y->coords[2],y->coords[3]);
	Box z00(z->coords[0],z->coords[1],z->coords[2]-(z->coords[2]-z->coords[0]+1)/2,z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z01(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1],z->coords[2],z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z10(z->coords[0],z->coords[1]+(z->coords[3]-z->coords[1]+1)/2,z->coords[2]-(z->coords[2]-z->coords[0]+1)/2,z->coords[3]);
	Box z11(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1]+(z->coords[3]-z->coords[1]+1)/2,z->coords[2],z->coords[3]);

	if(callStackDepth > recursionDepth+4)
	{


	C(&x00, xData, &y00, yData, &z10, zData, &z00, zData, callStackDepth+1);
	C(&x01, xData, &y00, yData, &z11, zData, &z01, zData, callStackDepth+1);
	C(&x10, xData, &y10, yData, &z10, zData, &z00, zData, callStackDepth+1);	C(&x11, xData, &y10, yData, &z11, zData, &z01, zData, callStackDepth+1);


	C(&x00, xData, &y01, yData, &a00, aData, &z10, zData, callStackDepth+1);
	C(&x01, xData, &y01, yData, &a01, aData, &z11, zData, callStackDepth+1);
	C(&x10, xData, &y11, yData, &a00, aData, &z10, zData, callStackDepth+1);	C(&x11, xData, &y11, yData, &a01, aData, &z11, zData, callStackDepth+1);


	return;

	}


#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x00, xData, &y00, yData, &z10, zData, &z00, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x01, xData, &y00, yData, &z11, zData, &z01, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x10, xData, &y10, yData, &z10, zData, &z00, zData, callStackDepth+1);	C(&x11, xData, &y10, yData, &z11, zData, &z01, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_sync;
#endif


#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x00, xData, &y01, yData, &a00, aData, &z10, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x01, xData, &y01, yData, &a01, aData, &z11, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_spawn
#endif
	C(&x10, xData, &y11, yData, &a00, aData, &z10, zData, callStackDepth+1);	C(&x11, xData, &y11, yData, &a01, aData, &z11, zData, callStackDepth+1);
#ifdef PARALLEL
	cilk_sync;
#endif


	return;
}

void A_unroll(Box* x, CELL_TYPE* xData, int parentTileIDx, int callStackDepth)
{
	if(callStackDepth == recursionDepth)
	{
		int writeTileID = GetTileID2D(parentTileIDx);
		bool localUpdate = false;
		FunctionCall* fnCall= NULL;
		tileUpdateLog[writeTileID]=fnCounter;
		if(IsLocalOwner(writeTileID))
		{
			CELL_TYPE* nullData=NULL;
			fnCall=new FunctionCall();
			fnCall->functionName = 'A';
			Parameter p;
			Box* b0=new Box(*x);
			p.data = b0;
			p.tile = xData;
			fnCall->params.push_back(p);
			std::tuple<_paramListFunctionA> t = std::make_tuple(b0, xData, callStackDepth+1);
			DeferredCall<_paramListFunctionA>* defdCall = new DeferredCall<_paramListFunctionA>();
			defdCall->params=t;
			defdCall->fptr=A;
			fnCall->fnCall = defdCall;
			fnCall->ID = fnCounter;
			if(fnCalls[writeTileID].size() > 0)
			{
				fnCall->numInPorts +=1;
				fnCall->wawSource = (fnCalls[writeTileID].back())->ID;
				FunctionCall* lastFunctionToUpdate = fnCalls[writeTileID].back();
				lastFunctionToUpdate->outPortOwners.insert(procRank);
			}
			fnCalls[writeTileID].push_back(fnCall);
			localUpdate = true;
		}
		fnCounter++;
		return;
	}

	Box x00(x->coords[0],x->coords[1],x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x01(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1],x->coords[2],x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x11(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2],x->coords[3]);

	A_unroll(&x00, xData, parentTileIDx*4+0, callStackDepth+1);	A_unroll(&x11, xData, parentTileIDx*4+3, callStackDepth+1);
	B_unroll(&x01, xData, parentTileIDx*4+1, &x00, xData, parentTileIDx*4+0, &x11, xData, parentTileIDx*4+3, &x11, xData, parentTileIDx*4+3, callStackDepth+1);

	return;
}

void B_unroll(Box* x, CELL_TYPE* xData, int parentTileIDx, Box* y, CELL_TYPE* yData, int parentTileIDy, Box* z, CELL_TYPE* zData, int parentTileIDz, Box* a, CELL_TYPE* aData, int parentTileIDa, int callStackDepth)
{
	if(callStackDepth == recursionDepth)
	{
		int writeTileID = GetTileID2D(parentTileIDx);
		bool localUpdate = false;
		FunctionCall* fnCall= NULL;
		tileUpdateLog[writeTileID]=fnCounter;
		if(IsLocalOwner(writeTileID))
		{
			CELL_TYPE* nullData=NULL;
			fnCall=new FunctionCall();
			fnCall->functionName = 'B';
			Parameter p;
			Box* b0=new Box(*x);
			p.data = b0;
			p.tile = xData;
			fnCall->params.push_back(p);
			Box* b3=new Box(*y);
			p.data = b3;
			p.tile = yData;
			fnCall->params.push_back(p);
			Box* b6=new Box(*z);
			p.data = b6;
			p.tile = zData;
			fnCall->params.push_back(p);
			Box* b9=new Box(*a);
			p.data = b9;
			p.tile = aData;
			fnCall->params.push_back(p);
			std::tuple<_paramListFunctionB> t = std::make_tuple(b0, xData, b3, yData, b6, zData, b9, aData, callStackDepth+1);
			DeferredCall<_paramListFunctionB>* defdCall = new DeferredCall<_paramListFunctionB>();
			defdCall->params=t;
			defdCall->fptr=B;
			fnCall->fnCall = defdCall;
			fnCall->ID = fnCounter;
			if(fnCalls[writeTileID].size() > 0)
			{
				fnCall->numInPorts +=1;
				fnCall->wawSource = (fnCalls[writeTileID].back())->ID;
				FunctionCall* lastFunctionToUpdate = fnCalls[writeTileID].back();
				lastFunctionToUpdate->outPortOwners.insert(procRank);
			}
			fnCalls[writeTileID].push_back(fnCall);
			localUpdate = true;
		}
		fnCounter++;
		int readTileIDs[3];
		CELL_TYPE* tiles[3];
		readTileIDs[0]=GetTileID2D(parentTileIDy);
		readTileIDs[1]=GetTileID2D(parentTileIDz);
		readTileIDs[2]=GetTileID2D(parentTileIDa);
		tiles[0]=yData;
		tiles[1]=zData;
		tiles[2]=aData;
		for(int i=0;i<3;i++)
		{
			int readTileID=readTileIDs[i];
			CELL_TYPE* curTile=tiles[i];
			if((curTile==NULL) && (readTileID != writeTileID))
			{
#ifdef TASK_AGGREGATION
				if(fnCalls[readTileID].size() > 0)
					(fnCalls[readTileID].back())->isReadBeforeNextWrite = true;
#endif
				if(localUpdate)
				{
					fnCall->numInPorts +=1;
					fnCall->params[i+1].portID = tileUpdateLog[readTileID];
				}
				if(fnCalls[readTileID].size() > 0)
				{
					FunctionCall* lastFunctionToUpdate = fnCalls[readTileID].back();
					lastFunctionToUpdate->outPortOwners.insert(GetOwner(writeTileID));
				}
			}
			else if((curTile!=NULL) && localUpdate)
				fnCall->params[i+1].portID = readTileID;
		}
		return;
	}

	Box a00(a->coords[0],a->coords[1],a->coords[2]-(a->coords[2]-a->coords[0]+1)/2,a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box a01(a->coords[0]+(a->coords[2]-a->coords[0]+1)/2,a->coords[1],a->coords[2],a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box x00(x->coords[0],x->coords[1],x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x01(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1],x->coords[2],x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x10(x->coords[0],x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]);
	Box x11(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2],x->coords[3]);
	Box y00(y->coords[0],y->coords[1],y->coords[2]-(y->coords[2]-y->coords[0]+1)/2,y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y01(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1],y->coords[2],y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y11(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1]+(y->coords[3]-y->coords[1]+1)/2,y->coords[2],y->coords[3]);
	Box z00(z->coords[0],z->coords[1],z->coords[2]-(z->coords[2]-z->coords[0]+1)/2,z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z01(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1],z->coords[2],z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z11(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1]+(z->coords[3]-z->coords[1]+1)/2,z->coords[2],z->coords[3]);

	B_unroll(&x10, xData, parentTileIDx*4+2, &y11, yData, parentTileIDy*4+3, &z00, zData, parentTileIDz*4+0, &a00, aData, parentTileIDa*4+0, callStackDepth+1);
	C_unroll(&x00, xData, parentTileIDx*4+0, &y01, yData, parentTileIDy*4+1, &a00, aData, parentTileIDa*4+0, &x10, xData, parentTileIDx*4+2, callStackDepth+1);	C_unroll(&x11, xData, parentTileIDx*4+3, &x10, xData, parentTileIDx*4+2, &z11, zData, parentTileIDz*4+3, &z01, zData, parentTileIDz*4+1, callStackDepth+1);
	B_unroll(&x00, xData, parentTileIDx*4+0, &y00, yData, parentTileIDy*4+0, &z00, zData, parentTileIDz*4+0, &x10, xData, parentTileIDx*4+2, callStackDepth+1);	B_unroll(&x11, xData, parentTileIDx*4+3, &y11, yData, parentTileIDy*4+3, &z11, zData, parentTileIDz*4+3, &a01, aData, parentTileIDa*4+1, callStackDepth+1);
	C_unroll(&x01, xData, parentTileIDx*4+1, &y01, yData, parentTileIDy*4+1, &a01, aData, parentTileIDa*4+1, &x11, xData, parentTileIDx*4+3, callStackDepth+1);
	C_unroll(&x01, xData, parentTileIDx*4+1, &x00, xData, parentTileIDx*4+0, &z11, zData, parentTileIDz*4+3, &z01, zData, parentTileIDz*4+1, callStackDepth+1);
	B_unroll(&x01, xData, parentTileIDx*4+1, &y00, yData, parentTileIDy*4+0, &z11, zData, parentTileIDz*4+3, &x11, xData, parentTileIDx*4+3, callStackDepth+1);

	return;
}

void C_unroll(Box* x, CELL_TYPE* xData, int parentTileIDx, Box* y, CELL_TYPE* yData, int parentTileIDy, Box* a, CELL_TYPE* aData, int parentTileIDa, Box* z, CELL_TYPE* zData, int parentTileIDz, int callStackDepth)
{
	if(callStackDepth == recursionDepth)
	{
		int writeTileID = GetTileID2D(parentTileIDx);
		bool localUpdate = false;
		FunctionCall* fnCall= NULL;
		tileUpdateLog[writeTileID]=fnCounter;
		if(IsLocalOwner(writeTileID))
		{
			CELL_TYPE* nullData=NULL;
			fnCall=new FunctionCall();
			fnCall->functionName = 'C';
			Parameter p;
			Box* b0=new Box(*x);
			p.data = b0;
			p.tile = xData;
			fnCall->params.push_back(p);
			Box* b3=new Box(*y);
			p.data = b3;
			p.tile = yData;
			fnCall->params.push_back(p);
			Box* b6=new Box(*a);
			p.data = b6;
			p.tile = aData;
			fnCall->params.push_back(p);
			Box* b9=new Box(*z);
			p.data = b9;
			p.tile = zData;
			fnCall->params.push_back(p);
			std::tuple<_paramListFunctionC> t = std::make_tuple(b0, xData, b3, yData, b6, aData, b9, zData, callStackDepth+1);
			DeferredCall<_paramListFunctionC>* defdCall = new DeferredCall<_paramListFunctionC>();
			defdCall->params=t;
			defdCall->fptr=C;
			fnCall->fnCall = defdCall;
			fnCall->ID = fnCounter;
			if(fnCalls[writeTileID].size() > 0)
			{
				fnCall->numInPorts +=1;
				fnCall->wawSource = (fnCalls[writeTileID].back())->ID;
				FunctionCall* lastFunctionToUpdate = fnCalls[writeTileID].back();
				lastFunctionToUpdate->outPortOwners.insert(procRank);
			}
			fnCalls[writeTileID].push_back(fnCall);
			localUpdate = true;
		}
		fnCounter++;
		int readTileIDs[3];
		CELL_TYPE* tiles[3];
		readTileIDs[0]=GetTileID2D(parentTileIDy);
		readTileIDs[1]=GetTileID2D(parentTileIDa);
		readTileIDs[2]=GetTileID2D(parentTileIDz);
		tiles[0]=yData;
		tiles[1]=aData;
		tiles[2]=zData;
		for(int i=0;i<3;i++)
		{
			int readTileID=readTileIDs[i];
			CELL_TYPE* curTile=tiles[i];
			if((curTile==NULL) && (readTileID != writeTileID))
			{
#ifdef TASK_AGGREGATION
				if(fnCalls[readTileID].size() > 0)
					(fnCalls[readTileID].back())->isReadBeforeNextWrite = true;
#endif
				if(localUpdate)
				{
					fnCall->numInPorts +=1;
					fnCall->params[i+1].portID = tileUpdateLog[readTileID];
				}
				if(fnCalls[readTileID].size() > 0)
				{
					FunctionCall* lastFunctionToUpdate = fnCalls[readTileID].back();
					lastFunctionToUpdate->outPortOwners.insert(GetOwner(writeTileID));
				}
			}
			else if((curTile!=NULL) && localUpdate)
				fnCall->params[i+1].portID = readTileID;
		}
		return;
	}

	Box a00(a->coords[0],a->coords[1],a->coords[2]-(a->coords[2]-a->coords[0]+1)/2,a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box a01(a->coords[0]+(a->coords[2]-a->coords[0]+1)/2,a->coords[1],a->coords[2],a->coords[3]-(a->coords[3]-a->coords[1]+1)/2);
	Box x00(x->coords[0],x->coords[1],x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x01(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1],x->coords[2],x->coords[3]-(x->coords[3]-x->coords[1]+1)/2);
	Box x10(x->coords[0],x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2]-(x->coords[2]-x->coords[0]+1)/2,x->coords[3]);
	Box x11(x->coords[0]+(x->coords[2]-x->coords[0]+1)/2,x->coords[1]+(x->coords[3]-x->coords[1]+1)/2,x->coords[2],x->coords[3]);
	Box y00(y->coords[0],y->coords[1],y->coords[2]-(y->coords[2]-y->coords[0]+1)/2,y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y01(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1],y->coords[2],y->coords[3]-(y->coords[3]-y->coords[1]+1)/2);
	Box y10(y->coords[0],y->coords[1]+(y->coords[3]-y->coords[1]+1)/2,y->coords[2]-(y->coords[2]-y->coords[0]+1)/2,y->coords[3]);
	Box y11(y->coords[0]+(y->coords[2]-y->coords[0]+1)/2,y->coords[1]+(y->coords[3]-y->coords[1]+1)/2,y->coords[2],y->coords[3]);
	Box z00(z->coords[0],z->coords[1],z->coords[2]-(z->coords[2]-z->coords[0]+1)/2,z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z01(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1],z->coords[2],z->coords[3]-(z->coords[3]-z->coords[1]+1)/2);
	Box z10(z->coords[0],z->coords[1]+(z->coords[3]-z->coords[1]+1)/2,z->coords[2]-(z->coords[2]-z->coords[0]+1)/2,z->coords[3]);
	Box z11(z->coords[0]+(z->coords[2]-z->coords[0]+1)/2,z->coords[1]+(z->coords[3]-z->coords[1]+1)/2,z->coords[2],z->coords[3]);

	C_unroll(&x00, xData, parentTileIDx*4+0, &y00, yData, parentTileIDy*4+0, &z10, zData, parentTileIDz*4+2, &z00, zData, parentTileIDz*4+0, callStackDepth+1);	C_unroll(&x01, xData, parentTileIDx*4+1, &y00, yData, parentTileIDy*4+0, &z11, zData, parentTileIDz*4+3, &z01, zData, parentTileIDz*4+1, callStackDepth+1);	C_unroll(&x10, xData, parentTileIDx*4+2, &y10, yData, parentTileIDy*4+2, &z10, zData, parentTileIDz*4+2, &z00, zData, parentTileIDz*4+0, callStackDepth+1);	C_unroll(&x11, xData, parentTileIDx*4+3, &y10, yData, parentTileIDy*4+2, &z11, zData, parentTileIDz*4+3, &z01, zData, parentTileIDz*4+1, callStackDepth+1);
	C_unroll(&x00, xData, parentTileIDx*4+0, &y01, yData, parentTileIDy*4+1, &a00, aData, parentTileIDa*4+0, &z10, zData, parentTileIDz*4+2, callStackDepth+1);	C_unroll(&x01, xData, parentTileIDx*4+1, &y01, yData, parentTileIDy*4+1, &a01, aData, parentTileIDa*4+1, &z11, zData, parentTileIDz*4+3, callStackDepth+1);	C_unroll(&x10, xData, parentTileIDx*4+2, &y11, yData, parentTileIDy*4+3, &a00, aData, parentTileIDa*4+0, &z10, zData, parentTileIDz*4+2, callStackDepth+1);	C_unroll(&x11, xData, parentTileIDx*4+3, &y11, yData, parentTileIDy*4+3, &a01, aData, parentTileIDa*4+1, &z11, zData, parentTileIDz*4+3, callStackDepth+1);

	return;
}

void Unroll()
{
	int parentTileIDx=0;
	Box *x = new Box(0,0,inputSizex[0],inputSizex[1]);
	if(0 == recursionDepth)
	{
		int writeTileID = GetTileID2D(parentTileIDx);
		bool localUpdate = false;
		FunctionCall* fnCall= NULL;
		tileUpdateLog[writeTileID]=fnCounter;
		if(IsLocalOwner(writeTileID))
		{
			CELL_TYPE* nullData=NULL;
			fnCall=new FunctionCall();
			fnCall->functionName = 'A';
			Parameter p;
			Box* b0=new Box(*x);
			p.data = b0;
			p.tile = nullData;
			fnCall->params.push_back(p);
			std::tuple<_paramListFunctionA> t = std::make_tuple(b0, nullData, 0+1);
			DeferredCall<_paramListFunctionA>* defdCall = new DeferredCall<_paramListFunctionA>();
			defdCall->params=t;
			defdCall->fptr=A;
			fnCall->fnCall = defdCall;
			fnCall->ID = fnCounter;
			if(fnCalls[writeTileID].size() > 0)
			{
				fnCall->numInPorts +=1;
				fnCall->wawSource = (fnCalls[writeTileID].back())->ID;
				FunctionCall* lastFunctionToUpdate = fnCalls[writeTileID].back();
				lastFunctionToUpdate->outPortOwners.insert(procRank);
			}
			fnCalls[writeTileID].push_back(fnCall);
			localUpdate = true;
		}
		fnCounter++;
		return;
	}
	A_unroll(x, NULL , 0, 0);
	delete x;
}

void ExecuteFunction(FunctionCall* fn, vector<CELL_TYPE*> dataRegions)
{
	switch(fn->functionName)
	{
		case 'A':
			{
				assert(dataRegions.size() == 1);
				DeferredCall<_paramListFunctionA>* df = (reinterpret_cast<DeferredCall<_paramListFunctionA>* >(fn->fnCall));
				std::get<1>(df->params) = dataRegions[0];
				(reinterpret_cast<DeferredCall<_paramListFunctionA>* >(fn->fnCall))->Run();
			}
			break;
		case 'B':
			{
				assert(dataRegions.size() == 4);
				DeferredCall<_paramListFunctionB>* df = (reinterpret_cast<DeferredCall<_paramListFunctionB>* >(fn->fnCall));
				std::get<1>(df->params) = dataRegions[0];
				std::get<3>(df->params) = dataRegions[1];
				std::get<5>(df->params) = dataRegions[2];
				std::get<7>(df->params) = dataRegions[3];
				(reinterpret_cast<DeferredCall<_paramListFunctionB>* >(fn->fnCall))->Run();
			}
			break;
		case 'C':
			{
				assert(dataRegions.size() == 4);
				DeferredCall<_paramListFunctionC>* df = (reinterpret_cast<DeferredCall<_paramListFunctionC>* >(fn->fnCall));
				std::get<1>(df->params) = dataRegions[0];
				std::get<3>(df->params) = dataRegions[1];
				std::get<5>(df->params) = dataRegions[2];
				std::get<7>(df->params) = dataRegions[3];
				(reinterpret_cast<DeferredCall<_paramListFunctionC>* >(fn->fnCall))->Run();
			}
			break;
		default: break;
	}
}

