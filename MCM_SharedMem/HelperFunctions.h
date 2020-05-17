#pragma once
#include<stdio.h>
#include<vector>
#include<set>
#include<iostream>
#include<cstring>
#include<cstdlib>
#include<fstream>
#include<algorithm>
#include<climits>
#include<float.h>
#include<boost/any.hpp>
#include<tuple>
#include<cassert>
#include<math.h>
#include<map>
#include<cstdarg>
#include<sstream>
#ifdef PARALLEL
#include<cilk/cilk.h>
#include<cilk/cilk_api.h>
#endif
#include"GlobalTypes.h"

using namespace std;

#define TASK_AGGREGATION

typedef enum ComputeGrid{
COMPUTE_FULL,
COMPUTE_UTM,
COMPUTE_SPARSE
}ComputeGrid;

typedef enum PartitioningScheme{
PARTITIONING_BLOCKED_VERTICAL,
PARTITIONING_BLOCKCYCLIC_HORIZONTAL,
PARTITIONING_BLOCKED_HORIZONTAL,
PARTITIONING_BLOCKCYCLIC_VERTICAL,
PARTITIONING_TILED
}PartitioningScheme;
extern int dataPartitioningScheme;
typedef struct Parameter
{
	int portID;	
	void* data;	
	CELL_TYPE* tile;
	bool operator==(const int k) const
	{
		return ((this->portID == k)?true:false);
	}
	Parameter():portID(-1),tile(NULL){}
}Parameter;

typedef struct FunctionCall
{
	int ID;
	char functionName;
	void* fnCall;
	vector<Parameter> params;
	set<int> outPortOwners;
	int numInPorts;
	int wawSource;
#ifdef TASK_AGGREGATION
	bool isReadBeforeNextWrite;
	FunctionCall():numInPorts(0),wawSource(-1), isReadBeforeNextWrite(false){}
#else
	FunctionCall():numInPorts(0),wawSource(-1){}
#endif
	void GetSignature(char* signature);
}FunctionCall;

class Task
{
	public:
	vector<FunctionCall*> functionCalls;
	vector<pair<bool,CELL_TYPE*> > recvBuffers;
	set<int> outPortOwners;
	int numInPorts;
	CELL_TYPE* updatedRegion;
	int updateLen;
	Task* next;
#ifdef TASK_AGGREGATION
	Task* aggrSafeTask;
	Task():numInPorts(0),next(NULL),updatedRegion(NULL),updateLen(0),aggrSafeTask(NULL){}
#else
	Task():numInPorts(0),next(NULL),updatedRegion(NULL),updateLen(0){}
#endif
};


/*typedef struct Box
{
	int  coords[DIMENSION*2]; //topleft and bottom right (x and y coordinates)
	Box(int coord0...)
	{
		va_list argList;
		coords[0]=coord0;
		va_start(argList,coord0);
		for(int argNum=1;argNum<2*DIMENSION;argNum++)
			coords[argNum]=va_arg(argList,int); 
		va_end(argList);
	}
	Box(){}
	void PrintStr(char *str)const 
	{
		ostringstream ostrstrm;
		for(int i=0;i<2*DIMENSION;i++)
		{
			ostrstrm<<coords[i];
			if(i!= (2*DIMENSION-1))
				ostrstrm<<string(" ");
		}
		//cout<<ostrstrm.str();
		strcpy(str,ostrstrm.str().c_str());
		//sprintf(str,"%d %d %d %d",coords[0],coords[1],coords[2],coords[3],coords[4],coords[5],coords[6],coords[7]);
	}
	Box& operator=(const Box& rhs)
	{
		for(int i = 0;i<DIMENSION*2;i++)
			this->coords[i]=rhs.coords[i];
		return *this;
	}
	Box(const Box& rhs)
	{
		for(int i = 0;i<DIMENSION*2;i++)
			this->coords[i]=rhs.coords[i];
	}
	bool operator==(const Box& rhs)
	{
		bool flag = true;
		for(int i = 0;i<DIMENSION*2;i++)
		{
			if(this->coords[i]!=rhs.coords[i])
			{
				flag = false;
				break;
			}
		}
		return flag;
	}
	long int GetBoxSize()
	{
		long int len = 1;
		for(int i=0;i<DIMENSION;i++)
		{
			len *= (coords[i+DIMENSION]-coords[i]+1);
		}	
		return len;
	}
}Box;*/

void ParseSystemParams(int argc, char* argv[]);
void GenerateTasks(int argc, char** argv);
void ExecuteTasks();
int  PrintGlobalMax();
void GetAggregate(long int partialSum);
bool IsLocalOwner(int tileID);
int GetOwner(int tileID);
int GetTileID1D(int zmortTileID);
int GetTileID2D(int zmortTileID);
int GetTileID3D(unsigned long zmortTileID);
int GetTileID(const Box* tile);
CELL_TYPE* GetTileOfDPTable(int i, int j);
#ifdef DEBUG
CELL_TYPE* GetLocalTile(int i);
CELL_TYPE* GetTileOfDPTable(int zMortonID);
void DebugPrintTileDetails();
void DebugPrintAllCells();
int GetZMortonTileIDFromCellCoords(int i, int j);
#endif



int ReadInput(int argc, char** argv);
void InitializeDPTable(Box* b, CELL_TYPE* data);
void PrintResults();
void Unroll();
void ExecuteFunction(FunctionCall* fn, vector<CELL_TYPE*> dataRegions);

template<int ...> struct Sequence{};
template<int N, int ...S> struct Generate_Seq : Generate_Seq<N-1, N-1, S...> {};
template<int ...S> struct Generate_Seq<0, S...>{ typedef Sequence<S...> type; }; //base case of recursive class template definition
template <typename ...Args> class DeferredCall {
public:
	std::tuple<Args...> params;
	void (*fptr)(Args...);
	void Run()
	{
		_Run(typename Generate_Seq<sizeof...(Args)>::type());
		return;
	}

	template<int ...S>
	void _Run(Sequence<S...>)
	{
		fptr(std::get<S>(params)...);
		return;
	}
};

extern int procRank;
extern int totalProcs;
extern long gPreprocessTime;
extern long gTaskExecutionTime;
extern int computeGrid;

//#define GRAPHVIZ_OUTPUT
#ifdef GRAPHVIZ_OUTPUT
void UpdateWorkSpanAnalysis(const int curFunctionID, const float curFunctionExecTime);
#endif

//#define IDLE_TIME
//#define PAPI
#ifdef PAPI
#include"papi.h"
extern int papiStat, eventSet;
extern long long values[2];
void handle_error(int retval);
#endif

extern int lNodewiseMsgLog[64], gNodewiseMsgLog[64];
//#define REMOTENODEFIRST
