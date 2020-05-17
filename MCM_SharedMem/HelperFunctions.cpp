#include<cassert>
#include<typeinfo>
#include<cstdint>

#include<sys/time.h>

#include "HelperFunctions.h"
#define MSG_DATA 0
#define MSG_EXIT 1

// Macro defined for debugging 
#define DEBUG 1

#ifdef GRAPHVIZ_OUTPUT
#include <algorithm>
#include <boost/graph/properties.hpp>
#include <boost/graph/directed_graph.hpp> // A subclass to provide reasonable arguments to adjacency_list for a typical directed graph
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
struct Vertex
{
    string name;
    int ID;
};

struct Edge
{
    string name;
};

typedef boost::directed_graph<Vertex, Edge> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_desc;
typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
vector<int> CheckIfAllPathsAreEnumerated(const vector<int>& pathLength);

vector<vector<int> > criticalPaths;
vector<double> longestPathExecTimes;
double T1=0., Tinf=0.;
#endif

#ifdef IDLE_TIME
#include<sys/time.h>
//FILE* idleTimeLog;
double idleTime;
int numMsgsRecvd=0;
#endif

int numMsgsSent=0;
#ifdef DEBUG
int numMsgsSaved=0;
#endif

int dataPartitioningScheme;
int cellsPerRowTile;
int tilesPerRowDPTable;
int tilesPerRowBlock;
int tilesPerColumnBlock;
//the following 3 are used in block cyclic partitioning.
float cycleLength;
int strideLength;
int strideOffset;

int numTilesToProcess;
int startingTileID;
int firstTileIDOfSmallerSizedChunks;
int avgNumTilesToProcess;
int remainingTiles;

#ifdef PARALLEL
char numCilkWorkers[4];
#endif
int recursionDepth;
map<long long,int> portReferenceCount;
extern map<int,vector<FunctionCall*> > fnCalls;
map<int,CELL_TYPE*> dpTileTable;
Task* tasks=NULL, *endTask=NULL; //head and tail of a linked list of tasks.
int numTasks; //number of elements in the list of tasks.


//local function declarations
int DetermineUnfoldingLevel(int numProcesses);
void GetTaskList();
void DistributeWork_DataParallel();
void CreateDataRegionsAndTasks();
void DoBackgroundWork();
void DoLocalComputation(Task& task);
bool ProcessAsyncRequest(int& numProcsReadyToExit);
bool UpdateTaskDependencies(CELL_TYPE* recvBuffer, bool localData=false, CELL_TYPE skipTask = -1);
void SendData(CELL_TYPE* data, int dataLen, int destProcID);
#ifdef DEBUG
unsigned int EncodeMorton(unsigned int x, unsigned int y);
#endif


/* This function sets system parameters from command line inputs 
 * recursion unrolling depth, partitioning scheme, number of threads (if compiled with PARALLEL mode), cycle_length (if using BLOCKCYCLIC_XX partitioning schemes) can be set. */
void ParseSystemParams(int argc, char* argv[])
{
	//setting default values.
	// Dead code, since we're not attempting a parallel execution 
#ifdef PARALLEL
	strcpy(numCilkWorkers,"1"); 
#endif
	cycleLength = 1; 
	recursionDepth = DetermineUnfoldingLevel(totalProcs); 
	dataPartitioningScheme = PARTITIONING_BLOCKCYCLIC_VERTICAL;
	
	//checking command line arguments to set system parameters.
	for (int i = 1; i < argc; i++) {
		char* arg = argv[i];
		if (strcmp(arg, "-recursion_depth") == 0) {
			recursionDepth = atoi(argv[++i]);
			if(recursionDepth < 0) {
				fprintf(stderr, "error: invalid recursion depth.\n");
				exit(0);
			}
		} else if (strcmp(arg, "-partition") == 0) {
			dataPartitioningScheme = atoi(argv[++i]);
			if((dataPartitioningScheme < 0) || (dataPartitioningScheme > PARTITIONING_TILED)){
				fprintf(stderr, "error: invalid data partitioning scheme.\n");
				exit(0);
			}
		} else if (strcmp(arg, "-cycle_length") == 0) {
			cycleLength = atof(argv[++i]);
			if(cycleLength <= 0) {
				fprintf(stderr, "error: invalid cycle_length.\n");
				exit(0);
			}
		}
#ifdef PARALLEL
		else if (strcmp(arg, "-t") == 0) {
			strcpy(numCilkWorkers,argv[++i]);
			if(atoi(numCilkWorkers) <= 0) {
				fprintf(stderr, "error: invalid number of cilk workers.\n");
				exit(0);
			}
			int nworkers=0;
			nworkers= __cilkrts_get_nworkers();
			if(procRank==0)
				printf("number of total Cilk worker threads available:%d\n", nworkers);
			if (0!= __cilkrts_set_param("nworkers","1")) {
			    printf("Failed to set worker count %s\n","1");
			    exit(0);
			}
		}
#endif

	}
/*#ifdef DEBUG
	printf("Recursiondepth:%d partitioningScheme:%d cycleLength:%f\n",recursionDepth,dataPartitioningScheme,cycleLength);
#endif*/
}


/* This function creates tasks, identifies inter-task dependencies, and partitions the tasks among different processes. 
 * The D2P inspector executes a process-specific code here. It is assumed that 'recursionDepth', the number of levels of recursion to unfold, is known before calling this function.
 * */ 
void GenerateTasks(int argc, char** argv)
{
 	int dpTableRowSize = ReadInput(argc, argv); //get the number of cells in a row of grid 
	struct timeval startTime, endTime;
	
	gettimeofday(&startTime,0);
	//A hierarchical decomposition of the data (DP table) creates tiles. First get the number of DP table cells per row in a tile of DP table.
	/* cellsPerRowTile is used to determine the tile ID given a box (or absolute coordinates in a grid). This computation is incorrect when there is uneven decomposition.
 	* hence, not used in most benchmarks. This code is still retained for comparing the inspector overhead with different approaches. */  	
	//dpTableRowSize = (number of tiles per row of the grid * number of cells per row of tile) 
	// = 2^recursion_depth * cellsPerRowTile
	//=>cellsPerRowTile = dpTableRowSize/2^recursionDepth
	Unroll(); //gets the list of function calls and the arguments with which they are called --- Populates fnCalls. fnCalls stores the list of method invocations.
	GetTaskList(); //creates a list of tasks from fnCalls.
	gettimeofday(&endTime,0);
	gPreprocessTime = (endTime.tv_sec-startTime.tv_sec)*1000000+(endTime.tv_usec-startTime.tv_usec);
}


void ExecuteTasks()
{
	bool done=false, readyToExit=false;
	int numProcsReadyToExit = 0;
	int exitMsg=0;
#ifdef PARALLEL
	__cilkrts_end_cilk();
	if (0!= __cilkrts_set_param("nworkers",numCilkWorkers)) {
            printf("Failed to set worker count %s\n",numCilkWorkers);
            return;
        }
#endif
#ifdef IDLE_TIME
	/*ostringstream logFileName;
	logFileName<<"proc"<<procRank<<"Log.txt";
	idleTimeLog = fopen(logFileName.str().c_str(),"w");*/
#endif

	struct timeval startTime, endTime;
	gettimeofday(&startTime,0);
	while(!done)
	{
		//Do background work.
		DoBackgroundWork();
		
		//If all local computation is done, notify all other processes and prepare for termination.
		if(numTasks == 0)
		{
			done = true;
		}

	}
	gettimeofday(&endTime,0);
	gTaskExecutionTime = (endTime.tv_sec-startTime.tv_sec)*1000000+(endTime.tv_usec-startTime.tv_usec);

	PrintResults();
#ifdef GRAPHVIZ_OUTPUT
	vector<double>::iterator largestIter = max_element(longestPathExecTimes.begin(),longestPathExecTimes.end());
	printf("T1:%lf CriticPathExecTime:%lf MaxSpeedup:%lf\n",T1, *largestIter, T1/(*largestIter));
#endif

#ifdef IDLE_TIME
	//fclose(idleTimeLog);
#endif
		
}


int DetermineUnfoldingLevel(int numProcesses)
{
	int unFoldingLevel = 0;
	int expTasksMax = numProcesses * numProcesses;
	int totalNodesAtCurrentLevel = 1;
	int numFnTypes = sizeof(fnCallHierarchySummary)/sizeof(fnCallHierarchySummary[0]);
	unsigned int* fnTypeCountsAtCurLevel = (unsigned int *)(calloc(numFnTypes,sizeof(int)));
	unsigned int* fnTypeCountsAtPrevLevel = (unsigned int *)(calloc(numFnTypes,sizeof(int)));
	//get the number of recursive function calls made at the top level (e.g. only function 'A' is called.
	fnTypeCountsAtCurLevel[0]=1; //says that function 'A' is called 1 time initially.
	unsigned int* numChildren = (unsigned int *)(calloc(numFnTypes,sizeof(int)));
	//get the total number of child nodes per recursive function call. A child node is a call to any recursive function.
	for(int i=0;i<numFnTypes;i++)
		for(int j=0;j<numFnTypes;j++)
			numChildren[i] += fnCallHierarchySummary[i][j]; 
	while(true)
	{
		if(totalNodesAtCurrentLevel >= expTasksMax)
			break;
		unFoldingLevel++;
		totalNodesAtCurrentLevel = 0;
		//total number of calls at current level = (number of calls to function 'X' * total number of children of 'X') summed over all 'X'
		for(int i=0;i<numFnTypes;i++)
			totalNodesAtCurrentLevel += (fnTypeCountsAtCurLevel[i]*numChildren[i]);
		memcpy(fnTypeCountsAtPrevLevel,fnTypeCountsAtCurLevel,sizeof(int)*numFnTypes);
		memset(fnTypeCountsAtCurLevel,0,sizeof(int)*numFnTypes);
		for(int i=0;i<numFnTypes;i++)
			for(int j=0;j<numFnTypes;j++)
				fnTypeCountsAtCurLevel[i] += fnTypeCountsAtPrevLevel[j] * fnCallHierarchySummary[j][i]; 
	}

	free(fnTypeCountsAtCurLevel);
	free(fnTypeCountsAtPrevLevel);
	free(numChildren);
	return unFoldingLevel;
}

#ifdef GRAPHVIZ_OUTPUT
vector<int> CheckIfAllPathsAreEnumerated(const vector<int>& pathLength)
{
	vector<int> critPath;
	vector<vector<int> >::iterator it = criticalPaths.begin();
	for(;it!=criticalPaths.end();it++)
	{
		if(pathLength[(*it)[0]] != 0)
		{
			critPath = *it;
			criticalPaths.erase(it);
			break;
		}
	}
	return critPath;
}	

void UpdateWorkSpanAnalysis(const int curFunctionID, const float execTime)
{
	T1 += execTime;
	
	vector<vector<int> >::iterator it = criticalPaths.begin();
	for(int i=0;i<criticalPaths.size();i++)
	{
		vector<int>::iterator longestPathIter = find((criticalPaths[i]).begin(), (criticalPaths[i]).end(), curFunctionID);
		if(longestPathIter != (criticalPaths[i]).end())
			longestPathExecTimes[i] += execTime;
	}
	return;
}

#endif


/* This function is a wrapper for creating tasks, data regions and initializing them */ 
void GetTaskList()
{
	/*struct timeval startTime, endTime;
	gettimeofday(&startTime,0);*/

#ifdef DEBUG
	int gNumLocalTiles, numLocalTiles=0, gNumLeaves, numLeaves = 0;
	map<int, vector<FunctionCall*> >::iterator fnIter = fnCalls.begin();
	while(fnIter != fnCalls.end())
	{
		numLeaves+=fnIter->second.size();
		if(fnIter->second.size() > 0)
			numLocalTiles++;
		fnIter++;
	}
	printf("%d: Number of tiles owned:%d\n",procRank,numLocalTiles);
			printf("Total number of leaves:%d tasks:%d\n",numLeaves,numLocalTiles);

	//for visualizing the data partitioning scheme
	/*int* rankTileMapping = new int[tilesPerRowDPTable*tilesPerRowDPTable];
	for(int i=0;i<tilesPerRowDPTable*tilesPerRowDPTable;i++)
		rankTileMapping[i]=INT_MAX;*/
#endif

	/*gettimeofday(&endTime,0);
	long gTaskDepDetTime, taskDepDetTime = (endTime.tv_sec-startTime.tv_sec)*1000000+(endTime.tv_usec-startTime.tv_usec);
	MPI_Allreduce(&taskDepDetTime,&gTaskDepDetTime, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD); 

	gettimeofday(&startTime,0);*/
#ifdef DEBUG
	/*map<int,vector<FunctionCall*> >::iterator fnIter = fnCalls.begin();
	while(fnIter != fnCalls.end())
	{

		Box* writeTile = (Box*)(((fnIter->second)[0])->params[0].data);
		int tileID = GetTileID2D(writeTile);
		tileID  = (writeTile->coords[1]/cellsPerRowTile)+(writeTile->coords[0]/cellsPerRowTile)*tilesPerRowDPTable;
		assert(tileID < tilesPerRowDPTable*tilesPerRowDPTable);
		assert(rankTileMapping[tileID] == INT_MAX);
		rankTileMapping[tileID] = procID;
		fnIter++;
	}*/
#endif
	/*gettimeofday(&endTime,0);
	long gPartitionTime, partitionTime = (endTime.tv_sec-startTime.tv_sec)*1000000+(endTime.tv_usec-startTime.tv_usec);
	MPI_Allreduce(&partitionTime,&gPartitionTime, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD); 
	gettimeofday(&startTime,0);*/

	if(computeGrid == COMPUTE_SPARSE)
		DistributeWork_DataParallel();
	
	CreateDataRegionsAndTasks();
	/*Task* it=tasks;
	while(it)
	{
		assert(it->functionCalls.size()==1);
		it->outPortOwners = it->functionCalls[0]->outPortOwners;
		it->numInPorts = it->functionCalls[0]->numInPorts;
		it= it->next;
	}*/
	/*gettimeofday(&endTime,0);
	long gDPTableInitTime, dpTableInitTime = (endTime.tv_sec-startTime.tv_sec)*1000000+(endTime.tv_usec-startTime.tv_usec);
	MPI_Allreduce(&dpTableInitTime,&gDPTableInitTime, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD); 
	if(procRank == 0)
	{
		printf("TaskDepDetermination time: %f (s) process-specific partition time: %f (s) DPTable Init time: %f (s)\n",gTaskDepDetTime/(float)1000000,gPartitionTime/(float)1000000,gDPTableInitTime/(float)1000000);
	}*/

#ifdef DEBUG
	/*if(procRank == 0)
	{
		for(int i=0;i<tilesPerRowDPTable;i++)
		{
			for(int j=0;j<tilesPerRowDPTable;j++)
			{
				printf("%d",rankTileMapping[i*tilesPerRowDPTable+j]);
				if(j!=tilesPerRowDPTable -1)
					printf(" ");		
			}
			printf("\n");
		}
	}
	delete [] rankTileMapping;*/
#endif

}

void DistributeWork_DataParallel()
{
	/*fnCalls contains the list of all tasks (leaves of the data dependency graph).
	Compute the starting index and the number of tasks that each process must process according to blocked distribution.*/
	int totalTasks=fnCalls.size();
	int avgNumTasksToProcess  = totalTasks / totalProcs;
	int* totalNumTasksToProcess = new int[totalProcs];
	int remainingTasks = totalTasks % totalProcs;
	int* numTasksToSkip = new int[totalProcs];
	for(int i=0;i<totalProcs;i++)
	{
		numTasksToSkip[i] = avgNumTasksToProcess * i; 
		totalNumTasksToProcess[i] = avgNumTasksToProcess; 
		if(i < remainingTasks)
			(totalNumTasksToProcess[i])++; 
		if(i > remainingTasks)
			numTasksToSkip[i] +=remainingTasks;
		else
			numTasksToSkip[i] += i;
	}
	
	/* Get process-specific pointers to the starting index and the ending index so that only the tasks assigned to the current process are stored in fnCalls. 
 	 * The pointers are used to delete the remaining tasks from fnCalls. */
	int endKey = -1;
	map<int,vector<FunctionCall*> >::iterator startIndxIter=fnCalls.end(), iter = fnCalls.begin();
	for(int i=0,j=0;iter!=fnCalls.end();iter++,i++)
	{

		if(i<numTasksToSkip[procRank])
			continue;

		j++;
		if((i!=0) && (j==1))
			startIndxIter = iter;

		if(j==totalNumTasksToProcess[procRank])
		{
			if(i != (fnCalls.size()-1))
			{
				iter++;
				endKey  = iter->first; 
			}
			break;
		}
	}
		
	/* Flatten the list of task IDs from fnCalls so that it is easier to compute the owner of a task given a task ID (used next for computing data dependencies /outPortOwner). */
	vector<set<int> > flatFnIds(totalProcs);
	iter = fnCalls.begin();
	for(int i=0,j=0;iter!=fnCalls.end();iter++,i++)
	{
		if(i>= (numTasksToSkip[j]+totalNumTasksToProcess[j]))
			j++;
		assert(j < totalProcs);
		for(int k=0;k<iter->second.size();k++)
		{
			FunctionCall* fnCall = iter->second[k];
			for(int l=1;l<fnCall->params.size();l++)
			{
				int readFnID = fnCall->params[l].portID;
				if(readFnID == -1)
					continue;
				flatFnIds[j].insert(readFnID);
			}
		}
	}
	delete [] numTasksToSkip;
	delete [] totalNumTasksToProcess;
	// delete the tasks that are not executed by the current process.
	if(startIndxIter != fnCalls.end())
	{
		fnCalls.erase(fnCalls.begin(),startIndxIter);
	}
	if(endKey != -1)
	{
		map<int,vector<FunctionCall*> >::iterator endIndxIter= fnCalls.find(endKey);
		assert(endIndxIter != fnCalls.end());
		fnCalls.erase(endIndxIter,fnCalls.end());
	}
	
#ifdef DEBUG
	int numLocalTiles = 0, numLeaves = 0;
	map<int,vector<FunctionCall*> >::iterator fnIter = fnCalls.begin();
	while(fnIter != fnCalls.end())
	{
		numLeaves+=fnIter->second.size();
		assert(fnIter->second.size() > 0);
		if(fnIter->second.size() > 0)
			numLocalTiles++;
		for(int k=0;k<fnIter->second.size();k++)
		{
			FunctionCall* fnCall = fnIter->second[k];
			int readParam1 = fnCall->params[1].portID;
			int readParam2 = fnCall->params[2].portID;
			printf("%d:%c(%d %d)\n",fnCall->ID,fnCall->functionName, readParam1, readParam2);
		}
		fnIter++;
	}
	
	printf("%d number of leaves:%d tasks:%d\n",procRank,numLeaves,numLocalTiles);
#endif

	//compute the process ID that reads the output of a task.
	iter = fnCalls.begin();
	for(;iter!=fnCalls.end();iter++)
	{
		for(int i=0;i<iter->second.size();i++)
		{
			FunctionCall* fnCall = (iter->second)[i];
			for(int k=0;k<flatFnIds.size();k++)
			{
				if(std::find(flatFnIds[k].begin(),flatFnIds[k].end(),fnCall->ID) != flatFnIds[k].end())
				{
					fnCall->outPortOwners.insert(k);
				}
			}
			if(i!=(iter->second.size()-1))
				fnCall->outPortOwners.insert(procRank);
		}
	}
}

/* This function first creates data regions (tiles) and initializes them. It then creates tasks and attaches the data regions and recursive method invocations to tasks. 
 * A work item is a recursive method invocation. Work items are created and distributed among processes before this function is called (during a call to Unroll). 
 * The data regions are created based on the per-process assignment of work items. A data region stores some metadata such as the bounding box of that region, and the ID of the task computing that region.
 * Tasks are then created. The work items and data regions created earlier are attached to tasks. This is done so that multiple work items can be attached to tasks for the purpose of coalescing 
 * those work items. 
 * Inputs: Task details (fnCalls): e.g. A(a,b,c) B(c,d,e). Here, it is assumed that the first parameter of the function represents the Tile being computed.
 * Output: tasks (global variable).
 */   
void CreateDataRegionsAndTasks()
{
	map<int,vector<FunctionCall*> >::iterator localIter = fnCalls.begin();
	while(localIter != fnCalls.end())
	{
		if(localIter->second.size() > 0)
		{
			Box* tmpRegion = (Box*)(localIter->second[0]->params[0].data);
			if((tmpRegion->GetBoxSize()+METADATASPACE) > INT_MAX)
			{
				printf("ERROR: memory requested too big to keep track of.\n");
				exit(0);
			}
			CELL_TYPE* tile = new CELL_TYPE[tmpRegion->GetBoxSize()+METADATASPACE];

			if(!tile)
			{
				printf("ERROR in initializing DP Table\n");
				exit(0);
			}
			//add metadata info
			for(int i=0;i<DIMENSION;i++)
			{
				tile[i]=tmpRegion->coords[i];
				tile[i+DIMENSION]=tmpRegion->coords[i+DIMENSION]-tmpRegion->coords[i]+1;
			}
			tile[2*DIMENSION]=localIter->first;
#ifdef DEBUG
			/*int tileID = GetZMortonTileIDFromCellCoords(tmpRegion->coords[0],tmpRegion->coords[1]);
			printf("Rank:%d (%d %d) ZMortonTileID:%d computed tileID:%d\n",procRank,tmpRegion->coords[0],tmpRegion->coords[1],localIter->first,tileID);*/
			//assert(tileID == localIter->first);
#endif
			InitializeDPTable(tmpRegion, tile+METADATASPACE);
			vector<FunctionCall*>& tileUpdatingFnList =localIter->second;
			vector<FunctionCall*>::iterator tileUpdatingFn = tileUpdatingFnList.begin();
			for(int i=0;i<tileUpdatingFnList.size();i++)
			{
				FunctionCall* tileUpdatingFn = tileUpdatingFnList[i];
				Task* t = new Task();
				t->functionCalls.push_back(tileUpdatingFn);
				t->updatedRegion = tile;
				
				int len = tmpRegion->GetBoxSize()+METADATASPACE;
				t->updateLen = len;
				//TODO: the below two fields of a task are duplicates. cleanup.
				t->outPortOwners = t->functionCalls[0]->outPortOwners;
				t->numInPorts = t->functionCalls[0]->numInPorts; 
	#ifdef TASK_AGGREGATION
				if(!(tileUpdatingFn)->isReadBeforeNextWrite)
				{
					if(i!=0)
						endTask->aggrSafeTask = t;
				}
	#endif
				if(endTask)
					endTask->next=t;
				else
					tasks=t;
				endTask=t;
				numTasks++;
			}
		}
		localIter++;
	}
}

void DoBackgroundWork()
{
	if(numTasks == 0)
		return;
	bool localWorkAvailable =true;

	while(localWorkAvailable)
	{
		localWorkAvailable = false;
		Task *prev=NULL, *it = tasks;
		while(it)
		{
			if(it->recvBuffers.size() == it->functionCalls[0]->numInPorts)
			{
				//printf("Executing task.\n");
				DoLocalComputation(*it);
				//printf("Finished executing task.\n");
				//printf("executed function %d\n",it->functionCalls[0]->ID);
				{
					bool tmpLocalWorkAvailable = UpdateTaskDependencies(it->updatedRegion, true, it->functionCalls[0]->ID);
					if(!localWorkAvailable)
						localWorkAvailable = tmpLocalWorkAvailable;
				}

#ifdef TASK_AGGREGATION
				int numCombinedTasks=0;
				Task* startOfChain = it, *curTaskExecuted = it;
				while(true)
				{
					Task* targetTask=curTaskExecuted->aggrSafeTask;
					if(targetTask && (targetTask->recvBuffers.size()==targetTask->functionCalls[0]->numInPorts))
					{
						numCombinedTasks++;
						DoLocalComputation(*targetTask);
						if(targetTask->outPortOwners.find(procRank) != targetTask->outPortOwners.end())
						{
							bool tmpLocalWorkAvailable = UpdateTaskDependencies(targetTask->updatedRegion, true, targetTask->functionCalls[0]->ID);
							if(!localWorkAvailable)
								localWorkAvailable = tmpLocalWorkAvailable;
						}
						curTaskExecuted = targetTask;
					}
					else
					{
						it = curTaskExecuted;
						break;
					}
				}
#ifdef DEBUG
				numMsgsSaved += (numCombinedTasks>0)?(numCombinedTasks-1):0;
#endif
#endif
#ifdef TASK_AGGREGATION
				it=startOfChain;
				while(numCombinedTasks>=0)
#endif	
				{
					assert(it);
					//deleting space reserved for function call arguments received from other processes.
					vector<pair<bool,CELL_TYPE*> >::iterator recvBufferIter = (it->recvBuffers).begin();
					for(;recvBufferIter!=(it->recvBuffers).end();recvBufferIter++)
					{
						if(!recvBufferIter->first) 
						{
							long long addr = reinterpret_cast<long long>(recvBufferIter->second);
							portReferenceCount[addr]--;
							if(!portReferenceCount[addr])
							{
								delete [] (recvBufferIter->second);
							}
						}
					}
			
					//deleting space reserved for function call arguments created during parsing.
					for(int j=0;j<it->functionCalls.size();j++)
					{
						FunctionCall* functionCall = (it->functionCalls)[j];
						for(int k=0;k<functionCall->params.size();k++)
						{
							Box* arg = (Box*)((functionCall->params[k]).data);
							delete arg;
						}
					}
					
					//deleting it;
					if(it==tasks)
					{
						if(tasks==endTask)
							endTask=NULL;
						tasks = it->next;	
					}
					else if(it==endTask)
					{
						assert(prev);
						prev->next=NULL;
						endTask = prev;
					}
					else
						prev->next=it->next;
					Task* next = it->next;

					dpTileTable[it->updatedRegion[2*DIMENSION]] = it->updatedRegion;
#ifdef TASK_AGGREGATION
					numCombinedTasks--;
					if(numCombinedTasks>=0)
						next = it->aggrSafeTask;
#endif
					
					delete it;
					it = next;
					numTasks--;
					//printf("%d numTasks:%d\n",procRank,numTasks);
				}
			}
			else
			{
				prev=it;
				it = it->next;
			}
		}
	}
	

	return;
}

void DoLocalComputation(Task& task)
{
	for(int i=0;i<task.functionCalls.size();i++)
	{
		FunctionCall* functionCall = (task.functionCalls)[i];
		Box* region = (Box*)(functionCall->params[0].data);
		vector<CELL_TYPE*> dataRegions;
		//printf("Loading data for %c(%d %d %d %d)\n",functionCall->functionName,region->coords[0],region->coords[1],region->coords[2],region->coords[3]);
		for(int j=0;j<functionCall->params.size();j++)
		{
			if(j==0)
				dataRegions.push_back(task.updatedRegion);
			else
			{
				int expectedDataSource = functionCall->params[j].portID;
				if(expectedDataSource == -1)
					dataRegions.push_back(task.updatedRegion);
				else
				{
					CELL_TYPE* data = (CELL_TYPE*)(functionCall->params[j].tile);
					//assert(data);
					dataRegions.push_back(data);
					if(data == NULL)
					{
					region = (Box*)(functionCall->params[j].data);
					printf("Function:%d ERROR: Data from fn (%d) not received.\n",functionCall->ID,expectedDataSource);
					exit(0);
					}
				}	
			}
		}
	
		/*if(task.aggrSafeTask)
			printf("executing function %d: %c(%d %d %d %d) hasAggrSafeTask:true\n",functionCall->ID, functionCall->functionName,dataRegions[0][0],dataRegions[0][1],dataRegions[0][0]+dataRegions[0][2]-1,dataRegions[0][1]+dataRegions[0][3]-1);
		else
			printf("executing function %d: %c(%d %d %d %d) hasAggrSafeTask:false\n",functionCall->ID, functionCall->functionName,dataRegions[0][0],dataRegions[0][1],dataRegions[0][0]+dataRegions[0][2]-1,dataRegions[0][1]+dataRegions[0][3]-1);*/
		/*if(functionCall->ID==8)
			printf("Debug break\n");*/
		ExecuteFunction(functionCall, dataRegions);
		task.updatedRegion[METADATASPACE-1]=functionCall->ID;
	}
	return;
}

/* Function to process incoming messages. 
 * Input: number of proceses from whom data is yet to be received.
 * Output: updated number of processes from whom data is yet to be received.
 * Comments:This function prioritizes handling incoming messages over doing background work (as indicated by the do-while loop). When there is no background work or waiting for producers to produce data,
 * a blocking call is executed. Otherwise, a non-blocking call is executed. Pending number of producers is updated as and when data is received from other processes (producers)*/    

bool UpdateTaskDependencies(CELL_TYPE* recvBuffer, bool localData, CELL_TYPE skipTask)
{
	int referenceCount = 0;
	bool localWorkAvailable = false;
	Task* it=tasks;
	/*vector<int> readFns;
	vector<int> newReadyList;*/
	while(it) 
	{
		if(it->functionCalls[0]->ID!=skipTask)
		{
			vector<FunctionCall*>& functionCalls = it->functionCalls; 	
			for(int j=0;j<functionCalls.size();j++)
			{
				FunctionCall* curFnCall = functionCalls[j];
				if(curFnCall->wawSource == recvBuffer[METADATASPACE-1])
				{
					curFnCall->numInPorts--;
					if(curFnCall->numInPorts == it->recvBuffers.size())
					{
						localWorkAvailable=true;
						//newReadyList.push_back(curFnCall->ID);
					}
					//readFns.push_back(curFnCall->ID);
					break;
				}	
				for(int k=1;k<curFnCall->params.size();k++)
				{
					int expectedDataSource = curFnCall->params[k].portID;
					if(expectedDataSource != -1)
					{
						if(recvBuffer[METADATASPACE-1]==expectedDataSource)
						{
							Box* region = (Box*)(curFnCall->params[0].data);
							/*if(curFnCall->ID==19)
							{
									printf("%d updating recvbuffer parameter %d from function %d\n",curFnCall->ID,k,recvBuffer[4]);
							}*/
							curFnCall->params[k].tile = recvBuffer;
							it->recvBuffers.push_back(make_pair(localData, recvBuffer));
							if(it->recvBuffers.size() == it->functionCalls[0]->numInPorts)
							{
							localWorkAvailable = true;
							//newReadyList.push_back(curFnCall->ID);
							}
							if(!localData)
								referenceCount++;
							//readFns.push_back(curFnCall->ID);
						}
					}
				}
			}
		}
		it = it->next;
	}

	/*if(recvBuffer[4]==7)
	{printf("number of read references to ");
	printf("%d",recvBuffer[4]);
	printf(":%d\n",readFns.size());


	printf("Functions reading: ");
	for(int i=0;i<readFns.size();i++)
		printf("%d ",readFns[i]);
	printf("\n"); 
	printf("Executed function: %d. ",recvBuffer[4]);
	printf("New ready list: ");
	for(int i=0;i<newReadyList.size();i++)
		printf("%d ",newReadyList[i]);
	printf("\n");
	}*/
	
	if(!localData)
	{
		long long addr = reinterpret_cast<long long>(recvBuffer);
		pair<map<long long, int>::iterator,bool> ret = portReferenceCount.insert(make_pair(addr,referenceCount));
		if(!ret.second)
		{	
			if(portReferenceCount[addr] != 0)
				printf("Warning: attempting to delete a DP table tile that is being referenced\n");
			portReferenceCount[addr]=referenceCount;
		}
	}
	return localWorkAvailable;
}



int PrintGlobalMax()
{
	struct
	{
		CELL_TYPE value;
		int rank;
	}localMax, gMax;

	if(typeid(int) == typeid(CELL_TYPE)){
		localMax.value = -INT_MAX;
	}
	else{
		localMax.value =-FLT_MAX;
	}
	localMax.rank = 0;

	CELL_TYPE* tileContainingMax=NULL;
	map<int,CELL_TYPE*>::iterator it = dpTileTable.begin();
	while(it!=dpTileTable.end())
	{
		CELL_TYPE* updatedRegion = it->second;
		int boxLen = 1;
		bool rankAndTileNotUpdated = true;
		
		for(int i=0;i<DIMENSION;i++)
			boxLen *= updatedRegion[i+DIMENSION];
		
		for(int i=METADATASPACE;i<METADATASPACE+boxLen;i++)
		{
			if(updatedRegion[i] > localMax.value)
			{
				localMax.value = updatedRegion[i];
				if(rankAndTileNotUpdated)
				{
					tileContainingMax = updatedRegion;
					rankAndTileNotUpdated = false;
					localMax.rank = procRank;
				}
			}
		}
		it++;
	}

		cout<<"Max value:"<<localMax.value<<endl;
	

	return gMax.rank;
}


bool IsLocalOwner(int tileID)
{
	return true;
}

int GetOwner(int tileID)
{
	return procRank;
}

int GetTileID1D(int zmortTileID)
{
	int tileID = zmortTileID;
	if ((dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_HORIZONTAL) || (dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_VERTICAL))
	{
		int oldTileID = tileID;
		tileID = ((tileID % strideOffset)/strideLength)*avgNumTilesToProcess + tileID % strideLength + (tileID/strideOffset)*strideLength;
	}
	return tileID;
}

unsigned int DecodeMortonHelper(unsigned int x)
{
  x &= 0x55555555;                  
  x = (x ^ (x >>  1)) & 0x33333333; 
  x = (x ^ (x >>  2)) & 0x0f0f0f0f; 
  x = (x ^ (x >>  4)) & 0x00ff00ff; 
  x = (x ^ (x >>  8)) & 0x0000ffff; 
  return x;
}

unsigned int DecodeMortonX(unsigned int code)
{
  return DecodeMortonHelper(code >> 0);
}

unsigned int DecodeMortonY(unsigned int code)
{
  return DecodeMortonHelper(code >> 1);
}

int GetTileID2D(int zmortTileID)
{
	int xCoord= DecodeMortonX(zmortTileID);
	int yCoord = DecodeMortonY(zmortTileID);
	int tileID;
	if(dataPartitioningScheme == PARTITIONING_BLOCKED_HORIZONTAL)
	{
		tileID = yCoord*tilesPerRowDPTable + xCoord; //rowmajor numbering of tiles.
		if(computeGrid == COMPUTE_UTM)
			tileID -= (yCoord * (yCoord+1)/2);
	}
	else if(dataPartitioningScheme == PARTITIONING_BLOCKED_VERTICAL)
	{
		tileID = xCoord*tilesPerRowDPTable + yCoord; //colmajor numbering of tiles.
		if(computeGrid == COMPUTE_UTM)
			tileID -= (xCoord * (tilesPerRowDPTable-(xCoord-1)-1+tilesPerRowDPTable-1))/2;
	}
	/*else if (dataPartitioningScheme == PARTITIONING_TILED)
	{
		//printf("tilesPerRowBlock:%d tilesPerColumnBlock:%d\n",1<<tilesPerRowBlock, 1<<tilesPerColumnBlock);
		int blockIDX = tile->coords[1]/(cellsPerRowTile *(tilesPerRowBlock));
		int blockIDY = tile->coords[0]/(cellsPerRowTile *(tilesPerColumnBlock));
		int numBlocksPerRow = tilesPerRowDPTable/(tilesPerRowBlock);
		int blockID = blockIDY*numBlocksPerRow+blockIDX;
		//printf("X:%d Y:%d blockID:%d\n",blockIDX,blockIDY, blockID);
		int firstTileOfBlock = blockID*tilesPerRowBlock*tilesPerColumnBlock;
		//assert((firstTileOfBlock == 0)||(firstTileOfBlock == 512));
		int xOffset= (blockIDX>0)?(tile->coords[1])/cellsPerRowTile-(blockIDX*tilesPerRowBlock):(tile->coords[1])/cellsPerRowTile;
		int yOffset = (blockIDY>0)?(tile->coords[0])/cellsPerRowTile - (blockIDY*tilesPerColumnBlock):(tile->coords[0])/cellsPerRowTile;
		//printf("xOffset:%d blockIDX:%d blockIDY:%d firstTileID:%d coords[1]:%d coords[0]:%d\n",xOffset, blockIDX, blockIDY, firstTileOfBlock,tile->coords[1],tile->coords[0]);
		tileID  = firstTileOfBlock + yOffset*tilesPerRowBlock+xOffset;
		assert(tileID < (tilesPerRowDPTable*tilesPerRowDPTable));
	}*/
	else if ((dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_HORIZONTAL) || (dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_VERTICAL))
	{
		if(computeGrid == COMPUTE_UTM)
		{
			assert(xCoord >= yCoord);
			int k = xCoord-yCoord;
			int i = 0, j=0;
			while(k > 0)
			{
				i+= tilesPerRowDPTable - j;
				j++;
				k--; 
			}
			i+= yCoord;
			tileID = i;
		}
		else
		{
			if((dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_HORIZONTAL))
				tileID = yCoord*tilesPerRowDPTable + xCoord; //rowmajor numbering of tiles.
			else
				tileID = xCoord*tilesPerRowDPTable + yCoord; //colmajor numbering of tiles.
			int oldTileID = tileID;
			tileID = ((tileID % strideOffset)/strideLength)*avgNumTilesToProcess + tileID % strideLength + (tileID/strideOffset)*strideLength;
			//printf("x: %d y: %d tileID: %d newTileID: %d\n",xCoord, yCoord, oldTileID, tileID); 
		}
	}
#ifdef DEBUG
	//printf("x: %d y: %d tileID: %d \n",tile->coords[0], tile->coords[1], tileID); 
#endif
	return tileID;
}

unsigned int DecodeMortonHelper3D(unsigned long x)
{
  x &= 0x0000249249249249;                  
  x = (x ^ (x >>  2)) & 0x00000c30c30c30c3; 
  x = (x ^ (x >>  4)) & 0x000000f00f00f00f; 
  x = (x ^ (x >>  8)) & 0x00000000ff0000ff; 
  x = (x ^ (x >>  16)) & 0x000000000000ffff; 
  return x;
}

unsigned int DecodeMortonX3D(unsigned long code)
{
  return DecodeMortonHelper3D(code >> 0);
}

unsigned int DecodeMortonY3D(unsigned long code)
{
  return DecodeMortonHelper3D(code >> 1);
}

unsigned int DecodeMortonZ3D(unsigned long code)
{
  return DecodeMortonHelper3D(code >> 2);
}

int GetTileID3D(unsigned long zmortTileID)
{
	int xCoord= DecodeMortonX3D(zmortTileID);
	int yCoord = DecodeMortonY3D(zmortTileID);
	int zCoord = DecodeMortonZ3D(zmortTileID);
	int tileID=-1;
	if(dataPartitioningScheme == PARTITIONING_BLOCKED_HORIZONTAL)
	{
		tileID =zCoord*tilesPerRowDPTable*tilesPerRowDPTable+yCoord*tilesPerRowDPTable + xCoord; //rowmajor numbering of tiles.
	}
	else if(dataPartitioningScheme == PARTITIONING_BLOCKED_VERTICAL)
	{
		tileID = xCoord*tilesPerRowDPTable*tilesPerRowDPTable + yCoord*tilesPerRowDPTable + zCoord; //colmajor numbering of tiles.
	}
	else if ((dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_HORIZONTAL) || (dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_VERTICAL))
	{
		if((dataPartitioningScheme == PARTITIONING_BLOCKCYCLIC_HORIZONTAL))
			tileID = zCoord*tilesPerRowDPTable*tilesPerRowDPTable+yCoord*tilesPerRowDPTable + xCoord; //rowmajor numbering of tiles.
		else
			tileID = zCoord*tilesPerRowDPTable*tilesPerRowDPTable+xCoord*tilesPerRowDPTable + yCoord; //colmajor numbering of tiles.
		int oldTileID = tileID;
		tileID = ((tileID % strideOffset)/strideLength)*avgNumTilesToProcess + tileID % strideLength + (tileID/strideOffset)*strideLength;
		//if(procRank ==0)printf("x: %d y: %d z:%d zMortTileID:%d oldTileID: %d newTileID: %d\n",xCoord, yCoord, zCoord, zmortTileID, oldTileID, tileID); 
	}
#ifdef DEBUG
	//printf("x: %d y: %d z: %d tileID: %d \n",tile->coords[0], tile->coords[1], tile->coords[2], tileID); 
#endif
	assert(tileID >= 0);
	return tileID;
}

int GetTileID(const Box* tile)
{
	int xCoord= (tile->coords[0])/cellsPerRowTile;
	int yCoord = (tile->coords[1])/cellsPerRowTile;
	int uCoord= (tile->coords[2])/cellsPerRowTile;
	int vCoord = (tile->coords[3])/cellsPerRowTile;
	int tileID;
	long int tmp = vCoord*tilesPerRowDPTable*tilesPerRowDPTable*tilesPerRowDPTable+uCoord*tilesPerRowDPTable*tilesPerRowDPTable+yCoord*tilesPerRowDPTable + xCoord; //rowmajor numbering of tiles.
	assert(tmp <= INT_MAX);
	tileID = tmp;
	return tileID;
}

CELL_TYPE* GetTileOfDPTable(int i, int j)
{
	CELL_TYPE* ret = NULL;
	map<int, CELL_TYPE*>::iterator it = dpTileTable.begin();
	while(it!=dpTileTable.end())
	{
		CELL_TYPE* curTile = it->second;
		if((i<=curTile[0]+curTile[2]-1) && (i>=curTile[0]) && (j<=curTile[1]+curTile[3]-1) && (j>=curTile[1]))
		{
			ret = curTile;			
			break;
		} 
		it++;
	}
	return ret;
}

#ifdef DEBUG
unsigned int EncodeMortonHelper(unsigned int x)
{
  x &= 0x0000ffff;                  
  x = (x ^ (x <<  8)) & 0x00ff00ff; 
  x = (x ^ (x <<  4)) & 0x0f0f0f0f; 
  x = (x ^ (x <<  2)) & 0x33333333; 
  x = (x ^ (x <<  1)) & 0x55555555; 
  return x;
}

unsigned int EncodeMorton(unsigned int x, unsigned int y)
{
        return (EncodeMortonHelper(y) << 1) + EncodeMortonHelper(x);
}

//reversing from (i,j) - (j,i) is intentional. The DP table is read and written in D2P as (row, col). In ZMorton order, (i,j) refers to ith column and jth row. 
int GetZMortonTileIDFromCellCoords(int i, int j)
{
	int x=i/cellsPerRowTile;
	int y=j/cellsPerRowTile;
	return EncodeMorton(x,y);
}

CELL_TYPE* GetTileOfDPTable(int tileID)
{
	CELL_TYPE* ret = NULL;
	if(dpTileTable.find(tileID)!=dpTileTable.end())
		ret = dpTileTable[tileID];
	return ret;
}

void DebugPrintTileDetails()
{
	map<int, CELL_TYPE*>::iterator it = dpTileTable.begin();
	while(it!=dpTileTable.end())
	{
		CELL_TYPE* curTile = it->second;
		ostringstream ostrm;
		ostrm<<procRank<<":";
		for(int i=0;i<DIMENSION;i++)
		{
			ostrm<<curTile[i];
			if(i!=DIMENSION-1)
				ostrm<<" ";
		}
		printf("%s\n",ostrm.str().c_str());
		it++;
	}
}

CELL_TYPE* GetLocalTile(int i)
{
	CELL_TYPE* ret = NULL;
	map<int,CELL_TYPE*>::iterator it = dpTileTable.begin();
	while(it!=dpTileTable.end())
	{
		CELL_TYPE* updatedRegion = it->second;
		if((i>=updatedRegion[0]) && (i<=updatedRegion[0]+updatedRegion[1]-1))
		{
			assert(ret == NULL);
			ret = updatedRegion;
			break;
		}
		it++;
	}
	return ret;
}


CELL_TYPE* GetLocalTile(int i, int j)
{
	CELL_TYPE* ret = NULL;
	map<int,CELL_TYPE*>::iterator it = dpTileTable.begin();
	while(it!=dpTileTable.end())
	{
		CELL_TYPE* updatedRegion = it->second;
		if((i>=updatedRegion[0]) && (j>=updatedRegion[1]) && (i<=updatedRegion[0]+updatedRegion[2]-1) && (j<=updatedRegion[1]+updatedRegion[3]-1))
		{
			assert(ret == NULL);
			ret = updatedRegion;
			break;
		}
		it++;
	}
	return ret;
}

CELL_TYPE* GetLocalTile(int i, int j, int k)
{
	CELL_TYPE* ret = NULL;
	map<int,CELL_TYPE*>::iterator it = dpTileTable.begin();
	while(it!=dpTileTable.end())
	{
		CELL_TYPE* updatedRegion = it->second;
		if((i>=updatedRegion[0]) && (j>=updatedRegion[1]) && (k>=updatedRegion[2]) && (i<=updatedRegion[0]+updatedRegion[3]-1) && (j<=updatedRegion[1]+updatedRegion[4]-1) && (k<=updatedRegion[2]+updatedRegion[5]-1))
		{
			assert(ret == NULL);
			ret = updatedRegion;
			break;
		}
		it++;
	}
	return ret;

}



#endif


