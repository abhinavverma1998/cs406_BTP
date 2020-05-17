#ifdef MPI
#include<mpi.h>
#endif
#include<sys/time.h>
#include "HelperFunctions.h"

int procRank;
int totalProcs=1;
//int argcCopy;
//char** argvCopy;
long gPreprocessTime;
long gTaskExecutionTime;
//void MakeCLACopy(int argc, char** argv);


#ifdef IDLE_TIME
#include<limits>
#include<math.h>
#define MAX_NUM_NODES 8
#define CORES_PER_NODE 16
extern double idleTime;
extern int numMsgsRecvd;
double GetGeoMean(const double* data, int numElems);
double GetVariance(const double* data, int numElems);
double GetStddev(const double* data, int numElems);
#endif

int lNodewiseMsgLog[64], gNodewiseMsgLog[64];
#ifdef PAPI
#include"papi.h"
int papiStat, eventSet=PAPI_NULL;	
long long values[2] = {(long long)0, (long long)0};
inline void handle_error(int retval)
{
	printf("PAPI error %d: %s\n",retval, PAPI_strerror(retval));
	exit(1);
}
#endif


int main(int argc, char* argv[])
{

	// DC
	if((argc==1) || strcmp(argv[1],"-h")==0)
	{
		printf("./<exe> <app-specific-arguments-parsed-by-ReadInput> <system-parameters>\n");
		printf("*****system-parameters****\n");
		printf("-t <>\t -- number of Cilk workers to use in parallel execution.-DPARALLEL must be specified during compilation.\n");
		printf("-recursion_depth <>\t -- unfolds the recursion a specified number of levels to generate tasks.\n");
		printf("-partition <>\t -- sets the data partitioning scheme: ");
		printf("PARTITIONING_BLOCKED_VERTICAL - 0, PARTITIONING_BLOCKCYCLIC_HORIZONTAL - 1, PARTITIONING_BLOCKED_VERTICAL - 2, PARTITIONING_BLOCKCYCLIC_VERTICAL - 3, PARTITIONING_TILED - 4\n");
		printf("-cycle_length <>\t -- sets the cycle length in case of BLOCKCYCLIC schemes. ignored in all other cases.\n");
		printf("*****App-specific parameters***\n");
	}

	//Get parameters.
	ParseSystemParams(argc, argv);
	//Create tasks
	GenerateTasks(argc, argv);

	/*for(int i=0;i<64;i++)
		lNodewiseMsgLog[i]=0;*/

	//execute tasks
	ExecuteTasks();

	if(procRank == 0)
	{
		printf("Preprocess time: %f seconds\n",gPreprocessTime/(float)1000000);
		printf("Task completion time: %f seconds\n",gTaskExecutionTime/(float)1000000);
	}
// #ifdef IDLE_TIME
// 	double avgTimeForDataReception = idleTime/numMsgsRecvd;
// 	if(numMsgsRecvd == 0)
// 		avgTimeForDataReception = 0;
// 	printf("%d: NumMsgsRecvd:%d avg time elapsed before receiving:%f\n",procRank,numMsgsRecvd,avgTimeForDataReception);
	
// 	double geomeanAvg, gAvgTimeForDataReception;
// 	/*double* arrAvgTimeForDataReception=nullptr;
// 	if(procRank == 0)
// 		arrAvgTimeForDataReception= new double[totalProcs];
// 	MPI_Gather(&avgTimeForDataReception,1, MPI_DOUBLE,arrAvgTimeForDataReception,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
// 	if(procRank == 0)
// 	{
// 		//printf("%f %f %f %f\n",arrAvgTimeForDataReception[0],arrAvgTimeForDataReception[1],arrAvgTimeForDataReception[2],arrAvgTimeForDataReception[3]);
// 		geomeanAvg = GetGeoMean(arrAvgTimeForDataReception,totalProcs);
// 	}
// 	MPI_Reduce(&avgTimeForDataReception,&gAvgTimeForDataReception, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
// 	if(procRank == 0)
// 		printf("Avg time elapsed between two messages (all processes):%f geomean: %f\n",gAvgTimeForDataReception/totalProcs, geomeanAvg); */

// 	if(procRank == 0)
// 		printf("Avg time elapsed between two messages (all processes):%f\n",gAvgTimeForDataReception/totalProcs); 
// 	double overhead = idleTime/(gTaskExecutionTime/(float)1000000);
// 	//printf("%d: idleTime:%f overhead:%f\n",procRank,idleTime,overhead);
// 	double perNodeOverhead[MAX_NUM_NODES], gPerNodeOverhead[MAX_NUM_NODES];
// 	int numProcsMappedToNodes[MAX_NUM_NODES], gNumProcsMappedToNodes[MAX_NUM_NODES];
// 	for(int i=0;i<MAX_NUM_NODES;i++)
// 	{
// 		perNodeOverhead[i]=0;
// 		numProcsMappedToNodes[i]=0;
// 	}
	
// 	numProcsMappedToNodes[(procRank/CORES_PER_NODE) % MAX_NUM_NODES]=1;

// 	int numNodesUsed = (totalProcs/CORES_PER_NODE >= MAX_NUM_NODES) ? MAX_NUM_NODES: ceil(totalProcs/(float)(CORES_PER_NODE));
// 	if(procRank == 0)
// 	{
// 		for(int i=0;i<MAX_NUM_NODES;i++)
// 		{
// 			if(gNumProcsMappedToNodes[i] != 0)
// 				gPerNodeOverhead[i] /= gNumProcsMappedToNodes[i];
// 			else
// 				gPerNodeOverhead[i] = 0;
// 		}

// 		for(int i=0;i<MAX_NUM_NODES;i++)
// 		{
// 			if(gPerNodeOverhead[i] != 0)
// 				printf("Wabash%d: NumProcsMapped: %d perNodeOverhead (avg): %f\n",i,gNumProcsMappedToNodes[i], gPerNodeOverhead[i]);
// 		}

// 		double mean=0;
// 		for(int i=0;i<numNodesUsed;i++)
// 			mean += gPerNodeOverhead[i];
// 		mean /= numNodesUsed;
// 		double var = GetVariance(gPerNodeOverhead, numNodesUsed);
// 		double stddev = sqrt(var);
// 		double cov = stddev / mean;
// 		printf("Number of nodes utilized:%d idle_time: Mean=%f stddev=%f CV(stddev/Mean)=%f\n",numNodesUsed,mean,stddev,cov);
// 	}
// #endif
	
	//PrintCellCost(procRank);
	return 0;
}

/*void MakeCLACopy(int argc, char** argv)
{
	argcCopy = argc;
	argvCopy = new char*[argcCopy];
	for(unsigned int i=0;i<argcCopy;i++)
	{
		int argLen = strlen(argv[i]);
		argvCopy[i] = new char[argLen+1];
		strcpy(argvCopy[i],argv[i]);
	}
}*/

// Utility functions 
#ifdef IDLE_TIME
double GetGeoMean(const double* data, int numElems)
{
    double partialSignificandProd = 1.0;
    long long exp = 0;
    double invN = 1.0 / numElems;

    for (int i=0;i<numElems; i++)
    {
        int curExp;
        double curSignificand = frexp(data[i],&curExp);
        partialSignificandProd *= curSignificand;
        exp += curExp;
    }

    return pow( std::numeric_limits<double>::radix,exp * invN) * pow(partialSignificandProd,invN);

}

double GetVariance(const double* data, int numElems)
{
    // Compute mean (average of elements)
    double sum = 0;
    for (int i = 0; i < numElems; i++)
        sum += data[i];

    double mean = sum / numElems;
 
    // Compute sum squared 
    // differences with mean.
    double sqDiff = 0;
    for (int i = 0; i < numElems; i++) 
        sqDiff += (data[i] - mean) * 
                  (data[i] - mean);
    return sqDiff / numElems;
}

double GetStddev(const double* data, int numElems)
{
	double var = GetVariance(data, numElems);
	return sqrt(var);
}

#endif
