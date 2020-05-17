#include<stdio.h>
#include<vector>
#include<set>
#include<iostream>
#include<cstring>
#include<cstdlib>
#include<fstream>
#include<algorithm>
#include<climits>
#include<ctime>
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


typedef struct Box
{   int coords[4];
    Box(){}
    Box(int a, int b, int c, int d)
    {
        coords[0]=a;
        coords[1]=b;
        coords[2]=c;
        coords[3]=d;
    }
    Box(const Box&b)
    {
        coords[0]=b.coords[0];coords[1]=b.coords[1];coords[2]=b.coords[2];coords[3]=b.coords[3];
    }
    bool operator==(const Box& rhs)
    {
        bool flag = true;
        if((this->coords[0]!=rhs.coords[0]) || (this->coords[1]!=rhs.coords[1]) || (this->coords[2]!=rhs.coords[2]) || (this->coords[3]!=rhs.coords[3]))
            flag=false;
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
    char* PrintStr(char *str)const
    {
        sprintf(str,"%d%d%d%d",coords[0],coords[1],coords[2],coords[3]);
    }
}Box;

// typedef struct Box
// {
//     int  coords[DIMENSION*2]; //topleft and bottom right (x and y coordinates)
//     Box(int coord0...)
//     {
//         va_list argList;
//         coords[0]=coord0;
//         va_start(argList,coord0);
//         for(int argNum=1;argNum<2*DIMENSION;argNum++)
//             coords[argNum]=va_arg(argList,int);
//         va_end(argList);
//     }
//     Box(){}
//     void PrintStr(char *str)const
//     {
//         std::ostringstream ostrstrm;
//         for(int i=0;i<2*DIMENSION;i++)
//         {
//             ostrstrm<<coords[i];
//             if(i!= (2*DIMENSION-1))
//                 ostrstrm<<std::string(" ");
//         }
//         //cout<<ostrstrm.str();
//         strcpy(str,ostrstrm.str().c_str());
//         //sprintf(str,"%d %d %d %d",coords[0],coords[1],coords[2],coords[3],coords[4],coords[5],coords[6],coords[7]);
//     }
//     Box& operator=(const Box& rhs)
//     {
//         for(int i = 0;i<DIMENSION*2;i++)
//             this->coords[i]=rhs.coords[i];
//         return *this;
//     }
//     Box(const Box& rhs)
//     {
//         for(int i = 0;i<DIMENSION*2;i++)
//             this->coords[i]=rhs.coords[i];
//     }
//     bool operator==(const Box& rhs)
//     {
//         bool flag = true;
//         for(int i = 0;i<DIMENSION*2;i++)
//         {
//             if(this->coords[i]!=rhs.coords[i])
//             {
//                 flag = false;
//                 break;
//             }
//         }
//         return flag;
//     }
//     long int GetBoxSize()
//     {
//         long int len = 1;
//         for(int i=0;i<DIMENSION;i++)
//         {
//             len *= (coords[i+DIMENSION]-coords[i]+1);
//         }
//         return len;
//     }
// }Box;


using namespace std;

int getNearestDepth(Box* b) {
    // Simple recursive function
    if (b->coords[3] <= 0) {
        return 0;
    }
    Box *b_new = new Box(0, 0, b->coords[2]/2, b->coords[3]/2);

    return (1 + getNearestDepth(b_new));
}

int main() {
    cout<<"Hello, entering this program\n";
    Box* b = new Box(0, 0, 256, 256);
    int depth = getNearestDepth(b);
	return 0;
}
