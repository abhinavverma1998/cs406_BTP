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
#include"GlobalTypes.h"

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
    Box* b = new Box(0, 0, 10240000, 10240000);
    int depth = getNearestDepth(b);
	return 0;
}
