/* The code in this file implements the D2P CodeGen. Start with the 'main' function. 
 * The comments in this file are a tribute to Joe Armstrong. "The Mess We're In", https://www.youtube.com/watch?v=lKXe3HUG2l4 */ 
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<sstream>
#include<string.h>
#include<cassert>
using namespace std;

int gridDimension=-1; //dimension of the grid computed.
//type definition and enumeration constants that indicate whether the grid is computed fuly or partially (triangular matrix, sparse).
enum ComputeGrid{
COMPUTE_FULL,
COMPUTE_UTM,
COMPUTE_SPARSE
};

int computeGrid = COMPUTE_FULL;
struct Invocation
{
	string name; //name of the function being called.
	vector<string> arguments; //arguments passed.
	void ExtractName(const string& line, int startPos);
};

struct RecMethod
{
	string name; //stores the name of the recursive method.
	vector<string> parameters; //stores method parameters.
	/*bodyStr initially contains a list of strings representing method body. populated after a method definition begins and before another method definition begins. 
 	* Call to GenerateFunctionBody clears the list and creates a single entry containing a C++ code for the method body. */  
	vector<string> bodyStr; 
	string signature, signature_unroll; //method signatures for inspector and executor versions of the method.
	vector<vector<Invocation> > calls; //method invocations per line of the recursive method body.
	void PopulateCalls();
	void GenerateFunctionBody();
	string GenerateStopCondWhileUnrolling(const vector<string>& fnArgs);
	string GenerateFunctionBodyForUnroll();
	string GenerateStopCond(const string& writeParam);
	string GenerateTopLevelFunctionInvocation();
	string GenerateParamListMacro() const;
};

/* This function basically splits an input string (tup) into a list of tokens based on the delimiter passed. 
 * Assumption: the input string is a tuple in the format '<' and '>'. There can be multiple tuples in the input as well.*/
void ProcessTuple(const string& tup, vector<vector<string> >& tokens, char delimiter)
{
	int pos = tup.find_first_of('<');	
	int endpos = tup.find_first_of('>');	

	assert(endpos > pos);
	while(endpos > pos)
	{
		vector<string> tokensOfaTuple;
		string tupCleaned = tup.substr(pos+1, endpos-pos-1);
		int startPosTok=0, endPosTok = tupCleaned.find_first_of(delimiter);	
		while(endPosTok != string::npos)
		{
			string tok = tupCleaned.substr(startPosTok, (endPosTok-startPosTok));
			tokensOfaTuple.push_back(tok);		
			startPosTok = endPosTok+1;
			endPosTok = tupCleaned.find_first_of(delimiter,startPosTok);	
		}
		string tok = tupCleaned.substr(startPosTok);
		tokensOfaTuple.push_back(tok);		
		pos = tup.find_first_of('<',pos+1);	
		endpos = tup.find_first_of('>',endpos+1);	
		tokens.push_back(tokensOfaTuple);
	}
	
}  

/* This function parses the input spec and returns a list of calls made within a recursive method body.
 * It takes as input a line of input spec (string) and returns a list of calls (Invocation objects) made in that line.
 */
void GetCalls(const string& line, vector<Invocation>& calls)
{
	//skip blank lines.
	int blnkLine = line.find_first_not_of(' ');	
	if(blnkLine == string::npos)
		return;
	int anchorPos = 0;
	int pos = line.find_first_of('(');	
	int endpos = line.find_first_of(')');	

	//Function calls are assumed to be in () format and hence, scan a line for '(' and ')'. There may be multiple calls in a line (indicating parallel invocations). 
	assert(endpos > pos);
	while((endpos != string::npos) && (pos != string::npos))
	{
		assert(pos > 0);
		//extract the name of the invocation/ function call
		Invocation inv;
		inv.ExtractName(line, anchorPos);
		string tupCleaned = line.substr(pos+1, endpos-pos-1);
		vector<vector<string> > arguments;
		map<string,int> uniqArgNames;
		int order=0;
		/*extract arguments. The spec can have combined data regions as arguments (<x,x,u><x v x>) to function calls. In such cases extract unique arguments from combined data regions.*/ 
		ProcessTuple(tupCleaned, arguments, ' '); 
		if(arguments.size() > 1)
		{
			//combined region tuples. extract unique arguments.
			for(int i=0;i<arguments.size();i++)
			{
				for(int j=0;j<arguments[i].size();j++)
				{
					if(uniqArgNames.find(arguments[i][j]) == uniqArgNames.end())
						uniqArgNames[arguments[i][j]] = order++;
				}
			}
			
			inv.arguments.resize(uniqArgNames.size());
			map<string,int>::iterator iter = uniqArgNames.begin();
			for(;iter!=uniqArgNames.end();iter++)
				(inv.arguments)[iter->second]=iter->first;
		}
		else//not a combined data region. Hence pass the arguments as is (i.e. set the arguments field of invocation structure).
			inv.arguments = arguments[0];
		//check if any other call exists on the same line.
		anchorPos = endpos+1;
		pos = line.find_first_of('(',pos+1);	
		endpos = line.find_first_of(')',endpos+1);	
		calls.push_back(inv);
	}
}


/* This function parses the raw string in the input spec and sets the 'name' field of a RecMethod structure. 
 * Assumption: a recursive method name is in this format: X_arbitrarystring(), where X is a single character that names the method. */
void GetMethodNameAndParams(const string& line, RecMethod* method)
{
	int pos = line.find("_");
	method->name = line.substr(0,pos);
	
	map<string,int> uniqParamNames;
	vector<vector<string> > parameters;
	ProcessTuple(line, parameters, ','); 
	int order=0;
	//replace duplicate parameters with uniq parameters.
	if(parameters.size() > 1)
	{
		//combined region tuples.
		for(int i=0;i<parameters.size();i++)
		{
			for(int j=0;j<parameters[i].size();j++)
			{
				if(uniqParamNames.find(parameters[i][j]) == uniqParamNames.end())
					uniqParamNames[parameters[i][j]] = order++;
			}
		}
	}
	else
	{
		uniqParamNames[parameters[0][0]] = order++;
		for(int i=1;i<parameters[0].size();i++)
		{
			assert(parameters[0][i].size() == 1);
			if(uniqParamNames.find(parameters[0][i]) == uniqParamNames.end())
				uniqParamNames[parameters[0][i]] = order++;
			else
			{
				char c = parameters[0][0][0]+i;	
				if(c > parameters[0][0][0]+25)
					c = parameters[0][0][0]+(c-(parameters[0][0][0]+25)-1);

				assert(uniqParamNames.find(string(1,c)) == uniqParamNames.end());
				uniqParamNames[string(1,c)]=order++;
			}
		}
	}
	method->parameters.resize(uniqParamNames.size());
	map<string,int>::iterator iter = uniqParamNames.begin();
	for(;iter!=uniqParamNames.end();iter++)
		(method->parameters)[iter->second]=iter->first;
}

/* This function generates the serial code from a parallel code consisting of cilk_spawn and cilk_sync keywords.
 * The generated serial code is exactly the same as the parallel code except for the absence of Cilk keywords.
 * So, this function basically does a search-and-delete of the Cilk keywords. */ 
string GenerateSerialCodeInFunctionBody(string fnBody)
{
	string serialCode("");
	std::string queryStr1("#ifdef PARALLEL\n\tcilk_spawn\n#endif\n");
	std::string queryStr2("#ifdef PARALLEL\n\tcilk_sync;\n#endif\n");
	int queryStr1Len = queryStr1.size();
	int queryStr2Len = queryStr2.size();
	bool parallelTasksFound=false;
	while(true)
	{
		size_t res1 = fnBody.find(queryStr1,0);
		if(res1==std::string::npos)
			break;
		else
			fnBody.erase(res1,queryStr1Len);
		parallelTasksFound = true;
	}
	
	while(true)
	{
		size_t res2 = fnBody.find(queryStr2,0);
		if(res2==std::string::npos)
			break;
		else
			fnBody.erase(res2,queryStr2Len);
	}

	/*if the method contains parallel code, then return the serial code enclosed in a conditional that tells when to execute the serial code.
	* by default, the executor stops spawning sub-tasks after 4 levels of executing a recursive implementation of a task. */ 
	if(parallelTasksFound)
	{
		serialCode = string("\n\tif(callStackDepth > recursionDepth+4)\n\t{\n");
		serialCode += fnBody;
		serialCode += string("\n\t}\n");
	}
	return serialCode;
}


/* This function generates function code for accessing cells within data regions. Given a multi-dimensional grid, its cells are accessed using integer coordinates.
 * The generated function code takes as input a list of integer coordinates and the data region. It computes the correct index (assuming a row-major ordering of the cells in memory)
 * and returns a pointer to the cell at that index.*/ 
string FetchCellFnCallCode()
{
	ostringstream ret("");
	vector<string> argList;
	//Construct function parameter list: first, generate variable names to hold integer coordinates. Each dimension needs a variable. Start the variable name with 'i'. 	
	char startVar = 'i';
	for (int i=0;i<gridDimension;i++)
	{
		argList.push_back(string("int ")+string(1,startVar)+string(", "));
		startVar++;
	}
	//next, the dataRegion accessed.
	argList.push_back(string("CELL_TYPE* data"));
	
	string sig("inline CELL_TYPE* GetDPTableCell(");
	for(int i=0;i<argList.size();i++)
		sig+=argList[i];
	sig +=string(")");
	
	/*skip the first few cells of the data region as they hold the meta data (e.g. for 2D grids, the first 5 cells contain: 
 	* <bounding_box_start_x, bounding_box_start_y, side_length_x, side_length_y, taskID>.) Compute the index based on the <x,y> parameter values passed as input.
 	* See the generated function code for GetDPTableCell in GlobalTypes.h for details.*/
	ret <<sig<<"\n";
	ret << "{\n";
	ret << "\tCELL_TYPE* cell = data+METADATASPACE;\n";
	ret << "\tint side = data[DIMENSION];\n";

	startVar='i';
	string lastVarName;
	for(int i=0;i<gridDimension;i++)
	{
		string varName;
		varName = string(1,startVar)+string("Offset");
		lastVarName = varName;
		ret << "\tint "<<varName <<"="<< startVar <<" - data["<<gridDimension-1-i<<"];\n";
		if(i!=gridDimension-1)
		{
			ret <<"\tcell += ("<<varName<<"*"<<gridDimension-1-i<<"* side);\n"; 
		}
		startVar++;
	}
	ret <<"\tcell += "<<lastVarName<<";\n";
	ret <<"\treturn cell;\n}\n\n"; 
	return ret.str();
}

/* This function converts a binary string into its integer equivalent number and returns the number.*/
int ConvertStringToBinary(const char* str)
{
	int len=strlen(str);
	int i=len-1, num=0, bit_position=0;
	do
	{
		num += ((1<<bit_position) * (str[i]-'0'));
		bit_position++;
		i--;
	}while(i>=0);
	return num;
}

/* This function generates code for defining and initializing variables passed as arguments to recursive method invocations. 
 * The input to this function is a list of unique variables with names that identify a data region (x01, y10 etc.).
 * Output is a string of variable definitions initialized appropriately. For example if x identifies data region (0,0,3,3) (x is a box with values <0,0,3,3>), then x01 identifies <2,2,3,3>.
 * This function generates code that derives <2,2,3,3> for initializing in x01 starting from <0,0,3,3> (and assuming that x01 refers to top-right quadrant)*/
string GenerateVarDefinitions(map<char, set<string> >& varsPendingDefinition)
{
	string declarations("");
	map<char, set<string> >::iterator itArg = varsPendingDefinition.begin();
	for(;itArg!=varsPendingDefinition.end();itArg++)
	{
	
		set<string>::iterator itOrth = itArg->second.begin();
		for(;itOrth!=itArg->second.end();itOrth++)
		{
			string var = *itOrth;
			string orthStr = itOrth->substr(1);
			ostringstream box;
			box<<"Box "<<var<<"(";
			int orthant = ConvertStringToBinary(orthStr.c_str());//stoi(itOrth->c_str(),NULL,2); 
			vector<string> offsets(gridDimension);
			//preparing offsets for topleft coordinates
			for(int dim=0;dim<gridDimension;dim++)
			{
				ostringstream side;
				side<<"("<<itArg->first<<"->coords["<<dim+gridDimension<<"]-"<<itArg->first<<"->coords["<<dim<<"]+1)/2";
				offsets[dim] = (orthant & 1)?side.str():""; 
				orthant >>= 1;
			}
			//setting topleft coordinates
			ostringstream topleft;
			for(int dim=0;dim<gridDimension;dim++)
			{
				topleft<<itArg->first<<"->coords["<<dim<<"]";
				if(offsets[dim].compare(""))
					topleft<<"+"<<offsets[dim];
				topleft<<",";
			}
			//preparing offsets for bottomright coordinates
			for(int dim=0;dim<gridDimension;dim++)
			{
				ostringstream side;
				side<<"("<<itArg->first<<"->coords["<<dim+gridDimension<<"]-"<<itArg->first<<"->coords["<<dim<<"]+1)/2";
				if(!(offsets[dim].compare("")))
					offsets[dim]="-"+side.str();
				else
					offsets[dim]="";
			}
			//setting bottomright coordinates
			ostringstream bottomright;
			for(int dim=0;dim<gridDimension;dim++)
			{
				bottomright<<itArg->first<<"->coords["<<gridDimension+dim<<"]";
				if(offsets[dim].compare(""))
					bottomright<<offsets[dim];
				if(dim != gridDimension-1)
					bottomright<<",";
			}
				
			box<<topleft.str()<<bottomright.str()<<");";

			declarations += "\t"+box.str();
			declarations +="\n";
		}
	}
	return declarations;
}

/*This function generates the executor code for executing the tasks. The tasks are created earlier by the inspector. 
 * Since tasks are basically recursive method invocations, the list of recursive method invocations is passed as an input.
 * The ExecuteFunction (generated function definition) itself takes as input a function pointer---which is of type from one among the methods in the list---, 
 * and the data regions (tiles) that the function reads/writes. This function returns the required function declaration and definition for executing tasks.*/
pair<string, string> GenerateDelayedFunctionExecution(const vector<RecMethod*>& recMethods)
{
	//generate function declaration for the header file. ExecuteFunction is called from HelperFunctions.cpp.
	string hdr(""),cpp("");
	hdr += "void ExecuteFunction(FunctionCall* fn, vector<CELL_TYPE*> dataRegions);\n";
	
	//generate function definition.
	cpp += "void ExecuteFunction(FunctionCall* fn, vector<CELL_TYPE*> dataRegions)\n{\n";
	/*depending upon the type of a recursive method, dereference the function pointer and pass the data regions as arguments:
	 *The remaining required arguments for executing a method (callStackDepth and the bounding boxes of the data regions passed) are already populated by the 
	 *inspector during unrolling. For example, if a task represents calling a recursive method A(x, xData, y, yData, callStackDepth)---here x and y are bounding 
	 *boxes (say <0,0,8,8> <9,9,17,17>) of the data regions pointed to by xData and yData, and callStackDepth (say 4) is the depth of the recursion at which the inspector encountered 
	 *this method invocation with the mentioned argument values.---, then the inspector populates the values of x, y, and the callStackDepth during unrolling. dataRegions (xData and yData) 
	 *are created (and hence known) only later when the executor begins to execute the tasks (via call to ExecuteFunction). Hence, inspector fills nullValues in place of xData, and yData. 
	 *When the executor calls ExecuteFunction after the inspector has already created the tasks, it has the information about dataRegions, which are passed as argument values 
	 *to the method invocation. Here it is assumed that the 1st item in the list of data regions (vector<CELL_TYPE*> dataRegions) has the bounding box x, the second item has the bounding 
	 *box y and so on.*/  
	cpp += "\tswitch(fn->functionName)\n\t{\n";
	for(unsigned int i=0;i<recMethods.size();i++)
	{
		cpp += "\t\tcase '"+recMethods[i]->name+"':\n";
		cpp += "\t\t\t{\n";
		ostringstream cppStr;
		cppStr <<"\t\t\t\tassert(dataRegions.size() == "<<recMethods[i]->parameters.size()<<");\n"; 
		cppStr <<"\t\t\t\tDeferredCall<_paramListFunction"<<recMethods[i]->name<<">* df = (reinterpret_cast<DeferredCall<_paramListFunction"<<recMethods[i]->name<<">* >(fn->fnCall));\n";
		for(unsigned int j=0;j<recMethods[i]->parameters.size();j++)
			cppStr<<"\t\t\t\tstd::get<"<<(j*2+1)<<">(df->params) = dataRegions["<<j<<"];\n";
		cppStr<<"\t\t\t\t(reinterpret_cast<DeferredCall<_paramListFunction"<<recMethods[i]->name<<">* >(fn->fnCall))->Run();\n";
		cpp += cppStr.str();
		cpp += "\t\t\t}\n\t\t\tbreak;\n";
	}
	cpp +="\t\tdefault: break;\n";	
	cpp += "\t}\n}\n\n";

	return make_pair(hdr,cpp);
}


/*This function is a wrapper for generating the call to unroll (by calling the top-level method 'A_unroll'). It also is a wrapper for generating the executor code 
 * that executes the tasks. This function returns the required function declarations and definitions that kick-start the unrolling and execution of tasks.*/
pair<string, string> GenerateTopLevelCalls(const vector<RecMethod*>& recMethods)
{
	//generates code for unrolling recursion and creating tasks.
	string hdrDecl("void Unroll();\n");
	string topLevelCallWrapper("void Unroll()\n{\n");
	string topLevelFunctionCall = recMethods[0]->GenerateTopLevelFunctionInvocation();
	topLevelCallWrapper += topLevelFunctionCall;
	topLevelCallWrapper += "}\n\n";
	
	//generates code for executing tasks.
	pair<string,string> defdFnExecCode = GenerateDelayedFunctionExecution(recMethods);
	topLevelCallWrapper += defdFnExecCode.second;
	hdrDecl += defdFnExecCode.first;
	return make_pair(hdrDecl,topLevelCallWrapper);
}

/* This function creates a summary of calls made by the recursive methods. When the recursion is unrolled starting from the top-level recursive method, the resulting tree structure has method 
 * invocations as vertices. Given the depth of recursion to unfold, the summary information helps in computing the number of leaf method invocations without actually unrolling. In other words, 
 * given a number of leaf method invocations (tasks), we can compute the amount of unrolling required. Hence, for a given run with a specific number of processes, we can create tasks as needed based on
 * the unroll depth. This function computes (and returns) the summary data structure's declaration and definition. The summary data structure, fnCallHierarchySummary, contains statistics about the number 
 * of invocations made to a specific recursive method within a method body. E.g. if method A calls method A 2 times, B 5 times, C 4 times then this structure would have an entry {2,5,4}. 
 * This information is used by HelperFunctions.cpp in computing the number of tasks created. */
pair<string,string> CreateFunctionCallHierarchySummary(const vector<RecMethod*>& listOfRecMethods)
{
	string hdrDecl(""), cppDefn("");
	
	//the below code computes <method, calls<method, number_of_times> >
	map<char, map<char, int> > fnCallHierarchySummary;
	for(unsigned int h=0;h<listOfRecMethods.size();h++)
	{
		char fnName = listOfRecMethods[h]->name[0];
		for(unsigned int i=0;i<listOfRecMethods[h]->calls.size();i++)
		{
			for(unsigned int j=0;j<listOfRecMethods[h]->calls[i].size();j++)
			{
				char fnInvoked = listOfRecMethods[h]->calls[i][j].name[0];
				map<char, int> invocations;
				invocations.insert(make_pair(fnInvoked,1));
				if(fnCallHierarchySummary.find(fnName) == fnCallHierarchySummary.end())
					fnCallHierarchySummary.insert(make_pair(fnName,invocations));
				else
				{
					if(fnCallHierarchySummary[fnName].find(fnInvoked) == fnCallHierarchySummary[fnName].end())
						fnCallHierarchySummary[fnName].insert(make_pair(fnInvoked, 1));
					else
						(fnCallHierarchySummary[fnName])[fnInvoked] +=1;
				}
			}
		}
	}

	//generate header declaration.
	ostringstream ostrHdr;
	ostrHdr<<"extern int fnCallHierarchySummary["<<fnCallHierarchySummary.size()<<"]["<<fnCallHierarchySummary.size()<<"];\n";
	hdrDecl += ostrHdr.str();

	//generate definition.
	ostringstream ostr;
	ostr<<"int fnCallHierarchySummary["<<fnCallHierarchySummary.size()<<"]["<<fnCallHierarchySummary.size()<<"]={";
	map<char, map<char, int> >::iterator callTypeIter = fnCallHierarchySummary.begin();
	for(int i=0;callTypeIter!=fnCallHierarchySummary.end();callTypeIter++,i++)
	{
		ostr<<"{";
		map<char, int>::iterator childrenTypeIter = (callTypeIter->second).begin();
		map<char, map<char, int> >::iterator callTypeIter2 = fnCallHierarchySummary.begin();
		unsigned int j=0;
		for(j=0;callTypeIter2!=fnCallHierarchySummary.end();callTypeIter2++,j++)
		{
			if(childrenTypeIter->first == callTypeIter2->first)
			{
				ostr<<childrenTypeIter->second;
				childrenTypeIter++;
			}
			else
				ostr<<"0";
			if(j!=(fnCallHierarchySummary.size()-1))
				ostr<<",";
		}
		ostr<<"}";
		if(i != (fnCallHierarchySummary.size()-1))
			ostr<<",";	
	}
	ostr<<"};\n";
	cppDefn += ostr.str();
	return make_pair(hdrDecl,cppDefn);
}

/* This function generates the definition of the type for the bounding box. 
 * A bounding box contains top-left and bottom-right coordinates and hence, depends on the grid dimensions.*/
string CreateBoxDefinition()
{
	ostringstream defn;
	defn<<"\ntypedef struct Box\n{";
	defn<<"\tint coords["<<gridDimension*2<<"];\n"; //field for storing topleft and bottom right coordinates
	//create default constructor
	defn<<"\tBox(){}\n\tBox("; 
	for(int i=0;i<2*gridDimension;i++)
	{
		defn<<"int "<<char('a'+i);
		if(i!=(2*gridDimension-1))
			defn<<", ";
	}
	defn<<")\n\t{\n\t\t";
	for(int i=0;i<2*gridDimension;i++)
	{
		defn<<"coords["<<i<<"]="<<char('a'+i)<<";";
	}
	defn<<"\n\t}\n";
	//create copy constructor	
	defn<<"\tBox(const Box&b)\n\t{\n\t\t";
	for(int i=0;i<2*gridDimension;i++)
		defn<<"coords["<<i<<"]=b.coords["<<i<<"];";
	defn<<"\n\t}\n";
	//overload assignment operator
	defn<<"\tbool operator==(const Box& rhs)\n\t{\n";
	defn<<"\t\tbool flag = true;\n\t\tif(";
	for(int i=0;i<2*gridDimension;i++)
	{
		defn<<"(this->coords["<<i<<"]!=rhs.coords["<<i<<"])";
		if(i!=(2*gridDimension-1))
			defn<<" || ";
	}
	defn<<")\n\t\t\tflag=false;\n\t\treturn flag;\n\t}\n";
	//define a method that returns the bounding box size (number of cells in the grid represented by the bounding box).
	defn<<"\tlong int GetBoxSize()\n\t{\n\t\tlong int len = 1;\n";
	defn<<"\t\tfor(int i=0;i<DIMENSION;i++)\n\t\t{\n";
	defn<<"\t\t\tlen *= (coords[i+DIMENSION]-coords[i]+1);\n\t\t}\n\t\treturn len;\n\t}\n";
	
	//define a helper method (used while debuggind only) that prints the coordinates.
	defn<<"\tchar* PrintStr(char *str)const\n\t{\n";
	defn<<"\t\tsprintf(str,\"";
	for(int i=0;i<2*gridDimension;i++)
		defn<<"%d";
	defn<<"\\n\",";
	for(int i=0;i<2*gridDimension;i++)
	{
		defn<<"coords["<<i<<"]";
		if(i!=(2*gridDimension-1))
			defn<<",";
	}
	defn<<");\n\t}\n";
	defn<<"}Box;\n";
	return defn.str();
}

/* This function parses the raw string in the input spec and sets the 'name' field of the Invocation structure. 
 * The invocation structure contains details about the call made within a method body.
 * Assumption: a recursive method call is in this format: X(...args....), where X is a single character that names the call. */
void Invocation::ExtractName(const string& line, int startPos)
{
	int nameEnd = line.find_first_of('(', startPos);
	assert(nameEnd != string::npos);
	//find the first blankspace before '(' if it exists.
	int nameStart=-1;
	int i=nameEnd-1;
	//skip any trailing blankspaces.
	for(i=nameEnd-1;i>=0;i--)
	{
		if(line[i]==' ')
			break;
	}
	nameStart=i+1;
	name = line.substr(nameStart, nameEnd-nameStart);
}

/* This function parses the raw string in the input spec and sets the 'calls' field of RecMethod structure. 
 * calls field contains a list of calls to other methods (or self) done within a recursive method body.*/
void RecMethod::PopulateCalls()
{
	calls.resize(bodyStr.size());
	for(int i=0;i<bodyStr.size();i++)
	{
		GetCalls(bodyStr[i], calls[i]);
	}
}

/* This function generates the 'Unroll' method, which is used by the inspector to create tasks.
 * Tasks are necessary during multi-process runs.
 * Unrolling begins with a call to the top-level recursive method 'A' (A_unroll). Hence, invoking this method with appropriate method arguments is necessary.
 * This function constructs the method arguments. Also, this function generates code to handle the case when no unrolling is specified. 
 * Single-process runs do not require unrolling (when input sizes are small enough to fit in the memory). 
 * Also no unrolling is required in multi-process runs when the program follows a data-parallel execution model (resulting code is not optimized in this case). */ 
string RecMethod::GenerateTopLevelFunctionInvocation()
{
	ostringstream topLevelCallArg;
	/* The below code creates a bounding box with (0, 0, inputSize, inputSize) as values. This means that an entire grid is to be computed. 
	Assumption: entire grid (0,0, inputSize, inputSize) is computed. In some scenarios (MCM) when the 0th row and column are used as initializer cells and not really computed, 
	an end-user needs to explicitly modify the arguments to A_unroll. E.g. MCM computes all the cells starting from 1st row and column. Therefore, in calling A_unroll(Box* x...), 
	end-user should provide the bounding box info used as an argument to A_unroll as x(1,1,inputSize,inputSize). */

	for(unsigned int i=0;i<parameters.size();i++)
	{
		//Create top-left coordinates of bounding box
		topLevelCallArg<<"\tint parentTileID"<<parameters[i]<<"=0;\n\tBox *"<<parameters[i]<<" = new Box(";
		for(unsigned int j=0;j<gridDimension;j++)
				topLevelCallArg<<0<<",";

		//Create bottom-right coordinates of bounding box
		for(unsigned int j=0;j<gridDimension;j++)
		{
			topLevelCallArg<<"inputSize"<<parameters[i]<<"["<<j<<"]";
			if(j!=gridDimension-1)
				topLevelCallArg<<",";
			else
				topLevelCallArg<<");\n";
		}
	
	}
	string topLevelCall = topLevelCallArg.str();
	
	/*When recursion depth is 0, no unrolling is necessary. Create a computation graph to be executed by the executor.
 	* The executor needs to know the method to call and the values of arguments to pass to the method. */
	//Create arguments to be used with A_unroll 
	vector<string> argList;
	for(unsigned int i=0;i<parameters.size();i++)	
	{
		argList.push_back(parameters[i]);
		if(i==0)
			argList.push_back("nullData");
		else
		{
			//more than one parameter to the top-level function indicates that the second parameter onwards are read-only parameters supplied as system inputs.
			ostringstream tmp;
			tmp<<parameters[i]<<"Data";
			argList.push_back(tmp.str());
		}
		argList.push_back("0"); //this is dummy argument added so that the number of arguments passed to GenerateStopCondWhileUnrolling is consistent among the two call sites to it.
	}
	argList.push_back("0");
	/*The below code does the equivalent of creating a vertex in the computation graph. 
	 The vertex here is the call A_unroll(x, nullptr, 0, 0), where x is the bounding box, nullptr is the data region of x (initially empty), 0 is the tile ID and callStack depth.*/
	topLevelCall += "\t"+GenerateStopCondWhileUnrolling(argList);
	topLevelCall += "\t"+name+"_unroll(";
	for(unsigned int i=0;i<argList.size()-1;i+=3)
	{
		topLevelCall += argList[i];
		if(i==0)
			topLevelCall+=", nullptr ";
		else
			topLevelCall+=", "+argList[i+1];
		topLevelCall += ", "+argList[i+2]+", ";
	}
	
	topLevelCall +="0);\n";
	for(unsigned int i=0;i<parameters.size();i++)
		topLevelCall+="\tdelete "+parameters[i]+";\n";
	return topLevelCall;
}


/* Function for naming recursive method types (executor versions). 
 * The defined names are a convenient (short) when creating pointers to method types used by the inspecor and executor codes.*/ 
string RecMethod::GenerateParamListMacro() const
{
	string paramList;
	paramList +="#define _paramListFunction"+name+" ";
	for(unsigned int j=0;j<parameters.size();j++)	
		paramList+="Box*, CELL_TYPE*, ";
	paramList+="int\n";
	return paramList;
}

/* Now the tricky bit.... */
string RecMethod::GenerateStopCondWhileUnrolling(const vector<string>& fnArgs)
{
	string stopCond("");
	stopCond += "if("+fnArgs[fnArgs.size()-1]+" == recursionDepth)\n";
	stopCond+="\t{\n";

	//assuming that the first parameter is the write-parameter.
	ostringstream tmpStr;
	if(computeGrid != COMPUTE_SPARSE)
		tmpStr<<"\t\tint writeTileID = GetTileID"<<gridDimension<<"D(parentTileID"<<fnArgs[0]<<");\n";
	else
		tmpStr<<"\t\tint writeTileID = GetTileID("<<fnArgs[0]<<");\n";
		
	stopCond+=tmpStr.str();
	stopCond+="\t\tbool localUpdate = false;\n";
	stopCond+="\t\tFunctionCall* fnCall= nullptr;\n";
	stopCond+="\t\ttileUpdateLog[writeTileID]=fnCounter;\n";
	if(computeGrid != COMPUTE_SPARSE)
		stopCond+="\t\tif(IsLocalOwner(writeTileID))\n\t\t{\n";
	else
		stopCond+="\t\tif(true)\n\t\t{\n";
	stopCond+="\t\t\tCELL_TYPE* nullData=nullptr;\n";
	stopCond+="\t\t\tfnCall=new FunctionCall();\n";
	stopCond+="\t\t\tfnCall->functionName = '"+name+"';\n";
	stopCond+="\t\t\tParameter p;\n";

	string argList("");
	for(unsigned int j=0;j<fnArgs.size()-1;j+=3)	
	{
		ostringstream stopCondStr;
		stopCondStr<<"\t\t\tBox* b"<<j<<"=new Box(*"<<fnArgs[j]<<");\n"; 
		stopCondStr<<"\t\t\tp.data = b"<<j<<";\n";
		stopCondStr<<"\t\t\tp.tile = "<<fnArgs[j+1]<<";\n";
		ostringstream tmpStr;
		tmpStr<<"b"<<j;
		tmpStr<<", "<< fnArgs[j+1]<<", "; 
		argList += tmpStr.str();
		stopCond += stopCondStr.str();
		stopCond +="\t\t\tfnCall->params.push_back(p);\n";
	}
	argList += fnArgs[fnArgs.size()-1]+"+1";

	stopCond +="\t\t\tstd::tuple<_paramListFunction"+name+"> t = std::make_tuple("+argList+");\n";
  	stopCond +="\t\t\tDeferredCall<_paramListFunction"+name+">* defdCall = new DeferredCall<_paramListFunction"+name+">();\n";
	stopCond +="\t\t\tdefdCall->params=t;\n\t\t\tdefdCall->fptr="+name+";\n";
	stopCond +="\t\t\tfnCall->fnCall = defdCall;\n";
	stopCond +="\t\t\tfnCall->ID = fnCounter;\n";
	stopCond +="\t\t\tif(fnCalls[writeTileID].size() > 0)\n\t\t\t{\n";
	stopCond +="\t\t\t\tfnCall->numInPorts +=1;\n";
	stopCond +="\t\t\t\tfnCall->wawSource = (fnCalls[writeTileID].back())->ID;\n"; 
	if(computeGrid != COMPUTE_SPARSE)
	{
		stopCond +="\t\t\t\tFunctionCall* lastFunctionToUpdate = fnCalls[writeTileID].back();\n";
		stopCond +="\t\t\t\tlastFunctionToUpdate->outPortOwners.insert(procRank);\n";
	}
	stopCond +="\t\t\t}\n\t\t\tfnCalls[writeTileID].push_back(fnCall);\n";
	stopCond +="\t\t\tlocalUpdate = true;\n\t\t}\n";
	stopCond +="\t\tfnCounter++;\n";

	if((fnArgs.size()-1)/3 > 1)
	{	
		ostringstream tmpStr;
		tmpStr <<"\t\tint readTileIDs["<<((fnArgs.size()-1)/3) -1<<"];\n";
		tmpStr <<"\t\tCELL_TYPE* tiles["<<((fnArgs.size()-1)/3) -1<<"];\n";
		for(int i=1;i<(fnArgs.size()-1)/3;i++)
		{
			if(computeGrid != COMPUTE_SPARSE)
				tmpStr <<"\t\treadTileIDs["<<i-1<<"]=GetTileID"<<gridDimension<<"D(parentTileID"<<fnArgs[3*i]<<");\n";	
			else
				tmpStr <<"\t\treadTileIDs["<<i-1<<"]=GetTileID("<<fnArgs[3*i]<<");\n";	
		}
		for(int i=1;i<(fnArgs.size()-1)/3;i++)
			tmpStr <<"\t\ttiles["<<i-1<<"]="<<fnArgs[3*i+1]<<";\n";	

		tmpStr<<"\t\tfor(int i=0;i<"<<((fnArgs.size()-1)/3)-1<<";i++)\n\t\t{\n";
		tmpStr<<"\t\t\tint readTileID=readTileIDs[i];\n";
		tmpStr<<"\t\t\tCELL_TYPE* curTile=tiles[i];\n";
		tmpStr<<"\t\t\tif((curTile==nullptr) && (readTileID != writeTileID))\n\t\t\t{\n";
		tmpStr<<"#ifdef TASK_AGGREGATION\n";
		tmpStr<<"\t\t\t\tif(fnCalls[readTileID].size() > 0)\n";	
		tmpStr<<"\t\t\t\t\t(fnCalls[readTileID].back())->isReadBeforeNextWrite = true;\n";
		tmpStr<<"#endif\n";			
		tmpStr<<"\t\t\t\tif(localUpdate)\n\t\t\t\t{\n";
		tmpStr<<"\t\t\t\t\tfnCall->numInPorts +=1;\n";
		tmpStr<<"\t\t\t\t\tfnCall->params[i+1].portID = tileUpdateLog[readTileID];\n\t\t\t\t}\n";
		if(computeGrid != COMPUTE_SPARSE)
		{
			tmpStr<<"\t\t\t\tif(fnCalls[readTileID].size() > 0)\n\t\t\t\t{\n";
			tmpStr<<"\t\t\t\t\tFunctionCall* lastFunctionToUpdate = fnCalls[readTileID].back();\n";
			tmpStr<<"\t\t\t\t\tlastFunctionToUpdate->outPortOwners.insert(GetOwner(writeTileID));\n\t\t\t\t}\n";
		}
		tmpStr<<"\t\t\t}\n";
		tmpStr<<"\t\t\telse if((curTile!=nullptr) && localUpdate)\n";
		tmpStr<<"\t\t\t\tfnCall->params[i+1].portID = readTileID;\n";
		tmpStr<<"\t\t}\n";
		stopCond += tmpStr.str();
	}
	
	stopCond +="\t\treturn;\n\t}\n";

	return stopCond;

}

/* This function is very similar to GenerateFunctionBody in structure except that the resulting code is executed by the D2P inspecor.
 * The inspector executes a slightly different version of recursive method bodies compared to those executed by the executor.
 * If executor executes a method A(Box* x, CELL_TYPE* xData int callStackDepth), then inspector executes A_unroll(Box* x, int parentTileID, int callStackDepth).
 * To efficiently create tasks, and facilitate their dependency computation and partitioning, inspector identifies each tile by its Z-number in a Z-Morton ordering.
 * E.g. For the top-level method 'A' computing the entire grid, the only tile ID it knows is 0. When this grid is decomposed into smaller parts due to recursive invocation of 'A'
 * we get 4 tiles/quadrants (in case of a 2D grid), which are given tile numbers 0-3. When a further decompsition of the quadrants into sub-quadrants happens, sub-quadrants in tile 0 
 * are given tile numbers 0-3 (older tile number 0 representing a bigger quadrant becomes invalid). This way of generating tile IDs follows a Z-Order numbering of cells in a 2-D grid.
 * Deriving tile numbers for the sub-quadrants from the tile ID of the parent quadrant is done as follows: parentTileID*K+offset where, K is maximum number of parts/tiles possible 
 * (4 for a 2D grid) and offset is the quadrant identifier (top-left=0, top-right=1, bottom-left=2, bottom-right=3).*/   
string RecMethod::GenerateFunctionBodyForUnroll()
{
	int K = 1 << gridDimension;
	//function signature
	ostringstream sig;
	vector<string> argList;
	sig<<"void "<<name<<"_unroll(";
	for(unsigned int j=0;j<parameters.size();j++)	
	{
		sig<<"Box* "<<parameters[j]<<", CELL_TYPE* "<<parameters[j]<<"Data, int parentTileID"<<parameters[j]<<", ";
		argList.push_back(parameters[j]);
		argList.push_back(parameters[j]+"Data");
		argList.push_back("parentTileID"+parameters[j]);
	}
	argList.push_back("callStackDepth");
	sig<<"int callStackDepth)";
	signature_unroll = sig.str();
	
	//Stop condition for the inspector code.
	string stopCond = "\t"+GenerateStopCondWhileUnrolling(argList);
	//function body
	map<char, set<string> > varsPendingDefinition;
	string ret("\n");
	for(unsigned int i=0;i<calls.size();i++)
	{
		bool parallelFlag = false;
		for(unsigned int j=0;j<calls[i].size();j++)
		{
			if(calls[i].size() > 1)
				parallelFlag=true;
			ret += "\t"+calls[i][j].name+"_unroll(";
			int numArgs = calls[i][j].arguments.size();
			for(int k=0;k<numArgs;k++)
			{
				string var="&"+calls[i][j].arguments[k]+", "+calls[i][j].arguments[k][0]+"Data";
				ostringstream  childZMortonCode;
				childZMortonCode<<"parentTileID"<<calls[i][j].arguments[k][0]<<"*"<<K<<"+"<<ConvertStringToBinary((calls[i][j].arguments[k].substr(1)).c_str());
				string varData = childZMortonCode.str();
				set<string> uniqArgs;
				uniqArgs.insert(calls[i][j].arguments[k]);
				pair<map<char,set<string> >::iterator, bool> status = varsPendingDefinition.insert(make_pair(calls[i][j].arguments[k][0],uniqArgs));
				if(!status.second)
					((status.first)->second).insert(calls[i][j].arguments[k]);
				ret += var;
				ret += ", "+varData;
				if(k != numArgs-1)
					ret += ", ";
			}
			ret += ", callStackDepth+1);";
		}
		ret += "\n";
	}
	ret+="\treturn;\n";

	//variable declarations corresponding to the variables used in function body.
	string declarations = GenerateVarDefinitions(varsPendingDefinition);
	return signature_unroll+"\n{\n"+stopCond+"\n"+declarations+ret+"}\n\n";
}

/*Generates the stop condition for a recursive method body.
*The recursion terminates when the hierarchical decomposition results in parts that cannot be broken down further:
*by default, terminates when a single cell of a grid is about to be computed */
string RecMethod::GenerateStopCond(const string& writeParam)
{
	ostringstream condStrm;
	string cond("\tif(");
	for(unsigned int i=0;i<gridDimension;i++)	
	{
		condStrm <<"("<<writeParam<<"->coords["<<i<<"]=="<<writeParam<<"->coords["<<i+gridDimension<<"])";
		if(i != gridDimension-1)
			condStrm<<" && ";
	}
	cond += condStrm.str();
	cond +=")\n\t{\n";
	cond +="\t\t//Write the code for terminating case here.\n";
	cond +="\t\treturn;\n";
	cond +="\t}\n";
	return cond;
	
}

/* This function translates the data stored in a string format into C++ code. 
 * The code generated is the body of a recursive method invocation. The body contains:
 * 1) terminating case 2) variable definitions 3) sequential code for serial execution 4) parallel version of the sequential code.
 * The resulting code is executed by the D2P executor. */   
void RecMethod::GenerateFunctionBody()
{
	//function signature
	ostringstream sig;
	sig<<"void "<<name<<"(";
	for(unsigned int j=0;j<parameters.size();j++)	
	{
		sig<<"Box* "<<parameters[j]<<", CELL_TYPE* "<<parameters[j]<<"Data, ";
	}
	/*for any method invocation, D2P codegen inserts an additional argument callStackDepth that is used for multiple purposes:
	* i) in the inspector code to stop unrolling the recursion as indicated by the variable recursionDepth. 
	* ii) in the executor code to stop executing the parallel version of the method body and switch to execution of the serial code to avoid
	* creating (spawning) too many fine-grained tasks.*/
	sig<<"int callStackDepth)";
	signature = sig.str();
	
	//function body
	map<char, set<string> > varsPendingDefinition;
	string ret("\n");
	for(unsigned int i=0;i<calls.size();i++)
	{
		bool parallelFlag = false;
		/*assumes that parallelism is specified in the input spec: multiple method invocations specified on a single line of input spec are capable of being executed in parallel.
 		* TODO: extract parallelism without this assumption.*/
		for(unsigned int j=0;j<calls[i].size();j++)
		{
			if(calls[i].size() > 1)
				parallelFlag=true;
			if(parallelFlag && (j!=calls[i].size()-1))
				ret += "\n#ifdef PARALLEL\n\tcilk_spawn\n#endif\n";
			ret += "\t"+calls[i][j].name+"(";
			int numArgs = calls[i][j].arguments.size();
			for(int k=0;k<numArgs;k++)
			{
				/* extract grid dimension only the first time when a method argument is encountered. gridDimension is a global variable used by multiple methods of CodeGen.
 				*Assumption: argument identifies the unique data region / tile that it is writing to. E.g. x00 specifies top-left quadrant in a 2D grid, x111 specifies
				the 8th octant in a 3D grid etc. Also assumed that argument length (excluding region identifiers) is 1. */ 
				if(gridDimension == -1)
					gridDimension =  (calls[i][j].arguments[k]).length()-1; 
				//arguments[k] represents the bounding box of a tile. Recursive methods are always invoked with pointers to these boxes to avoid copying.
				string var="&"+calls[i][j].arguments[k];
				//this argument represents the actual data region. It is a pointer to the tile. E.g. xData is the tile with bounding box x00. 
				string varData = string(1,calls[i][j].arguments[k][0])+"Data";
				/*Among all the references to variables made in method invocations, unique variable names are collected and defined at the beginning of the method body.
				varsPendingDefintion stores such unique variables for the purpose of defining them later.*/
				set<string> uniqArgs;
				uniqArgs.insert(calls[i][j].arguments[k]);
				pair<map<char,set<string> >::iterator, bool> status = varsPendingDefinition.insert(make_pair(calls[i][j].arguments[k][0],uniqArgs));
				if(!status.second)
					((status.first)->second).insert(calls[i][j].arguments[k]);
				ret += var;
				ret += ", "+varData;
				if(k != numArgs-1)
					ret += ", ";
			}
			ret += ", callStackDepth+1);";
			if(parallelFlag && (j==calls[i].size()-1))
				ret +="\n#ifdef PARALLEL\n\tcilk_sync;\n#endif\n";
		}
		ret += "\n";
	}
	ret+="\treturn;\n";
	
	/*A global variable computeGrid is set based on the number of quadrants (or octants) computed by the top-level recursive method named 'A' during a decomposition step.
	*Since we assume that there is only one top-level recursive method in the spec, the method's body can be inspected to correctly determine if an algorithm computes 
	*a triangular matrix (_UTM), full grid (_FULL) or some orthants of the grid.*/
	assert(name.length() == 1);
	if(name[0] == 'A') 
	{
		if(varsPendingDefinition[parameters[0][0]].size() == (1<<gridDimension)-1)
			computeGrid = COMPUTE_UTM;
		else if(varsPendingDefinition[parameters[0][0]].size() == (1<<gridDimension))
			computeGrid = COMPUTE_FULL;
		else
			computeGrid = COMPUTE_SPARSE;
	}
	
	//A stop condition with an empty terminating case (by default, the recursion stops when a single cell is to be computed)
	string stopCond = GenerateStopCond(parameters[0]);
	//variable declarations corresponding to the variables used in function body.
	string declarations = GenerateVarDefinitions(varsPendingDefinition);
	string serialCode = GenerateSerialCodeInFunctionBody(ret);
	//bodyStr field in RecMethod structure is set.
	bodyStr.push_back(signature+"\n{\n"+stopCond+"\n"+declarations+serialCode+ret+"}\n\n");
}

/* This function parses an input file and produces two files RecursiveFunctions.cpp and GlobalTypes.h.
 * GlobalTypes.h contains common data and function declarations used in RecursiveFunctions.cpp and HelperFunctions.cpp.
 * An end-user is required to add to GlobalyTypes.h the declarations of any variables or functions used in Recursion Terminating
 * cases of methods in RecursiveFunctions.cpp. The terminating conditions are part of the reference implementations
 * provided in xx_Aux.cpp, which contains problem-specific auxiliary methods such as reading input, initializing DPtable, Printing 
 * score, and the terminating case. Also, when applicable, an end-user is required to change the data type of the cell (of a DP table) 
 * in GlobalTypes.h.*/
int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		printf("Usage: ./exe <spec>\n");
		exit(0);
	}

	vector<RecMethod*> listOfRecMethods;
	RecMethod* newMethod=NULL;
	
	ifstream inFile(argv[1], ifstream::in);
	while(!inFile.eof())
	{
		string line;
		getline(inFile, line);
		//assume that a line that has an underscore is actually the first line of a recursive method definition.
		size_t pos = line.find("_");
		if(pos!=string::npos)
		{
			//if the line has an underscore, create a new RecMethod structure, get its name and parameters, and store it in a list.
			newMethod = new RecMethod();
			listOfRecMethods.push_back(newMethod);	
			GetMethodNameAndParams(line,newMethod);
		}
		else
		{
			//otherwise a line represents a recursive method body. store the raw format in the RecMethod structure last created.
			if(newMethod)
				newMethod->bodyStr.push_back(line);
		}
	}

	/*for each recursive method found while parsing, extract the list of invocations made to other (self) recursive methods.
	generate code containing the empty terminating case, variable definitions, and serial and parallel version of the body. 
	store the generated code in RecMethod structure (in the field bodyStr).*/
	for(int i=0;i<listOfRecMethods.size();i++)
	{
		(listOfRecMethods[i])->PopulateCalls();
		(listOfRecMethods[i])->bodyStr.clear();
		(listOfRecMethods[i])->GenerateFunctionBody();
		//cout<<body<<endl;
		
		
	}

	//create a handle to the file produced as an output. This file contains the D2P inspector and executor codes.
	string directory(".");
	filebuf fbcpp;
	string cppFile(directory+string("/RecursiveFunctions.cpp"));
	fbcpp.open(cppFile.c_str(),ios::out);
	ostream oscpp(&fbcpp);
	
	//Create call hierarchy summary.
	pair<string,string> fnCallHierarchySummary = CreateFunctionCallHierarchySummary(listOfRecMethods);
	oscpp<<"#include \"HelperFunctions.h\"\n";
	oscpp<<fnCallHierarchySummary.second;
	oscpp<<"map<int, vector<FunctionCall*> > fnCalls;\n"; //as an intermediate storage space to determine dependencies in inspector.
	oscpp<<"extern int recursionDepth;\n"; //for how much to unroll and controlling spawn overhead in executor. 
	oscpp<<"int fnCounter=0;\n"; //to identify leaf tasks (task ID).
	oscpp<<"map<int,int> tileUpdateLog;\n";

	//create a handle to the file produced as an output. This file contains common function and data structure declarations.
	filebuf globalTypesHeaderBuf;
	string globalTypesFileName("GlobalTypes.h");
	globalTypesHeaderBuf.open(globalTypesFileName.c_str(),ios::out);
	ostream globalTypesHeader(&globalTypesHeaderBuf);
	globalTypesHeader<<"#pragma once\n";
	globalTypesHeader<<"#define DIMENSION "<<gridDimension<<"\n"; //DIMENSION of the grid extracted from the spec.
	/*When a grid is decomposed, the smaller parts are stored as tiles. METADATASPACE reserves storage within each tile to identify the part of the grid it represents.
 	* For the problems encountered so far, it contains (top-left, bottom-right) coordinates of the grid. The extra +1 space is reserved to mark the task ID computing that grid.*/ 
	globalTypesHeader<<"#define METADATASPACE (2*DIMENSION+1)\n"; 
	globalTypesHeader<<"typedef int CELL_TYPE;\n"; //Write the user defined DP Table Cell type into GlobalTypes.h. Default:int
	globalTypesHeader<<fnCallHierarchySummary.first;
	globalTypesHeader<<FetchCellFnCallCode(); //Function definition for indexing into the tile.
	globalTypesHeader<<CreateBoxDefinition();
	globalTypesHeaderBuf.close();


	//Below code formats the previously parsed and stored data in RecMethod structures into an acceptable C++ syntax.
	string fnDeclarations("\n");
	string macroDefns("\n");;
	string recFunctionDefns("");
	//code for recursive method definitions and declarations (executor code). 
	for(unsigned int i=0;i<listOfRecMethods.size();i++)
	{
		fnDeclarations +=listOfRecMethods[i]->signature+";\n";
		macroDefns += listOfRecMethods[i]->GenerateParamListMacro();
		recFunctionDefns +=listOfRecMethods[i]->bodyStr[0];
	}
	
	//code for unroll methods (inspector code). 
	for(unsigned int i=0;i<listOfRecMethods.size();i++)
	{
		recFunctionDefns += listOfRecMethods[i]->GenerateFunctionBodyForUnroll();
		fnDeclarations += listOfRecMethods[i]->signature_unroll+";\n";
	}
	pair<string,string> topLevel = GenerateTopLevelCalls(listOfRecMethods);
	if(computeGrid == COMPUTE_UTM)
		oscpp<<"int computeGrid = COMPUTE_UTM;\n";
	else if(computeGrid == COMPUTE_FULL)
		oscpp<<"int computeGrid = COMPUTE_FULL;\n";
	else
		oscpp<<"int computeGrid = COMPUTE_SPARSE;\n";
	
	for(unsigned int i=0;i<listOfRecMethods[0]->parameters.size();i++)
	{
		oscpp<<"extern long int inputSize"<<listOfRecMethods[0]->parameters[i]<<"["<<gridDimension<<"];\n";
		if(i > 0)
			oscpp<<"extern CELL_TYPE* "<<listOfRecMethods[0]->parameters[i]<<"Data;\n";
	}

	oscpp<<macroDefns;
	oscpp<<fnDeclarations;
	oscpp<<recFunctionDefns;
	oscpp<<topLevel.second;
	fbcpp.close();

	inFile.close();
	return 0;
}
