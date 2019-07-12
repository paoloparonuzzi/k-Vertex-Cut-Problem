#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip> // for precision
#include "ilcplex/ilocplex.h"

using namespace std;

#include "Compact.h"
#include "StdPath_TerminalModel.h"
#include "StdPath_withLongPath.h"
#include "LP_Model.h"
#include "BilevelModel.h"
#include "BilevelModel_withLeaf.h"
#include "SmartModel.h"
#include "Graph_v4.h"
#include "global_functions.h"
#include "global_variables.h"
//#define perform_check_info

int main(int argc, char** argv) {

	const char* fileName;
	istname=new char[2000];
	if (argc == 10) {
		strcpy(istname, argv[1]);
		algorithm=atoi(argv[2]);
		option=atoi(argv[3]);
		k=atoi(argv[4]);
		tolerance=atof(argv[5]);
		timeLimit=atof(argv[6]);
		userFreq=atoi(argv[7]);
		useWeight=atoi(argv[8]);
		fileName = argv[9]; 
	}
	else {cout << "ERROR NUMBER STANDARD PARAMETERS" << endl; exit(2);}

	cout << "\n\nINSTANCE: ->\t" <<  istname << endl ;
	cout << "algorithm: ->\t" <<  algorithm << endl ;
	cout << "k: ->\t" <<  k << endl ;
	cout << "timeLimit: ->\t" <<  timeLimit << endl ;

	int *heads=NULL;
	int *tails=NULL;
	double *weights_arcs=NULL;
	double *weights_nodes=NULL;
	
	////////////////////////////////DIMACS////////////////////////////////
	//the heads and tails cannot be dimensioned inside the functions and in this
	//point we still don't know the dimensions
	heads=new int[10000000];
	tails=new int[10000000];
	weights_arcs=new double[10000000];
	weights_nodes=new double[10000000];

	ReadDIMACSFile(istname,&vertex_number,&edge_number,tails,heads,false);

	for(int i=0;i<edge_number;i++){weights_arcs[i]=1.0;}
	for(int i=0;i<vertex_number;i++){weights_nodes[i]=1.0;}

	cout << "GRAPH BUILDING\n";
	graphFF G = buildGraphFF(vertex_number,edge_number,heads,tails,weights_nodes,weights_arcs,1);
	cout << "DONE\n";

	/////////////////////////////////////////////////
	cout << "COMPLEMENTARY GRAPH BUILDING\n";
	graphFF G_bar = buildComplementaryGraphFF_undirected(G,1);
	cout << "DONE\n";
	/////////////////////////////////////////////////

	cout << "NODES\t" << G->n << endl;
	cout << "ARCS\t" << G->m << endl;

#ifdef 	perform_check_info
	checks_info(G,istname,G->n,G->m,G->T,G->H,false);
	exit(-1);
#endif

	//read vertex weights
	vertexWeight=new double[vertex_number];
	string nameFile(istname);
	string weightFileString = nameFile+".w";
	string line;
	const char* weightFile = weightFileString.c_str();
	ifstream file(weightFileString.c_str());
	if(file.is_open()){
		for(int i=0; i<vertex_number; i++){
			file >> vertexWeight[i];
			std::getline(file, line);
		}
	}

	int FSerr = CheckFS(G);
	printf("Errors FS: %1d\n", FSerr);

	int BSerr = CheckBS(G);
	printf("Errors BS: %1d\n", BSerr);

	delete[] heads;
	delete[] tails;
	delete[] weights_nodes;
	delete[] weights_arcs;
	
	//number of vertices that can be fixed in the separator	
	int *bound;
	bound = new int[vertex_number];
	preFixed = new bool[vertex_number]();
	
	for(int i=0; i<vertex_number; i++)
		bound[i]=maxStableSet(G,i);//prop2	
	
	int v=0;
	int totFixed=0;

	while(v<vertex_number){
		if(k>=bound[v]+2 && !preFixed[v]){
			preFixed[v]=true;
			totFixed++;
			for(int i=0; i<vertex_number; i++){
				if(!preFixed[i])
					bound[i]=maxStableSet(G,i);//prop2
			}
			v=0;			
		}
		else v++;
	}
	delete [] bound;

	///////////////////////////////////////////////////////////////////////

	
	output.open(fileName, std::ios_base::app);

	if(output.is_open()){
		output << istname << "\t" << G->n << "\t" << G->m << "\t" << algorithm << "\t";
	}
	if(output.is_open()){
		if(algorithm==1)
			output << "COMP\t";
			
		if(algorithm==2)
			output << "REP\t";
			
		if(algorithm==22)
			output << "REP_lp\t";
			
		if(algorithm==4)
			output << "NAT\t";
			
		if(algorithm==44)
			output << "NAT_s\t";
			
		if(algorithm==7)
			output << "HYB\t";
	}
	if(output.is_open()){
		output << k << "\t" << tolerance << "\t" << userFreq << "\t" << totFixed << "\t";
	}
	
	////////////////////////////////////////////////////////////////////////
	if(algorithm==1)
	{
		cout << "\n---->Compact Formulation\n\n";
		compact_model_load(G,k,q,b);
		compact_model_solve(istname,G,k,q,b);
		compact_model_free();
	}
	////////////////////////////////////////////////////////////////////////	
	if(algorithm!=1 && option==1)
	{
		cout << "\n---->Linear Relaxation\n\n";
		LP_model(G,k);
	}
	////////////////////////////////////////////////////////////////////////	
	else if(algorithm==2 && option!=1)
	{
		cout << "\n---->Std Path - Representative Formulation\n\n";
		StdPath_TerminalModel(G,k);	
	}
	////////////////////////////////////////////////////////////////////////	
	else if(algorithm==22 && option!=1)
	{
		cout << "\n---->Std Path - Representative Formulation-with Long\n\n";
		StdPath_withLongPath(G,k);	
	}
	////////////////////////////////////////////////////////////////////////	
	else if(algorithm==4 && option!=1)
	{
		cout << "\n---->Bilevel Formulation\n\n";
		BilevelModel(G,k);
	}
	////////////////////////////////////////////////////////////////////////
	else if(algorithm==44 && option!=1)
	{
		cout << "\n---->Bilevel Formulation with Leaves\n\n";
		BilevelModel_withLeaf(G,k);
	}		
	////////////////////////////////////////////////////////////////////////
	else if(algorithm==7 && option!=1)
	{
		cout << "\n---->Mixed Smart Formulation\n\n";
		SmartModel(G,k);
	}
	/////////////////////////////////////////////////////////////////////////
	deleteGraphFF(G);
	deleteGraphFF(G_bar);
	delete [] preFixed;
	delete [] vertexWeight;
	delete [] istname;
	cout << "Done." << endl;

	if(output.is_open()){
		output << endl;
		output.close();
	}
	//cin.get();
	////////////////////////////////////////////////////////////////////////
	return 1;
}
