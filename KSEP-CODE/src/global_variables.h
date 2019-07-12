#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <string>
#include <iomanip> // for precision
#include "ilcplex/ilocplex.h"

using namespace std;

//////////////////////////////////////////////////
extern int counter_clique;
extern int **clique;
extern int *sizeClique;

extern int disj_counter_clique;
extern int **disj_clique;
extern int *disj_sizeClique;

extern int *parent;
extern double *wweight;
extern int *iindex;
extern int *vindex;
extern int *statusDFS;
extern double *vertexWeight;

extern bool *visited;
extern double *distDijkstra;
extern int **rGraph;

extern bool *preFixed;

extern double no_path;
extern int userCuts;
extern int lazyCuts;
extern int longPathCuts;

extern IloExprArray capConstrArray;
//////////////////////////////////////////////////


/////////////////////////
extern int algorithm;
extern int option;
extern int k;
extern int q;
extern int b;
extern double timeLimit;
extern double tolerance;
extern int userFreq;
extern int useWeight;
extern int seed;
extern char* istname;
/////////////////////////

/////////////////////////////CUTTING STOCK/////////////////////////////////////
extern int itemTypeNumber;
extern int itemNumber;
extern int capacity;
extern int *itemProfit,*itemWeight;
extern int *itemCardinality;
extern double BIN_COST;

///////////////////////////////////////
//BEST FIT
extern int binNumber;
extern int *vector_item_bin;
///////////////////////////////////////


/////////////////////////////////////CPLEX/////////////////////////////////////

extern CPXENVptr env_COMPACT;
extern CPXLPptr lp_COMPACT;


extern CPXENVptr env_MASTER;
extern CPXLPptr lp_MASTER;
extern CPXENVptr env_SLAVE;
extern CPXLPptr lp_SLAVE;

extern int status_BP_FormB_model;
extern double time_total_FormB;


extern double BP_FormB_best_bd;
extern double BP_FormB_optValue;

extern CPXENVptr env_kp; 
extern CPXLPptr lp_kp; 


extern int status,ccnt,rcnt,nzcnt;
extern int* rmatbeg, *rmatind,*cmatbeg, *cmatind;
extern double* rmatval,*cmatval,*x,*pi,*obj, *lb, *ub,*rhs;
extern char *c_type,* sense;

extern IloNumArray userSolx;
extern IloNumArray userSolz;
extern IloNumArray lazySolx;
extern IloNumArray lazySolz;

extern IloIntArray goodCliqueU;
extern IloIntArray goodCliqueV;

extern IloIntArray longPath;
extern IloIntArray degreeTree;
extern IloNumArray treeCoeff;
///////////////////////////////
//OUTPUT
extern int _nodes;
extern int _var;
extern int _cons;
extern int  _stat_val;
extern double _cr_val;
extern double TIME_LIMIT;
extern double root_bound;
///////////////////////////////




/////////////////////////////////INFO/////////////////////////////
extern int nodecount;
extern int numcols;
extern int numrows;
extern int counter_infeasible_master;
extern int counter_together_branching;
extern int counter_separate_branching;
extern int counter_SeratorIN_branching;
extern int counter_SeratorOUT_branching;
extern int lpstat;

extern double objval;
extern double bestobjval;
extern double computation_time;
extern double BP_lp;
extern double BP_lp_time;
extern double BP_pricer_time;
extern double	BP_heur_time;
extern double	BP_opt_time;
extern bool STOP_TIME_LIMIT;
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
extern int edge_number;
extern int vertex_number;
extern int *fr;
extern int *to;
extern double *ca;
extern int source;
extern int sink;
extern double minmax_value;
extern int ncut;
extern double *flow;
extern short int *cut;

extern int edgeNumberNuc;
extern int nodeNumberNuc;
extern int *frNuc;
extern int *toNuc;
extern double *caNuc;
extern int sourceNuc;
extern int sinkNuc;
extern double minmax_valueNuc;
extern int ncutNuc;
extern double *flowNuc;
extern short int *cutNuc;
extern ofstream output;
//////////////////////////////////////////

#endif
