

#include "global_variables.h"

//////////////////////////////////////////////////
int counter_clique;
int **clique;
int *sizeClique;

int disj_counter_clique;
int **disj_clique;
int *disj_sizeClique;

int *parent;
double *wweight;
int *iindex;
int *vindex;
double *vertexWeight;

int *statusDFS;
bool *visited;
double *distDijkstra;
int **rGraph;

bool *preFixed;

double no_path=9999.0;
int userCuts=0;
int lazyCuts=0;
int longPathCuts=0;

IloExprArray capConstrArray;
//////////////////////////////////////////////////


/////////////////////////
int algorithm;
int option;
int k;
int q=-1;
int b=-1;
double timeLimit;
double tolerance;
int userFreq;
int useWeight;
int seed;
char* istname;
/////////////////////////


/////////////////////////////CUTTING STOCK/////////////////////////////////////
int itemTypeNumber;
int itemNumber;
int capacity;
int *itemProfit,*itemWeight;
int *itemCardinality;
double BIN_COST;

///////////////////////////////////////
//BEST FIT
int binNumber;
int *vector_item_bin;
///////////////////////////////////////


/////////////////////////////////////CPLEX/////////////////////////////////////
CPXENVptr env_COMPACT;
CPXLPptr lp_COMPACT;

CPXENVptr env_MASTER;
CPXLPptr lp_MASTER;
CPXENVptr env_SLAVE;
CPXLPptr lp_SLAVE;
int status_BP_FormB_modl;
double time_total_FormB;

double BP_FormB_best_bd;
double BP_FormB_optValue;

CPXENVptr env_kp; 
CPXLPptr lp_kp; 

int status,ccnt,rcnt,nzcnt; 
int* rmatbeg, *rmatind,*cmatbeg, *cmatind; 
double* rmatval,*cmatval,*x,*pi,*obj, *lb, *ub,*rhs;
char *c_type,* sense;

IloNumArray userSolx;
IloNumArray userSolz;
IloNumArray lazySolx;
IloNumArray lazySolz;

IloIntArray goodCliqueU;
IloIntArray goodCliqueV;

IloIntArray longPath;
IloIntArray degreeTree;
IloNumArray treeCoeff;
///////////////////////////////
//OUTPUT
int _nodes;
int _var;
int _cons;
int  _stat_val;
double _cr_val;
double TIME_LIMIT;
double root_bound;
///////////////////////////////



/////////////////////////////////INFO/////////////////////////////
int nodecount;
int numcols;
int numrows;
int counter_infeasible_master;
int counter_together_branching;
int counter_separate_branching;
int counter_SeratorIN_branching;
int counter_SeratorOUT_branching;
int lpstat;

double objval;
double bestobjval;
double computation_time;
double BP_lp;
double BP_lp_time;
double BP_pricer_time;
double	BP_heur_time;
double	BP_opt_time;
bool STOP_TIME_LIMIT=false;
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
int edge_number;
int vertex_number;
int *fr;
int *to;
double *ca;
int source;
int sink;
double minmax_value;
int ncut;
double *flow;
short int *cut;

int edgeNumberNuc;
int nodeNumberNuc;
int *frNuc;
int *toNuc;
double *caNuc;
int sourceNuc;
int sinkNuc;
double minmax_valueNuc;
int ncutNuc;
double *flowNuc;
short int *cutNuc;
ofstream output;
//////////////////////////////////////////
