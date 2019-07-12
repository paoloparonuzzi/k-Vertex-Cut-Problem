#ifndef FUNCION_local_HEADER
#define FUNCION_local_HEADER


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
#include <iomanip> // for precision
#include <queue> 

#include "global_variables.h"
#include "Graph_v4.h"


using namespace std;

/***************************************************************************/
int    goldberg(int, int, int, int, int*, int*, double*, double*, double*, int*, short int *);
/***************************************************************************/

/***************************************************************************/
void SORT_NON_INCR_INT(int *item,int *score,int n);
/***************************************************************************/

/***************************************************************************/
void kp_load_cplex(int n_item, double C, double *weights);
/***************************************************************************/

/***************************************************************************/
double kp_solve_cplex(int n_item,double *profits, double *solution);
/***************************************************************************/

// A utility function to find the subset of an element i
/***************************************************************************/
int find(int i);
/***************************************************************************/
 
// A utility function to do union of two subsets 
/***************************************************************************/
void Union(int x, int y);
/***************************************************************************/

// Clique algorithm for separation of clique path constraint
/***************************************************************************/
double findGoodClique(int uIndex, int vIndex, graphFF G, IloNumArray zStar, IloIntArray goodClique);
/***************************************************************************/

// Clique algorithm for separation of clique path constraint (second disjoint clique)
/***************************************************************************/
double findDisjClique(int vIndex, graphFF G, IloNumArray zStar, IloIntArray goodClique, IloIntArray firstClique);
/***************************************************************************/

// A utility function to find the vertex with minimum distance value, from the set of vertices not yet included in shortest path tree
/***************************************************************************/
int minDistance(int V);
/***************************************************************************/

// Function that implements Dijkstra's single source shortest path algorithm for a graph represented using adjacency matrix representation
/***************************************************************************/
double dijkstra(double** graph, int V, int* parent, int src, int dest);
/***************************************************************************/

/***************************************************************************/
int findNextEdge(graphFF G, double** weight, int *head, int *tail);
/***************************************************************************/

/***************************************************************************/
double findLongPath(graphFF G, double** weight, int V, int firstEdge);
/***************************************************************************/

/***************************************************************************/
int maxStableSet(graphFF G, int v);
/***************************************************************************/

#endif
