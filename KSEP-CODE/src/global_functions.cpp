

#include "global_functions.h"


/***************************************************************************/
void SORT_NON_INCR_INT(int *item,int *score,int n)
/***************************************************************************/
{
	int salto,i,j,tempItem;
	int tempScore;

	for (salto = n/2; salto > 0; salto /=2)
		for (i = salto; i < n; i++)
			for (j = i-salto; j >= 0; j-=salto)
			{
				if (score[j] >= score[j+salto]) break;
				tempScore = score[j];
				score[j]=score[j+salto];
				score[j+salto]=tempScore;
				tempItem = item[j];
				item[j]=item[j+salto];
				item[j+salto]=tempItem;
			}
}




/***************************************************************************/
void kp_load_cplex(int n_item, double C, double *weights)
/***************************************************************************/
{

	int i;



	env_kp=CPXopenCPLEX(&status);
	if(status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	lp_kp=CPXcreateprob(env_kp,&status,"KP");
	if(status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}


	CPXchgobjsen(env_kp,lp_kp,CPX_MAX);


	//variables
	ccnt=n_item;
	obj=(double*) calloc(ccnt,sizeof(double));
	lb=(double*) calloc(ccnt,sizeof(double));
	ub=(double*) calloc(ccnt,sizeof(double));
	char *_ctype=(char*) calloc(ccnt,sizeof(char));

	for(i=0; i<ccnt; i++)
	{
		obj[i]=(double)0.0;
		lb[i]=0.0;
		ub[i]=1.0;
		_ctype[i]='B';
	}

	status=CPXnewcols(env_kp,lp_kp,ccnt,obj,lb,ub,_ctype,NULL);
	if(status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	//constraints
	rcnt=1;
	nzcnt=n_item;
	rhs=(double*) calloc(rcnt,sizeof(double));
	sense=(char*) calloc(rcnt,sizeof(double));

	rhs[0]=(double)C;		
	sense[0]='L';


	rmatbeg=(int*) calloc(rcnt,sizeof(int));
	rmatind=(int*) calloc(nzcnt,sizeof(int));
	rmatval=(double*) calloc(nzcnt,sizeof(double));



	for(i=0; i<n_item; i++)
	{
		rmatval[i]=(double)weights[i];
		rmatind[i]=i;
	}

	rmatbeg[0]=0;

	status=CPXaddrows(env_kp,lp_kp,0,rcnt,nzcnt,rhs,sense,rmatbeg,rmatind,rmatval,NULL,NULL);
	if(status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(rmatbeg);
	free(rmatval);
	free(rmatind);
	free(rhs);
	free(sense);


	// writing the problem in a .lp format for control
	//	status=CPXwriteprob(env_kp,lp_kp,"KP.lp",NULL);
	//	if(status!=0) {printf("error in CPXwriteprob\n");	exit(-1);}
	//	cin.get();

	free(obj);
	free(lb);
	free(ub);
	free(_ctype);

}


/***************************************************************************/
double kp_solve_cplex(int n_item,double *profits, double *solution)
/***************************************************************************/
{

	int i;

#ifdef CPLEX_OUTPUT
	CPXsetintparam (env_kp, CPX_PARAM_SCRIND, CPX_ON);
#endif

	// * Set relative tolerance *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPAGAP, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPAGAP\n");
	}

	// * Set relative tolerance *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPGAP, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPGAP\n");
	}


	// * Set mip tolerances integrality *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPINT, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPINTP\n");
	}

	// * Set mip tolerances integrality *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPINT, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPINTP\n");
	}

	// * Set Feasibility tolerance *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPRHS, 1e-9);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPRHS\n");
	}

	// * Emphasizes precision in numerically unstable or difficult problems *
	//status = CPXsetdblparam (env_kp, CPX_PARAM_NUMERICALEMPHASIS, 1);
	//if (status)
	//{
	//	fprintf (stderr, "error for CPX_PARAM_NUMERICALEMPHASIS\n");
	//}


	int *ind = (int*) malloc(sizeof(int) * n_item);
	double *d = (double*) malloc(sizeof(double)* n_item);
	for (i = 0; i < n_item; i++)
	{
		d[i] = profits[i];
		ind[i] = i;
	}

	status = CPXchgobj(env_kp, lp_kp,n_item, ind, d);
	if (status != 0)
	{
		printf("error in CPXchgobj\n");
		exit(-1);
	}

	free(ind);
	free(d);



#ifdef writer_KP
	status=CPXwriteprob(env_kp,lp_kp,"kp.lp",NULL);
	if(status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);

	}
#endif



	status=CPXmipopt(env_kp,lp_kp);
	if(status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}


	//getting the solution
	double *_x=(double*) calloc(n_item,sizeof(double));
	double _objval_p=0;

	status=CPXgetmipx(env_kp,lp_kp,_x,0,n_item-1);
	if(status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}
	status=CPXgetmipobjval(env_kp,lp_kp,&_objval_p);
	if(status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

#ifdef print_ist_sol	
	cout << "\n\nOBJ_KP01 ->\t " << _objval_p << endl;
#endif

	for (i = 0; i < n_item; i++)
	{
		solution[i]=(int)(_x[i]+0.5);
	}

	free(_x);

#ifdef print_ist_sol
	cout << "\n\nSolution\n";
	for (i = 0; i < n_item; i++)
	{
		cout << solution[i] << "";
	}
	cout << endl;
	cout << endl;
#endif


	return _objval_p;

}

// A utility function to find the subset of an element i
/***************************************************************************/
int find(int i)
/***************************************************************************/
{
    if (parent[i] == -1)
        return i;
    return find(parent[i]);
}
 
// A utility function to do union of two subsets
/***************************************************************************/ 
void Union(int x, int y)
/***************************************************************************/
{
    int xset = find(x);
    int yset = find(y);
    parent[xset] = yset;
}

// A utility function to find the vertex with minimum distance
// value, from the set of vertices not yet included in shortest
// path tree
/***************************************************************************/
int minDistance(int V)
/***************************************************************************/
{
    // Initialize min value
	double min = no_path*vertex_number;
	int min_index;
 
    for (int v = 0; v < V; v++){
        if (visited[v] == false && distDijkstra[v] <= min){
            min = distDijkstra[v]; 
			min_index = v;
		}
	}
    return min_index;
}

// Function that implements Dijkstra's single source shortest path
// algorithm for a graph represented using adjacency matrix
// representation
/***************************************************************************/
double dijkstra(double** graph, int V, int* parent, int src, int dest)
/***************************************************************************/
{

	// visited[i] will true if vertex i is included / in shortest
	// path tree or shortest distance from src to i is finalized
	// Parent array to store shortest path tree
    
	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < V; i++)
	{
		parent[i] = -1;
		distDijkstra[i] = no_path;
		visited[i] = false;
	}
 
	// Distance of source vertex from itself is always 0
	distDijkstra[src] = 0.0;
 
	// Find shortest path for all vertices
	for (int count = 0; count < V-1; count++)
	{
		// Pick the minimum distance vertex from the set of
		// vertices not yet processed. u is always equal to src
		// in first iteration.
		int u = minDistance(V);
 
		// Mark the picked vertex as processed
		visited[u] = true;
 
		// Update distance value of the adjacent vertices of the
		// picked vertex.
		for (int v = 0; v < V; v++)
 
			// Update distance[v] only if is not in visited, there is
			// an edge from u to v, and total weight of path from
			// src to v through u is smaller than current value of
			// distance[v]
			if (!visited[v] && graph[u][v] < no_path &&
				distDijkstra[u] + graph[u][v] < distDijkstra[v])
			{
				parent[v]  = u;
				distDijkstra[v] = distDijkstra[u] + graph[u][v];
			}  
	}
	return distDijkstra[dest];
}

/***************************************************************************/
int findNextEdge(graphFF G, double** weight, int *head, int *tail)
/***************************************************************************/
{
	double maxWeight=-1;
	int next=-1;
	int which=0;

	for(int i=0; i<vertex_number; i++){

		if( (G->AMatrix[*head][i]==1 || G->AMatrix[i][*head]==1) && i != *tail){
			if(weight[*head][i] > maxWeight && longPath[i]==0){
				next=i;
				maxWeight=weight[*head][i];
				which=1;
			}
		}
		
		if( (G->AMatrix[*tail][i]==1 || G->AMatrix[i][*tail]==1) && i != *head){
			if(weight[*tail][i] > maxWeight && longPath[i]==0){
				next=i;
				maxWeight=weight[*tail][i];
				which=2;
			}
		}		
	}

	if(which > 0){
		longPath[next]=1;
		if(which==1)
			*head=next;
		if(which==2)
			*tail=next;
	}

	return which;
}

/***************************************************************************/
double findLongPath(graphFF G, double** weight, int V, int firstEdge)
/***************************************************************************/
{
	double longPathWeight=0;

	// Initialize all stpSet[] as false but firstEdge
	for (int i = 0; i < V; i++)
		longPath[i] = 0;
	longPath[firstEdge]=1;

	//starting vertices
	int head = G->H[firstEdge];
	int tail = G->T[firstEdge];

	int edgeInserted=1;

	while(edgeInserted>0){
		edgeInserted = findNextEdge(G, weight, &head, &tail);
	}

	return longPathWeight;
}

/***************************************************************************/
int maxStableSet(graphFF G, int v)
/***************************************************************************/
{
	IloEnv env;
	int output;
	try{
		IloModel model(env);
		IloIntVarArray x(env, vertex_number, 0, 1);

		//Objective Function: maximize the number of vertices in the clique
		IloObjective obj(env, IloSum(x), IloObjective::Maximize);
		model.add(obj);
		
		//v, his neighborhood and all fixed vertex can not belong to the stable set		
		for(int i=0; i<vertex_number; i++){
			if(i==v || G->AMatrix[i][v] == 1 || G->AMatrix[v][i] == 1 || preFixed[i])
				x[i].setUB(0);
		}
		
		//Constraints: for each edge (i, j), at most one vertex can be inside
		for(int i=0; i<edge_number; i++){
			int head = G->H[i];
			int tail = G->T[i];
			model.add(x[head]+x[tail]<=1);
		}

		//Optimize
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.solve();

		output = cplex.getObjValue()+0.5;
			
		cplex.end();
		x.end();
		model.end();
				
	}

	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	env.end();
	return output;
}

