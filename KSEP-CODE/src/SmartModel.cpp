#include "SmartModel.h"
#include <algorithm>

ILOSTLBEGIN

ILOLAZYCONSTRAINTCALLBACK2(LazyMixedSmartCallback, graphFF, G, IloIntVarArray, x)
{
	int i;
	IloEnv env = getEnv();
	getValues(lazySolx,x);
	
/////////////////////////////////////////// Bilevel Cuts //////////////////////////////////////////////////////////////////
	for(i=0;i<edge_number;i++){
		wweight[i]=1-lazySolx[G->H[i]]-lazySolx[G->T[i]];
		iindex[i]=i;
	}

	int randomSeed = rand() % 10000 + 1;
	shuffle(iindex, edge_number, randomSeed);

	int src, dest, a, b, next, sizeTree, rhsValue;
	next=0;
	sizeTree=0;
	double lhsValue;
	IloExpr tree(env);

	for(i=0; i<vertex_number; i++){
		parent[i]=-1;
		degreeTree[i]=-1;
	}
	
	int current_edge;
	while(sizeTree<vertex_number-1 && next < edge_number)
	{		
		current_edge=iindex[next];
		if(wweight[current_edge]<0.00001){
			next++;
			continue;
		}
		
		src=G->T[current_edge];
		dest=G->H[current_edge];
		a = find(src);
		b = find(dest);

		if(a != b){
			Union(a, b);
			sizeTree++;
			degreeTree[src]++;
			degreeTree[dest]++;
		}
		next++;
	}
	/////////////////////////////Leaves/////////////////////////////////////////
	next=0;
	while(sizeTree<vertex_number-1 && next < edge_number)
	{		
		current_edge=iindex[next];
		if(wweight[current_edge]>0.1 || wweight[current_edge]<-0.1){
			next++;
			continue;
		}
		
		src=G->T[current_edge];
		dest=G->H[current_edge];
		a = find(src);
		b = find(dest);

		if(a != b){
			Union(a, b);
			sizeTree++;
			degreeTree[src]++;
			degreeTree[dest]++;
		}
		next++;
	}
	////////////////////////////////////////////////////////////////////////////////
		
	rhsValue=k-vertex_number+sizeTree;
	IloInt zeroDegree=0;
	if(rhsValue >= 0){	
		i=0;
		while(zeroDegree>-1 && i<vertex_number){
			zeroDegree=degreeTree[i];
			i++;
		}
	}
	
	lhsValue=IloScalProd(degreeTree, lazySolx);	
				
	//case in which downlifting is not possible
	if((rhsValue<0 || zeroDegree<0) && lhsValue<rhsValue-tolerance)
	{
		add(IloScalProd(degreeTree, x) >= rhsValue).end();
		lazyCuts++;
	}
	//case in which downlifting is possible (rhs is non-negative and all vertices are in the found subgraph)
	else if(rhsValue>=0 && zeroDegree>=0)
	{
		lhsValue=0;
		for(i=0; i<vertex_number; i++)
			lhsValue=lhsValue+min(rhsValue,(int) degreeTree[i])*lazySolx[i];
		if(lhsValue<rhsValue-tolerance){
			for(i=0; i<vertex_number; i++)
				tree=tree+min(rhsValue,(int) degreeTree[i])*x[i];
			add(tree >= rhsValue).end();
			tree.clear();
			lazyCuts++;
			return;
		}
	}
	tree.end();
	return;
}

ILOUSERCUTCALLBACK2(UserMixedSmartCallback, graphFF, G, IloIntVarArray, x)
{
	IloInt nNodes = getNnodes();
	if(nNodes % userFreq != 0) return;
	
	IloEnv env = getEnv();
	int i;	
	getValues(userSolx,x);
	
///////////////////////// Bilevel Cuts /////////////////////////////////////////

	for(i=0;i<edge_number;i++){
		wweight[i]=1-userSolx[G->H[i]]-userSolx[G->T[i]];
		iindex[i]=i;
	}
	SORT_NON_INCR(iindex,wweight,edge_number);	
	
	int src, dest, a, b, next, sizeTree;
	next=0;
	sizeTree=0;

	for(i=0; i<vertex_number; i++){
		parent[i]=-1;
		degreeTree[i]=-1;
	}
	
	while(wweight[next]>=0 && sizeTree<vertex_number-1 && next < edge_number)
	{
		src=G->T[iindex[next]];
		dest=G->H[iindex[next]];
		a = find(src);
		b = find(dest);

		if(a != b){
			Union(a, b);
			sizeTree++;
			degreeTree[src]++;
			degreeTree[dest]++;
		}
		next++;
	}
	
	int rhsValue=k-vertex_number+sizeTree;
	IloInt zeroDegree=0;
	if(rhsValue >= 0){	
		i=0;
		while(zeroDegree>-1 && i<vertex_number){
			zeroDegree=degreeTree[i];
			i++;
		}
	}
	
	double lhsValue=IloScalProd(degreeTree, userSolx);	
			
	//case in which downlifting is not possible
	if((rhsValue<0 || zeroDegree<0) && lhsValue<rhsValue-tolerance)
	{
		add(IloScalProd(degreeTree, x) >= rhsValue).end();
		userCuts++;
	}
	//case in which downlifting is possible
	else if(rhsValue>=0 && zeroDegree>=0)
	{
		lhsValue=0;
		for(i=0; i<vertex_number; i++)
			lhsValue=lhsValue+min(rhsValue,(int) degreeTree[i])*userSolx[i];
		if(lhsValue<rhsValue-tolerance){
			IloExpr tree(env);			
			for(i=0; i<vertex_number; i++)
				tree=tree+min(rhsValue,(int) degreeTree[i])*x[i];
			add(tree >= rhsValue).end();
			tree.end();
			userCuts++;
			return;
		}
	}		
	return;
}

void SmartModel(graphFF G, int k)
{
	int i, j;

//bilevel formulation
	parent=new int[vertex_number];
	wweight = new double[edge_number];
	iindex = new int[edge_number];
	
//cplex
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		degreeTree = IloIntArray(env, vertex_number);		

		userSolx = IloNumArray(env, vertex_number);
		lazySolx = IloNumArray(env, vertex_number);

		IloIntVarArray x(env, vertex_number, 0, 1);
		IloIntVarArray z(env, vertex_number, 0, 1);

		//Objective Function: minimize the number of (weighted) vertices in the separator
		IloNumArray iloVertexWeight(env, vertex_number);
		if(useWeight>0){
			for(i=0; i<vertex_number; i++)
				iloVertexWeight[i]=vertexWeight[i];			
		}
		else{
			for(i=0; i<vertex_number; i++)
				iloVertexWeight[i]=1.0;
		}
		IloObjective obj(env, IloScalProd(iloVertexWeight, x), IloObjective::Minimize);
		model.add(obj);
		
		//////////////////////////////////////////////////////////////////////
		//Vertex fixed in the separator thanks to preprocessing
		for(i=0; i<vertex_number; i++){
			if(preFixed[i])
				x[i].setLB(1);
		}
		//////////////////////////////////////////////////////////////////////	
		
		//Constraints: the number of the representative vertices must be k
		model.add(IloSum(z) == k);

		/*Constraints: adjacent vertices can't be both representatives*/
		for(i=0; i<edge_number; i++)
			model.add(z[G->H[i]] + z[G->T[i]] <= 1.0);

		//a vertex can not be interdicted and representative at the same time
		for(i=0; i<vertex_number; i++)
			model.add(z[i]+x[i] <= 1.0);	
		
		//Degree constraint
		IloExpr degreeConstr(env);
		for(i=0; i<vertex_number; i++){
			degreeConstr = z[i] - (G->DT[i]-1)*x[i];
			for(j=0; j<vertex_number; j++){
				if(G->AMatrix[i][j] == 1 || G->AMatrix[j][i] == 1)
					degreeConstr += z[j];
			}
			model.add(degreeConstr <= 1.0);
			degreeConstr.clear();
		}
		degreeConstr.end();

		IloCplex::Callback cut;
		cut = cplex.use(UserMixedSmartCallback(env, G, x));
		
		IloCplex::Callback sec;	
		sec = cplex.use(LazyMixedSmartCallback(env, G, x));

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		
		clock_t time1 = clock(); 
		cplex.solve();
		clock_t time2 = clock();

		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;
		
		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){		
			IloNum bound = cplex.getBestObjValue();			
			IloInt nNode = cplex.getNnodes();		
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();		
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode << "\t" << userCuts << "\t" << lazyCuts;
			}
		}			
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();		
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode << "\t" << userCuts << "\t" << lazyCuts;
			}		
		}
		
		cut.end();
		sec.end();
		obj.end();
		z.end();
		x.end();		
		userSolx.end();
		lazySolx.end();
		degreeTree.end();
		iloVertexWeight.end();	
	}

	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();

	delete [] parent;
	delete [] wweight;
	delete [] iindex;
}
