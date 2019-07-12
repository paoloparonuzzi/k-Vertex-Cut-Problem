#include "BilevelModel.h"

ILOLAZYCONSTRAINTCALLBACK2(LazySubgraphCallback, graphFF, G, IloIntVarArray, x)
{

	int i;
	IloEnv env = getEnv();
	getValues(lazySolx,x);
	
	for(i=0;i<edge_number;i++){
		wweight[i]=1-lazySolx[G->H[i]]-lazySolx[G->T[i]];
		iindex[i]=i;
	}

	int randomSeed = rand() % 10000 + 1;
	shuffle(iindex, edge_number, randomSeed);

	int src, dest, a, b, next, sizeTree, rhsValue;
	IloInt zeroDegree;
	next=0;
	sizeTree=0;
	IloExpr tree(env);

	for(i=0; i<vertex_number; i++){
		parent[i]=-1;
		degreeTree[i]=-1;
	}
	
	int current_edge;
	while(sizeTree<vertex_number-1 && next < edge_number)
	{		
		current_edge=iindex[next];
		if(wweight[current_edge]<0.1){
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
			//cout << "(" << src+1 << ", " << dest+1 << "),  ";
		}
		next++;
	}

	rhsValue=k-vertex_number+sizeTree;
	zeroDegree=0;
	if(rhsValue >= 0){	
		i=0;
		while(zeroDegree>-1 && i<vertex_number){
			zeroDegree=degreeTree[i];
			i++;
		}
	}

	double lhsValue=IloScalProd(degreeTree, lazySolx);	
		
	//case in which downlifting is not possible
	if((rhsValue<0 || zeroDegree<0) && lhsValue<rhsValue-tolerance)
	{
		add(IloScalProd(degreeTree, x) >= rhsValue).end();
		lazyCuts++;
	}
	//case in which downlifting is possible (rhs is non-negative and all vertices are in the found subgraph)
	else if(rhsValue>=0 && zeroDegree>=0)
	{
		//cout << "lift" << endl;
		lhsValue=0;
		for(i=0; i<vertex_number; i++)
			lhsValue=lhsValue+min(rhsValue,(int) degreeTree[i])*lazySolx[i];
		if(lhsValue<rhsValue-tolerance){			
			for(i=0; i<vertex_number; i++)
				tree=tree+min(rhsValue,(int) degreeTree[i])*x[i];
			add(tree >= rhsValue).end();
			tree.clear();
			lazyCuts++;
		}
	}
	tree.end();
	return;
}

ILOUSERCUTCALLBACK2(UserSubgraphCallback, graphFF, G, IloIntVarArray, x)
{
	IloInt nNodes = getNnodes();
	if(nNodes % userFreq != 0) return;

	int i;
	IloEnv env = getEnv();
	getValues(userSolx,x);
	
	for(i=0;i<edge_number;i++){
		wweight[i]=1-userSolx[G->H[i]]-userSolx[G->T[i]];
		iindex[i]=i;
	}
	SORT_NON_INCR(iindex,wweight,edge_number);	
	
	int src, dest, a, b, next, sizeTree;
	next=0;
	sizeTree=0;
	IloExpr tree(env);

	for(i=0; i<vertex_number; i++){
		parent[i]=-1;
		degreeTree[i]=-1;
	}
	
	while(wweight[next]>0 && sizeTree<vertex_number-1 && next < edge_number)
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
			//cout << "(" << src+1 << ", " << dest+1 << "),  ";
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
			for(i=0; i<vertex_number; i++)
				tree=tree+min(rhsValue,(int) degreeTree[i])*x[i];
			add(tree >= rhsValue).end();
			tree.clear();
			userCuts++;
		}
	}

	tree.end();
	return;
}

void BilevelModel(graphFF G, int k)
{
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	parent = new int[vertex_number];
	wweight = new double[edge_number];
	iindex = new int[edge_number];

	try{
		
		IloIntVarArray x(env, vertex_number, 0, 1);
		userSolx = IloNumArray(env, vertex_number);
		lazySolx = IloNumArray(env, vertex_number);
		degreeTree = IloIntArray(env, vertex_number);

		//Objective Function: minimize the number of (weighted) vertices in the separator
		IloNumArray iloVertexWeight(env, vertex_number);
		if(useWeight>0){
			for(int i=0; i<vertex_number; i++)
				iloVertexWeight[i]=vertexWeight[i];			
		}
		else{
			for(int i=0; i<vertex_number; i++)
				iloVertexWeight[i]=1.0;
		}
		IloObjective obj(env, IloScalProd(iloVertexWeight, x), IloObjective::Minimize);
		model.add(obj);
		
		//////////////////////////////////////////////////////////////////////
		//Vertex fixed in the separator thanks to preprocessing
		for(int i=0; i<vertex_number; i++){
			if(preFixed[i])
				x[i].setLB(1);
		}
		//////////////////////////////////////////////////////////////////////	
						
		IloCplex::Callback sec;
		sec = cplex.use(LazySubgraphCallback(env, G, x));
		
		IloCplex::Callback cut;
		cut = cplex.use(UserSubgraphCallback(env, G, x));

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		
		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
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
		
		//cout << cplex.getStatus() << endl;
		//cplex.exportModel("model.lp");

		sec.end();
		cut.end();
		obj.end();
		x.end();
		userSolx.end();
		lazySolx.end();
		degreeTree.end();
		iloVertexWeight.end();	
		
	}

	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;

		IloInt nNode = cplex.getNnodes();
		IloNum bound = cplex.getBestObjValue();
		IloInt status = cplex.getStatus();
		IloInt nVar = cplex.getNcols();
		IloInt nConst = cplex.getNrows();

		if(output.is_open())
		{
			output.precision(10);
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t" 
				<< time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode << "\t" << userCuts << "\t" << lazyCuts;
		}
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
