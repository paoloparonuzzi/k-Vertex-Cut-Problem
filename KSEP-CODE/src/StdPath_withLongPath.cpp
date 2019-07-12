#include "StdPath_withLongPath.h"
#include <algorithm>

ILOSTLBEGIN

ILOLAZYCONSTRAINTCALLBACK6(LazyLongPathCallback, graphFF, G, IloIntVarArray, x, IloIntVarArray, z, double**, weight, double**, dist, int**, pred)
{
	int predIndex, i, j;
	double weightZ;
	IloEnv env = getEnv();

	getValues(lazySolx,x);
	getValues(lazySolz,z);	
	
	for(i=0; i<vertex_number; i++){
		for(j=0; j<vertex_number; j++){
			weight[i][j]=no_path;
		}
	}

	for(i=0;i<edge_number;i++){
		weight[G->T[i]][G->H[i]]=(lazySolx[G->T[i]]+lazySolx[G->H[i]])/2.0+0.000001;
		weight[G->H[i]][G->T[i]]=(lazySolx[G->H[i]]+lazySolx[G->T[i]])/2.0+0.000001;
	}

	IloExpr path(env);

	Floyd_Warshall(vertex_number,weight,dist,pred);

	for(i=0; i<vertex_number; i++){
		for(j=i+1; j<vertex_number; j++){
			
			if(G->AMatrix[i][j] == 1 || G->AMatrix[j][i] == 1 || dist[i][j]>=no_path)
				continue;

			predIndex = pred[i][j];
			weightZ=0;
			while(predIndex != i){
				weightZ=weightZ+lazySolz[predIndex];
				predIndex = pred[i][predIndex];
			}
	
			if(dist[i][j]-lazySolx[i]/2-lazySolx[j]/2<lazySolz[i]+lazySolz[j]+weightZ-1-tolerance){
				path = -z[i] -z[j] +1;
				predIndex = pred[i][j];
				while(predIndex != i){
					path = path +x[predIndex] -z[predIndex];
					predIndex = pred[i][predIndex];
				}
				add(path >= 0).end();
				lazyCuts++;
				path.clear();
			}
		}
	}

	path.end();
	return;
}

ILOUSERCUTCALLBACK6(UserLongPathCallback, graphFF, G, IloIntVarArray, x, IloIntVarArray, z, double**, weight, double**, dist, int**, pred)
{
	IloInt nNodes = getNnodes();
	if(nNodes % userFreq != 0) return;
	
	IloEnv env = getEnv();
	int predIndex, i, j;
	double weightZ;		
	
	getValues(userSolx,x);
	getValues(userSolz,z);

//Long path cut
		
	for(i=0; i<vertex_number; i++){
		for(j=0; j<vertex_number; j++){
			weight[i][j]=-no_path;
		}
	}
	
	double maxWeight=-1;
	int firstEdge=-1;
	for(i=0; i<edge_number; i++){
		weight[G->T[i]][G->H[i]]=(userSolz[G->T[i]]+userSolz[G->H[i]])/2.0-(userSolx[G->T[i]]+userSolx[G->H[i]])/2.0;
		weight[G->H[i]][G->T[i]]=(userSolz[G->H[i]]+userSolz[G->T[i]])/2.0-(userSolx[G->H[i]]+userSolx[G->T[i]])/2.0;
		
		if(weight[G->H[i]][G->T[i]]>maxWeight){
			maxWeight=weight[G->H[i]][G->T[i]];
			firstEdge=i;
		}
	}

	//starting vertices
	int head = G->H[firstEdge];
	int tail = G->T[firstEdge];
	
	// Initialize all longPath[] as 0 but firstEdge
	for (i=0; i<vertex_number; i++)
		longPath[i]=0;
	longPath[head]=1;
	longPath[tail]=1;
	
	int edgeInserted=1;
	while(edgeInserted>0){
		edgeInserted = findNextEdge(G, weight, &head, &tail);
	}
					
	if(IloScalProd(longPath, userSolx) -userSolx[head] -userSolx[tail] < IloScalProd(longPath, userSolz) -1 -tolerance){
		add(IloScalProd(longPath, x) -x[head] -x[tail] >= IloScalProd(longPath, z) -1).end();
		longPathCuts++;
		return;		
	}

///////////////////////////////STD PATH////////////////////////////////////////////////
//weights for standard path separation	
	for(i=0; i<vertex_number; i++){
		for(j=0; j<vertex_number; j++){
			weight[i][j]=no_path;
		}
	}
	
	for(i=0;i<edge_number;i++){
		weight[G->T[i]][G->H[i]]=(userSolx[G->T[i]]+userSolx[G->H[i]])/2.0+0.000001;
		weight[G->H[i]][G->T[i]]=(userSolx[G->H[i]]+userSolx[G->T[i]])/2.0+0.000001;
	}

	IloExpr path(env);

	Floyd_Warshall(vertex_number,weight,dist,pred);

	for(i=0; i<vertex_number; i++){
		for(j=i+1; j<vertex_number; j++){
		
			if(G->AMatrix[i][j] == 1 || G->AMatrix[j][i] == 1 || dist[i][j]>=no_path)
				continue;

			predIndex = pred[i][j];
			weightZ=0;
			while(predIndex != i){
				weightZ=weightZ+userSolz[predIndex];
				predIndex = pred[i][predIndex];
			}
	
			if(dist[i][j]-userSolx[i]/2-userSolx[j]/2<userSolz[i]+userSolz[j]+weightZ-1-tolerance){
				path = -z[i] -z[j] +1;
				predIndex = pred[i][j];
				while(predIndex != i){
					path = path +x[predIndex] -z[predIndex];
					predIndex = pred[i][predIndex];
				}
				add(path >= 0).end();
				userCuts++;
				path.clear();
			}
		}
	}
	
	path.end();		
	return;
}

void StdPath_withLongPath(graphFF G, int k)
{
	int i, j;

	double **userWeight=new double*[vertex_number];
	double **userDist=new double*[vertex_number];
	int **userPred=new int*[vertex_number];

	double **lazyWeight=new double*[vertex_number];
	double **lazyDist=new double*[vertex_number];
	int **lazyPred=new int*[vertex_number];

	for(i=0;i<vertex_number;i++){	
		userDist[i]=new double[vertex_number];
		userWeight[i]=new double[vertex_number];
		userPred[i]=new int[vertex_number];
		
		lazyDist[i]=new double[vertex_number];
		lazyWeight[i]=new double[vertex_number];
		lazyPred[i]=new int[vertex_number];
	}

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		userSolx = IloNumArray(env, vertex_number);
		userSolz = IloNumArray(env, vertex_number);
		lazySolx = IloNumArray(env, vertex_number);
		lazySolz = IloNumArray(env, vertex_number);	

		longPath = IloIntArray(env, vertex_number);

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
		model.add(IloSum(z) >= k);

		//a verttex can not be interdicted and representative at the same time
		for(i=0; i<vertex_number; i++)
			model.add(z[i]+x[i] <= 1);	
		
		//Constraints: adjacent vertices can't be both representative
		for(i=0; i<edge_number; i++)
			model.add(z[G->H[i]] + z[G->T[i]] <= 1); 

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
		cut = cplex.use(UserLongPathCallback(env, G, x, z, userWeight, userDist, userPred));		
		IloCplex::Callback sec;	
		sec = cplex.use(LazyLongPathCallback(env, G, x, z, lazyWeight, lazyDist, lazyPred));
		
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::TreLim, 8000);
		
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
				<< nConst << "\t" << nNode << "\t" << userCuts+longPathCuts << "\t" << lazyCuts;
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
				<< nConst << "\t" << nNode << "\t" << userCuts+longPathCuts << "\t" << lazyCuts;
			}		
		}

		cut.end();
		sec.end();
		obj.end();
		z.end();
		x.end();		
		userSolx.end();
		userSolz.end();
		lazySolx.end();
		lazySolz.end();	
		longPath.end();
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
			output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode << "\t" << userCuts+longPathCuts << "\t" << lazyCuts;
		}

	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();

	for(i=0; i<vertex_number; i++){
		delete [] userDist[i];
		delete [] userWeight[i];
		delete [] userPred[i];
		delete [] lazyDist[i];
		delete [] lazyWeight[i];
		delete [] lazyPred[i];
	}
	delete [] userDist;
	delete [] userWeight;
	delete [] userPred;
	delete [] lazyDist;
	delete [] lazyWeight;
	delete [] lazyPred;		
}
