#include "LP_Model.h"
#include <algorithm>

ILOSTLBEGIN

void LP_model(graphFF G, int k)
{
	int i, j, predIndex;

//representative formulation
	visited=new bool[vertex_number];
	parent=new int[vertex_number];
	distDijkstra=new double[vertex_number];

	double **weight=new double*[vertex_number];
	double **dist=new double*[vertex_number];
	int **pred=new int*[vertex_number];

	for(i=0;i<vertex_number;i++){	
		dist[i]=new double[vertex_number];
		pred[i]=new int[vertex_number];
	}

	for(i=0;i<vertex_number;i++){
		weight[i]=new double[vertex_number];
	}

//bilevel formulation
	wweight = new double[edge_number];
	iindex = new int[edge_number];
	
//cplex
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{

		longPath = IloIntArray(env, vertex_number);
		degreeTree = IloIntArray(env, vertex_number);

		userSolx = IloNumArray(env, vertex_number);
		userSolz = IloNumArray(env, vertex_number);

		IloNumVarArray x(env, vertex_number, 0.0, 1.0);
		IloNumVarArray z(env, vertex_number, 0.0, 1.0);

		//Objective Function: minimize the number of vertices in the separator
		IloObjective obj(env, IloSum(x), IloObjective::Minimize);
		model.add(obj);

		//////////////////////////////////////////////////////////////////////
		//Vertex fixed in the separator thanks to preprocessing
		for(i=0; i<vertex_number; i++){
			if(preFixed[i])
				x[i].setLB(1);
		}
		//////////////////////////////////////////////////////////////////////			
		
		if(algorithm != 4 && algorithm != 44){  
			//Constraints: the number of the representative vertices must be k
			model.add(IloSum(z) == k);

			//Constraints: adjacent vertices can't be both representatives
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
		}

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::IloCplex::Param::Preprocessing::Presolve, 0);

		//Path cuts
		IloInt status;
		bool cutFound = true;
		IloExpr path(env);
		double weightZ=0;
		int last=-1;
		
		clock_t time1 = clock();
		
		clock_t time_LP1;
		clock_t time_LP2;
		double time_LP=0.0;
		
		int countCut=0;
		
		while(cutFound){
			time_LP1 = clock();
			cplex.solve();
			time_LP2 = clock();
			time_LP += (double) (time_LP2 - time_LP1)/CLOCKS_PER_SEC;
			cutFound = false;

//Bilevel separation cut (no Leaf)
			if(algorithm == 4){
				time_LP1 = clock();
				cplex.solve();
				time_LP2 = clock();
				time_LP += (double) (time_LP2 - time_LP1)/CLOCKS_PER_SEC;

				for(i=0; i<vertex_number; i++)
					userSolx[i]=cplex.getValue(x[i]);
		
				status = cplex.getStatus();
				if(status == 3){
					cout << "Infeasible" << endl;
					break;
				}

				for(i=0;i<edge_number;i++){
					wweight[i]=1-userSolx[G->H[i]]-userSolx[G->T[i]];
					iindex[i]=i;
				}
				SORT_NON_INCR(iindex,wweight,edge_number);	
	
				int src, dest, a, b, next, sizeTree;
				next=0;
				sizeTree=0;
				
				//initialize also degreeTree to -1 because is useful for later when I use IloScalProd with var x
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
					model.add(IloScalProd(degreeTree, x) >= rhsValue);
					cutFound=true;
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
						model.add(tree >= rhsValue);
						tree.end();
						cutFound=true;
					}
				}
			}
			
//Bilevel separation cut (with Leaf)
			if(algorithm == 44 || algorithm == 5 || algorithm == 6 || algorithm == 7){
				time_LP1 = clock();
				cplex.solve();
				time_LP2 = clock();
				time_LP += (double) (time_LP2 - time_LP1)/CLOCKS_PER_SEC;

				for(i=0; i<vertex_number; i++)
					userSolx[i]=cplex.getValue(x[i]);
		
				status = cplex.getStatus();
				if(status == 3){
					cout << "Infeasible" << endl;
					break;
				}

				for(i=0;i<edge_number;i++){
					wweight[i]=1-userSolx[G->H[i]]-userSolx[G->T[i]];
					iindex[i]=i;
				}
				SORT_NON_INCR(iindex,wweight,edge_number);	
	
				int src, dest, a, b, next, sizeTree;
				next=0;
				sizeTree=0;
				
				//initialize also degreeTree to -1 because is useful for later when I use IloScalProd with var x
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
					model.add(IloScalProd(degreeTree, x) >= rhsValue);
					cutFound=true;
					countCut++;
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
						model.add(tree >= rhsValue);
						tree.end();
						cutFound=true;
						countCut++;
					}
				}
			}			
			
//Long path cut
			if(algorithm == 22){		
				time_LP1 = clock();
				cplex.solve();
				time_LP2 = clock();
				time_LP += (double) (time_LP2 - time_LP1)/CLOCKS_PER_SEC;

				for(i=0; i<vertex_number; i++){
					userSolx[i]=cplex.getValue(x[i]);
					userSolz[i]=cplex.getValue(z[i]);
				}

		
				status = cplex.getStatus();
				if(status == 3){
					cout << "Infeasible" << endl;
					break;
				}

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
					cutFound=true;
					model.add(IloScalProd(longPath, x) -x[head] -x[tail] >= IloScalProd(longPath, z)-1);
				}
			}


//Standard path cut
			if(algorithm == 2 || algorithm == 22){
				time_LP1 = clock();
				cplex.solve();
				time_LP2 = clock();
				time_LP += (double) (time_LP2 - time_LP1)/CLOCKS_PER_SEC;
				
				for(i=0; i<vertex_number; i++){
					userSolx[i]=cplex.getValue(x[i]);
					userSolz[i]=cplex.getValue(z[i]);
				}
		
				status = cplex.getStatus();
				if(status == 3){
					cout << "Infeasible" << endl;
					break;
				}

				for(i=0; i<vertex_number; i++){
					for(j=0; j<vertex_number; j++){
						weight[i][j]=no_path;
					}
				}
	
				for(i=0; i<edge_number; i++){
					weight[G->T[i]][G->H[i]]=(userSolx[G->T[i]]+userSolx[G->H[i]])/2.0+0.000001;
					weight[G->H[i]][G->T[i]]=(userSolx[G->H[i]]+userSolx[G->T[i]])/2.0+0.000001;
				}

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
							//cout << path << endl;
							model.add(path >= 0);
							cutFound=true;
							path.clear();
						}
					}
				}
			}

			clock_t checkT = clock();
			double Tcheck = (double) (checkT - time1)/CLOCKS_PER_SEC;
			if(Tcheck > timeLimit)
				break;
		}

		time_LP1 = clock();
		cplex.solve();
		time_LP2 = clock();
		time_LP += (double) (time_LP2 - time_LP1)/CLOCKS_PER_SEC;

		clock_t time2 = clock();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;
		
		IloNum objective = cplex.getObjValue();
		cout << "LP obj: " << cplex.getObjValue() << endl;
		cout << "LP status: " << cplex.getStatus() << endl;
		cout << "number of found cut: " << countCut << endl;

		if(output.is_open())
			output << objective <<  "\t" << time << "\t";

		for(int i=0; i<vertex_number; i++)
			cout << "x" << i << ": " << cplex.getValue(x[i]) << "\t";
		cout << endl;
		
		obj.end();
		z.end();
		x.end();				

		longPath.end();
		degreeTree.end();

		userSolx.end();
		userSolz.end();
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
			output << "null" <<  "\t" << time << "\t";
		}

	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();

	for(i=0; i<vertex_number; i++){
		delete [] dist[i];
		delete [] pred[i];
	}
	delete [] dist;
	delete [] pred;
	
	for(i=0;i<vertex_number;i++){
		delete [] weight[i];
	}
	delete [] weight;

	delete [] visited;
	delete [] parent;
	delete [] distDijkstra;

	delete [] wweight;
	delete [] iindex;
}
