#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "epanet2.h" 
#include "gurobi_c.h"

//September 06, 2013
//L1-Approximation LP for leak detection in a water distribution network 
//using EPANET, and Gurobi Optimizer
//
//

//Changing delta, emitterCoeff, iterations and numOfLeaks manipulates the 
//full program
double delta = 1, emitterCoeff = 1, binaryLeakLimit = 2.0;
int numOfLeaks = 2, iterations = 50;
char inputFile[50] = "hanoi-1.inp";
char reportFile[50] = "hanoi.rpt";
char directoryString[50] = "L1_MIP/Current";

int totalNodeCount;
int *leakNodes;
double totalDemand, bigM = 999999.99;
double *baseCasePressureMatrix, *observedPressure, *coefficients, *b, *bhat,
	*realLeakValues, *singleRunErrors, *leakDemands, *leakMagnitudes, 
	*modelError, **largePressureMatrix, **largeA, **Ahat,  **I, 
	*objectiveValues; 
	
FILE *ptr_file;

void initializeArrays();
void populateMatricies(int);
void randomizeLeaks(double, int, int);
void printLeakInfo(int);
void analyzeBaseCase(int);
void oneLeak(int, double, int, int);
void nLeaks(int, int);
double calculateError(int, double[]);
int writeSummaryFile(int, int, double, double[]);
int writeRawResults(int, int, double[]);
int writeLeakFile(int);
int writeErrorFile();

int main(int argc, char *argv[]) 
{
	GRBenv *env = NULL;
	GRBmodel *model = NULL;	
	int  i, j, k, numNodes, storage;
	double errorSum;
	
	//Randomize the leak locations, commented out will use the same seeding 
	//for each run
	//srand(time(NULL));
	
	i = j = k = numNodes = 0;
	errorSum = 0.0;
	
	//Open EPANET & Input file
	ENopen(inputFile,reportFile,"");
	
	// Get the number of nodes
	ENgetcount(EN_NODECOUNT, &numNodes);
	ENgetcount(EN_TANKCOUNT, &storage);
	totalNodeCount = numNodes - storage;
	
	int       error = 0;
	double    sol[(totalNodeCount * 3)];
	int       ind[(totalNodeCount * 3)];
	double    val[(totalNodeCount * 3)];
	double    obj[(totalNodeCount * 3)];
	char      vtype[(totalNodeCount * 3)];	
	int       optimstatus;
	double    objval;

	
	baseCasePressureMatrix = (double *) calloc(totalNodeCount, sizeof(double));
	observedPressure = (double *) calloc(totalNodeCount, sizeof(double));
	coefficients = (double *) calloc((totalNodeCount * 2), sizeof(double));
	realLeakValues = (double *) calloc(totalNodeCount, sizeof(double));
	singleRunErrors = (double *) calloc(totalNodeCount, sizeof(double));
	b = (double *) calloc(totalNodeCount, sizeof(double));
	bhat = (double *) calloc((totalNodeCount * 2), sizeof(double));
	leakDemands = (double *) calloc(numOfLeaks, sizeof(double));
	leakNodes = (int *) calloc(numOfLeaks,sizeof(int));
	leakMagnitudes = (double *) calloc(numOfLeaks,sizeof(double));
	modelError = (double *) calloc(iterations, sizeof(double));
	objectiveValues = (double *) calloc(iterations, sizeof(double));
	
	largePressureMatrix = (double **) malloc(totalNodeCount * sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
		largePressureMatrix[i] = malloc(totalNodeCount * sizeof(double));
	
	largeA = (double **) malloc(totalNodeCount * sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
		largeA[i] = malloc(totalNodeCount * sizeof(double));
	
	I = (double **) malloc(totalNodeCount * sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
		I[i] = malloc(totalNodeCount * sizeof(double));
	
	Ahat = (double **) malloc( (totalNodeCount * 2) * sizeof(double *) );
	for(i = 0; i < (totalNodeCount * 2); i++)
		Ahat[i] = malloc( (totalNodeCount * 2) * sizeof(double) );
	
	/* Create environment */
 	error = GRBloadenv(&env, "L1_MIP.log");
 	if (error) goto QUIT;
		
	//Create observation	
	for (k = 0; k < iterations; k++)
	{		
		initializeArrays();
		
		randomizeLeaks(emitterCoeff, totalNodeCount, numOfLeaks);
		
		analyzeBaseCase(totalNodeCount);
		
		nLeaks(numOfLeaks, totalNodeCount);
		
		printLeakInfo(numOfLeaks);
		
		populateMatricies(totalNodeCount);
		
		// Create an empty model 		
 		error = GRBnewmodel(env, &model, "L1MIP", 0, NULL, NULL, NULL, NULL, 
 			NULL);
 		if (error) goto QUIT;
 		 	
 		// Add variables 
 		for (i = 0; i < (totalNodeCount * 2); i++)
 		{
 			obj[i] = coefficients[i]; 			
 			vtype[i] = GRB_CONTINUOUS; 			
 		}
 		
 		for (i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
 		{
 			obj[i] = 0.0;
 			vtype[i] = GRB_BINARY;
 		}
 		 				
		error = GRBaddvars(model, (totalNodeCount * 3), 0, NULL, NULL, NULL, obj,
			NULL, NULL, vtype, NULL);
		if (error) goto QUIT;
		
		// Integrate new variables		
		error = GRBupdatemodel(model);
		if (error) goto QUIT;
				
		// First constraint: Ax <= b						
		for (i = 0; i < (totalNodeCount * 2); i++)
		{
			for (j = 0; j < (totalNodeCount * 2); j++)
			{
				ind[j] = j;
				val[j] = Ahat[i][j];			
			}								
			error = GRBaddconstr(model, (totalNodeCount * 2), ind, val, 
				GRB_LESS_EQUAL, bhat[i],NULL);			
			if (error) goto QUIT;
		}
		
		//Leak magnitude - (binary * bigM) <= 0
		for (i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
		{		
			ind[0] = (i - (totalNodeCount * 2)); 	ind[1] = i; 
			val[0] = 1.0; 		val[1] = -bigM ;
									
			error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL,0.0,NULL);
			if (error) goto QUIT;
		}
		
		// Limit sum of binaries to number of leaks searching for...		
		for (i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
		{		
			ind[i-(totalNodeCount * 2)] = i;
			val[i-(totalNodeCount * 2)] = 1.0;
		}								
		error = GRBaddconstr(model, totalNodeCount, ind, val, GRB_LESS_EQUAL,
			binaryLeakLimit,NULL);
		if (error) goto QUIT;	

		error = GRBoptimize(model);
		if (error) goto QUIT;
		
		// Write model to 'L1Approx.lp'		
		error = GRBwrite(model, "L1_MIP.lp");
		if (error) goto QUIT;
		
		error = GRBwrite(model, "L1_MIP.sol");
		if (error) goto QUIT;
		
		// Capture solution information		
		error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
		if (error) goto QUIT;
		
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
		if (error) goto QUIT;
		
		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
			(totalNodeCount * 3), sol);
		if (error) goto QUIT;
		
		printf("\nOptimization complete\n");
		if (optimstatus == GRB_OPTIMAL)
		{
			printf("Optimal objective: %.4e\n", objval);
			//for(i = 0; i < (totalNodeCount * 2); i++)
			//{		  	
				//printf("  sol[%d] = %f \n", (i+1), sol[i]);
			//}
			//for(i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
			//{		  	
				//printf("  sol[%d] = %1.0f \t binary for sol[%d] \n", (i+1), 
				//	sol[i], (i-(totalNodeCount * 2) + 1));
			//}
		} else if (optimstatus == GRB_INF_OR_UNBD) 
		{
			printf("Model is infeasible or unbounded\n");
		} else 
		{
		  printf("Optimization was stopped early\n");
		}
		
		objectiveValues[k] = objval;
		modelError[k] = calculateError(totalNodeCount, sol);		
		
		writeSummaryFile(k, optimstatus, objval, sol);
		writeRawResults(k, optimstatus, sol);
		writeLeakFile(k);
		
		/* Free model */
		GRBfreemodel(model);		
	}
	
	ENclose();
	
	//writeErrorFile();
	
	free(leakNodes);
	free(leakMagnitudes);
	free(leakDemands);
	free(modelError);
	free(objectiveValues);
	free(baseCasePressureMatrix);
	free(observedPressure);
	free(coefficients);
	free(b);
	free(bhat);
	free(realLeakValues);
	free(singleRunErrors);
	
	for(i = 0; i < totalNodeCount; i++)
		free((void *)largePressureMatrix[i]);
	free((void *)largePressureMatrix);
	
	for(i = 0; i < totalNodeCount; i++)
		free((void *)largeA[i]);
	free((void *)largeA);
	
	for(i = 0; i < totalNodeCount; i++)
		free((void *)I[i]);
	free((void *)I);
	
	for(i = 0; i < (totalNodeCount * 2); i++)
		free((void *)Ahat[i]);
	free((void *)Ahat);
	
	
		QUIT:

		/* Error reporting */
		
		if (error) {
		  printf("ERROR: %s\n", GRBgeterrormsg(env));
		  exit(1);
		}
		
		/* Free model */
		
		GRBfreemodel(model);
		
		/* Free environment */
		
		GRBfreeenv(env);
			
	return 0;
}

//FUNCTION
//Initialze various arrays to be populated during simulation
void initializeArrays()
{
	int i, j;	
	i = j = 0;
	
	//Array initialization	
	for (i = 0; i < totalNodeCount; i++)
	{
		observedPressure[i] = 0;
		baseCasePressureMatrix[i] = 0;
		b[i] = 0;
		realLeakValues[i] = 0.0;
		singleRunErrors[i] = 0.0;	
	}
	
	for (i = 0; i < (totalNodeCount * 2); i++)
	{
		bhat[i] = 0;		
	}
	
	for (i = 0; i < (totalNodeCount * 2); i++)
	{
		coefficients[i] = 0;		
	}
	
	for (i = 0; i < totalNodeCount; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{
			largePressureMatrix[i][j] = 0;
			largeA[i][j] = 0;		
		}
	}
	
	for (i = 0; i < (totalNodeCount * 2); i++)
	{
		for (j = 0; j < (totalNodeCount * 2); j++)
		{
			Ahat[i][j] = 0;		
		}
	}
	
	for (i = 0; i < totalNodeCount; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{
			I[i][j] = 0;
		}
	}	

	//Create Identity Matrix
	for(i = 0; i < totalNodeCount; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{
			if (i==j)
				I[i][j] = 1;
		}
	}
	
	//Create c-transpose
	for (i = 0; i < totalNodeCount; i++)
	{
		coefficients[i] = 0.0;
	}	
	for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
	{
		coefficients[i] = 1.0;
	}
	
}

//FUNCTION
//Populate array values for the L1 Approximation
//Also calls single leak simulations for each node in the network
void populateMatricies(int numNodes)
{
	int i, j, temp;
	
	i = j = 0;
	
	//Update b matrix
	for (i = 0; i < numNodes; i++)
	{
		b[i] = (baseCasePressureMatrix[i] - observedPressure[i]);	
		//printf("b[%d] = %f\n",i,b[i]);
	}
	
	//Create b-hat
	for (i = 0; i < numNodes; i++)
	{
		bhat[i] = b[i];
	}
	for (i = numNodes; i < (numNodes * 2); i++)
	{
		bhat[i] = -b[i-numNodes];
	}
	
	
	for(i = 1; i <= numNodes; i++)
	{		
		oneLeak(i, delta, numNodes, i-1);		
	}
	
	//Update A matrix		
	for(i = 0; i < numNodes; i++)
	{		
		for(j = 0; j < numNodes; j++)
		{
			largeA[i][j] = (baseCasePressureMatrix[i] - 
				largePressureMatrix[i][j]) / delta;			
		}			
	}	
	
	//Create A-hat
	for(i = 0; i < numNodes; i++)
	{
		for(j = 0; j < numNodes; j++)
		{
			Ahat[i][j] =  largeA[i][j];
		}
	}
	for(i = numNodes; i < (numNodes * 2); i++)
	{
		for(j = 0; j < numNodes; j++)
		{
			Ahat[i][j] = -largeA[i-numNodes][j];
		}
	}
	for(i = 0; i < numNodes; i++)
	{
		for(j = numNodes; j < (numNodes * 2); j++)
		{
			Ahat[i][j] = -I[i][j-numNodes];
		}
	}
	for(i = numNodes; i < (numNodes * 2); i++)
	{
		for(j = numNodes; j < (numNodes * 2); j++)
		{
			Ahat[i][j] = -I[i-numNodes][j-numNodes];
		}
	}
	
	//Keep track of the emitter coefficient at every network node (most should
	//	be zero)
	for (i = 0; i < numOfLeaks; i++)
	{
		temp = (leakNodes[i]-1);		
		realLeakValues[temp] = leakMagnitudes[i];
	}
}

void randomizeLeaks(double emitterCoeff, int numNodes, int numOfLeaks)
{	
	int i, j;
	
	i = j = 0;
	
	for (i = 0; i < numOfLeaks; i++)
		{
			leakDemands[i] = 0;
			leakNodes[i] = 0;
			leakMagnitudes[i] = 0.0;
		}

		for (i = 0; i < numOfLeaks; i++)
		{			
			leakNodes[i] = (int)(rand()%numNodes)+1;
			/*
			if (i > 0)
			{	
				for (j = 0; j < i; j++)
				{
					if (leakNodes[i] == leakNodes[j])
					{
						do
						{						
							leakNodes[i] = (int)(rand()%numNodes)+1;
						}while (leakNodes[i] == leakNodes[j]);
					}
				}
			}
			*/
			/*
			if (leakNodes[i] == 32)
			{
				do
				{			
					leakNodes[i] = (int)(rand()%numNodes);//+1;
				}while(leakNodes[i] == 32);
			}
			*/
			
			j = i;
			
			while(j >= 1)
			{
				if (leakNodes[i] == leakNodes[j-1])
				{
					do
					{			
						leakNodes[i] = (int)(rand()%numNodes);//+1;
					}while(leakNodes[i] == leakNodes[j-1]);										
				}
				j--;
			}
		}			
		
		for(i = 0; i < numOfLeaks; i++)
		{
			leakMagnitudes[i] = drand48() * 10;
		}		

}


//FUNCTION
//Print the location and magnitude of leaks
void printLeakInfo(int numOfLeaks)
{
	int i;
	char name[10];
	i = 0;
	
	for(i = 0; i < numOfLeaks; i++)
	{			
		ENgetnodeid(leakNodes[i], name);
		printf("\n leak info Node: %d: %d \t Node ID: %s \t Magnitude: %f \n",
			i, leakNodes[i], name, leakMagnitudes[i]);
	}
}

//FUNCTION
//Runs the hydraulic analysis on the base case scenario
void analyzeBaseCase(int nodeCount)
{		
	long t, tstep, hydraulicTimeStep, duration;
	float pressure;
	int i;	
	//char name[20];
	
	i = 0;
	pressure = 0.0;
	
	ENgettimeparam( EN_HYDSTEP, &hydraulicTimeStep );
	ENgettimeparam( EN_DURATION, &duration );
	
	//Open and initialize the hydraulic solver
	ENopenH();  
	ENinitH(0);  

	//Run the hydraulic solver one hydraulic time step at a time
	do 
	{  		
		ENrunH(&t);		
		// Retrieve hydraulic results for time t
		for (i=1; i <= nodeCount; i++)
		{
			ENgetnodevalue(i, EN_PRESSURE, &pressure);
			//ENgetnodeid(i, name);
			baseCasePressureMatrix[i-1] = pressure;		
		}		
		ENnextH(&tstep);  	
	} while (tstep > 0); 
	
	//Close the hydraulic solver
	ENcloseH();  
}

//FUNCTION
//Place a single leak at the index location in the network and run hydraulic analysis
//Determines how many pressure violations occur in the network by leak location
void oneLeak(int index, double emitterCoeff, int nodeCount, int columnNumber) 
{	
	int i;
	long t, tstep, hydraulicTimeStep;
	float pressure;
	
	i = 0;
	pressure = 0;
	
	ENgettimeparam(EN_HYDSTEP, &hydraulicTimeStep);
	
	//Create the leak
	ENsetnodevalue(index, EN_EMITTER, emitterCoeff);
	
	ENopenH();  
	ENinitH(0);

	//Run the hydraulic analysis
	do {  	
		ENrunH(&t);		
		if (t%hydraulicTimeStep == 0)
		{
			for (i = 1; i <= nodeCount; i++)
			{			
				ENgetnodevalue(i, EN_PRESSURE, &pressure);
				largePressureMatrix[i-1][columnNumber] = pressure;			
            }
         }
		ENnextH(&tstep); 
	} while (tstep > 0); 
	
	//Close the hydraulic solver
	ENcloseH();
	
	//"Fix" the leak
	ENsetnodevalue(index, EN_EMITTER, 0.0);
}

//FUNCTION
//Generalized multi-leak simulator
void nLeaks(int leakCount, int nodeCount) 
{
	long t, tstep, hydraulicTimeStep, duration;	
	float pressure, baseDemand, demand;
	int i;
	//char name[20];

	i = 0;
	totalDemand = pressure = baseDemand = demand = 0.0;
	
	ENgettimeparam(EN_HYDSTEP, &hydraulicTimeStep);
	
	ENgettimeparam( EN_DURATION, &duration );
	
	//Create the leaks
	for (i = 0; i < leakCount; i++)
	{
		//printf("\n\n Inside nLeaks: %d \n\n", leakNodes[i]);
		ENsetnodevalue(leakNodes[i], EN_EMITTER, leakMagnitudes[i]);
	}
	
	ENopenH();  
	ENinitH(0);
	
	//Run the hydraulic analysis
	do 
	{  	
		ENrunH(&t);
		
		for (i = 1; i <= nodeCount; i++)
		{			
			ENgetnodevalue(i, EN_PRESSURE, &pressure);						
			ENgetnodevalue(i, EN_DEMAND, &demand);												
			observedPressure[i-1] = (double)pressure;			
			totalDemand += demand;	
		}
		
		for (i = 0; i < leakCount; i++)
		{
			ENgetnodevalue(leakNodes[i], EN_BASEDEMAND, &baseDemand);					
			ENgetnodevalue(leakNodes[i], EN_DEMAND, &demand);			
			leakDemands[i] = (demand - baseDemand);
		}
		
		ENnextH(&tstep); 		
	} while (tstep > 0); 
	
	//Close the hydraulic solver
	ENcloseH();
	
	//"Fix" the leak
	for (i = 0; i < leakCount; i++)
	{
		ENsetnodevalue(leakNodes[i], EN_EMITTER, 0);
	}	
}


//FUNCTION
//Sum model error
double calculateError(int numNodes, double solution[])
{
	double errorSum;
	int i;
	
	i = 0;
	errorSum = 0.0;
	
	for (i = 0; i < numNodes; i++)
	{
		singleRunErrors[i] = fabs(realLeakValues[i] - solution[i]);
		//printf("node %d model error %f \n", (i+1), modelError[i]);
		errorSum += singleRunErrors[i];
	}
	//printf("Total Error: %f \n", errorSum);
	
	return errorSum;
}


//FUNCTION
//Create an output file for each simulation/optimization run
int writeSummaryFile(int k, int optimstatus, double objval, double sol[])
{	
	char sequentialFile[100], buffer[10], name[10];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, "/home/andrew/Ubuntu One/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Summary_");
	sprintf(buffer,"%d",k);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	for (i = 0; i < numOfLeaks; i++)
	{
		ENgetnodeid(leakNodes[i], name);
		fprintf(ptr_file, "Leak %d:, Node %d, Node ID:, %s, Magnitude:, %2.2f  \n",
			i, leakNodes[i], name, leakMagnitudes[i] );										
	}
	
	fprintf(ptr_file, "Delta:,%2.2f \n",delta);
	fprintf(ptr_file, "Total Demand: %f \n", totalDemand);
	fprintf(ptr_file, "Run #, %d, Model Error:, %f \n", (k + 1), modelError[k]);
	
	for (i = 0; i < numOfLeaks; i++)
	{
		fprintf(ptr_file, "Leak %d demand:, %f, Demand Fraction:, %f, %% \n",
			i, leakDemands[i], ((leakDemands[i]/totalDemand)* 100));
	}
	
	fprintf(ptr_file, "\nOptimization complete\n");
	if (optimstatus == GRB_OPTIMAL) 
	{
		fprintf(ptr_file, "Optimal objective:, %.4e\n", objval);
		for(i = 0; i < (totalNodeCount * 2); i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f \n", (i+1), sol[i]);
		}
		for(i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f, binary for, sol[%d] \n",
				(i+1), sol[i]), (i - (totalNodeCount * 2) + 1);
		}
	} else if (optimstatus == GRB_INF_OR_UNBD) 
	{
		fprintf(ptr_file, "Model is infeasible or unbounded\n");
	} else 
	{
		fprintf(ptr_file, "Optimization was stopped early\n");
	}
	
	fclose(ptr_file);
	return 0;
}


//FUNCTION
//Create an output file for each simulation/optimization run
int writeRawResults(int k, int optimstatus, double sol[])
{	
	char sequentialFile[100], buffer[10];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, "/home/andrew/Ubuntu One/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Run_");
	sprintf(buffer,"%d",k);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;	
	
	if (optimstatus == GRB_OPTIMAL) 
	{	
		for(i = 0; i < ((totalNodeCount * 2) - 1); i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f \n", (i+1), sol[i]);
		}
		for(i = ((totalNodeCount * 2) - 1); i < (totalNodeCount * 2); i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f", (i+1), sol[i]);
		}
	} else if (optimstatus == GRB_INF_OR_UNBD) 
	{
		fprintf(ptr_file, "Model is infeasible or unbounded\n");
	} else 
	{
		fprintf(ptr_file, "Optimization was stopped early\n");
	}
	
	fclose(ptr_file);
	return 0;
}


//FUNCTION
//Create an output file for each set of iterations
int writeErrorFile()
{	
	char sequentialFile[100], buffer[10];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, "/home/andrew/Ubuntu One/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Error.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Objective_Value, Model_Error\n");
	
	for (i = 0; i < (iterations - 1); i++)
	{
		fprintf(ptr_file, "%f, %f\n", objectiveValues[i], modelError[i]);										
	}
	for (i = (iterations - 1); i < iterations; i++)
	{
		fprintf(ptr_file, "%f, %f", objectiveValues[i], modelError[i]);										
	}
	
	fclose(ptr_file);
	return 0;
}

//FUNCTION
//Print the location and magnitude of leaks to file
int writeLeakFile(int k)
{
	int i;
	char sequentialFile[100], buffer[10], name[10];
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, "/home/andrew/Ubuntu One/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Leaks_");
	sprintf(buffer,"%d",k);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;	
	
	for(i = 0; i < (numOfLeaks - 1); i++)
	{			
		ENgetnodeid(leakNodes[i], name);
		fprintf(ptr_file,"leak %d, %d, Node ID, %s, Magnitude, %f \n",
			i, leakNodes[i], name, leakMagnitudes[i]);
	}
	for(i = (numOfLeaks - 1); i < numOfLeaks; i++)
	{			
		ENgetnodeid(leakNodes[i], name);
		fprintf(ptr_file,"leak %d, %d, Node ID, %s, Magnitude, %f",
			i, leakNodes[i], name, leakMagnitudes[i]);
	}
	
	fclose(ptr_file);
	return 0;	
}