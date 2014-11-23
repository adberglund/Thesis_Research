#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include "epanet2.h" 
#include "gurobi_c.h"

//September 10, 2013
//L1-Approximation (L1 calculates absolute error, in this case, between
//	simulated and "observed" scenarios and is being used for linear 
//	approximation of water distrubition network operations) using EPANET
//	and Guribo Optimizer

//The // bracketed variables currently serve as the adjustable parameters
//	for number of leaks and number of simulations  
//
//
int numOfLeaks = 2, iterations = 50;
double delta = 1, minLeakSize = 1.0, maxLeakSize = 10.0,
	binaryLeakLimit = 2.0, minLeakThreshold = 0.5;
char inputFile[50] = "Net3.inp";//"hanoi-1.inp"; //
char reportFile[50] = "Net3.rpt";//"hanoi.rpt"; //
char directoryString[50] = "L1_Iterative/";

//
//

int totalNodeCount, EPANETsimCounter;
int *leakNodes, *MIPStartSolution;
double totalDemand, averageDelta, averagePreviousDelta, bigM = 9999999999.99,
	totalTime, timePerIteration;
double *baseCasePressureMatrix, *observedPressure, *coefficients, *b, *bhat,
	*realLeakValues, *singleRunErrors, *leakDemands, *leakMagnitudes, 
	*LPmodelError, *MIPmodelError, *deltas, *previousDeltas, *leakGuesses,
	**largePressureMatrix, **largeA, **Ahat,  **I, *LPSolutions, *MIPSolutions, *tempSolutions,
	*LPobjectiveValues, *MIPobjectiveValues; 
char globalDirName[100];
clock_t startTime, endTime, iterationStartTime, iterationEndTime;


FILE *ptr_file;


void initializeArrays();
void populateMatricies(int);
void populateBMatrix(int);
void randomizeLeaks(int, int);
void printLeakInfo(int);
void analyzeBaseCase(int);
void oneLeak(int, double, int, int);
void nLeaks(int, int);
void findHighestMagnitudes(double *);
void forgeMIPStartSolution(double []);
double calculateError(int, double[]);
int writeSummaryFile(int, int, double, double[]);
int writeRawResults(int, int, double[]);
int writeInterimResults(int, int, int, double [], char *);
int writeLeakFile(int);
int writeErrorFile();
int setOutputDirectory();

int main(int argc, char *argv[]) 
{
	startTime = clock();
	GRBenv *env = NULL;
	GRBmodel *model = NULL;
	int  i, j, k, l, numNodes, storage, counter, directoryCode;
	double previousObjectiveValue;
	
	//Randomize the leak locations, commented out will use the same seeding 
	//for each run
	//srand(time(NULL));
	
	i = j = k = l = numNodes = counter = EPANETsimCounter = 0;
	averageDelta = averagePreviousDelta = previousObjectiveValue = 0.0;
	
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
	
	leakNodes = (int *) calloc(numOfLeaks, sizeof(int));
	MIPStartSolution = (int *) calloc(totalNodeCount, sizeof(int));
	
	baseCasePressureMatrix = (double *) calloc(totalNodeCount, sizeof(double));
	observedPressure = (double *) calloc(totalNodeCount, sizeof(double));
	coefficients = (double *) calloc((totalNodeCount * 2), sizeof(double));
	b = (double *) calloc(totalNodeCount, sizeof(double));
	bhat = (double *) calloc((totalNodeCount * 2), sizeof(double));
	realLeakValues = (double *) calloc(totalNodeCount, sizeof(double));
	singleRunErrors = (double *) calloc(totalNodeCount, sizeof(double));
	leakDemands = (double *) calloc(numOfLeaks, sizeof(double));
	leakMagnitudes = (double *) calloc(numOfLeaks, sizeof(double));
	LPmodelError = (double *) calloc(iterations, sizeof(double));
	MIPmodelError = (double *) calloc(iterations, sizeof(double));
	LPobjectiveValues = (double *) calloc(iterations, sizeof(double));
	MIPobjectiveValues = (double *) calloc(iterations, sizeof(double));
	
	LPSolutions = (double *) calloc(totalNodeCount * 2, sizeof(double));
	MIPSolutions = (double *) calloc(totalNodeCount * 2, sizeof(double));
	tempSolutions = (double *) calloc(totalNodeCount * 2, sizeof(double));
	
	deltas = (double *) calloc(totalNodeCount, sizeof(double));
	previousDeltas = (double *) calloc(totalNodeCount, sizeof(double));
	leakGuesses = (double *) calloc(binaryLeakLimit, sizeof(double));
	
	largePressureMatrix = (double **) calloc(totalNodeCount, sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
	{
		largePressureMatrix[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	
	largeA = (double **) calloc(totalNodeCount, sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
	{
		largeA[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	
	I = (double **) calloc(totalNodeCount, sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
	{
		I[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	
	Ahat = (double **) calloc( (totalNodeCount * 2), sizeof(double *) );
	for(i = 0; i < (totalNodeCount * 2); i++)
	{
		Ahat[i] = (double *) calloc( (totalNodeCount * 2), sizeof(double) );
	}
		 
	// Create environment 
 	error = GRBloadenv(&env, "L1_Iterative.log");
 	if (error) goto QUIT;
 	
 	directoryCode = setOutputDirectory();
 	
	//Create observation	
	for (k = 0; k < iterations; k++)
	{
		iterationStartTime = clock();
		
		EPANETsimCounter = 0;
		
		initializeArrays();
		
		randomizeLeaks(totalNodeCount, numOfLeaks);
							
		objval = 9999;
		counter = 0;
 		
		analyzeBaseCase(totalNodeCount);
		
		nLeaks(numOfLeaks, totalNodeCount);
		
		populateBMatrix(totalNodeCount);
		
		do{
		

			counter++;

			populateMatricies(totalNodeCount);		
		
			// Create an empty model 		
 			error = GRBnewmodel(env, &model, "L1Approx", 0, NULL, NULL, NULL, NULL, 
 				NULL);
 			if (error) goto QUIT;
 			 	
 			// Add variables 
 			for (i = 0; i < (totalNodeCount * 2); i++)
 			{
 				obj[i] = coefficients[i]; 			
 				vtype[i] = GRB_CONTINUOUS; 			
 			}
 			 				
			error = GRBaddvars(model, (totalNodeCount * 2), 0, NULL, NULL, NULL, 
				obj, NULL, NULL, vtype, NULL);
			if (error) goto QUIT;
			
			// Integrate new variables		
			error = GRBupdatemodel(model);
			if (error) goto QUIT;
			
			
			// First constraint: Ax <= b						
			for (i = 0; i < (totalNodeCount); i++)
			{		
				for (j = 0; j < (totalNodeCount); j++)
				{
					ind[j] = j;
					val[j] = Ahat[i][j];			
				}								
				ind[totalNodeCount] = j + i;
				val[totalNodeCount] = Ahat[i][j+i];
				error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
					GRB_LESS_EQUAL, bhat[i],NULL);			
				if (error) goto QUIT;
			}
			
			for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
			{
				
					for (j = 0; j < (totalNodeCount); j++)
					{
						ind[j] = j;
						val[j] = Ahat[i][j];			
					}								
					ind[totalNodeCount] = j + (i-totalNodeCount);
					val[totalNodeCount] = Ahat[i][j+(i-totalNodeCount)];
					error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
						GRB_LESS_EQUAL, bhat[i],NULL);			
					if (error) goto QUIT;
			}
			
			/*		
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
			*/
			
			error = GRBoptimize(model);
			if (error) goto QUIT;
							
			previousObjectiveValue = objval;
			
			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;
			
			error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
					(totalNodeCount * 2), sol);
				if (error) goto QUIT;
			
				
				
			if ((objval - previousObjectiveValue) < 0)
			{								
				writeInterimResults(k, counter, optimstatus, sol, "LP");
				
				for (i = 0; i < totalNodeCount * 2; i++)
				{
					LPSolutions[i] = sol[i];
				}
				
				binaryLeakLimit = 0.0;			
								
				for (i = 0; i < totalNodeCount; i++)
				{
					deltas[i] = 1.0;
					if (sol[i] >= minLeakThreshold)
						deltas[i] = sol[i];
					//printf("\t\tLP deltas[%d] = %f\n", i, deltas[i]);
					if (sol[i] > minLeakThreshold)
						binaryLeakLimit++;
				}
				
				LPobjectiveValues[k] = objval;
				
			}
			
			// Free model 
			GRBfreemodel(model);
		}while((objval - previousObjectiveValue) < 0);
		
		forgeMIPStartSolution(LPSolutions);
		objval = 999999;
		counter = 0;
		
		do
		{
			counter++;
										
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
 			 				
			error = GRBaddvars(model, (totalNodeCount * 3), 0, NULL, NULL, NULL,
				obj, NULL, NULL, vtype, NULL);
			if (error) goto QUIT;
			
			// Integrate new variables		
			error = GRBupdatemodel(model);
			if (error) goto QUIT;
			
			// First constraint: Ax <= b						
				for (i = 0; i < (totalNodeCount); i++)
				{

						for (j = 0; j < (totalNodeCount); j++)
						{
							ind[j] = j;
							val[j] = Ahat[i][j];			
						}								
						ind[totalNodeCount] = j + i;
						val[totalNodeCount] = Ahat[i][j+i];
						error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
							GRB_LESS_EQUAL, bhat[i],NULL);			
						if (error) goto QUIT;

				}
				
				for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
				{
					
						for (j = 0; j < (totalNodeCount); j++)
						{
							ind[j] = j;
							val[j] = Ahat[i][j];			
						}								
						ind[totalNodeCount] = j + (i-totalNodeCount);
						val[totalNodeCount] = Ahat[i][j+(i-totalNodeCount)];
						error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
							GRB_LESS_EQUAL, bhat[i],NULL);			
						if (error) goto QUIT;					
				}
			
			
			
			/*
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
			*/
			
			//Leak magnitude - (binary * bigM) <= 0
			for (i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
			{		
				ind[0] = (i - (totalNodeCount * 2)); 	ind[1] = i; 
				val[0] = 1.0; 		val[1] = -bigM ;
										
				error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL,0.0,
					NULL);
				if (error) goto QUIT;
			}
			
			// Limit sum of binaries to number of leaks searching for...		
			for (i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
			{		
				ind[i-(totalNodeCount * 2)] = i;
				val[i-(totalNodeCount * 2)] = 1.0;
			}								
			error = GRBaddconstr(model, totalNodeCount, ind, val, 
				GRB_LESS_EQUAL, binaryLeakLimit,NULL);
			if (error) goto QUIT;	
			
			for(i = 0; i < totalNodeCount; i++)
			{
				error = GRBsetdblattrelement(model, "Start", 
					i + (totalNodeCount * 2), MIPStartSolution[i]);
				if (error) goto QUIT;
			}
        	
			error = GRBoptimize(model);
			if (error) goto QUIT;	
			
			// Capture solution information		
			error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
			if (error) goto QUIT;
			
			previousObjectiveValue = objval;
			
			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;
			
			error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
				(totalNodeCount * 3), sol);
			if (error) goto QUIT;
			
			for (i = 0; i < totalNodeCount; i++)
			{				
				if  (sol[i] >= minLeakThreshold)
				{
					deltas[i] = sol[i];
				}
				else
					deltas[i] = 1.0;
				//printf("\t\tMIP deltas[%d] = %f\n", i, deltas[i]);
			}
			
			if ((objval - previousObjectiveValue) < 0)
			{
				writeInterimResults(k, counter, optimstatus, sol, "MIP");
				for (i = 0; i < totalNodeCount * 2; i++)
				{
					MIPSolutions[i] = sol[i];
				}				
				MIPobjectiveValues[k] = objval;
			}
			
			//for (i = totalNodeCount*2; i < totalNodeCount*3; i++)
			//{	
				//printf( "\n\t\t\tsol[%d] = %f", i, sol[i]);
			//}
			
			//findHighestMagnitudes(sol);	
			/*
			averagePreviousDelta = averageDelta;
			averageDelta = 0;
			
			for (i = 0; i < numOfLeaks; i++)
			{
				averageDelta += leakGuesses[i];
				//printf("leakGuess[%d] = %f \n", i, leakGuesses[i]);
			}
			
			averageDelta = averageDelta / numOfLeaks;
							
			for (i = 0; i < totalNodeCount; i++)
			{
				previousDeltas[i] = deltas[i];				
				deltas[i] = averageDelta;
			}
			*/
			//objectiveValues[k] = objval;
			
			
			forgeMIPStartSolution(sol);
			
			// Free model
			GRBfreemodel(model);
			
		}while((objval - previousObjectiveValue) < 0); 		
		
		LPmodelError[k] = calculateError(totalNodeCount, LPSolutions);
		MIPmodelError[k] = calculateError(totalNodeCount, MIPSolutions);
		/*
		
		//Set delta values for polishing step
		for (i = 0; i < totalNodeCount; i++)
		{			
			deltas[i] = sol[i];
		}
			
		//Polishing step that uses individual solution values for response 
		//	matrix generation in L1 approximation instead of averaging n
		//	highest magnitudes
		do
		{
			counter++;
			
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
 			 				
			error = GRBaddvars(model, (totalNodeCount * 3), 0, NULL, NULL, NULL,
				obj, NULL, NULL, vtype, NULL);
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
										
				error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL,0.0,
					NULL);
				if (error) goto QUIT;
			}
			
			// Limit sum of binaries to number of leaks searching for...		
			for (i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
			{		
				ind[i-(totalNodeCount * 2)] = i;
				val[i-(totalNodeCount * 2)] = 1.0;
			}								
			error = GRBaddconstr(model, totalNodeCount, ind, val, 
				GRB_LESS_EQUAL, binaryLeakLimit,NULL);
			if (error) goto QUIT;	
        	
			error = GRBoptimize(model);
			if (error) goto QUIT;
						
			// Capture solution information		
			error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
			if (error) goto QUIT;
			
			previousObjectiveValue = objval;
			
			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;
			
			error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
				(totalNodeCount * 3), sol);
			if (error) goto QUIT;
			
			for (i = totalNodeCount*2; i < totalNodeCount*3; i++)
			{	
				printf( "\n\t\t\tsol[%d] = %f", i, sol[i]);
			}
			
			for (i = 0; i < totalNodeCount; i++)
			{	
				if  (sol[i] > 0.01)
				{
					deltas[i] = sol[i];
				}
			}
			
			objectiveValues[k] = objval;
			modelError[k] = calculateError(totalNodeCount, sol);
			
			// Free model
			GRBfreemodel(model);
			
		}while((objval - previousObjectiveValue) < 0);		
		*/
		
		iterationEndTime = clock();
		timePerIteration = ((double)(iterationEndTime - iterationStartTime)) / CLOCKS_PER_SEC;
		
		printf("\nSolution Time: %.9f\n", timePerIteration);
		
		writeSummaryFile(k, optimstatus, objval, sol);
		writeRawResults(k, optimstatus, sol);
		writeLeakFile(k);
	}
	
	
	
	writeErrorFile();
	
	/*
	for (i = 0; i < totalNodeCount; i++) 
	{
		printf("realLeakValues[%d] = %f\n", i, realLeakValues[i]);
	}
	for (i = 0; i < numOfLeaks; i++)
	{
		printf("leakNodes[%d] = %d \t", i, leakNodes[i]);
		printf("leakMagnitudes[%d] = %f\n", i, leakMagnitudes[i]);
	}
	*/
	
	free(leakNodes);	
	free(MIPStartSolution);
	free(baseCasePressureMatrix);	
	free(observedPressure);	
	free(coefficients);	
	free(b);
	free(bhat);
	free(realLeakValues);	
	free(singleRunErrors);	
	free(leakDemands);	
	free(leakMagnitudes);	
	free(LPmodelError);
	free(MIPmodelError);	
	free(LPobjectiveValues);
	free(MIPobjectiveValues);
	free(LPSolutions);
	free(MIPSolutions);
	free(tempSolutions);
	free(deltas);
	free(previousDeltas);	
	free(leakGuesses);
	
	for(i = 0; i < totalNodeCount; i++)
	{
		free((void *)largePressureMatrix[i]);
	}
	free((void *)largePressureMatrix);
	
	for(i = 0; i < totalNodeCount; i++)
	{
		free((void *)largeA[i]);
	}
	free((void *)largeA);
	
	for(i = 0; i < totalNodeCount; i++)
	{
		free((void *)I[i]);
	}
	free((void *)I);
	
	for(i = 0; i < (totalNodeCount * 2); i++)
	{
		free((void *)Ahat[i]);
	}
	free((void *)Ahat);	
	
	
	
		QUIT:

		// Error reporting
		if (error) 
		{
			printf("ERROR: %s\n", GRBgeterrormsg(env));
			exit(1);
		}
		
		// Free model
		GRBfreemodel(model);
		
		// Free environment
		GRBfreeenv(env);
		
	ENclose();
	
	endTime = clock();
 	
 	totalTime = ((double)(endTime - startTime)) / CLOCKS_PER_SEC;
 	
 	printf("\nTotal Time Taken: %.9f\n\n", totalTime);
	
	return 0;
}

//FUNCTION
//Initialze various arrays to be populated during simulation
void initializeArrays()
{
	int i, j;	
	i = j = 0;
	delta = 1.0;
	
	//Array initialization	
	for (i = 0; i < totalNodeCount; i++)
	{
		observedPressure[i] = 0;
		baseCasePressureMatrix[i] = 0;
		b[i] = 0;
		realLeakValues[i] = 0.0;
		singleRunErrors[i] = 0.0;
		deltas[i] = delta;
		previousDeltas[i] = delta;		
		//printf("deltas[%d] = %f \n", i, deltas[i]);
	}
	
	for( i = 0; i < numOfLeaks; i++)
	{
		leakGuesses[i] = 0.0;
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


void populateBMatrix(int numNodes)
{
	int i, j;
	
	//printf("local numNodes variable = %d", numNodes);
	//getchar();
	i = j = 0;
	
	
	
	for (i = 0; i < numNodes; i++)
	{
		b[i] = 0;
	}
	
	
	for (i = 0; i < (numNodes * 2); i++)
	{
		bhat[i] = 0;
	}
	
	//Update b matrix
	

	for (j = 0; j < numNodes; j++)
	{
		b[j] = (baseCasePressureMatrix[j] - observedPressure[j]);	
		//printf("b[%d] = %f\n",i,b[i]);
	}
	//getchar();
	
	//Create b-hat
	
	for (i = 0; i < numNodes; i++)
	{
		bhat[i] = b[i];
		bhat[i + numNodes] = -b[i];
	}
	
	//for (i = numNodes; i < (numNodes * 2); i++)
	//{
		//bhat[i] = -b[i-numNodes];
		//printf("\t\t\t\tbhat[%d] = %f", i-numNodes, bhat[i-numNodes]);
		//printf("\tbhat[%d] = %f\n", i, bhat[i]);
	//}
	//getchar();
}



//FUNCTION
//Populate array values for the L1 Approximation
//Also calls single leak simulations for each node in the network
void populateMatricies(int numNodes)
{
	int i, j;
	
	//printf("local numNodes variable = %d", numNodes);
	//getchar();
	i = j = 0;
	
	
	for (i = 0; i < numNodes; i++)
	{
		b[i] = 0;		
	}
	
	for (i = 0; i < (numNodes * 2); i++)
	{
		bhat[i] = 0;		
	}
	
	//Update b matrix
	for (i = 0; i < numNodes; i++)
	{
		b[i] = (baseCasePressureMatrix[i] - observedPressure[i]);	
		//printf("b[%d] = %f\n",i,b[i]);
	}
	//getchar();
	//Create b-hat
	for (i = 0; i < numNodes; i++)
	{
		bhat[i] = b[i];
	}
	for (i = numNodes; i < (numNodes * 2); i++)
	{
		bhat[i] = -b[i-numNodes];
		//printf("\t\t\t\tbhat[%d] = %f", i-numNodes, bhat[i-numNodes]);
		//printf("\tbhat[%d] = %f\n", i, bhat[i]);
	}
	
	
	for(i = 1; i <= numNodes; i++)
	{		
		oneLeak(i, deltas[i-1], numNodes, i-1);		
	}
	
	//Update A matrix		
	for(i = 0; i < numNodes; i++)
	{		
		for(j = 0; j < numNodes; j++)
		{
			//THIS MAY BE WRONG!!!
			if(deltas[j] != 0)
			{
				largeA[i][j] = (baseCasePressureMatrix[i] - 
					largePressureMatrix[i][j]) / deltas[j];
			}
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
		//temp = (leakNodes[i]-1);		
		//realLeakValues[temp] = leakMagnitudes[i];
		realLeakValues[(leakNodes[i]-1)] = leakMagnitudes[i];
	}
	
}

void randomizeLeaks(int numNodes, int numOfLeaks)
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
			
			j = i;
			
			while(j >= 1)
			{
				if (leakNodes[i] == leakNodes[j-1])
				{
					do
					{			
						leakNodes[i] = (int)(rand()%numNodes)+1;
					}while(leakNodes[i] == leakNodes[j-1]);										
				}
				j--;
			}
		}			
		
		for(i = 0; i < numOfLeaks; i++)
		{
			leakMagnitudes[i] = drand48() * maxLeakSize;
			if(leakMagnitudes[i] < minLeakSize)
			{
				do
				{
					leakMagnitudes[i] = drand48() * maxLeakSize;
				}while(leakMagnitudes[i] < minLeakSize);
			}
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
	
	i = 0;
	pressure = 0.0;
	EPANETsimCounter++;
	
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
	EPANETsimCounter++;
	
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

	i = 0;
	totalDemand = pressure = baseDemand = demand = 0.0;
	EPANETsimCounter++;
	
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
//Find the n highest leak magnitudes in the current solution
void findHighestMagnitudes(double *solutions)
{
	int i,j,k, placeHolder;
	
	for (i = 0; i < binaryLeakLimit; i++)
	{
		leakGuesses[i] = 0.0;
	}
	
	for (i = 0; i < binaryLeakLimit; i++)
	{	
		for (j = 0; j < totalNodeCount; j++)
		{
			placeHolder = 0;				
			if(solutions[j] > leakGuesses[i])
			{
				for(k = 0; k < binaryLeakLimit; k++)
				{
					if(solutions[j] == leakGuesses[k])
						placeHolder++;
				}
				if (placeHolder == 0)
					leakGuesses[i] = solutions[j];
			}
		}
	}
}

void forgeMIPStartSolution(double sol[])
{
	int i;
	i = 0;
	
	for (i = 0; i < totalNodeCount; i++)
	{
		MIPStartSolution[i] = 0;
	}
	/*
	for (i = 0; i < binaryLeakLimit; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{
			if (leakGuesses[i] == sol[j])
				MIPStartSolution[currentPeriod][i] = 1;
		}
	}
	*/
	
	for (i = 0; i < totalNodeCount; i++)
	{
		if (sol[i] >= minLeakThreshold)
			MIPStartSolution[i] = 1;
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
		singleRunErrors[i] = 0;
		//printf("\t\t\t realLeakValues[%d] = %f \t solution[%d] = %f \n", 
			//i, realLeakValues[i], i, solution[i]);
	}
	
	for (i = 0; i < numNodes; i++)
	{
		singleRunErrors[i] = fabs(realLeakValues[i] - solution[i]);
		
		//printf("node %d model error %f \n", (i+1), singleRunErrors[i]);
		errorSum += singleRunErrors[i];
	}
	//printf("Total Error: %f \n", errorSum);
	//getchar();
	
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
	strcat(sequentialFile, globalDirName);
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
		fprintf(ptr_file, "Leak %d:, Node %d, Node ID:, %s, Magnitude:, %f\n",
			i, leakNodes[i], name, leakMagnitudes[i] );										
	}
	
	fprintf(ptr_file, "Delta:,%2.2f \n",delta);
	fprintf(ptr_file, "Total Demand: %f \n", totalDemand);
	fprintf(ptr_file, "Run #, %d, LP_Model Error:,%f,MIP_Model Error:,%f\n", (k + 1), LPmodelError[k],MIPmodelError[k]);
	
	for (i = 0; i < numOfLeaks; i++)
	{
		fprintf(ptr_file, "Leak %d demand:, %f, Demand Fraction:, %f, %% \n",
			i, leakDemands[i], ((leakDemands[i]/totalDemand)* 100));
	}
	fprintf(ptr_file, "Time To Solution:, %f, seconds\n", timePerIteration);
	
	fprintf(ptr_file, "\nOptimization complete\n");
	if (optimstatus == GRB_OPTIMAL) 
	{
		fprintf(ptr_file, "Optimal objective:, %.4e\n", objval);
		for (i = 0; i < totalNodeCount; i++)
		{
			ENgetnodeid((i+1), name);		  	
			fprintf(ptr_file, "sol[%d] =, %f, Node ID:, %s \n", 
				(i+1), sol[i], name);
		}
		for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
		{		  	
			fprintf(ptr_file, "sol[%d] =, %f, Error for sol[%d] \n",
				(i+1), sol[i], (i + 1 - totalNodeCount));
		}
		for(i = (totalNodeCount * 2); i < (totalNodeCount * 3); i++)
		{		  	
			fprintf(ptr_file, "sol[%d] =, %f, binary for, sol[%d] \n", 
				(i + 1), sol[i], (i + 1 - (totalNodeCount * 2)) );
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
	strcat(sequentialFile, globalDirName);
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

int writeInterimResults(int k, int count, int optimstatus, double sol[], char *method)
{	
	char sequentialFile[150], buffer[10], bufferDos[10];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/");
	strcat(sequentialFile, method);
	strcat(sequentialFile, "_Run_");
	sprintf(buffer,"%d",k);
	sprintf(bufferDos,"%d",count);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, "_");
	strcat(sequentialFile, bufferDos);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;	
	
	//if (optimstatus == GRB_OPTIMAL) 
	//{	
	for(i = 0; i < ((totalNodeCount * 2) - 1); i++)
	{		  	
		fprintf(ptr_file, "  sol[%d] =, %f \n", (i+1), sol[i]);
	}
	for(i = ((totalNodeCount * 2) - 1); i < (totalNodeCount * 2); i++)
	{		  	
		fprintf(ptr_file, "  sol[%d] =, %f", (i+1), sol[i]);
	}
	/*	
	} else if (optimstatus == GRB_INF_OR_UNBD) 
	{
		fprintf(ptr_file, "Model is infeasible or unbounded\n");
	} else 
	{
		fprintf(ptr_file, "Optimization was stopped early\n");
	}
	*/
	fclose(ptr_file);
	return 0;
}


//FUNCTION
//Create an output file for each set of iterations
int writeErrorFile()
{	
	char sequentialFile[100];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/Error.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "LP_Objective_Value,MIP_Objective_Value, LP_Model_Error,MIP_Model_Error\n");
	
	for (i = 0; i < (iterations - 1); i++)
	{
		fprintf(ptr_file, "%f,%f,%f,%f\n",LPobjectiveValues[i],MIPobjectiveValues[i],LPmodelError[i], MIPmodelError[i]);										
	}
	for (i = (iterations - 1); i < iterations; i++)
	{
		fprintf(ptr_file, "%f,%f,%f,%f",LPobjectiveValues[i],MIPobjectiveValues[i],LPmodelError[i], MIPmodelError[i]);										
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
	strcat(sequentialFile, globalDirName);
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

int setOutputDirectory()
{
	int status;
	char dirName[100], date[25];
	
	time_t seconds;
	struct tm *time_struct;
	
	time(&seconds);
	time_struct = localtime(&seconds);
	
	//Create summary CSV file for each set of leaks
	dirName[0] = '\0';
	strcat(dirName, "/home/andrew/Ubuntu One/Research/Thesis_Results/");
	strcat(dirName, directoryString);
	strftime(date, 50, "%Y_%m_%d", time_struct);
	strcat(dirName, date);
	status = mkdir(dirName, S_IRWXU | S_IRWXG | S_IRWXO);
	
	strcpy(globalDirName, dirName);
	
	//printf("Global dirname \t %s\n\n\n", globalDirName);
	
	if (status != 0)
	{
		printf("\nDirectory Creation Error\n");
		return status;
	}
	
	return 0;
}