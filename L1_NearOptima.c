#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include "epanet2.h" 
#include "gurobi_c.h"

#define SECONDS_PER_HOUR 3600
#define SECONDS_PER_DAY 86400
#define HOURS_PER_DAY 24
#define WARMUP_PERIOD (259200 + (0 * 3600))

//September 10, 2013
//L1-Approximation (L1 calculates absolute error, in this case, between
//	simulated and "observed" scenarios and is being used for linear 
//	approximation of water distrubition network operations) using EPANET
//	and Guribo Optimizer

//The // bracketed variables currently serve as the adjustable parameters
//	for number of leaks and number of simulations  
//
//
int numOfLeaks = 2, iterations = 1, numOfHours = 12, numOfNodesToIgnore = 8;	
double delta = 1, minLeakSize = 1.0, maxLeakSize = 10.0, minLeakThreshold =0.5,
	binaryLeakLimit = 0.0, sensorPercentOfTotalNodes = 1.0,
	numPeriodsPerSimulation;
char inputFile[50] = "Net3.inp";
char reportFile[50] = "Net3.rpt";
char directoryString[50] = "L1_HeatMap/";
char *nodesToIgnore[8] = {"10", "20", "40", "50", "60", "61", "601", "123"};

//
//

long simDuration;
int totalNodeCount, EPANETsimCounter, nearOptimaBinaryLimit;
int *leakNodes, *sensorNodes, **MIPStartSolution, *isConstraintWorthy;
double averageDelta, averagePreviousDelta, bigM = 999,
	totalTime, timePerIteration;
double ***baseCasePressureMatrix, **baseCaseDemand, ***observedPressure,
	**observedDemand, *coefficients, **b, **bhat,
	*realLeakValues, *singleRunErrors, *totalDemand,  **leakDemands, *leakMagnitudes, 
	**LPmodelError, **MIPmodelError, **deltas, *previousDeltas, *leakGuesses, 
	***largePressureMatrix, ***largeA, ***Ahat,  **I, **LPobjectiveValues, 
	**MIPobjectiveValues, **nearOptimaObjectiveValues,  **LPSolutions, 
	**MIPSolutions, **nearOptimaSolutions, *LPHeatMap, *MIPHeatMap, *LPWeightedHeatMap, 
	*MIPWeightedHeatMap; 
char globalDirName[100];
clock_t startTime, endTime, iterationStartTime, iterationEndTime;


FILE *ptr_file;


void initializeArrays(int, int);
void populateMatricies(int, int, int, int *);
void populateNearOptimaMatricies(int, int, int, int, int *);
void populateBMatrix(int, int);
void randomizeLeaks(int, int);
int setNumOfSensors(double);
void divinePressureSensorLocations(int);
int compare();
void printLeakInfo(int);
void analyzeBaseCase(int, int);
void oneLeak(int, double, int, int, int);
void nLeaks(int, int, int);
void findHighestMagnitudes(double *);
void calculateLeakDemand(int);
void forgeMIPStartSolution(int, double[]);
double calculateError(int, int, double **);
int makeHeatMap(int);
int writeIndividualSolutions(int, double **, char *);
int writeNearOptimaSolutions(int, int, double **);
int writeNearOptimaObjectives(int, double **);
int writeSummaryFile(int, int, int, double, double[]);
int writeRawResults(int, int, double[]);
int writeLeakFile(int);
int writePumpFile(int);
int writeErrorFile();
int writeAhat(int, char *);
int setOutputDirectory();

int main(int argc, char *argv[]) 
{
	startTime = clock();
	GRBenv *env = NULL;
	GRBmodel *model = NULL;
	int  i, j, k, l, m, numNodes, storage, counter, directoryCode, MIPCounter,
		numOfPressureSensors, periodCount;
	double previousObjectiveValue;
	
	//Randomize the leak locations, commented out will use the same seeding 
	//for each run
	//srand(time(NULL));
	
	i = j = k = l = numNodes = counter = EPANETsimCounter = MIPCounter = 0;
	periodCount = 0;
	averageDelta = averagePreviousDelta = previousObjectiveValue = 0.0;
	
	//Open EPANET & Input file
	ENopen(inputFile,reportFile,"");
	
	// Get the number of nodes
	ENgetcount(EN_NODECOUNT, &numNodes);
	ENgetcount(EN_TANKCOUNT, &storage);
	totalNodeCount = numNodes - storage;
	//totalNodeCount = totalNodeCount - numOfNodesToIgnore;
	
	ENgettimeparam(EN_DURATION, &simDuration);
	//printf("simDuration from EPANET = %ld\n", simDuration);
	simDuration = (simDuration - WARMUP_PERIOD) / SECONDS_PER_HOUR;
	//printf("simDuration = %ld\n", simDuration);
	
	numPeriodsPerSimulation = simDuration / numOfHours;
 	//printf("number of periods per simulation = %f", numPeriodsPerSimulation);
 	//getchar();
	
	int       error = 0;
	double    sol[(totalNodeCount * 3)];
	int       ind[(totalNodeCount * 3)];
	double    val[(totalNodeCount * 3)];
	double    obj[(totalNodeCount * 3)];
	char      vtype[(totalNodeCount * 3)];	
	int       optimstatus;
	double    objval = 999999;
	
	numOfPressureSensors = setNumOfSensors(sensorPercentOfTotalNodes);
	//printf("numOfPressureSensors = %d", numOfPressureSensors);	
	
	leakNodes = (int *) calloc(numOfLeaks, sizeof(int));
	sensorNodes = (int *) calloc(numOfPressureSensors, sizeof(int));
	isConstraintWorthy = (int *) calloc(totalNodeCount, sizeof(int));
	
	baseCasePressureMatrix = (double ***) calloc(numPeriodsPerSimulation, sizeof(double **));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		baseCasePressureMatrix[i] = (double **) calloc(numOfHours, sizeof(double *));
		for (j = 0; j < numOfHours; j++)
		{
			baseCasePressureMatrix[i][j] = (double *) calloc(totalNodeCount, sizeof(double));
		}
	}
	
	baseCaseDemand = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		baseCaseDemand[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	
	observedPressure = (double ***) calloc(numPeriodsPerSimulation, sizeof(double **));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		observedPressure[i] = (double **) calloc(numOfHours, sizeof(double *));
		for (j = 0; j < numOfHours; j++)
		{
			observedPressure[i][j] = (double *) calloc(totalNodeCount, sizeof(double));
		}
	}
	
	MIPStartSolution = (int **) calloc(numPeriodsPerSimulation, sizeof(int *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		MIPStartSolution[i] = (int *) calloc(totalNodeCount, sizeof(int));
	}
		
	observedDemand = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		observedDemand[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}	
	//observedDemand = (double *) calloc(totalNodeCount, sizeof(double));
	
	coefficients = (double *) calloc((totalNodeCount * 2), sizeof(double));
	
	b = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		b[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	bhat = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		bhat[i] = (double *) calloc((totalNodeCount * 2), sizeof(double));
	}
	
	realLeakValues = (double *) calloc(totalNodeCount, sizeof(double));
	singleRunErrors = (double *) calloc(totalNodeCount, sizeof(double));
	
	totalDemand = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	
	leakDemands = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		leakDemands[i] = (double *) calloc(numOfLeaks, sizeof(double));
	}
	
	leakMagnitudes = (double *) calloc(numOfLeaks, sizeof(double));
	LPmodelError = (double **) calloc(iterations, sizeof(double *));
	for (i = 0; i < iterations; i++)
	{
		LPmodelError[i] = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	}
	
	MIPmodelError = (double **) calloc(iterations, sizeof(double *));	
	for (i = 0; i < iterations; i++)
	{
		MIPmodelError[i] = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	}
	
	
	deltas = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i ++)
	{
		deltas[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	
	previousDeltas = (double *) calloc(totalNodeCount, sizeof(double));
	
	
	largePressureMatrix = (double ***) calloc(numOfHours, sizeof(double **));
	for(i = 0; i < numOfHours; i++)
	{
		largePressureMatrix[i] = (double **) calloc(totalNodeCount, sizeof(double *));
		for (j = 0; j < totalNodeCount; j++)
		{
			largePressureMatrix[i][j] = (double *) calloc(totalNodeCount, sizeof(double));
		}
	}
	
	largeA = (double ***) calloc(numPeriodsPerSimulation, sizeof(double **));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		largeA[i] = (double **) calloc(totalNodeCount, sizeof(double *));
		for(j = 0; j < totalNodeCount; j++)
		{
			largeA[i][j] = (double *) calloc(totalNodeCount, sizeof(double));
		}
	}
	
	I = (double **) calloc(totalNodeCount, sizeof(double *));
	for(i = 0; i < totalNodeCount; i++)
	{
		I[i] = (double *) calloc(totalNodeCount, sizeof(double));
	}
	
	Ahat = (double ***) calloc( numPeriodsPerSimulation, sizeof(double **));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		Ahat[i] = (double **) calloc( (totalNodeCount * 2), sizeof(double *) );
		for(j = 0; j < (totalNodeCount * 2); j++)
		{
			Ahat[i][j] = (double *) calloc( (totalNodeCount * 2), sizeof(double) );
		}
	}
	
	LPobjectiveValues = (double **) calloc(iterations, sizeof(double *));
	for (i = 0; i < iterations; i++)
	{
		LPobjectiveValues[i] = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	}
	
	MIPobjectiveValues = (double **) calloc(iterations, sizeof(double *));
	for (i = 0; i < iterations; i++)
	{
		MIPobjectiveValues[i] = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	}
	
	nearOptimaObjectiveValues = (double **) calloc((numPeriodsPerSimulation), sizeof(double *));
	for (i = 0; i < (numPeriodsPerSimulation); i++)
	{
		nearOptimaObjectiveValues[i] = (double *) calloc((numPeriodsPerSimulation), sizeof(double));
	}
	
	LPSolutions = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		LPSolutions[i] = (double *) calloc((totalNodeCount * 2), sizeof(double));
	}
	
	MIPSolutions = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		MIPSolutions[i] = (double *) calloc((totalNodeCount * 2), sizeof(double));
	}

	nearOptimaSolutions = (double **) calloc(numPeriodsPerSimulation, sizeof(double *));
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		nearOptimaSolutions[i] = (double *) calloc((totalNodeCount * 2), sizeof(double));
	}	
	
	LPHeatMap = (double *) calloc(totalNodeCount, sizeof(double));
	MIPHeatMap = (double *) calloc(totalNodeCount, sizeof(double));
	
	LPWeightedHeatMap = (double *) calloc(totalNodeCount, sizeof(double));
	MIPWeightedHeatMap = (double *) calloc(totalNodeCount, sizeof(double));
	
	// Create environment 
 	error = GRBloadenv(&env, "L1_Iterative.log");
 	if (error) goto QUIT;
 	
 	directoryCode = setOutputDirectory();
 	if (directoryCode != 0)
 		printf("\nDirectory Creation Error\n");
 	
 	divinePressureSensorLocations(numOfPressureSensors);
 	
 	//for (i = 0; i < numOfPressureSensors; i++)
 	//{
 		//printf("\n sensorNodes[%d] = %d", i, sensorNodes[i]);
 	//}
 	//getchar();
 	
	//Create observation
	for (k = 0; k < iterations; k++)
	{
		iterationStartTime = clock();
		
		EPANETsimCounter = 0;
		counter = 0;
		periodCount = 0;
		
		initializeArrays(numOfPressureSensors, k);
		
		randomizeLeaks(totalNodeCount, numOfLeaks);
		
		do{		
 		
			objval = 9999;
			
			analyzeBaseCase(numOfPressureSensors, periodCount);
			
			nLeaks(numOfLeaks, numOfPressureSensors, periodCount);										
			
			populateBMatrix(periodCount, totalNodeCount);
			
			printLeakInfo(numOfLeaks);
						
			//calculateLeakDemand(periodCount);
			
			//for (i = 0; i < totalNodeCount; i++)
			//{
				//printf("\n\tisConstraintWorthy[%d] = %d", i, isConstraintWorthy[i]);
			//}
			
				
			do{
				
				populateMatricies(periodCount, totalNodeCount, numOfPressureSensors, sensorNodes);
				//writeAhat(k, "LP");		
				//printf("\n\n\tSEG FAULT TEST 9999\n");
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
					if (isConstraintWorthy[i] == 1)
					{
						for (j = 0; j < (totalNodeCount); j++)
						{
							ind[j] = j;
							val[j] = Ahat[periodCount][i][j];			
						}								
						ind[totalNodeCount] = j + i;
						val[totalNodeCount] = Ahat[periodCount][i][j+i];
						error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
							GRB_LESS_EQUAL, bhat[periodCount][i],NULL);			
						if (error) goto QUIT;
					}
				}
				
				for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
				{
					if (isConstraintWorthy[i-totalNodeCount] == 1)
					{
						for (j = 0; j < (totalNodeCount); j++)
						{
							ind[j] = j;
							val[j] = Ahat[periodCount][i][j];			
						}								
						ind[totalNodeCount] = j + (i-totalNodeCount);
						val[totalNodeCount] = Ahat[periodCount][i][j+(i-totalNodeCount)];
						error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
							GRB_LESS_EQUAL, bhat[periodCount][i],NULL);			
						if (error) goto QUIT;
					}
				}
				
				error = GRBoptimize(model);
				if (error) goto QUIT;
				
				error = GRBwrite(model, "L1_Iterative_LP.lp");
				if (error) goto QUIT;
				
				error = GRBwrite(model, "L1_Iterative_LP.sol");
				if (error) goto QUIT;
				
				previousObjectiveValue = objval;
				
				error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
				if (error) goto QUIT;
				
				if ((objval - previousObjectiveValue) < 0)
				{
				
					error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
						(totalNodeCount * 2), sol);
					if (error) goto QUIT;
					
					for (i = 0; i < totalNodeCount * 2; i++)
					{
						LPSolutions[periodCount][i] = sol[i];
					}
        	    	
	    	    	
					binaryLeakLimit = 0.0;			
									
					for (i = 0; i < totalNodeCount; i++)
					{
						deltas[periodCount][i] = 1.0;
						if (sol[i] >= minLeakThreshold)
							deltas[periodCount][i] = sol[i];
						//printf("\t\tLP deltas[%d] = %f\n", i, deltas[i]);
						if (sol[i] > minLeakThreshold)
							binaryLeakLimit++;
					}
					
					
					
					leakGuesses = (double *) calloc(binaryLeakLimit, sizeof(double));
					
					for( i = 0; i < binaryLeakLimit; i++)
					{
						leakGuesses[i] = 0.0;
					}
			    	
					averageDelta = averagePreviousDelta = 0.0;
					
					findHighestMagnitudes(sol);
					
					for (i = 0; i < binaryLeakLimit; i++)
					{
						averageDelta += leakGuesses[i];
					}
					
					averageDelta = averageDelta / binaryLeakLimit;
					
					//printf("\n\n averageDelta = %f", averageDelta);
					
					for (i = 0; i < totalNodeCount; i++)
					{
						previousDeltas[i] = deltas[periodCount][i];
						averagePreviousDelta += previousDeltas[i];			
						//deltas[i] = averageDelta;
						//printf("\t\t\t This is the new delta: %f \n", deltas[i]);
					}		
					averagePreviousDelta = averagePreviousDelta / totalNodeCount;
					
					LPobjectiveValues[k][periodCount] = objval;
					
					forgeMIPStartSolution(periodCount, sol);
					free(leakGuesses);
				}
				//getchar();
				// Free model 
				GRBfreemodel(model);
				
				
				counter++;	
								
			}while((objval - previousObjectiveValue) < 0);
			
		
		
			//Reset the objval to a high value to give MIP a chance to iterate
			objval = 999999;
			MIPCounter = 0;
			//periodCount = 0;
			
			//Run MIP Polishing Step
			do
			{
				counter++;
				
				populateMatricies(periodCount, totalNodeCount, numOfPressureSensors, sensorNodes);
				//writeAhat(k, "Polish");
			
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
					if (isConstraintWorthy[i] == 1)
					{
						for (j = 0; j < (totalNodeCount); j++)
						{
							ind[j] = j;
							val[j] = Ahat[periodCount][i][j];			
						}								
						ind[totalNodeCount] = j + i;
						val[totalNodeCount] = Ahat[periodCount][i][j+i];
						error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
							GRB_LESS_EQUAL, bhat[periodCount][i],NULL);			
						if (error) goto QUIT;
					}
				}
				
				for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
				{
					if (isConstraintWorthy[i-totalNodeCount] == 1)
					{
						for (j = 0; j < (totalNodeCount); j++)
						{
							ind[j] = j;
							val[j] = Ahat[periodCount][i][j];			
						}								
						ind[totalNodeCount] = j + (i-totalNodeCount);
						val[totalNodeCount] = Ahat[periodCount][i][j+(i-totalNodeCount)];
						error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
							GRB_LESS_EQUAL, bhat[periodCount][i],NULL);			
						if (error) goto QUIT;
					}
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
					GRB_LESS_EQUAL, binaryLeakLimit, NULL);
				if (error) goto QUIT;	
        		
				
				//if(MIPCounter == 0)
				//{
				for(i = 0; i < totalNodeCount; i++)
				{
					error = GRBsetdblattrelement(model, "Start", 
						i + (totalNodeCount * 2), MIPStartSolution[periodCount][i]);
					if (error) goto QUIT;
				}
				//}
				
				error = GRBoptimize(model);
				if (error) goto QUIT;
				
				error = GRBwrite(model, "L1_Iterative_Polish.lp");
				if (error) goto QUIT;
			
				error = GRBwrite(model, "L1_Iterative_Polish.sol");
				if (error) goto QUIT;
							
				// Capture solution information		
				error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
				if (error) goto QUIT;
				
				printf("previous obj = %f \n objval = %f\n", previousObjectiveValue, objval);
				previousObjectiveValue = objval;
				
				error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
				if (error) goto QUIT;
				
				
				error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
					(totalNodeCount * 3), sol);
				if (error) goto QUIT;
				
				//for (i = 0; i < totalNodeCount*3; i++)
				//{	
					//printf( "\nsol[%d] = %f", i, sol[i]);
				//}
				//printf( "\n");
				//getchar();
				
				//for (i = 0; i < totalNodeCount; i++)
				//{
					//deltas[i] = 1.0;
				//}
				
				for (i = 0; i < totalNodeCount; i++)
				{	
					if  (sol[i] >= minLeakThreshold)
					{
						deltas[periodCount][i] = sol[i];
					}
					else
						deltas[periodCount][i] = 1.0;
					//printf("\t\tMIP deltas[%d] = %f\n", i, deltas[i]);
				}
				
				if ((objval - previousObjectiveValue) < 0)
				{
					for (i = 0; i < totalNodeCount * 2; i++)
					{
						MIPSolutions[periodCount][i] = sol[i];
					}				
					MIPobjectiveValues[k][periodCount] = objval;
				}				
				
				forgeMIPStartSolution(periodCount, sol);
				
				// Free model
				GRBfreemodel(model);
				MIPCounter++;
				
			}while((objval - previousObjectiveValue) < 0);		
			
			LPmodelError[k][periodCount] = calculateError(totalNodeCount, 
				periodCount, LPSolutions);
			MIPmodelError[k][periodCount] = calculateError(totalNodeCount, 
				periodCount, MIPSolutions);
			
			periodCount++;
			printf("\n\n\n\t\t\t Iteration: %d\t Sub-Period: %d\n", k, periodCount);
		
		}while(periodCount < numPeriodsPerSimulation);
		
		
		//Insert near optima exploration mheah
		for (l = 0; l < (numPeriodsPerSimulation); l++)
			{
				nearOptimaBinaryLimit = 0;
				for (i = 0; i < totalNodeCount; i++)
				{
					if(MIPSolutions[l][i] > 0);
					{
						nearOptimaBinaryLimit++;
					}
				}
				for (i = 0; i < totalNodeCount; i++)
				{
					deltas[l][i] = 1.0;
					if (MIPSolutions[l][i] >= minLeakThreshold)
					{
						deltas[l][i] = MIPSolutions[l][i];
					}
				}
				for (m = 0; m < (numPeriodsPerSimulation); m++)	
				{				
					
					populateNearOptimaMatricies(l, m, totalNodeCount, numOfPressureSensors, sensorNodes);
					//writeAhat(k, "Polish");
			    	
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
						if (isConstraintWorthy[i] == 1)
						{
							for (j = 0; j < (totalNodeCount); j++)
							{
								ind[j] = j;
								val[j] = Ahat[m][i][j];			
							}								
							ind[totalNodeCount] = j + i;
							val[totalNodeCount] = Ahat[m][i][j+i];
							error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
								GRB_LESS_EQUAL, bhat[m][i],NULL);			
							if (error) goto QUIT;
						}
					}
					
					for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
					{
						if (isConstraintWorthy[i-totalNodeCount] == 1)
						{
							for (j = 0; j < (totalNodeCount); j++)
							{
								ind[j] = j;
								val[j] = Ahat[m][i][j];			
							}								
							ind[totalNodeCount] = j + (i-totalNodeCount);
							val[totalNodeCount] = Ahat[m][i][j+(i-totalNodeCount)];
							error = GRBaddconstr(model, (totalNodeCount + 1), ind, val, 
								GRB_LESS_EQUAL, bhat[m][i],NULL);			
							if (error) goto QUIT;
						}
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
						GRB_LESS_EQUAL, nearOptimaBinaryLimit, NULL);
					if (error) goto QUIT;	
        			
					/*
					//if(MIPCounter == 0)
					//{
					for(i = 0; i < totalNodeCount; i++)
					{
						error = GRBsetdblattrelement(model, "Start", 
							i + (totalNodeCount * 2), MIPStartSolution[periodCount][i]);
						if (error) goto QUIT;
					}
					//}
					*/
					
					//Force leak magnitude solution via constraint to obtain objective
					//value at each time point
					for (i = 0; i < totalNodeCount; i++)
					{		
						ind[0] = i; 
						val[0] = 1.0;
												
						error = GRBaddconstr(model, 1, ind, val, GRB_EQUAL,MIPSolutions[l][i],
							NULL);
						if (error) goto QUIT;
					}
					
					error = GRBoptimize(model);
					if (error) goto QUIT;
					
					error = GRBwrite(model, "L1_NearOptima.lp");
					if (error) goto QUIT;
			    	
					error = GRBwrite(model, "L1_NearOptima.sol");
					if (error) goto QUIT;
								
					// Capture solution information		
					error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
					if (error) goto QUIT;
					
					//printf("previous obj = %f \n objval = %f\n", previousObjectiveValue, objval);
					//previousObjectiveValue = objval;
					
					error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
					if (error) goto QUIT;
					
					nearOptimaObjectiveValues[l][m] = objval;
					printf("\n\n\n\t\t\t Solution from time: %d\t Evaluated at time: %d\n", l, m);
					
					error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 
						(totalNodeCount * 3), sol);
					if (error) goto QUIT;
					printf("\n\n\n\t\t\tseg fault test No. 0\n");
					for (i = 0; i < totalNodeCount * 2; i++)
					{
						nearOptimaSolutions[m][i] = sol[i];
					}
					printf("\n\n\n\t\t\tseg fault test No. 0.5\n");
					writeNearOptimaSolutions(k, l, nearOptimaSolutions);						
					//forgeMIPStartSolution(periodCount, sol);
					
					// Free model
					GRBfreemodel(model);					
				
				}			
			}
			
			writeNearOptimaObjectives(k, nearOptimaObjectiveValues);
			
		
		
		iterationEndTime = clock();
		timePerIteration = ((double)(iterationEndTime - iterationStartTime))
			/ CLOCKS_PER_SEC;
				
		printf("\nSolution Time: %.9f\n", timePerIteration);				
		
		//makeHeatMap(k);
		//writeSummaryFile(k, optimstatus, numOfPressureSensors, objval, sol);
		//writeRawResults(k, optimstatus, sol);
		writeLeakFile(k);
		//writeIndividualSolutions(k, LPSolutions, "LP");
		writeIndividualSolutions(k, MIPSolutions, "MIP");
				
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
	free(sensorNodes);
	free(isConstraintWorthy);
	
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		for (j = 0; j < numOfHours; j++)
		{
			free((void *)baseCasePressureMatrix[i][j]);
		}
		free((void *)baseCasePressureMatrix[i]);
	}
	free((void *)baseCasePressureMatrix);
	
	/*
	for (i = 0; i < numOfHours; i++)
	{
		free((void *)baseCasePressureMatrix[i]);
	}
	//printf("\n\t\t\t\t\tERROR TESTERrrrrr\n");
	free((void *)baseCasePressureMatrix);
	*/
	
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)baseCaseDemand[i]);
	}
	free(baseCaseDemand);
	
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		for (j = 0; j < numOfHours; j++)
		{
			free((void *)observedPressure[i][j]);
		}
		free((void *)observedPressure[i]);
	}
	free((void *)observedPressure);
	
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)observedDemand[i]);
	}
	free(observedDemand);
	//free(observedDemand);
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)MIPStartSolution[i]);
	}
	free(MIPStartSolution);
	
	free(coefficients);	
	
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)b[i]);
	}
	free(b);
	
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)bhat[i]);
	}
	free(bhat);
	//free(b);
	//free(bhat);
	
	free(realLeakValues);	
	free(singleRunErrors);	
	
	free((void *)totalDemand);
	for(i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)leakDemands[i]);
	}
	free((void *)leakDemands);	
	free(leakMagnitudes);	
	
	for (i = 0; i < iterations; i++)
	{
		free((void *)LPmodelError[i]);
		free((void *)MIPmodelError[i]);
	}
	free((void *)LPmodelError);
	free((void *)MIPmodelError);
	
	for (i = 0; i < numPeriodsPerSimulation; i ++)
	{
		free((void *)deltas[i]);
	}
	free(deltas);
	
	free(previousDeltas);	
	
	
	for(i = 0; i < numOfHours; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{
			free((void *)largePressureMatrix[i][j]);
		}
		free((void *)largePressureMatrix[i]);
	}
	free((void *)largePressureMatrix);
	
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		for(j = 0; j < totalNodeCount; j++)
		{
			free((void *)largeA[i][j]);
		}
		free((void *)largeA[i]);
	}
	free((void *)largeA);
	
	for(i = 0; i < totalNodeCount; i++)
	{
		free((void *)I[i]);
	}
	free((void *)I);
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		for(j = 0; j < (totalNodeCount * 2); j++)
		{
			free((void *)Ahat[i][j]);
		}
		free((void *)Ahat[i]);
	}
	free((void *)Ahat);
	/*
	for(i = 0; i < (totalNodeCount * 2); i++)
	{
		free((void *)Ahat[i]);
	}
	free((void *)Ahat);
	*/
	
	for (i = 0; i < iterations; i++)
	{
		free((void *)LPobjectiveValues[i]);
	}
	free((void *)LPobjectiveValues);
	
	for (i = 0; i < iterations; i++)
	{
		free((void *)MIPobjectiveValues[i]);
	}
	free((void *)MIPobjectiveValues);
	
	for (i = 0; i < (numPeriodsPerSimulation); i++)
	{
		free((void *)nearOptimaObjectiveValues[i]);		 
	}
	free((void *)nearOptimaObjectiveValues);
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)LPSolutions[i]);
	}
	free((void *)LPSolutions);
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)MIPSolutions[i]);
	}
	free((void *)MIPSolutions);
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		free((void *)nearOptimaSolutions[i]);
	}
	free((void *)nearOptimaSolutions);
	
	free((void *)LPHeatMap);
	free((void *)MIPHeatMap);
	free((void *)LPWeightedHeatMap);
	free((void *)MIPWeightedHeatMap);
	
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
void initializeArrays(int numOfPressureSensors, int currentRun)
{
	int i, j, k;	
	i = j = k = 0;
	delta = 1.0;
	
	//Array initialization	
	for (i = 0; i < totalNodeCount; i++)
	{	
		realLeakValues[i] = 0.0;
		singleRunErrors[i] = 0.0;		
		previousDeltas[i] = delta;		
		LPHeatMap[i] = 0;
		MIPHeatMap[i] = 0;
		//printf("deltas[%d] = %f \n", i, deltas[i]);
	}
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{
			observedDemand[i][j] = 0;
			baseCaseDemand[i][j] = 0;
			b[i][j] = 0;
			deltas[i][j] = delta;
		}
	}
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		for (j = 0; j < (totalNodeCount * 2); j++)
		{
			bhat[i][j] = 0;
		}
	}
	
	for (i = 0; i < (numPeriodsPerSimulation); i++)
	{
		for (j = 0; j < (numPeriodsPerSimulation); j++)
		{
			nearOptimaObjectiveValues[i][j] = 0;
		}
	}
	
	

	for (i = 0; i < (totalNodeCount * 2); i++)
	{		
		coefficients[i] = 0;		
	}
	
	if (currentRun == 0)
	{
		for (i = 0; i < numPeriodsPerSimulation; i++)
		{
			for (j = 0; j < (totalNodeCount * 2); j++)
			{
				LPSolutions[i][j] = 0;
				MIPSolutions[i][j] = 0;						
				nearOptimaSolutions[i][j] = 0;	
			}
		}	
	
		for (i = 0; i < iterations; i++)
		{
			for (j = 0; j < numPeriodsPerSimulation; j++)
			{
				LPobjectiveValues[i][j] = 0;
				MIPobjectiveValues[i][j] = 0;			
			}
		}
	
	}
	
	for (k = 0; k < numOfHours; k++)
	{
		for (i = 0; i < totalNodeCount; i++)
		{
			for (j = 0; j < totalNodeCount; j++)
			{
				largePressureMatrix[k][i][j] = 0;
						
			}
		}
	}
		
	for (k = 0; k < numPeriodsPerSimulation; k++)
	{
		for (i = 0; i < (totalNodeCount * 2); i++)
		{
			for (j = 0; j < (totalNodeCount * 2); j++)
			{
				Ahat[k][i][j] = 0;		
			}
		}
	}
	
	for (k = 0; k < numPeriodsPerSimulation; k++)
	{
		for (i = 0; i < totalNodeCount; i++)
		{
			for (j = 0; j < totalNodeCount; j++)
			{
				largeA[k][i][j] = 0;		
			}
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
		coefficients[i + totalNodeCount] = 1.0;
	}	
	//for (i = totalNodeCount; i < (totalNodeCount * 2); i++)
	//{
		//coefficients[i] = 1.0;
	//}
	
}

void populateBMatrix(int currentPeriod, int numNodes)
{
	int i, j;
	
	//printf("local numNodes variable = %d", numNodes);
	//getchar();
	i = j = 0;
	
	
	
	for (i = 0; i < numNodes; i++)
	{
		b[currentPeriod][i] = 0;
	}
	
	
	for (i = 0; i < (numNodes * 2); i++)
	{
		bhat[currentPeriod][i] = 0;
	}

	
	//for (i = 0; i < numNodes; i++)
	//{
		//b[i] = 0;		
	//}
	
	//for (i = 0; i < (numNodes * 2); i++)
	//{
		//bhat[i] = 0;		
	//}
	
	//Update b matrix
	
		for (i = 0; i < numOfHours; i++)
		{
			for (j = 0; j < numNodes; j++)
			{
				b[currentPeriod][j] += (baseCasePressureMatrix[currentPeriod][i][j] - observedPressure[currentPeriod][i][j]);	
				//printf("b[%d] = %f\n",i,b[i]);
			}
		}
	
	
	for (i = 0; i < numNodes; i++)
	{
		b[currentPeriod][i] = b[currentPeriod][i] / numOfHours;		
		//printf("b[%d] = %f\n", j, b[j]);
	}

	//getchar();
	
	//Create b-hat
	
	for (i = 0; i < numNodes; i++)
	{
		bhat[currentPeriod][i] = b[currentPeriod][i];
		bhat[currentPeriod][i + numNodes] = -b[currentPeriod][i];
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
void populateMatricies(int currentPeriod, int numNodes, int numOfPressureSensors, int *sensorNodes)
{
	int i, j, k;
	
	//printf("local numNodes variable = %d", numNodes);
	//getchar();
	i = j = k = 0;
	
	/*
	for (i = 0; i < numNodes; i++)
	{
		b[i] = 0;		
	}
	
	for (i = 0; i < (numNodes * 2); i++)
	{
		bhat[i] = 0;		
	}
	
	//Update b matrix
	for (i = 0; i < numOfHours; i++)
	{
		for (j = 0; j < numNodes; j++)
		{
			b[j] += (baseCasePressureMatrix[i][j] - observedPressure[i][j]);	
			//printf("b[%d] = %f\n",i,b[i]);
		}
	}
	
	for (i = 0; i < numNodes; i++)
	{
		b[i] = b[i] / numOfHours;		
		//printf("b[%d] = %f\n", j, b[j]);
	}
	//getchar();
	
	//Create b-hat
	for (i = 0; i < numNodes; i++)
	{
		bhat[i] = b[i];
		bhat[i + numNodes] = -b[i];
	}
	for (i = numNodes; i < (numNodes * 2); i++)
	{
		//bhat[i] = -b[i-numNodes];
		printf("\t\t\t\tbhat[%d] = %f", i-numNodes, bhat[i-numNodes]);
		printf("\tbhat[%d] = %f\n", i, bhat[i]);
	}
	*/
	
	for (k = 0; k < numOfHours; k++)
	{
		for (i = 0; i < totalNodeCount; i++)
		{
			for (j = 0; j < totalNodeCount; j++)
			{
				largePressureMatrix[k][i][j] = 0;
				//largeA[i][j] = 0;
			}
		}
	}

	for (i = 0; i < totalNodeCount; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{	
			largeA[currentPeriod][i][j] = 0;
		}
	}

	for(i = 1; i <= numNodes; i++)
	{		
		oneLeak(i, deltas[currentPeriod][i-1], numOfPressureSensors, i-1, currentPeriod);		
	}
	
	
	//Update A matrix
	
	for (k = 0; k < numOfHours; k++)
	{
		for(i = 0; i < numNodes; i++)
		{		
			for(j = 0; j < numNodes; j++)
			{
				largeA[currentPeriod][i][j] += (baseCasePressureMatrix[currentPeriod][k][i] - 
					largePressureMatrix[k][i][j]) / deltas[currentPeriod][j]; // / delta;			
			}			
		}
	}
	
	
	for(i = 0; i < numNodes; i++)
	{		
		for(j = 0; j < numNodes; j++)
		{
			//if (deltas[j] != 0)
			//{
				largeA[currentPeriod][i][j] = largeA[currentPeriod][i][j] / (numOfHours); // / delta;
			//}
		}			
	}
	/*
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
	*/
	
	//Create A-hat
	for(i = 0; i < numNodes; i++)
	{
		for(j = 0; j < numNodes; j++)
		{
			Ahat[currentPeriod][i][j] =  largeA[currentPeriod][i][j];
		}
	}
	for(i = numNodes; i < (numNodes * 2); i++)
	{
		for(j = 0; j < numNodes; j++)
		{
			Ahat[currentPeriod][i][j] = -largeA[currentPeriod][i-numNodes][j];
		}
	}
	for(i = 0; i < numNodes; i++)
	{
		for(j = numNodes; j < (numNodes * 2); j++)
		{
			Ahat[currentPeriod][i][j] = -I[i][j-numNodes];
		}
	}
	for(i = numNodes; i < (numNodes * 2); i++)
	{
		for(j = numNodes; j < (numNodes * 2); j++)
		{
			Ahat[currentPeriod][i][j] = -I[i-numNodes][j-numNodes];
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


void populateNearOptimaMatricies(int solutionPeriod, int currentPeriod, 
	int numNodes, int numOfPressureSensors, int *sensorNodes)
{
	int i, j, k;
	i = j = k = 0;
	
	for (k = 0; k < numOfHours; k++)
	{
		for (i = 0; i < totalNodeCount; i++)
		{
			for (j = 0; j < totalNodeCount; j++)
			{
				largePressureMatrix[k][i][j] = 0;
				//largeA[i][j] = 0;
			}
		}
	}

	for (i = 0; i < totalNodeCount; i++)
	{
		for (j = 0; j < totalNodeCount; j++)
		{	
			largeA[currentPeriod][i][j] = 0;
		}
	}	
	
	for(i = 1; i <= numNodes; i++)
	{		
		oneLeak(i, deltas[solutionPeriod][i-1], numOfPressureSensors, i-1, currentPeriod);		
	}
	
	//Update A matrix
	
	for (k = 0; k < numOfHours; k++)
	{
		for(i = 0; i < numNodes; i++)
		{		
			for(j = 0; j < numNodes; j++)
			{
				largeA[currentPeriod][i][j] += (baseCasePressureMatrix[currentPeriod][k][i] - 
					largePressureMatrix[k][i][j]) / deltas[solutionPeriod][j]; // / delta;			
			}			
		}
	}
	
	for(i = 0; i < numNodes; i++)
	{		
		for(j = 0; j < numNodes; j++)
		{
			//if (deltas[j] != 0)
			//{
				largeA[currentPeriod][i][j] = largeA[currentPeriod][i][j] / (numOfHours); // / delta;
			//}
		}			
	}
	
	
	//Create A-hat
	for(i = 0; i < numNodes; i++)
	{
		for(j = 0; j < numNodes; j++)
		{
			Ahat[currentPeriod][i][j] =  largeA[currentPeriod][i][j];
		}
	}
	for(i = numNodes; i < (numNodes * 2); i++)
	{
		for(j = 0; j < numNodes; j++)
		{
			Ahat[currentPeriod][i][j] = -largeA[currentPeriod][i-numNodes][j];
		}
	}
	for(i = 0; i < numNodes; i++)
	{
		for(j = numNodes; j < (numNodes * 2); j++)
		{
			Ahat[currentPeriod][i][j] = -I[i][j-numNodes];
		}
	}
	for(i = numNodes; i < (numNodes * 2); i++)
	{
		for(j = numNodes; j < (numNodes * 2); j++)
		{
			Ahat[currentPeriod][i][j] = -I[i-numNodes][j-numNodes];
		}
	}

}


void randomizeLeaks(int numNodes, int numOfLeaks)
{	
	int i, j, k, compareResult;
	char name[20];
	
	i = j = k = 0;
	
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
						leakNodes[i] = (int)(rand()%numNodes);//+1;
					}while(leakNodes[i] == leakNodes[j-1]);										
				}
				j--;
			}
		}			

		for (i = 0; i < numOfLeaks; i++)
		{
			ENgetnodeid(leakNodes[i], name);
			for (j = 0; j < numOfNodesToIgnore; j++)
			{
				compareResult = strncmp(name, nodesToIgnore[j], 20);				
				if (compareResult == 0)
				{
					do
					{
						leakNodes[i] = (int)(rand()%numNodes)+1;
						ENgetnodeid(leakNodes[i], name);
						compareResult = strncmp(name, nodesToIgnore[j], 20);
					}while(compareResult == 0);
				}
			}
		}
		
		for(i = 0; i < numOfLeaks; i++)
		{
			leakMagnitudes[i] = drand48() * maxLeakSize;
			if (leakMagnitudes[i] < minLeakSize)
			{
				do
				{
					leakMagnitudes[i] = drand48() * maxLeakSize;
				}while (leakMagnitudes[i] < minLeakSize);
			}
		}		

}

int setNumOfSensors(double percentOfTotalNodes)
{
	int numOfSensors;
	
	numOfSensors = round(percentOfTotalNodes * totalNodeCount);
	
	return numOfSensors;
}


//FUNCTION
//Randomly determines the location of n number of pressure sensors in the network
void divinePressureSensorLocations(int numOfSensors)
{
	int i, j, k, compareResult;
	char name[20];
	
	i = j = k = 0;
	
	for (i = 0; i < numOfSensors; i++)
	{
		sensorNodes[i] = 0;
	}
	
	for (i = 0; i < totalNodeCount; i++)
	{
		isConstraintWorthy[i] = 0;
	}
	
	if (numOfSensors == totalNodeCount)
	{
		for (i = 0; i < totalNodeCount; i++)
		{
			sensorNodes[i] = i+1;
		}
	}
	
	else
	{	
		for (i = 0; i < numOfSensors; i++)
		{			
			sensorNodes[i] = (int)(rand()%totalNodeCount)+1;			
			j = i;
			
			while(j >= 1)
			{
				if (sensorNodes[i] == sensorNodes[j-1])
				{
					do
					{			
						sensorNodes[i] = (int)(rand()%totalNodeCount);//+1;
					}while(sensorNodes[i] == sensorNodes[j-1]);										
				}
				j--;
			}
		}			
    	
		for (i = 0; i < numOfSensors; i++)
		{
			ENgetnodeid(sensorNodes[i], name);
			for (j = 0; j < numOfNodesToIgnore; j++)
			{
				compareResult = strncmp(name, nodesToIgnore[j], 20);				
				if (compareResult == 0)
				{
					do
					{
						sensorNodes[i] = (int)(rand()%totalNodeCount)+1;
						ENgetnodeid(sensorNodes[i], name);
						compareResult = strncmp(name, nodesToIgnore[j], 20);
					}while(compareResult == 0);
				}
			}
		}				
	}
	qsort(sensorNodes, numOfSensors, sizeof(int), compare);
	
	for (i = 0; i < totalNodeCount; i++)
	{
		for (j = 0; j < numOfSensors; j++)
		{
			if (sensorNodes[j] == (i+1))
			{
				isConstraintWorthy[i] = 1;
				break;
			}
		}
		//printf("\n\tisConstraintWorthy[%d] = %d", i, isConstraintWorthy[i]);
	}
	//sortSensorLocations(numOfSensors);
	//getchar();
}

int compare(const void * a, const void * b)
{
	return ( *(int *)a - *(int *)b );
}
/*
void sortSensorLocations(int numOfSensors)
{
	int i, j, temp;
	int *tempSort;
	
	i = j = temp = 0;
	
	tempSort = (int *) calloc(numOfSensors, sizeof(int));
	
	temp = sensorNodes[0];
	tempSort[0] = temp;
	
	for (i = 0; i < numOfSensors - 1; i++)
	{
		for (j = 0; j < numOfSensors; j++)
		{
			if (sensorNodes[j] < tempSort[i])
			{
				temp = tempSort[i];
				tempSort[i] = sensorNodes[j];
				tempSort[i+1] = temp;
			}			
		}		 		
	}
	
	for (i = 0; i < numOfSensors; i++)
	{
		printf("\n\tsorted sensors[%d] = %d", i, tempSort[i]);
	}
	getchar();
	
	free(tempSort);
	
}
*/

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
void analyzeBaseCase(int sensorCount, int currentPeriod)
{		
	long t, tstep, hydraulicTimeStep, duration, warmUpExtension;
	float pressure, demand;
	int i, j, currentTime, compareResult;	
	char name[20];
	
	i = j = currentTime = 0;
	pressure = demand = 0.0;
	EPANETsimCounter++;
	warmUpExtension = (numOfHours * currentPeriod * 3600);
	
	for (i = 0; i < numOfHours; i++)
	{		
		for (j = 0; j < totalNodeCount; j++)
		{
			baseCasePressureMatrix[currentPeriod][i][j] = 0;
			baseCaseDemand[currentPeriod][j] = 0;
		}
	}
	
	//for (i=1; i <= sensorCount; i++)
	//{
		//baseCasePressureMatrix[i-1] = 0;
		
	//}
	
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
		if (t%hydraulicTimeStep == 0 && t >= (WARMUP_PERIOD + warmUpExtension)
			&& currentTime < numOfHours)
		{
			for (i=0; i < sensorCount; i++)
			{
				ENgetnodeid(sensorNodes[i], name);
				for (j = 0; j < numOfNodesToIgnore; j++)
				{					
					compareResult = strncmp(name, nodesToIgnore[j], 20);
					if (compareResult == 0)
					{								
						break;
					}
				}
				if (compareResult != 0)
				{
					ENgetnodevalue(sensorNodes[i], EN_PRESSURE, &pressure);						
					baseCasePressureMatrix[currentPeriod][currentTime][sensorNodes[i]-1] = pressure;	
				}
			}
			for (i = 1; i <= totalNodeCount; i++)
			{
				ENgetnodevalue(i, EN_DEMAND, &demand);
				baseCaseDemand[currentPeriod][i-1] += demand;				
			}
			currentTime++;
		}		
		ENnextH(&tstep);  	
	} while (tstep > 0); 
	
	//Close the hydraulic solver
	ENcloseH(); 
	
	for (i=1; i <= totalNodeCount; i++)
	{
		//baseCasePressureMatrix[i-1] = baseCasePressureMatrix[i-1] 
			// / numOfHours;
		baseCaseDemand[currentPeriod][i-1] = baseCaseDemand[currentPeriod][i-1] / numOfHours;
		//printf("\tnode %d base case demand = %f\n", i, baseCaseDemand[i-1]);
		
		//ENgetnodeid(i, name);
		//printf("base Case pressure node %s = %f\n", name, baseCasePressureMatrix[i-1]);
	}
	//getchar();
}

//FUNCTION
//Place a single leak at the index location in the network and run hydraulic analysis
//Determines how many pressure violations occur in the network by leak location
void oneLeak(int index, double emitterCoeff, int sensorCount, int columnNumber, int currentPeriod) 
{	
	int i, j, k, currentTime, compareResult, currentSensorNode;
	long t, tstep, hydraulicTimeStep, warmUpExtension;
	float pressure;
	char name[20];
	
	i = j = k = currentTime = currentSensorNode = 0;
	pressure = 0;
	EPANETsimCounter++;
	warmUpExtension = (numOfHours * currentPeriod * 3600);
	
	ENgettimeparam(EN_HYDSTEP, &hydraulicTimeStep);
	
	//Create the leak
	ENsetnodevalue(index, EN_EMITTER, emitterCoeff);
	
	ENopenH();  
	ENinitH(0);
	
	//Run the hydraulic analysis
	do {  	
		ENrunH(&t);		
		if (t%hydraulicTimeStep == 0 && t >= (WARMUP_PERIOD + warmUpExtension)
			&& currentTime < numOfHours)
		{		
			for (i = 0; i < sensorCount; i++)
			{	
				//printf("\n\nsensorNodes[%d] = %d, currentSensorNode = %d\n\n",i,sensorNodes[i],currentSensorNode);
				ENgetnodeid(sensorNodes[i], name);
				for (j = 0; j < numOfNodesToIgnore; j++)
				{
					compareResult = strncmp(name, nodesToIgnore[j], 20);
					if (compareResult == 0)
					{								
						break;
					}
				}
				
				if (compareResult != 0)
				{
					
					currentSensorNode = sensorNodes[i] - 1;
					
					ENgetnodevalue(sensorNodes[i], EN_PRESSURE, &pressure);
					
					largePressureMatrix[currentTime][currentSensorNode][columnNumber] = pressure;
					//printf("\n\n\tSEG FAULT TEST 9999\n");
				}
            }
            currentTime++;
         }
		ENnextH(&tstep); 
	} while (tstep > 0);
	
	//Close the hydraulic solver
	ENcloseH();
	
	//"Fix" the leak
	ENsetnodevalue(index, EN_EMITTER, 0.0);
	
	//for (i = 1; i <= sensorCount; i++)
	//{					
		//largePressureMatrix[i-1][columnNumber] = 
			//largePressureMatrix[i-1][columnNumber] / numOfHours;
		//ENgetnodeid(i, name);
		//printf("One leak Pressure @ node %s = %f\n", name, largePressureMatrix[i-1][columnNumber]);			
	//}
	//getchar();
}

//FUNCTION
//Generalized multi-leak simulator
void nLeaks(int leakCount, int sensorCount, int currentPeriod) 
{
	long t, tstep, hydraulicTimeStep, duration, warmUpExtension;	
	float pressure, baseDemand, demand;
	int i, j, currentTime, compareResult;
	char name[20];

	i = j = currentTime = 0;
	pressure = baseDemand = demand = 0.0;
	EPANETsimCounter++;
	warmUpExtension = (numOfHours * currentPeriod * 3600);
		
	for (i = 0; i < numOfHours; i++)
	{
		for (j = 0; j <= totalNodeCount; j++)
		{
			observedPressure[currentPeriod][i][j] = 0;
			observedDemand[currentPeriod][j] = 0;
		}
	}
	
	/*
	for (i=1; i <= sensorCount; i++)
	{
		observedPressure[i-1] = 0;
		observedDemand[i-1] = 0;
	}
	*/
	
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
		if (t%hydraulicTimeStep == 0 && t >= (WARMUP_PERIOD + warmUpExtension)
			&& currentTime < numOfHours)
		{
			for (i = 0; i < sensorCount; i++)
			{					
				ENgetnodeid(sensorNodes[i], name);
				for (j = 0; j < numOfNodesToIgnore; j++)
				{
					compareResult = strncmp(name, nodesToIgnore[j], 20);
					if (compareResult == 0)
					{								
						break;
					}	
				}	
				if (compareResult != 0)
				{
					ENgetnodevalue(sensorNodes[i], EN_PRESSURE, &pressure);																
					observedPressure[currentPeriod][currentTime][sensorNodes[i]-1] = pressure;								
					//printf("%f\t",demand);
				}
				
			}
			for (i = 1; i <= totalNodeCount; i++)
			{
				ENgetnodevalue(i, EN_DEMAND, &demand);
				observedDemand[currentPeriod][i-1] += demand;
			}
			
			//for (i = 0; i < leakCount; i++)
			//{
				//ENgetnodevalue(leakNodes[i], EN_BASEDEMAND, &baseDemand);					
				//ENgetnodevalue(leakNodes[i], EN_DEMAND, &demand);			
				//leakDemands[i] = (demand - baseDemand);
			//}
			currentTime++;
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
	
	for (i=1; i <= totalNodeCount; i++)
	{
		//observedPressure[i-1] = observedPressure[i-1] / numOfHours;
		observedDemand[currentPeriod][i-1] = observedDemand[currentPeriod][i-1] / numOfHours;
		//printf("\tnode %d observed demand = %f\n", i, observedDemand[i-1]);
		//ENgetnodeid(i, name);
		//printf("observed pressure @ node %s = %f\n", name, observedPressure[i-1]);
	}
	//getchar();
}

//FUNCTION
//Find the n highest leak magnitudes in the current solution
void findHighestMagnitudes(double *solutions)
{
	int i,j,k, placeHolder;
	
	for (i = 0; i < binaryLeakLimit; i++)
	{
		//leakGuesses[i] = 0.0;
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
				{
					leakGuesses[i] = solutions[j];				
				}
			}
		}
		//printf("\n\t\tleakGuesses[%d] = %f\t", i, leakGuesses[i]);
	}
	
}

void forgeMIPStartSolution(int currentPeriod, double sol[])
{
	int i;
	i = 0;
	
	for (i = 0; i < totalNodeCount; i++)
	{
		MIPStartSolution[currentPeriod][i] = 0;
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
			MIPStartSolution[currentPeriod][i] = 1;
	}
	
}

void calculateLeakDemand(int currentPeriod)
{
	int i, j;
	
	i = j = 0;
	

	for (i = 0; i < numOfLeaks; i++)
	{
		leakDemands[currentPeriod][i] = observedDemand[currentPeriod][leakNodes[i] - 1] - baseCaseDemand[currentPeriod][leakNodes[i] - 1];
		//printf("\t\tleakDemands[%d] = %f\n", leakNodes[i], leakDemands[i]);		
	}
		
	for (i = 0; i < totalNodeCount; i++)
	{
		totalDemand[currentPeriod] += observedDemand[currentPeriod][i];
	}
	
}

//FUNCTION
//Sum model error
double calculateError(int numNodes, int currentPeriod, double **solution)
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
		singleRunErrors[i] = fabs(realLeakValues[i] - solution[currentPeriod][i]);
		
		//printf("node %d model error %f \n", (i+1), singleRunErrors[i]);
		errorSum += singleRunErrors[i];
	}
	//printf("Total Error: %f \n", errorSum);
	//getchar();
	
	return errorSum;
}

//FUNCTION
//Creates the heat map results for a single 24 hour period combining objective
//weighted results from the specified number of sub-period solutions

//Nothing good can come from this function, silly and broken
int makeHeatMap(int currentRun)
{
	int i, j;
	double LPMax, MIPMax, LPWMax, MIPWMax;
	char sequentialFile[100], buffer[10], name[10];
	i = j = 0;	
	
	//LPMax = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	//MIPMax = (double *) calloc(numPeriodsPerSimulation, sizeof(double));
	
	//for (i = 0; i < numPeriodsPerSimulation; i++)
	//{
		//LPMax[i] = 0;
		//MIPMax[i] = 0;
	//}
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{		
		for (j = 0; j < totalNodeCount; j++)
		{			
			LPHeatMap[j] += LPSolutions[i][j];
			MIPHeatMap[j] += MIPSolutions[i][j];
			LPWeightedHeatMap[j] += LPSolutions[i][j] * LPobjectiveValues[currentRun][i];
			MIPWeightedHeatMap[j] += MIPSolutions[i][j] * MIPobjectiveValues[currentRun][i];
		}
	}
	
	for (j = 0; j < totalNodeCount; j++)
	{
		if (LPHeatMap[j] > LPMax)
			LPMax = LPHeatMap[j];
		if (MIPHeatMap[j] > MIPMax)
			MIPMax = MIPHeatMap[j];
		if (LPWeightedHeatMap[j] > LPWMax)
			LPWMax = LPWeightedHeatMap[j];
		if (MIPWeightedHeatMap[j] > MIPWMax)
			MIPWMax = MIPWeightedHeatMap[j];			
	}
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/HeatMap_");
	sprintf(buffer,"%d",currentRun);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;	
	fprintf(ptr_file, ",LP_HeatMap,LP_Weighted_HeatMap,MIP_HeatMap,MIP_Weighted_HeatMap\n");
	for(i = 0; i < (totalNodeCount); i++)
	{
		ENgetnodeid(i+1, name);		
		fprintf(ptr_file, "%s,", name);								
		fprintf(ptr_file, "%f,%f,%f,%f\n", (LPHeatMap[i]/LPMax),
			(LPWeightedHeatMap[i]/LPWMax), (MIPHeatMap[i]/MIPMax),
			(MIPWeightedHeatMap[i]/MIPWMax));								
	}
	
	fclose(ptr_file);
	/*
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/IndividualResults_");
	sprintf(buffer,"%d",currentRun);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, ",Objectives");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	fprintf(ptr_file, "\n,");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%f,", LPobjectiveValues[currentRun][i]);
	}
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%f,", MIPobjectiveValues[currentRun][i]);
	}	
	fprintf(ptr_file, "\n\n");
	
	fprintf(ptr_file, ",Solutions");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	fprintf(ptr_file, "\n");
	
	for (i = 0; i < totalNodeCount; i++)
	{
		ENgetnodeid(i+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPSolutions[j][i]);			
		}
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPSolutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
	}
	
	fclose(ptr_file);
	*/
	
	//free((void *) LPMax);
	//free((void *) MIPMax);
	
	return 0;	
}

int writeIndividualSolutions(int currentRun, double **solutions, char *method)
{	
	int i, j, timeBlock;
	char sequentialFile[150], buffer[10], name[10];
	
	i = j = timeBlock = 0;	
	
	//Create CSV file for each time period's solutions
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/R_");
	strcat(sequentialFile, method);
	strcat(sequentialFile, "_IndividualResults_");
	sprintf(buffer,"%d",currentRun);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	/*
	fprintf(ptr_file, ",LP_Objectives_");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	fprintf(ptr_file, "MIP_Objectives_");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	fprintf(ptr_file, "\n,");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%f,", LPobjectiveValues[currentRun][i]);
	}
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%f,", MIPobjectiveValues[currentRun][i]);
	}	
	fprintf(ptr_file, "\n\n");
	*/
	fprintf(ptr_file, ",");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d-", timeBlock);
		timeBlock = timeBlock + numOfHours;
		fprintf(ptr_file, "%d,", timeBlock);
	}
	fprintf(ptr_file, "\n");
	/*
	fprintf(ptr_file, "MIP_Solutions_");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d,", i);
	}
	fprintf(ptr_file, "\n");
	*/
	
	for (i = 0; i < totalNodeCount - 1; i++)
	{
		ENgetnodeid(i+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
		/*
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPSolutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
		*/
	}
	
	for (i = totalNodeCount - 1; i < totalNodeCount; i++)
	{
		ENgetnodeid(i+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}		
	}

	fclose(ptr_file);
	
	timeBlock = 0;
	
	//Create CSV file for each time period's solutions
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/P_");
	strcat(sequentialFile, method);
	strcat(sequentialFile, "_PressureError_");
	sprintf(buffer,"%d",currentRun);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, ",");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d-", timeBlock);
		timeBlock = timeBlock + numOfHours;
		fprintf(ptr_file, "%d,", timeBlock);
	}
	fprintf(ptr_file, "\n");
	
	for (i = totalNodeCount; i < (totalNodeCount * 2) - 1; i++)
	{
		ENgetnodeid((i-totalNodeCount)+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
	}
	
	for (i = (totalNodeCount * 2) - 1; i < (totalNodeCount * 2); i++)
	{
		ENgetnodeid((i-totalNodeCount)+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}		
	}

	fclose(ptr_file);
	

	//free((void *) LPMax);
	//free((void *) MIPMax);

	return 0;	

}

int writeNearOptimaSolutions(int currentRun, int currentSolution, double **solutions)
{	
	int i, j, timeBlock;
	char sequentialFile[200], buffer[10], name[10];
	
	i = j = timeBlock = 0;	
	
	//Create CSV file for each time period's solutions
	sequentialFile[0] = '\0';
	
	strcat(sequentialFile, globalDirName);
	
	strcat(sequentialFile, "/R_");	
	
	strcat(sequentialFile, "_NearOptimaResults_");
	
	sprintf(buffer,"%d",currentRun);		
	
	strcat(sequentialFile, buffer);
	
	strcat(sequentialFile, "_");	
	
	sprintf(buffer,"%d",currentSolution);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	

	
	fprintf(ptr_file, ",");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d-", timeBlock);
		timeBlock = timeBlock + numOfHours;
		fprintf(ptr_file, "%d,", timeBlock);
	}
	fprintf(ptr_file, "\n");
	
	
	for (i = 0; i < totalNodeCount - 1; i++)
	{
		ENgetnodeid(i+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
		/*
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPSolutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
		*/
	}
	printf("\n\n\n\t\t\tseg fault test No. 2\n");
	for (i = totalNodeCount - 1; i < totalNodeCount; i++)
	{
		ENgetnodeid(i+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}		
	}
	printf("\n\n\n\t\t\tseg fault test No. 3\n");
	fclose(ptr_file);
	
	timeBlock = 0;
	
	//Create CSV file for each time period's solutions
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/P_");	
	strcat(sequentialFile, "_PressureError_");
	sprintf(buffer,"%d",currentRun);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, "_");	
	strcat(sequentialFile, currentSolution);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, ",");
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "%d-", timeBlock);
		timeBlock = timeBlock + numOfHours;
		fprintf(ptr_file, "%d,", timeBlock);
	}
	fprintf(ptr_file, "\n");
	
	for (i = totalNodeCount; i < (totalNodeCount * 2) - 1; i++)
	{
		ENgetnodeid((i-totalNodeCount)+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}
		fprintf(ptr_file, "\n");
	}
	printf("\n\n\n\t\t\tseg fault test No. 4\n");
	for (i = (totalNodeCount * 2) - 1; i < (totalNodeCount * 2); i++)
	{
		ENgetnodeid((i-totalNodeCount)+1, name);		
		fprintf(ptr_file, "%s,", name);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", solutions[j][i]);			
		}		
	}
	printf("\n\n\n\t\t\tseg fault test No. 5\n");
	fclose(ptr_file);
	

	//free((void *) LPMax);
	//free((void *) MIPMax);

	return 0;	

}


int writeNearOptimaObjectives(int currentRun, double **objectives)
{	
	int i, j;
	char sequentialFile[150], buffer[10], name[10];
	
	i = j = 0;	
	
	//Create CSV file for each time period's solutions
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/E_NearOptimaObjectives_");	
	sprintf(buffer,"%d",currentRun);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, ",");
	for (i = 0; i < (numPeriodsPerSimulation); i++)
	{
		fprintf(ptr_file, "Objective at Hour %d,", i);
	}
	fprintf(ptr_file, "\n");
	
	for (i = 0; i < (numPeriodsPerSimulation); i++)
	{
		fprintf(ptr_file, "Solution from Hour %d,", i);
		for (j = 0; j < (numPeriodsPerSimulation ); j++)
		{	
			fprintf(ptr_file, "%f,", nearOptimaObjectiveValues[i][j]);
		}
		fprintf(ptr_file, "\n");
	}

	return 0;	

}
	

//FUNCTION
//Create an output file for each simulation/optimization run
int writeSummaryFile(int k, int optimstatus, int numOfSensors, double objval, double sol[])
{
	/*	
	char sequentialFile[150], buffer[10], name[10];
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
	
	fprintf(ptr_file, "No. Of Pressure Sensors:,%d \n",numOfSensors);
	fprintf(ptr_file, "Total Demand:, %f \n", totalDemand);
	fprintf(ptr_file, "Run #, %d, Model Error:, %f \n", (k + 1), modelError[k]);
	
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
	}c
	
	fclose(ptr_file);
	*/
	return 0;
}


//FUNCTION
//Create an output file for each simulation/optimization run
int writeRawResults(int k, int optimstatus, double sol[])
{	
	/*
	char sequentialFile[100], buffer[10], name[20];
	int i, j; //compareResult; 
	
	i = j = 0;	
	
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
		for(i = 0; i < (totalNodeCount); i++)
		{
			ENgetnodeid(i+1, name);
			//for (j = 0; j < numOfNodesToIgnore; j++)
			//{
				//printf("\nnodesToIgnore[%d]=%s,name=%sblob\n",j,nodesToIgnore[j],name);
				//compareResult = strncmp(name, nodesToIgnore[j], 20);
				//if (compareResult == 0)
				//{
					//printf("\n\n\n\t\t\t%sblob in raw\n\n\n", name);					
					//break;
				//}
			//}
			//if (compareResult != 0)
			//{
				//printf("in raw results name = %s \t and sol = %f\n", name, sol[i]);
				fprintf(ptr_file, "%s,", name);								
				fprintf(ptr_file, "%f,%f\n", sol[i], LPSolutions[i]); //((i+1) + (counter * totalNodeCount * 3)), sol[i + (counter * (totalNodeCount * 3))]);
			//}			
				
		}
	} else if (optimstatus == GRB_INF_OR_UNBD) 
	{
		fprintf(ptr_file, "Model is infeasible or unbounded\n");
	} else 
	{
		fprintf(ptr_file, "Optimization was stopped early\n");
	}
	
	fclose(ptr_file);
	*/
	return 0;
}


//FUNCTION
//Create an output file for each set of iterations
int writeErrorFile()
{	
	char sequentialFile[150];
	int i, j; 
	double averageLPObj, averageMIPObj, averageLPModelError, averageMIPModelError;
	i = j = 0;	
	averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
	
	/*
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/Error.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Run_#,,");
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "Time_%d,", i);
	}
	fprintf(ptr_file, "Average\n");
	
	//fprintf(ptr_file, "LP_Objective_Value,MIP_Objective_Value,Model_Error\n");
	
	for (i = 0; i < (iterations - 1); i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);
		fprintf(ptr_file, "LP_objective,");
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPobjectiveValues[i][j]);
			averageLPObj += LPobjectiveValues[i][j];
		}
		averageLPObj = averageLPObj / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n,LP_solution_error,", averageLPObj);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPmodelError[i][j]);
			averageLPModelError += LPmodelError[i][j];
		}
		averageLPModelError = averageLPModelError / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n,MIP_objective,", averageLPModelError);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPobjectiveValues[i][j]);	
			averageMIPObj += MIPobjectiveValues[i][j];
		}
		averageMIPObj = averageMIPObj / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n,MIP_solution_error,", averageMIPObj);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPmodelError[i][j]);
			averageMIPModelError += MIPmodelError[i][j];
		}
		averageMIPModelError = averageMIPModelError / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n", averageMIPModelError);		
	}
	for (i = (iterations - 1); i < iterations; i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);
		fprintf(ptr_file, "LP_objective,");
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPobjectiveValues[i][j]);
			averageLPObj += LPobjectiveValues[i][j];
		}
		fprintf(ptr_file, "%f\n,LP_solution_error,", averageLPObj);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPmodelError[i][j]);
			averageLPModelError += LPmodelError[i][j];
		}
		fprintf(ptr_file, "%f\n,MIP_objective,", averageLPModelError);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPobjectiveValues[i][j]);	
			averageMIPObj += MIPobjectiveValues[i][j];
		}
		fprintf(ptr_file, "%f\n,MIP_solution_error,", averageMIPObj);
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPmodelError[i][j]);
			averageMIPModelError += MIPmodelError[i][j];
		}
		fprintf(ptr_file, "%f", averageMIPModelError);						
	}
	
	fclose(ptr_file);
	*/
	
	//LP Objectives
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/E_LP_objectives.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Run_#,");
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "Time_%d,", i);
	}
	fprintf(ptr_file, "Average\n");
	
	//fprintf(ptr_file, "LP_Objective_Value,MIP_Objective_Value,Model_Error\n");
	
	for (i = 0; i < (iterations - 1); i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);		
		
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPobjectiveValues[i][j]);
			averageLPObj += LPobjectiveValues[i][j];
		}
		averageLPObj = averageLPObj / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n", averageLPObj);				
	}
	
	for (i = (iterations - 1); i < iterations; i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);		
		
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPobjectiveValues[i][j]);
			averageLPObj += LPobjectiveValues[i][j];
		}
		averageLPObj = averageLPObj / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f", averageLPObj);								
	}
	
	fclose(ptr_file);
	
	
	//LP Solution Error
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/E_LP_solution_error.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Run_#,");
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "Time_%d,", i);
	}
	fprintf(ptr_file, "Average\n");
		
	
	for (i = 0; i < (iterations - 1); i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);				
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPmodelError[i][j]);
			averageLPModelError += LPmodelError[i][j];
		}
		averageLPModelError = averageLPModelError / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n", averageLPModelError);				
	}
	
	for (i = (iterations - 1); i < iterations; i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);		
		
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", LPmodelError[i][j]);
			averageLPModelError += LPmodelError[i][j];
		}
		averageLPModelError = averageLPModelError / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f", averageLPModelError);								
	}
	
	fclose(ptr_file);
	
	
	//MIP Objectives
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/E_MIP_objectives.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Run_#,");
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "Time_%d,", i);
	}
	fprintf(ptr_file, "Average\n");		
	
	for (i = 0; i < (iterations - 1); i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);		
		
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPobjectiveValues[i][j]);	
			averageMIPObj += MIPobjectiveValues[i][j];
		}
		averageMIPObj = averageMIPObj / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n", averageMIPObj);				
	}
	
	for (i = (iterations - 1); i < iterations; i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);			
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPobjectiveValues[i][j]);	
			averageMIPObj += MIPobjectiveValues[i][j];
		}
		averageMIPObj = averageMIPObj / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f", averageMIPObj);								
	}
	
	fclose(ptr_file);
		
	
	//MIP Solution Error
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/E_MIP_solution_error.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Run_#,");
	
	for (i = 0; i < numPeriodsPerSimulation; i++)
	{
		fprintf(ptr_file, "Time_%d,", i);
	}
	fprintf(ptr_file, "Average\n");		
	
	for (i = 0; i < (iterations - 1); i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);
		
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPmodelError[i][j]);
			averageMIPModelError += MIPmodelError[i][j];
		}
		averageMIPModelError = averageMIPModelError / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f\n", averageMIPModelError);		
	}
	
	for (i = (iterations - 1); i < iterations; i++)
	{
		averageLPObj = averageMIPObj = averageLPModelError = averageMIPModelError = 0.0;
		fprintf(ptr_file, "%d,",i);	
		
		for (j = 0; j < numPeriodsPerSimulation; j++)
		{
			fprintf(ptr_file, "%f,", MIPmodelError[i][j]);
			averageMIPModelError += MIPmodelError[i][j];
		}
		averageMIPModelError = averageMIPModelError / numPeriodsPerSimulation;
		fprintf(ptr_file, "%f", averageMIPModelError);						
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

int writeAhat(int k, char *version)
{
	/*
	int i, j;
	char sequentialFile[100], buffer[10];
	i = j = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	//strcat(sequentialFile, "/home/andrew/Ubuntu One/Research/Thesis_Results/");
	//strcat(sequentialFile, directoryString);
	strcat(sequentialFile, globalDirName);
	strcat(sequentialFile, "/Ahat_");
	sprintf(buffer,"%d",k);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, version);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;	
	
	for(i = 0; i < (totalNodeCount * 2); i++)
	{			
		for (j = 0; j < (totalNodeCount * 2); j++)
		{
			fprintf(ptr_file,"%f,", Ahat[i][j]);				
		}
		fprintf(ptr_file,"\n");
	}
	
	
	fclose(ptr_file);
	*/
	return 0;	
}

int setOutputDirectory()
{
	int status;
	char dirName[100], date[25], buffer[10];
	
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
	sprintf(buffer,"%d",numOfHours);
	strcat(dirName, "/");
	strcat(dirName, buffer);
	strcat(dirName, "hourPeriods");
	status = mkdir(dirName, S_IRWXU | S_IRWXG | S_IRWXO);
	strcpy(globalDirName, dirName);
	
	//printf("Global dirname \t %s\n\n\n", globalDirName);
	
	//if (status != 0)
	//{
		//printf("\nDirectory Creation Error\n");
		//return status;
	//}
	
	return status;
}
