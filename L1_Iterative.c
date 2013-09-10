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
double delta = 1, emitterCoeff = 1, binaryLeakLimit = 2.0;
int numOfLeaks = 2, iterations =100;
char inputFile[50] = "hanoi-1.inp";
char reportFile[50] = "hanoi.rpt";
char directoryString[50] = "L1_Iterative/";
//
//

char globalDirName[100];
int totalNodeCount;
int *leakNodes;
double totalDemand, averageDelta, averagePreviousDelta, bigM = 999999.99;
double *baseCasePressureMatrix, *observedPressure, *coefficients, *b, *bhat,
	*realLeakValues, *singleRunErrors, *leakDemands, *leakMagnitudes, 
	*modelError, **largePressureMatrix, **largeA, **Ahat,  **I, 
	*objectiveValues, *deltas, *previousDeltas, *leakGuesses; 

FILE *ptr_file;

void initializeArrays();
void populateMatricies(int);
void randomizeLeaks(int, int);
void printLeakInfo(int);
void analyzeBaseCase(int);
void oneLeak(int, double, int, int);
void nLeaks(int, int);
double calculateError(int, double[]);
int writeSummaryFile(int, int, double, double[]);
int writeRawResults(int, int, double[]);
int writeLeakFile(int);
int writeErrorFile();
int setOutputDirectory();

int main(int argc, char *argv[]) 
{
	GRBenv *env = NULL;
	GRBmodel *model = NULL;
	//leakDemands = NULL;
	//leakNodes = NULL;
	//leakMagnitudes = NULL;
	//double *leakMagnitudes; //*leakSolutions,
	int  i, j, k, l, numNodes, storage, counter, placeHolder;//, *leakNodes;
	double errorSum, previousObjectiveValue;
	
	//Randomize the leak locations, commented out will use the same seeding 
	//for each run
	//srand(time(NULL));
	
	i = j = k = l = numNodes = counter = placeHolder =  0;
	averageDelta = averagePreviousDelta = previousObjectiveValue = 0.0;
	//leakSolutions = malloc(sizeof(double) * iterations);
	leakDemands = (double *) calloc(numOfLeaks, sizeof(double));
	leakNodes = (int *) calloc(numOfLeaks,sizeof(int));
	leakMagnitudes = (double *) calloc(numOfLeaks,sizeof(double));
	modelError = (double *) calloc(iterations, sizeof(double));
	objectiveValues = (double *) calloc(iterations, sizeof(double));
	leakGuesses = (double *) calloc(numOfLeaks, sizeof(double));
	
	//Open EPANET & Input file
	ENopen("hanoi-1.inp","hanoi.rpt","");
	
	// Get the number of nodes
	ENgetcount(EN_NODECOUNT, &numNodes);
	ENgetcount(EN_TANKCOUNT, &storage);	
		 
	/* Create environment */
 	error = GRBloadenv(&env, "L1_Iterative_MIP.log");
 	if (error) goto QUIT;
		
 	deltas = (double *) calloc((numNodes-storage), sizeof(double));
	previousDeltas = (double *) calloc((numNodes-storage), sizeof(double));

	//Create observation	
	for (k = 0; k < iterations; k++)
	{
	
		delta = 1.0;
	
		for( i = 0; i < (numNodes-storage); i++)
		{
			deltas[i] = delta;
			previousDeltas[i] = delta;		
			//printf("deltas[%d] = %f \n", i, deltas[i]);
		}
		for( i = 0; i < numOfLeaks; i++)
		{
			leakGuesses[i] = 0.0;
		}
		
		initializeArrays();
		randomizeLeaks(emitterCoeff, numNodes, storage, numOfLeaks);
 		
		//First Run LP to generate a better leak magnitude guess, represented
		//by delta
		
		analyzeBaseCase(numNodes);				 			
		nLeaks(numOfLeaks, numNodes);										
		populateMatricies(numNodes);		
		
		/* Create an empty model */ 		
 		error = GRBnewmodel(env, &model, "L1Approx", 0, NULL, NULL, NULL, NULL, 
 			NULL);
 		if (error) goto QUIT;
 		 	
 		/* Add variables */
 		for (i = 0; i < 62; i++)
 		{
 			obj[i] = coefficients[i]; 			
 			vtype[i] = GRB_CONTINUOUS; 			
 		}
 		 				
		error = GRBaddvars(model, 62, 0, NULL, NULL, NULL, obj, NULL, NULL, 
			vtype, NULL);
		if (error) goto QUIT;
		
		/* Integrate new variables */		
		error = GRBupdatemodel(model);
		if (error) goto QUIT;
				
		/* First constraint: Ax <= b */						
		for (i = 0; i < 62; i++)
		{
			for (j = 0; j < 62; j++)
			{
				ind[j] = j;
				val[j] = Ahat[i][j];			
			}						
			error = GRBaddconstr(model, 62, ind, val, GRB_LESS_EQUAL, 
				bhat[i],NULL);			
			if (error) goto QUIT;
		}
		
		error = GRBoptimize(model);
		if (error) goto QUIT;
		
		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 62, sol);
		if (error) goto QUIT;
		
    	
		
		averageDelta = averagePreviousDelta = 0.0;
		
		for (j = 0; j < numOfLeaks; j++)
		{
			leakGuesses[j] = 0.0;
		}
		for (j = 0; j < numOfLeaks; j++)
		{	
			for (i = 0; i < 31; i++)
			{
				placeHolder = 0;				
				if(sol[i] > leakGuesses[j])
				{
					for(l = 0; l < numOfLeaks; l++)
					{
						if(sol[i] == leakGuesses[l])
							placeHolder++;
					}
					if (placeHolder == 0)
						leakGuesses[j] = sol[i];
				}
			}
		}
		//printf("\n\t\t\t\t\t FIRST Run # %d \n\n", counter);
		for (i = 0; i < 31; i++)
		{
			//printf("%f \t", sol[i]);
		}
		for (i = 0; i < numOfLeaks; i++)
		{
			averageDelta += leakGuesses[i];
			//printf("leakGuess[%d] = %f \n", i, leakGuesses[i]);
		}
		averageDelta = averageDelta / numOfLeaks;
		//printf("\t\t\taverage Delta: %f\n\n", averageDelta);
		for (i = 0; i < 31; i++)
		{
			previousDeltas[i] = deltas[i];
			averagePreviousDelta += previousDeltas[i];			
			deltas[i] = averageDelta;
			//printf("\t\t\t This is the new delta: %f \n", deltas[i]);
		}		
		averagePreviousDelta = averagePreviousDelta / 31;
		
		/* Free model */
		GRBfreemodel(model);
		
		
		//Create observation	
		//for (k = 0; k < iterations; k++)
		do
		{
			counter++;
			initializeArrays();
			//randomizeLeaks(emitterCoeff, numNodes, numOfLeaks);				
			analyzeBaseCase(numNodes);				 			
			nLeaks(numOfLeaks, numNodes);							
			//printLeakInfo(numOfLeaks);		
			populateMatricies(numNodes);		
		
			/* Create an empty model */ 		
 			error = GRBnewmodel(env, &model, "L1MIP", 0, NULL, NULL, NULL, NULL, 
 				NULL);
 			if (error) goto QUIT;
 			 	
 			/* Add variables */
 			for (i = 0; i < 62; i++)
 			{
 				obj[i] = coefficients[i]; 			
 				vtype[i] = GRB_CONTINUOUS; 			
 			}
 			
 			for (i = 62; i < 93; i++){
 				obj[i] = 0.0;
 				vtype[i] = GRB_BINARY;
 			}
 			 				
			error = GRBaddvars(model, 93, 0, NULL, NULL, NULL, obj, NULL, NULL, 
				vtype, NULL);
			if (error) goto QUIT;
			
			/* Integrate new variables */		
			error = GRBupdatemodel(model);
			if (error) goto QUIT;
					
			/* First constraint: Ax <= b */						
			for (i = 0; i < 62; i++)
			{
				for (j = 0; j < 62; j++)
				{
					ind[j] = j;
					val[j] = Ahat[i][j];			
				}						
				error = GRBaddconstr(model, 62, ind, val, GRB_LESS_EQUAL, 
					bhat[i],NULL);			
				if (error) goto QUIT;
			}
			
			//Leak magnitude - (binary * bigM) <= 0
			for (i = 62; i < 93; i++)
			{		
				ind[0] = (i - 62); 	ind[1] = i; 
				val[0] = 1.0; 		val[1] = -bigM ;
										
				error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL,0.0,NULL);
				if (error) goto QUIT;
			}
			
			// Limit sum of binaries to number of leaks searching for...		
			for (i = 62; i < 93; i++)
			{		
				ind[i-62] = i;
				val[i-62] = 1.0;
			}								
			error = GRBaddconstr(model, 31, ind, val, GRB_LESS_EQUAL,
				binaryLeakLimit,NULL);
			if (error) goto QUIT;
    	
			
			error = GRBoptimize(model);
			if (error) goto QUIT;
			
			/* Write model to 'L1Approx.lp' */		
			//error = GRBwrite(model, "L1_MIP.lp");
			//if (error) goto QUIT;
			
			//error = GRBwrite(model, "L1_MIP.sol");
			//if (error) goto QUIT;
			
			/* Capture solution information */		
			error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
			if (error) goto QUIT;
			
			previousObjectiveValue = objval;
			
			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;
			
			error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 93, sol);
			if (error) goto QUIT;
			
			//averageDelta = averagePreviousDelta = 0.0;
			for (i = 0; i < 31; i++)
			{
				//printf("%f \t", sol[i]);			
			}
			
			for (j = 0; j < numOfLeaks; j++)
			{
				leakGuesses[j] = 0.0;
			}
			for (j = 0; j < numOfLeaks; j++)
			{	
				for (i = 0; i < 31; i++)
				{
					placeHolder = 0;				
					if(sol[i] > leakGuesses[j])
					{
						for(l = 0; l < numOfLeaks; l++)
						{
							if(sol[i] == leakGuesses[l])
								placeHolder++;
						}
						if (placeHolder == 0)
							leakGuesses[j] = sol[i];
					}
				}
			}
			
			//printf("\n\t\t\t\t\t Run # %d \n\n", counter);		
			
			averagePreviousDelta = averageDelta;
			averageDelta = 0;
			
			for (i = 0; i < numOfLeaks; i++)
			{
				averageDelta += leakGuesses[i];
				//printf("leakGuess[%d] = %f \n", i, leakGuesses[i]);
			}
			
			averageDelta = averageDelta / numOfLeaks;
			
			//printf("\t\t\tprevious averageDelta %f \t\t average Delta: %f\n\n", 
			//	averagePreviousDelta, averageDelta);
							
			for (i = 0; i < 31; i++)
			{
				previousDeltas[i] = deltas[i];
				//averagePreviousDelta += previousDeltas[i];
				
				//July 30, 2013
				//Alternating methods for determining next iteration's deltas
				//deltas[i] = deltas[i] + 1;
				deltas[i] = averageDelta;
				
				//averageDelta += deltas[i];
			}
			
			//averagePreviousDelta = averagePreviousDelta / 31;
			//averageDelta = averageDelta / 31;
			
			//for (i = 0; i < 31; i++)
			//{
				//previousDeltas[i] = deltas[i];
				//averagePreviousDelta += previousDeltas[i];			
				//deltas[i] = averageDelta;			
			//}		
			//averagePreviousDelta = averagePreviousDelta / 31;
			
			/*
			printf("\nOptimization complete\n");
			if (optimstatus == GRB_OPTIMAL)
			{
				printf("Optimal objective: %.4e\n", objval);
				for(i = 0; i < 62; i++)
				{		  	
					printf("  sol[%d] = %2.4f \n", (i+1), sol[i]);
				}
				for(i = 62; i < 93; i++)
				{		  	
					printf("  sol[%d] = %1.0f \t binary for sol[%d] \n", (i+1), sol[i], (i-61));
				}
			} else if (optimstatus == GRB_INF_OR_UNBD) 
			{
				printf("Model is infeasible or unbounded\n");
			} else 
			{
			  printf("Optimization was stopped early\n");
			}
			*/
			
			objectiveValues[k] = objval;
			modelError[k] = calculateError(numNodes, sol);
			
			/* Free model */
			GRBfreemodel(model);
			
		}while((objval - previousObjectiveValue) < 0); //while(fabs(averageDelta-averagePreviousDelta) > 0.10 );//while(counter < 20);		
		
		for (i = 0; i < 31; i++)
			{			
				deltas[i] = sol[i];
			}
		
		
		do
		{
			counter++;
			initializeArrays();
			//randomizeLeaks(emitterCoeff, numNodes, numOfLeaks);				
			analyzeBaseCase(numNodes);				 			
			nLeaks(numOfLeaks, numNodes);							
			//printLeakInfo(numOfLeaks);		
			populateMatricies(numNodes);		
		
			/* Create an empty model */ 		
 			error = GRBnewmodel(env, &model, "L1MIP", 0, NULL, NULL, NULL, NULL, 
 				NULL);
 			if (error) goto QUIT;
 			 	
 			/* Add variables */
 			for (i = 0; i < 62; i++)
 			{
 				obj[i] = coefficients[i]; 			
 				vtype[i] = GRB_CONTINUOUS; 			
 			}
 			
 			for (i = 62; i < 93; i++){
 				obj[i] = 0.0;
 				vtype[i] = GRB_BINARY;
 			}
 			 				
			error = GRBaddvars(model, 93, 0, NULL, NULL, NULL, obj, NULL, NULL, 
				vtype, NULL);
			if (error) goto QUIT;
			
			/* Integrate new variables */		
			error = GRBupdatemodel(model);
			if (error) goto QUIT;
					
			/* First constraint: Ax <= b */						
			for (i = 0; i < 62; i++)
			{
				for (j = 0; j < 62; j++)
				{
					ind[j] = j;
					val[j] = Ahat[i][j];			
				}						
				error = GRBaddconstr(model, 62, ind, val, GRB_LESS_EQUAL, 
					bhat[i],NULL);			
				if (error) goto QUIT;
			}
			
			//Leak magnitude - (binary * bigM) <= 0
			for (i = 62; i < 93; i++)
			{		
				ind[0] = (i - 62); 	ind[1] = i; 
				val[0] = 1.0; 		val[1] = -bigM ;
										
				error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL,0.0,NULL);
				if (error) goto QUIT;
			}
			
			// Limit sum of binaries to number of leaks searching for...		
			for (i = 62; i < 93; i++)
			{		
				ind[i-62] = i;
				val[i-62] = 1.0;
			}								
			error = GRBaddconstr(model, 31, ind, val, GRB_LESS_EQUAL,
				binaryLeakLimit,NULL);
			if (error) goto QUIT;
    	
			
			error = GRBoptimize(model);
			if (error) goto QUIT;
			
			/* Write model to 'L1Approx.lp' */		
			//error = GRBwrite(model, "L1_MIP.lp");
			//if (error) goto QUIT;
			
			//error = GRBwrite(model, "L1_MIP.sol");
			//if (error) goto QUIT;
			
			/* Capture solution information */		
			error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
			if (error) goto QUIT;
			
			previousObjectiveValue = objval;
			
			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
			if (error) goto QUIT;
			
			error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 93, sol);
			if (error) goto QUIT;
			
			//averageDelta = averagePreviousDelta = 0.0;
			for (i = 0; i < 31; i++)
			{
				//printf("%f \t", sol[i]);			
			}
			
			for (j = 0; j < numOfLeaks; j++)
			{
				leakGuesses[j] = 0.0;
			}
			for (j = 0; j < numOfLeaks; j++)
			{	
				for (i = 0; i < 31; i++)
				{
					placeHolder = 0;				
					if(sol[i] > leakGuesses[j])
					{
						for(l = 0; l < numOfLeaks; l++)
						{
							if(sol[i] == leakGuesses[l])
								placeHolder++;
						}
						if (placeHolder == 0)
							leakGuesses[j] = sol[i];
					}
				}
			}
			
			//printf("\n\t\t\t\t\t Run # %d \n\n", counter);		
			
			averagePreviousDelta = averageDelta;
			averageDelta = 0;
			
			for (i = 0; i < numOfLeaks; i++)
			{
				averageDelta += leakGuesses[i];
				//printf("leakGuess[%d] = %f \n", i, leakGuesses[i]);
			}
			
			averageDelta = averageDelta / numOfLeaks;
			
			//printf("\t\t\tprevious averageDelta %f \t\t average Delta: %f\n\n", 
			//	averagePreviousDelta, averageDelta);
							
			for (i = 0; i < 31; i++)
			{	
				if  (sol[i] > 0.01)
				{
					deltas[i] = sol[i];
				}
			}
			
			objectiveValues[k] = objval;
			modelError[k] = calculateError(numNodes, sol);
			
			/* Free model */
			GRBfreemodel(model);
			
		}while((objval - previousObjectiveValue) < 0); //while(fabs(averageDelta-averagePreviousDelta) > 0.10 );//while(counter < 20);		
		
	
		
		//writeSummaryFile(k);
		writeRawResults(k);
		writeLeakFile(k);
		/*
		printf("\nOptimization complete\n");
		if (optimstatus == GRB_OPTIMAL)
		{
			printf("Optimal objective: %.4e\n", objval);
			for(i = 0; i < 31; i++)
			{		  	
				//printf("  delta[%d] = %f \n", (i+1), deltas[i]);
			}
			for(i = 0; i < 62; i++)
			{		  	
				printf("  sol[%d] = %2.4f \n", (i+1), sol[i]);
			}
			for(i = 62; i < 93; i++)
			{		  	
				printf("  sol[%d] = %1.0f \t binary for sol[%d] \n", (i+1), sol[i], (i-61));
			}
		} else if (optimstatus == GRB_INF_OR_UNBD) 
		{
			printf("Model is infeasible or unbounded\n");
		} else 
		{
		  printf("Optimization was stopped early\n");
		}
		
		printLeakInfo(numOfLeaks);
		*/
	}
	
	
	writeErrorFile();
	
	free(leakNodes);
	
	free(leakMagnitudes);
	
	free(leakDemands);
		
	free(modelError);
	free(objectiveValues);
	free(deltas);
	//printf("\n \t\t\t\t Seg Fault Errrror Test \n");
	free(previousDeltas);
	free(leakGuesses);
	
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
		
	ENclose();
	return 0;
}

//FUNCTION
//Initialze various arrays to be populated during simulation
void initializeArrays()
{
	
	int i, j;
	
	i = j = 0;
	
	//Array initialization
	for (i = 0; i < 33; i++)
	{
		xhat[i] = 0;
	}
	for (i = 0; i < 31; i++)
	{
		observedPressure[i] = 0;
		baseCasePressureMatrix[i] = 0;
		b[i] = 0;
		realLeakValues[i] = 0.0;
		singleRunErrors[i] = 0.0;
	}
	for (i = 0; i < 62; i++)
	{
		bhat[i] = 0;		
	}
	for (i = 0; i < 62; i++)
	{
		coefficients[i] = 0;		
	}
	for (i = 0; i < 31; i++)
	{
		for (j = 0; j < 2; j++)
		{
			pressureMatrix[i][j] = 0;
			A[i][j] = 0;
		}
	}
	for (i = 0; i < 31; i++)
	{
		for (j = 0; j < 31; j++)
		{
			largePressureMatrix[i][j] = 0;
			largeA[i][j] = 0;		
		}
	}
	for (i = 0; i < 62; i++)
	{
		for (j = 0; j < 62; j++)
		{
			Ahat[i][j] = 0;		
		}
	}
	for (i = 0; i < 31; i++)
	{
		for (j = 0; j < 31; j++)
		{
			I[i][j] = 0;
		}
	}	

	//Create Identity Matrix
	for(i = 0; i< 31; i++)
	{
		for (j = 0; j<31; j++)
		{
			if (i==j)
				I[i][j] = 1;
		}
	}
	
	//Create c-transpose
	for (i = 0; i < 31; i++)
	{
		coefficients[i] = 0;
	}
	for (i = 31; i < 62; i++)
	{
		coefficients[i] = 1;
	}	
}

//FUNCTION
//Populate array values for the L1 Approximation
//Also calls single leak simulations for each node in the network
void populateMatricies(int numNodes)
{
	int i, j;
	
	i = j = 0;
	
	//Update b matrix
	for (i = 0; i < 31; i++)
	{
		b[i] = (baseCasePressureMatrix[i] - observedPressure[i]);
		//printf("b[%d] = %f\n",i,b[i]);
	}
	
	//Create b-hat
	for (i = 0; i < 31; i++)
	{
		bhat[i] = b[i];
	}
	for (i = 31; i < 62; i++)
	{
		bhat[i] = -b[i-31];
	}
	  
	for(i = 1; i <= 31; i++)
	{		
		oneLeak(i, deltas[i-1], numNodes, i-1);		
	}
	
	//Update A matrix		
	for(i = 0; i < 31; i++)
	{		
		for(j = 0; j < 31; j++)
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
	for(i = 0; i < 31; i++)
	{
		for(j = 0; j< 31; j++)
		{
			Ahat[i][j] =  largeA[i][j];
		}
	}
	for(i = 31; i < 62; i++)
	{
		for(j = 0; j< 31; j++)
		{
			Ahat[i][j] = -largeA[i-31][j];
		}
	}
	for(i = 0; i < 31; i++)
	{
		for(j = 31; j< 62; j++)
		{
			Ahat[i][j] = -I[i][j-31];
		}
	}
	for(i = 31; i < 62; i++)
	{
		for(j = 31; j< 62; j++)
		{
			Ahat[i][j] = -I[i-31][j-31];
		}
	}
	for (i = 0; i < numOfLeaks; i++)
	{
		realLeakValues[(leakNodes[i]-1)] = leakMagnitudes[i];		
	}
}

void randomizeLeaks(double emitterCoeff, int numNodes, int storage,
	int numOfLeaks)
{	
	int i, j, randLimit;
	
	i = j = 0;
	randLimit = numNodes - storage;
	
	for(i = 0; i < numOfLeaks; i++)
		{
			leakDemands[i] = 0;
			leakNodes[i] = 0;
			leakMagnitudes[i] = 0.0;
		}

		for(i = 0; i < numOfLeaks; i++)
		{			
			leakNodes[i] = (int)(rand()%randLimit)+1;			
			if (leakNodes[i] == 32)
			{
				do
				{			
					leakNodes[i] = (int)(rand()%randLimit)+1;
				}while(leakNodes[i] == 32);
			}			
			j = i;
			while(j >= 1)
			{
				if (leakNodes[i] == leakNodes[j-1])
				{
					do
					{			
						leakNodes[i] = (int)(rand()%randLimit)+1;
					}while(leakNodes[i] == leakNodes[j-1]);										
				}
				j--;
			}
			//printf("\n\n\n\n %d \n\n\n", leakNodes[i]);
		}			
		
		for(i = 0; i < numOfLeaks; i++)
		{
			//Each leak node is now set with a unique random emitter coefficient
			leakMagnitudes[i] = (double)(rand()%10 +1);
			//leakMagnitudes[i] = emitterCoeff;
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
		printf("\n\n leak info Node %d: %d \t Node ID: %s \t Magnitude: %f \n\n",
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
	char name[20];
	
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
		for (i=1; i < nodeCount; i++)
		{
			ENgetnodevalue(i, EN_PRESSURE, &pressure);
			ENgetnodeid(i, name);
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
			for (i = 1; i < nodeCount; i++)
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
		for (i = 1; i < nodeCount; i++)
		{			
			ENgetnodevalue(i, EN_PRESSURE, &pressure);			
			ENgetnodevalue(i, EN_DEMAND, &demand);			
			//ENgetnodeid(i,name);			
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
//Runs the hydraulic analysis on the base case scenario
void modAnalyzeBaseCase(int leak1, int leak2, double emitterCoeff1, double emitterCoeff2, int nodeCount){
	
	int i = 0;
	long t, tstep, hydraulicTimeStep, duration;
	float pressure = 0.0;	
	char name[20];
	
	ENgettimeparam( EN_HYDSTEP, &hydraulicTimeStep );
	ENgettimeparam( EN_DURATION, &duration );
	
	//Open and initialize the hydraulic solver
	ENopenH();  
	ENinitH(0);  

	ENsetnodevalue(leak1, EN_EMITTER, emitterCoeff1);
	ENsetnodevalue(leak2, EN_EMITTER, emitterCoeff2);

	//Run the hydraulic solver one hydraulic time step at a time
	do {  
		
		ENrunH(&t);
		
		// Retrieve hydraulic results for time t
		for (i=1; i < nodeCount; i++){
			ENgetnodevalue(i, EN_PRESSURE, &pressure);
			ENgetnodeid(i, name);
			baseCasePressureMatrix[i-1] = pressure;

			//Test print statements
			//printf("%i", i);
			//printf("\t");
			//printf("%s", name);
			//printf("\t");
			//printf("%f", pressure);			
			//printf("\t");
			//printf("%f", baseCasePressureMatrix[i-1]);
			//printf("\n");
		}
		
		ENnextH(&tstep);  
	
	} while (tstep > 0); 
	
	ENsetnodevalue(leak1, EN_EMITTER, 0);
	ENsetnodevalue(leak2, EN_EMITTER, 0);
	//Close the hydraulic solver
	ENcloseH();  

}

//FUNCTION
//Sum model error
double calculateError(int numNodes, double sol[])
{
	double errorSum;
	int i;
	
	i = 0;
	errorSum = 0.0;
	
	for (i = 0; i < 31; i++)
	{
		singleRunErrors[i] = fabs(realLeakValues[i] - sol[i]);
		//printf("node %d model error %f \n", (i+1), modelError[i]);
		errorSum += singleRunErrors[i];
	}
	//printf("Total Error: %f \n", errorSum);
	
	return errorSum;
}

//FUNCTION
//Create an output file for each simulation/optimization run
int writeSummaryFile(int k)
{	
	char sequentialFile[100], buffer[10], name[10];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, "/home/andrew/Documents/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Solution_");
	for (i = 0; i < numOfLeaks; i++)
	{
		sprintf(buffer,"%d",leakNodes[i]);
		strcat(sequentialFile, buffer);
		strcat(sequentialFile, "_");
	}		
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
	fprintf(ptr_file, "Run #, %d, Model Error:, %f \n", k, modelError[k]);
	
	for (i = 0; i < numOfLeaks; i++)
	{
		fprintf(ptr_file, "Leak %d demand:, %f, Demand Fraction:, %f, %% \n",
			i, leakDemands[i], ((leakDemands[i]/totalDemand)* 100));
	}
	
	fprintf(ptr_file, "\nOptimization complete\n");
	if (optimstatus == GRB_OPTIMAL) 
	{
		fprintf(ptr_file, "Optimal objective:, %.4e\n", objval);
		for(i = 0; i < 62; i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f \n", (i+1), sol[i]);
		}
		for(i = 62; i < 93; i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f, binary for, sol[%d] \n",
				(i+1), sol[i], (i-61));
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
int writeRawResults(int k)
{	
	char sequentialFile[100], buffer[10];
	int i; 
	
	i = 0;	
	
	//Create summary CSV file for each set of leaks
	sequentialFile[0] = '\0';
	strcat(sequentialFile, "/home/andrew/Documents/Research/Thesis_Results/");
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
		for(i = 0; i < 61; i++)
		{		  	
			fprintf(ptr_file, "  sol[%d] =, %f \n", (i+1), sol[i]);
		}
		for(i = 61; i < 62; i++)
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
	strcat(sequentialFile, "/home/andrew/Documents/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Error.csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;
	
	fprintf(ptr_file, "Objective_Value, Model_Error\n");
	
	for (i = 0; i < iterations-1; i++)
	{
		fprintf(ptr_file, "%f, %f\n", objectiveValues[i], modelError[i]);										
	}
	for (i = iterations-1; i < iterations; i++)
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
	strcat(sequentialFile, "/home/andrew/Documents/Research/Thesis_Results/");
	strcat(sequentialFile, directoryString);
	strcat(sequentialFile, "/Leaks_");
	sprintf(buffer,"%d",k);
	strcat(sequentialFile, buffer);
	strcat(sequentialFile, ".csv");
	
	ptr_file = fopen(sequentialFile, "w");
	if (!ptr_file)
		return 1;	
	
	for(i = 0; i < numOfLeaks-1; i++)
	{			
		ENgetnodeid(leakNodes[i], name);
		fprintf(ptr_file,"leak %d, %d, Node ID, %s, Magnitude, %f \n",
			i, leakNodes[i], name, leakMagnitudes[i]);
	}
	for(i = numOfLeaks-1; i < numOfLeaks; i++)
	{			
		ENgetnodeid(leakNodes[i], name);
		fprintf(ptr_file,"leak %d, %d, Node ID, %s, Magnitude, %f",
			i, leakNodes[i], name, leakMagnitudes[i]);
	}
	
	fclose(ptr_file);
	return 0;	
}