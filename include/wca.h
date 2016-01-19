#ifndef WCA_H
#define WCA_H

#include "opt.h"

#define MU 0.1

typedef struct _RainDropPopulation{
    int m; /* size of population */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best raindrop */
    double best_fitness; /* value of best fitness */
    double nsr; /* number of rivers */
    double dmax; /* learning rate */
    gsl_matrix *x; /* position */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}RainDropPopulation;

/* Allocation and deallocation */
RainDropPopulation *CreateRainDropPopulation(int m, int n); /*It allocates the search space */
void DestroyRainDropPopulation(RainDropPopulation **P); /* It deallocates the search space */
RainDropPopulation *ReadRainDropPopulationFromFile(char *fileName); /* It creates a search space specified in a file */
RainDropPopulation *CopyRainDropPopulation(RainDropPopulation *P); /* It copies an entire search space */

/* Auxiliary functions */
void CheckRainDropPopulationLimits(RainDropPopulation *P); /* it checks the limits of each decision variable */
void InitializeRainDropPopulation(RainDropPopulation *P); /* it initializes the search space */
void ShowRainDropPopulation(RainDropPopulation *P); /* It displays the search space's content */
void ShowRainDropPopulationInformation(RainDropPopulation *P); /* It displays the search space's main information */

/* Main algorithm */
double EvaluateRainDrop(RainDropPopulation *P, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates a raindrop */
void EvaluateRainDropPopulation(RainDropPopulation *P, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates the entire population */
void runWCA(RainDropPopulation *P, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Water Cycle Algorithm for function minimization */
void runERWCA(RainDropPopulation *P, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the ER-WCA for function minimization */

#endif
