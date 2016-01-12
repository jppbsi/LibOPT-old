#ifndef CS_H
#define CS_H

#include "opt.h"

#define BETA 1.5
#define MEAN 0
#define VARIANCE 1

typedef struct _NestPopulation{
    int m; /* number of nests */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best nest */
    double best_fitness; /* value of best fitness */
    double alpha; /* step size */
    double alpha_min; /* lower bound for alpha */
    double alpha_max; /* upper bound for alpha */
    double p; /* probability */
    double p_min; /* lower bound for probability */
    double p_max; /* upper bound for probability */
    gsl_matrix *x; /* position */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}NestPopulation;

/* Allocation and deallocation */
NestPopulation *CreateNestPopulation(int m, int n); /*It allocates the search space */
void DestroyNestPopulation(NestPopulation **P); /* It deallocates the search space */
NestPopulation *ReadNestPopulationFromFile(char *fileName); /* It creates a search space specified in a file */
NestPopulation *CopyNestPopulation(NestPopulation *P); /* It copies an entire search space */

/* Auxiliary functions */
void CheckNestPopulationLimits(NestPopulation *P); /* it checks the limits of each decision variable */
void InitializeNestPopulation(NestPopulation *P); /* it initializes the search space */
void ShowNestPopulation(NestPopulation *P); /* It displays the search space's content */
void ShowNestPopulationInformation(NestPopulation *P); /* It displays the search space's main information */
int NestLossParameter(int size, float probability);
void LevyFlightNest(NestPopulation *P, double alpha, double sigma);
void SortingNestPopulation(NestPopulation *P);
gsl_vector *AllocGaussianVector(int n); /* It computes the Box-Muller Transform generating a vector whose elements are positives and belong to a normal distribution */

/* Main algorithm */
double EvaluateNest(NestPopulation *P, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg);
void EvaluateNestPopulation(NestPopulation *P, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates all nests */
void runCS(NestPopulation *P, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Cuckoo Search for function minimization */
void runACS(NestPopulation *P, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Adaptative Cuckoo Search for function minimization */

#endif
