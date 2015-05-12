#ifndef GA_H
#define GA_H

#include "opt.h"

typedef struct _Population{
    int m; /* number of individuals */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best individual */
    double best_fitness; /* value of best fitness */
    double prob_crossover; /* crossover probability */
    double prob_mutation; /* mutation probability */
    gsl_matrix *gene; /* chromosome */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}Population;

/* Allocation and deallocation */
Population *CreatePopulation(int m, int n); /*It allocates the search space */
void DestroyPopulation(Population **P); /* It deallocates the search space */
Population *ReadPopulationFromFile(char *fileName); /* It creates a search space specified in a file */
Population *CopyPopulation(Population *P); /* It copies an entire search space */

/* Auxiliary functions */
void CheckPopulationLimits(Population *P); /* it checks the limits of each decision variable */
void InitializePopulation(Population *P); /* it initializes the search space */
void ShowPopulation(Population *P); /* It displays the search space's content */
void ShowPopulationInformation(Population *P); /* It displays the search space's main information */

/* Main algorithm */
void EvaluatePopulation(Population *P, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates all individuals */
void runGA(Population *P, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Genetic Algorithm for function minimization */

#endif
