#ifndef MBA_H
#define MBA_H

#include "opt.h"

typedef struct _MineField{
    int m; /* number of shrapnel pieces */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best individual */
    double best_fitness; /* value of best fitness */
    double alpha; /* reduction factor */
    double mi; /* exploration factor */
    gsl_matrix *xe; /* location of the exploded mine */
    gsl_matrix *x; /* location of the improved mine */
    gsl_matrix *d; /* distance of shrapnel pieces */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *x_fitness; /* fitness values of improved bomb */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}MineField;

/* Allocation and deallocation */
MineField *CreateMineField(int m, int n); /*It allocates the search space */
void DestroyMineField(MineField **P); /* It deallocates the search space */
MineField *ReadMineFieldFromFile(char *fileName); /* It creates a search space specified in a file */
MineField *CopyMineField(MineField *P); /* It copies an entire search space */

/* Auxiliary functions */
void CheckMineFieldLimits(MineField *P); /* it checks the limits of each decision variable */
void InitializeMineField(MineField *P); /* it initializes the search space */
void ShowMineField(MineField *P); /* It displays the search space's content */
void ShowMineFieldInformation(MineField *P); /* It displays the search space's main information */

/* Main algorithm */
double EvaluateMine(MineField *P, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg);
void EvaluateMineField(MineField *P, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates the entire MineField */
void runMBA(MineField *P, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Mine Blast Algorithm for function minimization */

#endif
