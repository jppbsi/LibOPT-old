/* The implementation of the Bat Algorithm is based on the paper "Bat algorithm: a novel approach for global engineering optimization",
available at http://www.emeraldinsight.com/doi/abs/10.1108/02644401211235834 */

#ifndef BA_H
#define BA_H

#include "opt.h"

typedef struct _Bats{
    int m; /* number of bats */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best bat */
    double best_fitness; /* value of best fitness */
    double f_min; /* minimum frequency */
    double f_max; /* maximum frequency */
    double A_min; /* minimum loudness */
    double A_max; /* maximum loudness */
    double mean_A; /* mean loudness */
    double alpha; /* used to update the loudness */
    double gamma; /* used to update the pulse rate */
    gsl_matrix *x; /* position */
    gsl_matrix *v; /* velocity */
    gsl_vector *f; /* frequency */
    gsl_vector *r0; /* original pulse rate */
    gsl_vector *r; /* pulse rate */
    gsl_vector *A; /* loudness */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
    gsl_vector *step; /* step size for each decision variable */
}Bats;

/* Allocation and deallocation */
Bats *CreateBats(int m, int n); /*It allocates the search space */
void DestroyBats(Bats **B); /* It deallocates the search space */
Bats *ReadBatsFromFile(char *fileName); /* Tt creates a search space specified in a file */
Bats *CopyBats(Bats *B); /* It copies an entire search space */

/* Auxiliary functions */
void CheckBatsLimits(Bats *H); /* it checks the limits of each decision variable */
void InitializeBats(Bats *B); /* it initializes the search space */
void ShowBats(Bats *H); /* It displays the search space's content */
void ShowBatsInformation(Bats *H); /* It displays the search space's main information */

/* Main algorithm */
inline void SetBatFrequency(Bats *B, int bat_id); /* It sets the frequency of each bat */
inline void UpdateBatVelocity(Bats *B, int bat_id); /* It updates the velocity of each bat */
inline void UpdateBatTemporaryPosition(Bats *B, int bat_id, gsl_vector *tmp); /* It updates the position of each bat */
inline void GenerateLocalSolutionNearBest(Bats *B, int best, gsl_vector *tmp); /* It generates a local solution near the best solution */
double EvaluateNewSolution(gsl_vector *tmp, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates a new solution */
void EvaluateBats(Bats *B, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates all bats */
void LocalSearchAndUpdateBest(Bats *B, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It performs the local search and updates the best bat */
void UpdateLoudness(Bats *B); /* It updates the loudness */
void UpdatePulseRate(Bats *B, int t); /* It updates the pulse rate */
void runBA(Bats *B, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Bat Algorithm for function minimization */

/* Hybrid approaches */

#endif

