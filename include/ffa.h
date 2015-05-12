#ifndef FFA_H
#define FFA_H

#include "opt.h"

typedef struct _FireflySwarm{
    int m; /* number of fireflies */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best firefly */
    double best_fitness; /* value of best fitness */
    double gamma; /* light absorption coefficient */
    double beta_0; /* attractiveness coefficient */
    double alpha; /* step size */
    gsl_matrix *x; /* position */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}FireflySwarm;

/* Allocation and deallocation */
FireflySwarm *CreateFireflySwarm(int m, int n); /*It allocates the search space */
void DestroyFireflySwarm(FireflySwarm **F); /* It deallocates the search space */
FireflySwarm *ReadFireflySwarmFromFile(char *fileName); /* It creates a search space specified in a file */
FireflySwarm *CopyFireflySwarm(FireflySwarm *F); /* It copies an entire search space */

/* Auxiliary functions */
void CheckFireflySwarmLimits(FireflySwarm *F); /* it checks the limits of each decision variable */
void InitializeFireflySwarm(FireflySwarm *F); /* it initializes the search space */
void ShowFireflySwarm(FireflySwarm *F); /* It displays the search space's content */
void ShowFireflySwarmInformation(FireflySwarm *F); /* It displays the search space's main information */

/* Main algorithm */
inline void UpdateFireflyPosition(FireflySwarm *F, int firefly_id); /* It updates the position of each firefly */
inline void UpdateBestFireflyPosition(FireflySwarm *F, int best_firefly_id); /* It updates the position of the best firefly */
void EvaluateFireflySwarm(FireflySwarm *F, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates all fireflies */
void runFFA(FireflySwarm *F, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Firefly Algorithm for function minimization */

#endif
