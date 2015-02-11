/* This implementation is based on the paper "Migration Birds Optimization: A new metaheuristic approach and its performance on quadratic assignment problem". Notice some parameters' names do not follow the paper's nomenclature, since we opted to standardize them in LibOPT. */

#ifndef MBO_H
#define MBO_H

#include "opt.h"

typedef struct _BirdFlock{
    int m; /* number of birds in the flock */
    int n; /* number of decision variables */
    int max_iterations; /* maximum number of iterations */
    int k; /* number of neighbours solutions to be considered */
    int X; /* number of neighbor solutions to be shared with the next solution */
    int M; /* number of tours, i.e., the number of iterations for the leader */
    int leader; /* index of leader bird  */
    int *left; /* indeces of the birds that are on the left of leader bird */
    int *right; /* indeces of the birds that are on the right of leader bird */
    gsl_matrix *x; /* possible solutions (birds) */
    gsl_vector *fitness; /* fitness value */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}BirdFlock;

/* Allocation and Deallocation  */
BirdFlock *CreateBirdFlock(int m, int n); /* It alocattes the bird flock */
void DestroyBirdFlock(BirdFlock **B); /* It desalocates the bird flock */
BirdFlock *ReadBirdFlockFromFile(char *fileName); /* It creates a bird flock specified in a file */

/* Auxiliary Functions */
void InitializeBirdFlock(BirdFlock *B); /* It initializes the flock of birds */
void ShowBirdFlock(BirdFlock *B); /* It displays the harmomy memory's content */
void ShowBirdFlockInformation(BirdFlock *B); /* It displays the harmomy memory's main information */
double EvaluateBird(BirdFlock *B, gsl_vector *x,  int bird_id, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates a bird solution */

/* Main Algoritm */
void ImproveLeaderSolution(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It improves the lead bird by evaluating its neighbours */

#endif
