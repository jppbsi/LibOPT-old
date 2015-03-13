/* This implementation is based on the paper "Migration Birds Optimization: A new metaheuristic approach and its performance on quadratic assignment problem". Notice some parameters' names do not follow the paper's nomenclature, since we opted to standardize them in LibOPT. */

#ifndef MBO_H
#define MBO_H

#include "opt.h"


//COLOCAR UMA BREVE DESCRICAO EXPLICANDO COMO SAO AS ESTRUTURAS .. LEADER SEMPRE 0, LEFT SIDE TEM ceil((B->M-1)/2) ....
//ARRUMAR PARAMETROS - DESCRICAO
//VERIFICAR SE O PROBLEMA TEM LEAK QUANDO TEMOS ERROS DO USUARIO

typedef struct _BirdFlock{
    int m; /* number of birds in the flock */
    int n; /* number of decision variables */
    int max_iterations; /* maximum number of iterations */
    int k; /* number of neighbours solutions to be considered */
    int X; /* number of neighbor solutions to be shared with the next solution */
    int M; /* number of tours, i.e., the number of iterations for the leader */
	int leftSide; /* a flag to know what bird will be changed */
    gsl_vector **left; /* indeces of the birds that are on the left of leader bird */
    gsl_vector **right; /* indeces of the birds that are on the right of leader bird */
    gsl_matrix *x; /* possible solutions (birds) */
    gsl_vector *fitness; /* fitness value */
    gsl_matrix *nb_left; /* neighbor set of the current iteration regarding left birds */
    gsl_vector *nb_fitness_left; /* neighbors's fitness of the current iteration regarding left birds */
    gsl_matrix *nb_right; /* neighbor set of the current iteration regarding right birds */
    gsl_vector *nb_fitness_right; /* neighbors's fitness of the current iteration regarding right birds */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}BirdFlock;

/* Allocation and Deallocation  */
BirdFlock *CreateBirdFlock(int m, int n, int k); /* It alocattes the bird flock */
void DestroyBirdFlock(BirdFlock **B); /* It desalocates the bird flock */
BirdFlock *ReadBirdFlockFromFile(char *fileName); /* It creates a bird flock specified in a file */

/* Auxiliary Functions */
void InitializeBirdFlock(BirdFlock *B); /* It initializes the flock of birds */
void ShowBirdFlock(BirdFlock *B); /* It displays the harmomy memory's content */
void ShowBirdFlockInformation(BirdFlock *B); /* It displays the harmomy memory's main information */
double EvaluateBird(BirdFlock *B, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates a bird solution */
void EvaluateBirdFlock(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...); /* It evaluates a bird flock */

/* Main Algoritm */
void ImproveLeaderSolution(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...); /* It improves the lead bird by evaluating its neighbours */
void ImproveOtherSolutions(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...); /* It improves the other birds by evaluating its neighbours */
void ReplaceLeader(BirdFlock *B); /* It replaces the leader bird by the next bird of V formation */

#endif
